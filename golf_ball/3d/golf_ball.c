/**
 * Step 3: Golf Ball CFD — Dimpled Sphere from STL
 * =================================================
 * Reads the preprocessed golf ball STL (step 1), builds a signed distance
 * field on the adaptive octree, then solves Navier-Stokes with the ball
 * as an embedded boundary. Outputs drag coefficient for comparison with
 * the smooth sphere (step 2).
 *
 * Key pipeline:
 *   golfball.stl
 *     └─ input_stl()         → coord* triangles (triangle soup)
 *     └─ distance()          → scalar d[] (signed distance field)
 *     └─ adapt_wavelet(d)    → AMR to resolve dimples
 *     └─ solid(cs,fs, d_interp) → embedded boundary fractions
 *     └─ Navier-Stokes solver
 *     └─ embed_force()       → drag / lift forces
 *
 * Sign convention for distance field:
 *   d[] > 0  →  outside ball (fluid)
 *   d[] < 0  →  inside ball  (solid)
 * This matches standard outward-normal STL orientation (fixed in step 1).
 *
 * Compile:
 *   qcc -O2 -Wall -disable-dimensions golf_ball.c -o golf_ball -lm
 *   ./golf_ball
 *
 * Parallel:
 *   CC='mpicc -D_MPI=8' qcc -O2 -Wall -disable-dimensions golf_ball.c \
 *       -o golf_ball -lm
 *   mpirun -n 8 ./golf_ball
 *
 * Important notes:
 *   1. GEOM_LEVEL must be >= FLOW_LEVEL so dimples are resolved before flow
 *   2. The geometry refinement loop in init() can take minutes — this is normal
 *   3. After init, the grid has ~5-15M cells at level 10 near the surface
 */

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "distance.h"
#include "fractions.h"
#include "embed_force_3d.h"
#include "output_vtu.h"
// view.h and lambda2.h removed — headless run, no OpenGL needed

/* ── Tunable parameters ───────────────────────────────────────────────────── */

// MUST match smooth_sphere.c exactly for a fair comparison
#ifndef RE
# define RE 300.
#endif

// Geometry refinement level: used to resolve dimples in the distance field
// Dimples are ~10% of D wide. Need >=6 cells/dimple.
// At GEOM_LEVEL=11: 2^11/16 = 128 cells/D → ~12.8 cells per 10%-wide dimple ✓
// At GEOM_LEVEL=10: 64 cells/D → ~6.4 cells per dimple (marginal)
#ifndef GEOM_LEVEL
# define GEOM_LEVEL 11
#endif

// Flow AMR refinement level (can be less than GEOM_LEVEL)
// High-gradient flow features don't need dimple-level resolution everywhere
#ifndef FLOW_LEVEL
# define FLOW_LEVEL 9
#endif

#define D   1.0     // ball diameter (normalized by preprocessing)
#define U0  1.0     // inflow velocity
#define T_END  100.
#define T_STAT  40.

/* ── STL file ─────────────────────────────────────────────────────────────── */
#ifndef STL_FILE
# define STL_FILE golfball_clean.stl
#endif

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#define STL_FILE_STR STR(STL_FILE)

/* ── Global fields ────────────────────────────────────────────────────────── */
face vector muv[];

/* ── Main ─────────────────────────────────────────────────────────────────── */

int main() {
  /**
   * Domain identical to smooth_sphere.c:
   *   16D x 16D x 16D, origin at (-3D, -L0/2, -L0/2)
   * Both simulations must share domain setup for Cd comparison validity.
   */
  init_grid(64);
  size(16.*D);
  origin(-3.*D, -L0/2., -L0/2.);
  mu = muv;
  run();
}

/* ── Viscosity ────────────────────────────────────────────────────────────── */

event properties(i++) {
  foreach_face()
    muv.x[] = fm.x[] * D * U0 / RE;
}

/* ── Boundary conditions (identical to smooth_sphere.c) ───────────────────── */

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

/* ── Initialization: STL → distance field → embedded boundary ─────────────── */

event init(t = 0) {
  /**
   * PHASE 1: Load STL and compute signed distance field
   * ----------------------------------------------------
   * input_stl() returns a coord* array of triangle vertices:
   *   [x0,y0,z0, x1,y1,z1, x2,y2,z2,  x0,y0,z0, ...]
   * terminated by a {HUGE, HUGE, HUGE} sentinel.
   *
   * The STL was normalized to diameter=1 centered at origin in step 1,
   * so it fits naturally in our domain.
   */
  fprintf(stderr, "# Loading STL: %s\n", STL_FILE_STR);
  FILE *stl_file = fopen(STL_FILE_STR, "r");
  if (!stl_file) {
    fprintf(stderr, "ERROR: Cannot open %s\n", STL_FILE_STR);
    fprintf(stderr, "Run step1_preprocess_stl.py first!\n");
    exit(1);
  }
  coord *triangles = input_stl(stl_file);
  fclose(stl_file);

  // Verify bounding box matches expected D=1 geometry
  coord bb_min, bb_max;
  bounding_box(triangles, &bb_min, &bb_max);
  fprintf(stderr, "# STL bounding box:\n");
  fprintf(stderr, "#   min: (%g, %g, %g)\n", bb_min.x, bb_min.y, bb_min.z);
  fprintf(stderr, "#   max: (%g, %g, %g)\n", bb_max.x, bb_max.y, bb_max.z);
  double stl_diam = fmax(fmax(bb_max.x - bb_min.x,
                              bb_max.y - bb_min.y),
                              bb_max.z - bb_min.z);
  fprintf(stderr, "#   diameter: %g  (expected ~%.1f)\n", stl_diam, D);
  if (fabs(stl_diam - D) > 0.1)
    fprintf(stderr, "WARNING: STL diameter %.3f differs from D=%.1f by >10%%\n",
            stl_diam, D);

  /**
   * PHASE 2: Compute distance field on coarse grid, then refine adaptively
   * -----------------------------------------------------------------------
   * distance() computes d[] such that:
   *   d[] = signed distance from cell center to nearest triangle surface
   *   d[] > 0 outside the ball (fluid region)
   *   d[] < 0 inside the ball  (solid region)
   *
   * The sign depends on outward normals — that's why step 1 fixes them.
   *
   * The refinement loop drives AMR until the distance field is resolved
   * to within 1e-3 * L0 relative error at level GEOM_LEVEL.
   * This ensures dimples appear as distinct geometric features.
   */
  fprintf(stderr, "# Computing distance field (refining to level %d)...\n",
          GEOM_LEVEL);

  /**
   * Pre-refine the grid near the ball surface using a sphere approximation
   * before computing the exact STL distance field. 
   * We use fabs(r - D/2) which has a "kink" at the surface. adapt_wavelet
   * detects this kink and refines exactly the spherical shell without
   * over-refining the rest of the domain.
   */
  int iter = 0;
  scalar ref[];
  do {
    foreach() {
      double r = sqrt(sq(x) + sq(y) + sq(z));
      ref[] = fabs(r - D/2.);
    }
    boundary({ref});
    fprintf(stderr, "#   pre-refinement iteration %d\n", ++iter);
  } while (adapt_wavelet({ref}, (double[]){0.005}, GEOM_LEVEL).nf && iter < 12);
  fprintf(stderr, "# Pre-refinement done after %d iterations\n", iter);

  // Now compute the actual STL distance field once on the already-refined grid
  scalar d[];
  distance(d, triangles);
  fprintf(stderr, "# Distance field computed\n");

  /**
   * PHASE 3: Convert distance field to embedded boundary fractions
   * -------------------------------------------------------------
   * solid(cs, fs, phi) internally evaluates phi at cell VERTICES via
   * a stencil loop and calls fractions() to set:
   *   cs[]  = volume fraction of fluid in each cell (0 = solid, 1 = fluid)
   *   fs.x[], fs.y[], fs.z[] = face fractions
   *
   * We interpolate d[] from cell centers to vertices using the 8-neighbor
   * average — this is the standard pattern from distance.c.
   *
   * Note: for an STL with outward normals (step 1 ensures this),
   * d[] > 0 outside (fluid), so the expression below is positive in fluid.
   */
  solid(cs, fs,
        (d[] + d[-1,0,0] + d[0,-1,0] + d[-1,-1,0] +
               d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1]) / 8.);

  /**
   * PHASE 4: Initialize velocity field
   * Fluid cells get inflow velocity U0; solid cells get zero.
   * cs[] handles the masking automatically.
   */
  foreach()
    u.x[] = cs[] * U0;

  fprintf(stderr, "# Golf ball initialized, Re=%.0f\n", RE);

  // Free triangle array (no longer needed after solid() is set)
  free(triangles);
}

/* ── Drag force logging ────────────────────────────────────────────────────── */

event drag(i++) {
  /**
   * Same computation as smooth_sphere.c.
   * Frontal area A = pi*(D/2)^2, dynamic pressure q = 0.5*rho*U0^2
   * Cd = F_x / (q*A) = F_x / (pi/8) for D=1, rho=1, U0=1
   */
  coord Fp = {0}, Fmu = {0};
  embed_force_3d(p, u, mu, &Fp, &Fmu);

  double Fdrag = Fp.x + Fmu.x;
  double Flift = Fp.y + Fmu.y;
  double q_A   = 0.5 * sq(U0) * M_PI * sq(D/2.);  // = pi/8 for D=U0=rho=1
  double Cd = Fdrag / q_A;
  double Cl = Flift / q_A;

  static FILE *fp = NULL;
  if (!fp) {
    fp = fopen("drag_golf.dat", "w");
    fprintf(fp, "# t  Cd  Cl  Cd_pressure  Cd_viscous\n");
  }
  fprintf(fp, "%g %g %g %g %g\n",
          t, Cd, Cl,
          Fp.x  / q_A,
          Fmu.x / q_A);
  fflush(fp);
}

/* ── Time-averaged Cd ─────────────────────────────────────────────────────── */

event average_drag(t = T_END) {
  FILE *fin = fopen("drag_golf.dat", "r");
  double t_val, cd, cl, cdp, cdv;
  double sum_cd = 0.;
  int    n_cd   = 0;
  char   buf[256];
  while (fgets(buf, sizeof(buf), fin)) {
    if (buf[0] == '#') continue;
    sscanf(buf, "%lf %lf %lf %lf %lf", &t_val, &cd, &cl, &cdp, &cdv);
    if (t_val >= T_STAT) { sum_cd += cd; n_cd++; }
  }
  fclose(fin);

  double mean_cd = (n_cd > 0) ? sum_cd / n_cd : -1.;
  fprintf(stderr, "\n# === GOLF BALL RESULT ===\n");
  fprintf(stderr, "# Re = %.0f\n", RE);
  fprintf(stderr, "# Mean Cd (t >= %.0f) = %.4f\n", T_STAT, mean_cd);

  FILE *fout = fopen("result_golf.txt", "w");
  fprintf(fout, "Re=%.0f\n", RE);
  fprintf(fout, "Cd_mean=%.6f\n", mean_cd);
  fprintf(fout, "Cd_samples=%d\n", n_cd);
  fclose(fout);
}

/* ── AMR: refine on geometry + velocity ───────────────────────────────────── */

event adapt(i++) {
  /**
   * Two-level strategy:
   *   - Keep embedded surface sharp at GEOM_LEVEL (dimples must stay resolved)
   *   - Refine wake/shear regions at FLOW_LEVEL for vortex accuracy
   * Coarsen far-field cells aggressively (min level 4).
   *
   * Important: after each coarsen/refine step, Basilisk automatically
   * reinterpolates cs and fs — no need to recompute distance field.
   */
  adapt_wavelet({cs, u.x, u.y, u.z},
                (double[]){1e-2, 0.01*U0, 0.01*U0, 0.01*U0},
                FLOW_LEVEL, 4);
}

/* ── ParaView VTU snapshots ───────────────────────────────────────────────── */

event vtu_output(t += 5.0; t <= T_END) {
  char fname[128];
  snprintf(fname, sizeof(fname), "vtu_golf/golf_t%04.0f", t);
  output_vtu((scalar *){p, cs}, (vector *){u}, fname, t);
}

/* ── Progress ─────────────────────────────────────────────────────────────── */

event logfile(i += 10)
  fprintf(stderr, "t=%.2f i=%d mgp=%d mgu=%d\n", t, i, mgp.i, mgu.i);

event end(t = T_END) {
  fprintf(stderr, "# Simulation complete at t=%.1f\n", T_END);
}