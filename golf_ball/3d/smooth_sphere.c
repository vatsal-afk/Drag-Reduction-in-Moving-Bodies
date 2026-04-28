/**
 * Step 2: Smooth Sphere CFD — Baseline for Drag Comparison
 * =========================================================
 * Solves incompressible Navier-Stokes around a smooth sphere using
 * Basilisk's embedded boundary method on an adaptive octree.
 *
 * This is the CONTROL case. The golf ball simulation (step3) uses
 * the exact same domain, Re, and flow setup. Differing only in geometry.
 *
 * Key outputs:
 *   drag_smooth.dat   — time series of Cd (pressure + viscous)
 *   smooth_sphere.mp4 — lambda2 vortex animation (if bview available)
 *
 * Compile & run:
 *   qcc -O2 -Wall -disable-dimensions smooth_sphere.c -o smooth_sphere -lm
 *   ./smooth_sphere
 *
 * Parallel (e.g. 8 cores):
 *   CC='mpicc -D_MPI=8' qcc -O2 -Wall -disable-dimensions smooth_sphere.c \
 *       -o smooth_sphere -lm
 *   mpirun -n 8 ./smooth_sphere
 *
 * Physical setup:
 *   D   = 1.0     (sphere diameter, normalized)
 *   U0  = 1.0     (inflow velocity)
 *   rho = 1.0     (fluid density)
 *   Re  = U0*D/nu (set via RE macro below, default 300)
 *   Cd  = F_drag / (0.5 * rho * U0^2 * pi * (D/2)^2)
 */

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "embed_force_3d.h"
#include "output_vtu.h"
// view.h and lambda2.h removed — headless run, no OpenGL needed

/* ── Tunable parameters ───────────────────────────────────────────────────── */

// Reynolds number: Re = U0*D/nu
// Re=300  → laminar wake, vortex shedding, fast to simulate
// Re=3700 → transitional (need higher maxlevel ~10)
// Re=1e5  → turbulent (use wall model, see wall-model.c example)
#ifndef RE
# define RE 300.
#endif

// Maximum AMR refinement level
// Level 8  → 256^3 equiv, fast, coarse dimple resolution
// Level 9  → 512^3 equiv, good for Re~300-1000
// Level 10 → 1024^3 equiv, needed for Re>3000
#ifndef MAXLEVEL
# define MAXLEVEL 8
#endif

// Sphere diameter (normalized)
#define D  1.0

// Inflow velocity
#define U0 1.0

// Simulation end time (non-dimensional: t * U0 / D)
#define T_END 100.

// When to start collecting drag statistics (let transient die out)
#define T_STAT 40.

/* ── Global fields ────────────────────────────────────────────────────────── */

// Dynamic viscosity field (face-centered, for embed compatibility)
face vector muv[];

/* ── Domain setup ─────────────────────────────────────────────────────────── */

int main() {
  /**
   * Domain: 16D x 16D x 16D box
   *   - Sphere center at origin (x=0, y=0, z=0)
   *   - Inflow at x = -3D  (3 diameters upstream)
   *   - Outflow at x = 13D (13 diameters downstream, room for wake)
   */
  init_grid(64);           // Start with 64^3 coarse grid
  size(16.*D);             // Domain edge length
  origin(-3.*D,            // x: start 3D upstream
         -L0/2.,           // y: centered
         -L0/2.);          // z: centered
  mu = muv;
  run();
}

/* ── Viscosity (set every timestep, handles AMR correctly) ─────────────────── */

event properties(i++) {
  foreach_face()
    muv.x[] = fm.x[] * D * U0 / RE;
  // fm.x[] is the face fraction from embed — zero inside solid
}

/* ── Boundary conditions ───────────────────────────────────────────────────── */

// LEFT face: uniform inflow
u.n[left] = dirichlet(U0);
p[left]   = neumann(0.);
pf[left]  = neumann(0.);

// RIGHT face: free outflow (zero pressure)
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// TOP, BOTTOM, FRONT, BACK: free-slip walls
// (default in Basilisk is no-penetration + free-slip, which is correct here)

// EMBEDDED boundary: no-slip on sphere surface
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

/* ── Initialization ────────────────────────────────────────────────────────── */

event init(t = 0) {
  /**
   * Geometry: unit sphere centered at origin.
   *   solid(cs, fs, phi) sets the embedded boundary where:
   *   phi > 0 → fluid (cs = 1)
   *   phi < 0 → solid (cs = 0)
   *   phi = 0 → interface (0 < cs < 1)
   *
   *   For a sphere: phi = R^2 - (x^2 + y^2 + z^2)
   *   i.e. positive inside sphere... wait — we want fluid OUTSIDE.
   *   So: phi = x^2 + y^2 + z^2 - R^2  → positive outside (fluid), negative inside (solid)
   *
   *   Note: Basilisk's solid() uses the convention that phi > 0 is FLUID.
   *   The expression below matches sphere.c from basilisk.fr exactly.
   */
  refine(x*x + y*y + z*z < sq(0.7*D) && level < MAXLEVEL);
  solid(cs, fs, sq(x) + sq(y) + sq(z) - sq(D/2.));

  // Set initial velocity: U0 everywhere outside the sphere
  foreach()
    u.x[] = cs[] * U0;

  fprintf(stderr, "# Smooth sphere initialized, Re=%.0f, maxlevel=%d\n",
          RE, MAXLEVEL);
}

/* ── Drag force logging ────────────────────────────────────────────────────── */

event drag(i++) {
  /**
   * embed_force() integrates pressure + viscous stress over the
   * embedded surface to give net force components.
   *
   * Fp  = pressure (form) drag
   * Fmu = viscous (skin friction) drag
   * Total drag in x-direction: Fp.x + Fmu.x
   *
   * Drag coefficient:
   *   Cd = F_drag / (q * A_frontal)
   *   q  = 0.5 * rho * U0^2 = 0.5  (rho=1, U0=1)
   *   A  = pi * (D/2)^2 = pi/4     (D=1)
   *   q*A = pi/8
   *   Cd = 8 * F_drag / pi
   */
  coord Fp = {0}, Fmu = {0};
  embed_force_3d(p, u, mu, &Fp, &Fmu);

  double Fdrag  = Fp.x + Fmu.x;
  double Flift  = Fp.y + Fmu.y;
  double Cd = Fdrag / (0.5 * sq(U0) * M_PI * sq(D/2.));
  double Cl = Flift / (0.5 * sq(U0) * M_PI * sq(D/2.));

  // Log to file: time, Cd, Cl, Fp.x (pressure drag), Fmu.x (viscous drag)
  static FILE *fp = NULL;
  if (!fp) {
    fp = fopen("drag_smooth.dat", "w");
    fprintf(fp, "# t  Cd  Cl  Cd_pressure  Cd_viscous\n");
  }
  fprintf(fp, "%g %g %g %g %g\n",
          t,
          Cd,
          Cl,
          Fp.x  / (0.5 * sq(U0) * M_PI * sq(D/2.)),
          Fmu.x / (0.5 * sq(U0) * M_PI * sq(D/2.)));
  fflush(fp);
}

/* ── Time-averaged Cd (written once at end) ────────────────────────────────── */

event average_drag(t = T_END) {
  // Re-read drag_smooth.dat and compute mean Cd from T_STAT onward
  FILE *fin = fopen("drag_smooth.dat", "r");
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
  fprintf(stderr, "\n# === SMOOTH SPHERE RESULT ===\n");
  fprintf(stderr, "# Re = %.0f\n", RE);
  fprintf(stderr, "# Mean Cd (t >= %.0f) = %.4f\n", T_STAT, mean_cd);
  fprintf(stderr, "# (Literature: Re=300 → Cd ≈ 0.65, Re=1e5 → Cd ≈ 0.47)\n");

  // Write summary
  FILE *fout = fopen("result_smooth.txt", "w");
  fprintf(fout, "Re=%.0f\n", RE);
  fprintf(fout, "Cd_mean=%.6f\n", mean_cd);
  fprintf(fout, "Cd_samples=%d\n", n_cd);
  fclose(fout);
}

/* ── AMR: refine on velocity gradients + geometry ─────────────────────────── */

event adapt(i++) {
  /**
   * Adapt on:
   *   cs  — geometry, keep embedded surface sharp (threshold 1e-2)
   *   u.x, u.y, u.z — velocity components (threshold 0.01 * U0 for turbulent wake)
   * Minimum refinement level 4 keeps coarse cells far from ball.
   */
  adapt_wavelet({cs, u.x, u.y, u.z},
                (double[]){1e-2, 0.01*U0, 0.01*U0, 0.01*U0},
                MAXLEVEL, 4);
}

/* ── ParaView VTU snapshots ───────────────────────────────────────────────── */

event vtu_output(t += 5.0; t <= T_END) {
  char fname[128];
  snprintf(fname, sizeof(fname), "vtu_smooth/smooth_t%04.0f", t);
  output_vtu((scalar *){p, cs}, (vector *){u}, fname, t);
}

/* ── Progress logging ────────────────────────────────────────────────────────  */

event logfile(i += 10)
  fprintf(stderr, "t=%.2f i=%d mgp=%d mgu=%d\n", t, i, mgp.i, mgu.i);

event end(t = T_END) {
  fprintf(stderr, "# Simulation complete at t=%.1f\n", T_END);
}