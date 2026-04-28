/**
 * 2D Cylinder: Smooth vs Bumpy (Corrected)
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "utils.h"

// Note: Ensure output_vtu.h is in your local directory
#include "output_vtu.h" 

double Re   = 11000.0;
double D    = 0.125;
double U0   = 1.0;
int maxLevel = 11;

face vector mu_v[];

/* Boundary conditions */
u.n[left]   = dirichlet(U0);
p[left]     = neumann(0.0);
pf[left]    = neumann(0.0);

u.n[right]  = neumann(0.0);
p[right]    = dirichlet(0.0);
pf[right]   = dirichlet(0.0);

/* Top and Bottom: Outflow/Symmetry to prevent blockage */
u.n[top]    = neumann(0.0);
u.n[bottom] = neumann(0.0);
u.t[top]    = neumann(0.0);
u.t[bottom] = neumann(0.0);

/* Embedded Boundary: No-slip */
u.n[embed] = dirichlet(0.0);
u.t[embed] = dirichlet(0.0);

int main(int argc, char * argv[]) {
  if (argc > 1) {
    Re = atof(argv[1]);
  }
  fprintf(stderr, "Running simulation with Re = %g\n", Re);

  L0 = 10.0;
  // Place cylinder at (0,0). Inlet is at -2.5 (20 diameters away)
  origin(-2.5, -L0/2.0); 
  N = 512;
  mu = mu_v;
  run();
}

event properties(i++) {
  foreach_face()
    mu_v.x[] = (fm.x[] * D * U0) / Re;
}

event init(t = 0) {
  double R = D/2.0;

#ifdef BUMPY
  double amp     = 0.05; // 5% roughness
  int    n_bumps = 24;
  /* Correct logic: Result must be > 0 for fluid, < 0 for solid */
  solid(cs, fs, sq(x) + sq(y) - sq(R * (1.0 + amp * cos(n_bumps * atan2(y, x)))));
#else
  solid(cs, fs, sq(x) + sq(y) - sq(R));
#endif

  foreach() {
    u.x[] = cs[] * U0;
    u.y[] = 0.0;
  }
}

event drag(i++) {
  coord Fp, Fmu;
  embed_force(p, u, mu, &Fp, &Fmu);

  // Reference for 2D Drag Coefficient: Cd = F / (0.5 * rho * U^2 * D)
  double ref = 0.5 * sq(U0) * D;
  double Cd  = (Fp.x + Fmu.x) / ref;

  static FILE *out = NULL;
  if (!out) {
    char dat_name[100];
#ifdef BUMPY
    sprintf(dat_name, "drag_bumpy_Re%g.dat", Re);
#else
    sprintf(dat_name, "drag_smooth_Re%g.dat", Re);
#endif
    out = fopen(dat_name, "w");
  }

  if (t > 0) {
    fprintf(out, "%g %g %g %g\n", t, Cd, Fp.x/ref, Fmu.x/ref);
    fflush(out);
  }
}

event movies(i += 10; t <= 50.0) {
  scalar omega[];
  vorticity(u, omega);

  char mp4_name[100];
#ifdef BUMPY
  sprintf(mp4_name, "vort_bumpy_Re%g.mp4", Re);
#else
  sprintf(mp4_name, "vort_smooth_Re%g.mp4", Re);
#endif

  output_ppm(omega, file = mp4_name,
    box = {{-1.0, -1.5}, {6.0, 1.5}},
    min = -5, max = 5,
    linear = true,
    n = 1024,
    mask = cs
  );
}

event vtk_output (t += 0.5) {
  char name[100];
#ifdef BUMPY
  sprintf(name, "bumpy_Re%g-%g.vtu", Re, t);
#else
  sprintf(name, "smooth_Re%g-%g.vtu", Re, t);
#endif

  scalar omega[];
  vorticity(u, omega);
  output_vtu((scalar *){cs, p, omega}, (vector *){u}, name, t);
}

event adapt(i++) {
  adapt_wavelet({cs, u}, (double[]){1e-4, 1e-3, 1e-3}, maxLevel, 5);
}

event logfile(i++) {
  if (i % 10 == 0)
    fprintf(stderr, "Step: %d, Time: %g\n", i, t);
}