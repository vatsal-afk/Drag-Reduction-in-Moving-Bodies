/**
 * Golf Ball Aerodynamics (2D Planar) — Basilisk C
 * =================================================
 * Compares smooth ball vs dimpled ball.
 *
 * Re is set to 1000 by default — the laminar NS solver has no turbulence
 * model, so running at Re=110,000 gives unphysical separation and Cd.
 * At Re=1000 the laminar solution is valid and Cd converges to a sensible
 * value (~0.4-0.6 for a circle). You can raise Re via the GUI slider to
 * explore higher Re, but Cd will drift from experiment above ~Re=5000.
 *
 * KEY FIX vs previous version:
 *   Drag now computed via embed_force() — Basilisk's built-in function
 *   that handles normal direction conventions and viscous contribution
 *   correctly. Manual embed_geometry integration had a sign error causing
 *   the ~100x under-prediction of Cd.
 *
 * Compile and run:
 *   Smooth:  qcc -O2 -Wall -o smooth  golfball.c -lm && ./smooth
 *   Dimpled: qcc -O2 -Wall -DDIMPLED -o dimpled golfball.c -lm && ./dimpled
 *
 * Outputs:
 *   vort_{smooth,dimpled}.mp4    — vorticity field
 *   tracer_{smooth,dimpled}.mp4  — passive scalar (top/bottom split)
 *   stderr                       — t, Cd, Cl, solver iterations
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"

// ── Physical parameters ──────────────────────────────────────────────────────

double Re       = 1000.0;  // Laminar NS is valid up to ~Re=5000 for a circle.
                           // Raise via GUI slider to observe wake transition.
double D        = 0.125;   // Ball diameter
double U0       = 1.0;     // Freestream velocity
int    maxLevel = 10;      // AMR max level (10 is sufficient; 11 is expensive)

// ── Dimple geometry ──────────────────────────────────────────────────────────
//
// 24 circular dimples cut evenly around the perimeter.
// Each dimple is a circle of radius DIMPLE_R whose centre is inset
// DIMPLE_DEPTH beneath the ball surface.
//
// Level-set CSG union: phi = max(phi_smooth, phi_dimple_k)
// max() → fluid if outside ball OR inside any dimple, solid otherwise.

#define N_DIMPLES    24
#define DIMPLE_R     (0.20 * D / 2.0)
#define DIMPLE_DEPTH (0.45 * DIMPLE_R)

face vector mu_v[];

// ── Boundary conditions ──────────────────────────────────────────────────────

u.n[left]  = dirichlet(U0);
u.n[right] = neumann(0.0);
p[left]    = neumann(0.0);
p[right]   = dirichlet(0.0);
pf[left]   = neumann(0.0);
pf[right]  = dirichlet(0.0);
u.n[embed] = dirichlet(0.0);
u.t[embed] = dirichlet(0.0);

// Tracer: f=1 for bottom half inflow, f=0 for top half.
// The mixing interface reveals the separation point and wake width.
scalar f[];
scalar *tracers = {f};
f[left]  = dirichlet(y < 0 ? 1.0 : 0.0);
f[right] = neumann(0.0);

// ── Level-set functions ──────────────────────────────────────────────────────

static double smooth_ball(double x, double y) {
    return sqrt(sq(x) + sq(y)) - D / 2.0;
}

static double dimpled_ball(double x, double y) {
    double R   = D / 2.0;
    double phi = smooth_ball(x, y);
    for (int k = 0; k < N_DIMPLES; k++) {
        double theta = k * 2.0 * M_PI / N_DIMPLES;
        double xc    = (R - DIMPLE_DEPTH) * cos(theta);
        double yc    = (R - DIMPLE_DEPTH) * sin(theta);
        phi = max(phi, DIMPLE_R - sqrt(sq(x - xc) + sq(y - yc)));
    }
    return phi;
}

// ── Viscosity ─────────────────────────────────────────────────────────────────

event properties(i++) {
    foreach_face()
        mu_v.x[] = fm.x[] * (D * U0) / Re;
}

// ── Initialisation ────────────────────────────────────────────────────────────

event init(t = 0) {
#ifdef DIMPLED
    solid(cs, fs, dimpled_ball(x, y));
#else
    solid(cs, fs, smooth_ball(x, y));
#endif
    foreach() {
        u.x[] = cs[] ? U0 : 0.0;
        u.y[] = 0.0;
        f[]   = cs[] ? (y < 0 ? 1.0 : 0.0) : 0.0;
    }
}

// ── Force and drag logging ────────────────────────────────────────────────────
//
// embed_force() computes total (pressure + viscous) force and torque on
// the embedded solid with correct sign conventions internally.
// Cd = Fx / (0.5 * rho * U0^2 * D)  with rho = 1
// Cl = Fy / (0.5 * rho * U0^2 * D)
// Reporting Cl alongside Cd reveals vortex shedding (Cl oscillates at
// the Strouhal frequency when shedding is established).

event logfile(i++) {
    coord Fp = {0}, Fv = {0};         // pressure force, viscous force
    embed_force(p, u, mu_v, &Fp, &Fv);
    double Cd = 2.0 * (Fp.x + Fv.x) / (sq(U0) * D);
    double Cl = 2.0 * (Fp.y + Fv.y) / (sq(U0) * D);
    fprintf(stderr,
            "i=%-6d  t=%.4f  Cd=%+.4f  Cl=%+.4f  mgp.i=%d  mgu.i=%d\n",
            i, t, Cd, Cl, mgp.i, mgu.i);
}

// ── Movie output ──────────────────────────────────────────────────────────────

event movies(i += 4; t <= 30.0) {
    scalar omega[], m[];
    vorticity(u, omega);
    foreach()
        m[] = cs[] - 0.5;

#ifdef DIMPLED
    output_ppm(omega, file = "vort_dimpled.mp4",
               box = {{-0.5, -0.5}, {7.5, 0.5}},
               min = -10, max = 10, linear = true, mask = m, n = 1500);
    output_ppm(f, file = "tracer_dimpled.mp4",
               box = {{-0.5, -0.5}, {7.5, 0.5}},
               min = 0, max = 1, linear = false, mask = m, n = 1500);
#else
    output_ppm(omega, file = "vort_smooth.mp4",
               box = {{-0.5, -0.5}, {7.5, 0.5}},
               min = -10, max = 10, linear = true, mask = m, n = 1500);
    output_ppm(f, file = "tracer_smooth.mp4",
               box = {{-0.5, -0.5}, {7.5, 0.5}},
               min = 0, max = 1, linear = false, mask = m, n = 1500);
#endif
}

// ── Adaptive mesh refinement ──────────────────────────────────────────────────

event adapt(i++) {
    adapt_wavelet({cs,   u.x,  u.y,  f},
                  (double[]){1e-3, 3e-2, 3e-2, 5e-2},
                  maxLevel, 4);
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    L0 = 8.0;
    origin(-0.5, -L0 / 2.0);
    N  = 256;
    mu = mu_v;

    display_control(Re,       100,  50000);
    display_control(maxLevel, 6,    12);

    run();
}sin