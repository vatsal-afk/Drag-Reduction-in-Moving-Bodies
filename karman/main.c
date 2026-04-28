#include "embed.h"
#include "navier-stokes/centered.h" // computes fluxes across cell faces
#include "tracer.h"

double Re = 160000.0; // Reynold's number
int maxLevel = 9;
face vector mu_v[]; // dynamic viscosity field (stored on cell face)

double D = 0.125, U0 = 2.0; // cylinder diameter, inflow velocity
// boundary conditions

// velocity vector field
u.n[left] = dirichlet(U0); // n is normal component
u.n[right] = neumann(0.0);

// pressure
p[left] = neumann(0.0);
p[right] = dirichlet(0.0);

// face pressure used in projection step
pf[left] = neumann(0.0);
pf[right] = dirichlet(0.0);

// cylinder is no-slip
u.n[embed] = dirichlet(0.0);
u.t[embed] = dirichlet(0.0);

// tracer
scalar f[];
scalar *tracers = {f};
// boundary conditions
f[left] = dirichlet(y < 0);
f[right] = neumann(0.0);

event properties (i++) {
	foreach_face() {
		mu_v.x[] = (fm.x[]*D*U0) / Re; // fm is fraction of each face open to fluid
	}	
}

event init (t = 0) {
	solid (cs, fs, sqrt(sq(x) + sq(y)) - D/2.0); // cs is cell fraction field, for solid cs = 0, liquid cs > 0
	
	foreach() {
    		u.x[] = cs[] ? U0 : 0.0; // velocity is 0 for solid cell
    	}
}

event logfile (i++) {
	fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}
  
event movies (i += 4; t <= 15.0) {
	scalar omega[], m[]; // vorticity, visualization mask fields
	vorticity (u, omega);
	foreach() {
		m[] = cs[] - 0.5; // mask so solid regions are hidden in plots
	}
    	
	output_ppm (omega, file = "vorticity_field.mp4", box = {{-0.5, -0.5},{7.5, 0.5}}, min = -10, max = 10, linear = true, mask = m);
	output_ppm (f, file = "tracer_field.mp4", box = {{-0.5, -0.5},{7.5, 0.5}}, linear = false, min = 0, max = 1, mask = m);
}

// refine the grid based on wavelet error estimation
event adapt (i++) {
	adapt_wavelet ({cs, u, f}, {1e-2, 3e-2, 3e-2, 3e-2}, maxLevel, 4);
}

int main() {
	L0 = 8.0 [1]; // length of computational box, [1] is the length dimension
	origin (-0.5, -L0/2.0); // allow flow to develop behind the cylinder
	N = 512; // base grid resolution
	mu = mu_v;
	
	// GUI sliders
	display_control (Re, 10, 1000);
  	display_control (maxLevel, 6, 12);
  
  	run(); 
}
