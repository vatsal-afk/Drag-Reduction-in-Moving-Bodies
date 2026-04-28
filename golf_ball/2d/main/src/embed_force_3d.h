/**
 * embed_force_3d.h — 3D-compatible embed_force replacement
 * =========================================================
 * Basilisk's built-in embed_force() (in embed.h) only supports 2D
 * (embed_interpolate asserts dimension==2, and the viscous stress
 * formula is 2D-only).
 *
 * This header provides embed_force_3d() that works in both 2D and 3D.
 *
 * Include AFTER embed.h:
 *   #include "embed.h"
 *   #include "embed_force_3d.h"
 *
 * Then call embed_force_3d(p, u, mu, &Fp, &Fmu) instead of embed_force().
 */

#ifndef EMBED_FORCE_3D_H
#define EMBED_FORCE_3D_H

/**
 * Interpolation at embedded boundary barycentre, extended to 3D.
 * 2D: bilinear.  3D: trilinear.  Fallback: gradient-based.
 */
static inline
double embed_interpolate_3d (Point point, scalar s, coord p)
{
#if dimension == 2
  int i = sign(p.x), j = sign(p.y);
  if (cs[i] && cs[0,j] && cs[i,j])
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) +
            (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
#else /* dimension == 3 */
  int i = sign(p.x), j = sign(p.y), k = sign(p.z);
  if (cs[i]     && cs[0,j]   && cs[i,j]   &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k])
    return
      (((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) +
        (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y))
       * (1. - fabs(p.z)) +
       ((s[0,0,k]*(1. - fabs(p.x)) + s[i,0,k]*fabs(p.x))*(1. - fabs(p.y)) +
        (s[0,j,k]*(1. - fabs(p.x)) + s[i,j,k]*fabs(p.x))*fabs(p.y))
       * fabs(p.z));
#endif
  /* Fallback: gradient-based interpolation for degenerate cells */
  double val = s[];
  foreach_dimension() {
    int ii = sign(p.x);
    if (cs[ii])
      val += fabs(p.x)*(s[ii] - s[]);
    else if (cs[-ii])
      val += fabs(p.x)*(s[] - s[-ii]);
  }
  return val;
}

/**
 * 3D-compatible embedded boundary force computation.
 * Pressure drag Fp and viscous drag Fmu over the embedded surface.
 */
void embed_force_3d (scalar p, vector u, face vector mu,
                     coord * Fp, coord * Fmu)
{
  coord Fps = {0}, Fmus = {0};
  foreach (reduction(+:Fps) reduction(+:Fmus))
    if (cs[] > 0. && cs[] < 1.) {

      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      /* Pressure force */
      double Fn = area * embed_interpolate_3d (point, p, b);
      foreach_dimension()
        Fps.x += Fn*n.x;

      /* Viscous force */
      if (constant(mu.x) != 0.) {
        double mua = 0., fa = 0.;
        foreach_dimension() {
          mua += mu.x[] + mu.x[1];
          fa  += fm.x[] + fm.x[1];
        }
        mua /= fa;

        coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
        foreach_dimension()
          Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) +
                              dudn.y*n.x*n.y);
#else /* dimension == 3 */
        /**
         * 3D viscous stress (symmetric deformation tensor D dotted with n):
         *   F_mu,x = -mu*area*[dudn.x*(nx^2+1) + dudn.y*nx*ny + dudn.z*nx*nz]
         * with foreach_dimension() rotating x->y->z automatically.
         */
        foreach_dimension()
          Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) +
                              dudn.y*n.x*n.y +
                              dudn.z*n.x*n.z);
#endif
      }
    }

  *Fp = Fps;
  *Fmu = Fmus;
}

#endif /* EMBED_FORCE_3D_H */
