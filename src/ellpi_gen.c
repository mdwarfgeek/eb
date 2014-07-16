#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "machconst.h"

/* Complete elliptic integral of the third kind PI(p, y) using the
   method of Bulirsch (1965), algorithm 4 with the modifications
   from paper II.  In terms of the standard Legendre variables k
   and n, y = 1-k^2 and n = p-1 (for - sign of the characteristic;
   the reverse Abramowitz & Stegun sign convention of + would
   need n = 1-p).  The routine takes sp = sqrt(p) to reduce
   calculations for the caller. 

   Parameter ranges:

   p is assumed to be positive.  The value can be as large as
   sqrt(FLT_MAX) or sqrt(DBL_MAX) for the single and double
   precision versions, respectively.  See machconst.h for macros
   with precomputed values of these quantities.

   y is in (0,1].  Arguments as small as machine epsilon can be
   tolerated without underflow occurring. */

DATATYPE FUNC (DATATYPE sp, DATATYPE y) {
  DATATYPE c, d, e, f, g, kc, mu;
#ifndef NDEBUG
  int iter;
#endif

  /* Convert argument */
  kc = TSQRT(y);

  c = 1.0;
  d = 1.0 / sp;

  e = kc;
  mu = 1.0;

  /* Main iteration loop.  Iteration limit with warning is implemented
     in debug mode to avoid hangs. */
#ifdef NDEBUG
  for(;;) {
#else
  for(iter = 0; iter < 1000; iter++) {
#endif
    /* Update sqrt(p), c, d */
    f   = c;
    c  += d / sp;
    g   = e / sp;
    d   = 2*(f*g + d);
    sp += g;

    /* Update mu, kc */
    g   = mu;
    mu += kc;
    
    if(TABS(g-kc) < TPREC*g)
      break;

    kc = 2*TSQRT(e);
    e = kc * mu;
  }

#ifndef NDEBUG
#define STR(f) #f
  if(iter >= 1000)
    fprintf(stderr,
            STR(FUNC) " failed to converge\n");
#endif

  /* Expression from paper II */
  return(M_PI_2 * (c*mu + d) / (mu*(mu+sp)));
}
