#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Do not try to compile this file directly */

/* ellppi: compute Carlson's R_J(0, y, 1, p)
   
   These specific values are all that is needed for the limb darkened
   occulations, and the calculation has been optimized for this case. */

DATATYPE FUNC (DATATYPE p, DATATYPE y) {
  DATATYPE powq, sumc;
  DATATYPE x, z, a, b;
  DATATYPE srx, sry, srz, srp, sra, srb;
  DATATYPE lc, mc, rmc, sc, poly;
  DATATYPE lj, mj, rmj, tol;
  DATATYPE sx, sy, sz, sp;
  DATATYPE tp, ts;
  DATATYPE ea, eb, ec, e2, e3;
#ifndef NDEBUG
  int iter;
#endif

  powq = 1.0;
  sumc = 0.0;

  /* First iteration, simplifies */
  sry = TSQRT(y);
  srp = TSQRT(p);

  lj = sry;
  sra = p*(sry+1);
  p += lj;
  srb = p*srp;
  p *= 0.25;

  x = 0.25*lj;
  y = 0.25*(y+lj);
  z = 0.25*(1+lj);

  /* Main iteration loop.  Iteration limit with warning is implemented
     in debug mode to avoid hangs. */
#ifdef NDEBUG
  for(;;) {
#else
  for(iter = 0; iter < 1000; iter++) {
#endif
    /* Compute R_C(a, b), first iteration simplifies */
    a = sra*sra;
    b = srb*srb;

    lc = 2*sra*srb + b;

    a = 0.25*(a+lc);
    b = 0.25*(b+lc);

    for(;;) {
      /* Compute mu and test for end of loop.  mu > 0 is guaranteed. */
      mc = (a+b+b)/3;
      sc = b-mc;

      if(TABS(sc) < CTOL*mc)
        break;

      /* Next iteration.  These square roots must be done separately
         rather than doing sqrt(a*b) to prevent overflow. */
      sra = TSQRT(a);
      srb = TSQRT(b);

      lc = 2*sra*srb + b;

      a = 0.25*(a+lc);
      b = 0.25*(b+lc);
    }

    rmc = 1.0/mc;
    sc *= rmc;

    /* Sum polynomial for R_C */
    poly = ((((9.0/22.0)*sc + (3.0/8.0))*sc + (1.0/7.0))*sc + (3.0/10.0))*sc*sc;

    /* Accumulate sum of R_Cs */
    sumc += powq * (poly+1) * TSQRT(rmc);
    powq *= 0.25;

    /* Compute mu and epsilon for R_J, and test for end of loop.
       Because 0 < y < 1 and x_0 = 0, z_0 = 1, we can guarantee
       x < y < z.  This means the test on y is not needed.  We
       also know mu > 0, so do not need to take its absolute
       value. */
    mj = 0.2*(x+y+z+p+p);

    tol = mj*JTOL;

    sx = x-mj;
    sy = y-mj;
    sz = z-mj;
    sp = p-mj;

    if(TABS(sx) < tol &&
       TABS(sy) < tol &&
       TABS(sz) < tol &&
       TABS(sp) < tol)
      break;

    /* Next iteration */
    srx = TSQRT(x);
    sry = TSQRT(y);
    srz = TSQRT(z);
    srp = TSQRT(p);

    tp = srx*sry;
    ts = srx+sry;

    lj = tp + ts*srz;
    sra = p*(ts+srz) + tp*srz;
    p += lj;
    srb = srp * p;
    p *= 0.25;

    x = 0.25*(x+lj);
    y = 0.25*(y+lj);
    z = 0.25*(z+lj);
  }

#ifndef NDEBUG
#define STR(f) #f
  if(iter >= 1000)
    fprintf(stderr,
            STR(FUNC) " failed to converge (x,y,z,p) = (%g,%g,%g,%g)\n",
            x, y, z, p);
#endif

  rmj = 1.0/mj;

  sx *= rmj;
  sy *= rmj;
  sz *= rmj;
  sp *= rmj;

  /* Polynomial summation, based on the method Carlson used in the
     SLATEC implementation, file drj.f. */
  tp = sy*sz;
  ts = sx*(sy+sz);

  ea = tp+ts;
  eb = sx*tp;
  ec = sp*sp;
  
  e2 = ea-3*ec;
  e3 = sp*(ea-ec);
  e3 = eb+e3+e3;

  poly = ((9.0/88.0)*e2 - (9.0/52.0)*e3 - (3.0/14.0)) * e2 + 
         (((3.0/26.0)*sp - (3.0/11.0))*sp + (1.0/6.0)) * eb + 
         ((-(3.0/22.0)*sp + (1.0/3.0))*ea - (1.0/3.0)*ec) * sp + 1.0;

  return(3*sumc + powq * poly * rmj*TSQRT(rmj));
}
