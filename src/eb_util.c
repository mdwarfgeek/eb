#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "eb.h"
#include "machconst.h"

double eb_phiperi (double esinw, double ecosw) {
  double esinwsq, ecoswsq, esq, num, phi;

  esinwsq = esinw*esinw;
  ecoswsq = ecosw*ecosw;
  esq = esinwsq + ecoswsq;

  /* Mean anomaly offset from inferior conjunction to periastron. */
  num = ecosw * sqrt(1.0 - esq);
  phi = -atan2(num, esq + esinw) + num / (1.0 + esinw);
  phi = remainder(phi, TWOPI);

  return(phi / TWOPI);
}

double eb_phisec (double esinw, double ecosw) {
  double esinwsq, ecoswsq, esq, num, phi;

  esinwsq = esinw*esinw;
  ecoswsq = ecosw*ecosw;
  esq = esinwsq + ecoswsq;

  /* Mean anomaly offset from inferior conjunction to superior
     conjunction (= approximate phase of secondary eclipse). */
  num = 2 * ecosw * sqrt(1.0 - esq);
  phi = atan2(-num, esq + ecoswsq - 1.0) + num / (1.0 - esinwsq);
  phi = fmod(phi, TWOPI);
  if(phi < 0)
    phi += TWOPI;

  return(phi / TWOPI);
}

static double eb_cont_f (double cvw,
                         double esinw, double ecosw,
                         double cosi, double targ,
                         double sgn) {
  double cvwsq, svwsq, denom;

  cvwsq = cvw*cvw;
  svwsq = 1.0 - cvwsq;

  denom = 1.0 + ecosw*cvw + sgn * sqrt(svwsq) * esinw;

  return((cvwsq + svwsq * cosi*cosi) / (denom*denom) - targ);
}

#define BISECT_MAXITER 100

static double eb_cont_bisect (double a, double b,
                              double fa, double fb,
                              double esinw, double ecosw,
                              double cosi, double targ,
                              double sgn) {
  int iter;
  double c, fc;

  for(iter = 0; iter < BISECT_MAXITER; iter++) {
    c = 0.5*(a+b);
    fc = eb_cont_f(c, esinw, ecosw, cosi, targ, sgn);
    if(fc == 0 || 0.5*(b-a) < DBL_EPSILON)
      break;

    if((fc < 0 && fa < 0) ||
       (fc >= 0 && fa >= 0)) {
      a = c;
      fa = fc;
    }
    else {
      b = c;
      fb = fc;
    }
  }

  if(iter > BISECT_MAXITER)
    fprintf(stderr, "eb_cont_bisect: iteration limit reached\n");

  return c;
}

static void eb_cont_find (double esinw, double ecosw,
                          double cosi, double targ,
                          double sgn,
                          double *cvwout) {
  double a, b, c, d, h, w;
  double fa, fb, fc, fd;
  int iter, niter;

  /* Interval bracketing minimum */
  a = -1;
  b =  1;

  fa = eb_cont_f(a, esinw, ecosw, cosi, targ, sgn);
  fb = eb_cont_f(b, esinw, ecosw, cosi, targ, sgn);

  /* Abort if endpoints are in eclipse.  This happens when the bodies
     are touching or inside each other.  The eclipse solutions are
     invalid here, but it's possible somebody might try to calculate
     the contact points without checking. */
  if(fa <= 0 || fb <= 0) {
    /* Return conjunction to signal error */
    cvwout[0] = 0;
    cvwout[1] = 0;
    return;
  }

  /* Calculate value at conjunction */
  c = 0;
  fc = eb_cont_f(c, esinw, ecosw, cosi, targ, sgn);

  /* Check if conjunction is in eclipse, this is the most common case. */
  if(fc > 0) {
    /* It is not.  Golden section search until we bracket a root. */
    h = b-a;

    /* 1/phi */
    w = (sqrt(5) - 1) / 2;

    /* Number of iterations to bracket down to sqrt(epsilon) */
    niter = ceil(log(SQRT_DBL_EPSILON/h)/log(w));

    /* c,d are probe points partitioning interval [a,b] in ratio w */
    c = a + (1-w)*h;
    d = a + w*h;
    
    fc = eb_cont_f(c, esinw, ecosw, cosi, targ, sgn);
    fd = eb_cont_f(d, esinw, ecosw, cosi, targ, sgn);

    /* Search, terminating when there's a root in [a,c] or [d,b] */
    for(iter = 0; fc > 0 && fd > 0 && iter < niter; iter++) {
      h *= w;
      
      if(fc < fd) {
        /* Minimum in [a,d], set new endpoint */
        b = d;
        fb = fd;

        /* We've already calculated the upper probe point */
        d = c;
        fd = fc;

        /* Calculate the new lower probe point */
        c = a + (1-w)*h;
        fc = eb_cont_f(c, esinw, ecosw, cosi, targ, sgn);
      }
      else {
        /* Minimum in [c,b], set new start point */
        a = c;
        fa = fc;

        /* We've already calculated the lower probe point */
        c = d;
        fc = fd;

        /* Calculate the new upper probe point */
        d = a + w*h;
        fd = eb_cont_f(d, esinw, ecosw, cosi, targ, sgn);
      }
    }

    /* Use whichever of c, d is the smaller, for simplicity */
    if(fd < fc) {
      c = d;
      fc = fd;
    }

    /* There's now a root in [a,c] or [c,b] or no root at all */
  }

  /* Is there a root? */
  if(fc >= 0) {
    /* Nope, means there's no eclipse so return conjunction */
    cvwout[0] = 0;
    cvwout[1] = 0;
    return;
  }

  /* Bisect roots */
  if(sgn > 0) {
    cvwout[1] = eb_cont_bisect(a, c, fa, fc, esinw, ecosw, cosi, targ, sgn);
    cvwout[0] = eb_cont_bisect(c, b, fc, fb, esinw, ecosw, cosi, targ, sgn);
  }
  else {
    cvwout[0] = eb_cont_bisect(a, c, fa, fc, esinw, ecosw, cosi, targ, sgn);
    cvwout[1] = eb_cont_bisect(c, b, fc, fb, esinw, ecosw, cosi, targ, sgn);
  }
}

void eb_phicont (double esinw, double ecosw, double cosi,
                 double d, double *phi) {
  double esq, ecc, omesq, roe, sinw, cosw;
  double dsec, dcec, dc;
  double sv, cv, dse, dce, dd;
  double dphi;
  double targ, cvw[4], svw;
  int i;

  double sgn[4] = { 1.0, 1.0, -1.0, -1.0 };

  esq = esinw*esinw + ecosw*ecosw;
  ecc = sqrt(esq);
  omesq = 1.0 - esq;
  roe = sqrt(omesq);

  targ = d / omesq;
  targ *= targ;

  if(ecc > 0) {
    cosw = ecosw / ecc;
    sinw = esinw / ecc;
  }
  else {
    cosw = 0.0;
    sinw = 1.0;
  }

  /* d* sin, cos eccentric anomaly at inferior conjunction */
  dsec = roe * cosw;
  dcec = ecc + sinw;
  dc = 1.0 + esinw;

  /* Primary */
  eb_cont_find(esinw, ecosw, cosi, targ, sgn[0], &(cvw[0]));

  /* Secondary */
  eb_cont_find(esinw, ecosw, cosi, targ, sgn[2], &(cvw[2]));

  /* Convert to phase */
  for(i = 0; i < 4; i++) {
    svw = sgn[i] * sqrt(1.0 - cvw[i]*cvw[i]);

    /* sin(v) and cos(v) */
    sv = svw * cosw - cvw[i] * sinw;
    cv = cvw[i] * cosw + svw * sinw;

    /* d*sin(E) and d*cos(E) */
    dse = sv * roe;
    dce = ecc + cv;
    dd = 1.0 + ecc*cv;

    /* M-MC = E-EC - e (sin E - sin EC)
       
       E-EC = atan2(dc * sin (E-EC),
                    dc * cos (E-EC))
       
       and use sin E and dc*sin(EC) from above for rest. */
    
    dphi = atan2(dse * dcec - dce * dsec,
                 dce * dcec + dse * dsec)
         - ecc * (dse*dc - dsec*dd)/(dd*dc);
    dphi = fmod(dphi, TWOPI);
    if(dphi < 0)
      dphi += TWOPI;

    phi[i] = dphi / TWOPI;
  }
}

void eb_getvder (double *v, double gamma, double ktot, double *vder) {
  double esq, roe, sini, qpo, tmp, omega;
  double phi[4], dpp, dps;

  vder[EB_PAR_I] = acos(v[EB_PAR_COSI]) * 180.0/M_PI;
  vder[EB_PAR_R1A] = v[EB_PAR_RASUM] / (1.0 + v[EB_PAR_RR]);
  vder[EB_PAR_R2A] = v[EB_PAR_RR] * vder[EB_PAR_R1A];

  esq = v[EB_PAR_ECOSW]*v[EB_PAR_ECOSW] + v[EB_PAR_ESINW]*v[EB_PAR_ESINW];

  vder[EB_PAR_E] = sqrt(esq);
  vder[EB_PAR_OMEGA] = atan2(v[EB_PAR_ESINW], v[EB_PAR_ECOSW]) * 180.0/M_PI;
  if(vder[EB_PAR_OMEGA] < 0.0)  /* wrap to conventional range [0,360) */
    vder[EB_PAR_OMEGA] += 360.0;
  
  roe = sqrt(1.0-esq);
  sini = sqrt(1.0-v[EB_PAR_COSI]*v[EB_PAR_COSI]);
  qpo = 1.0+v[EB_PAR_Q];

  /* Orbital angular frequency. */ 
  omega = 2.0*M_PI / (v[EB_PAR_P]*EB_DAY);

  /* Average apsidal motion is included in the period, so take that off. */
  omega -= v[EB_PAR_DWDT] / EB_DAY;

  /* Correct to system Barycenter.  The factor (1.0+gamma/c) accounts
   * for the Doppler shift due to the systemic motion relative to the
   * solar system barycenter, which is the rest frame the period was
   * calculated in.
   */
  omega += omega * gamma*1000/EB_LIGHT;

  tmp = ktot*1000 * roe / sini;

  vder[EB_PAR_A] = tmp / (EB_RSUN*omega);
  vder[EB_PAR_MTOT] = tmp*tmp*tmp / (EB_GMSUN*omega);
  vder[EB_PAR_M1] = vder[EB_PAR_MTOT] / qpo;
  vder[EB_PAR_M2] = v[EB_PAR_Q] * vder[EB_PAR_M1];
  vder[EB_PAR_RTOT] = v[EB_PAR_RASUM] * vder[EB_PAR_A];
  vder[EB_PAR_R1] = vder[EB_PAR_R1A] * vder[EB_PAR_A];
  vder[EB_PAR_R2] = vder[EB_PAR_R2A] * vder[EB_PAR_A];

  /* factor of 100 converts to cgs */
  vder[EB_PAR_LOGG1] = log10(100 * EB_GMSUN * vder[EB_PAR_M1] /
                        (EB_RSUN*EB_RSUN*vder[EB_PAR_R1]*vder[EB_PAR_R1]));
  vder[EB_PAR_LOGG2] = log10(100 * EB_GMSUN * vder[EB_PAR_M2] /
                        (EB_RSUN*EB_RSUN*vder[EB_PAR_R2]*vder[EB_PAR_R2]));

  /* I think this is only for circular orbits */
  vder[EB_PAR_VSYNC1] = 1.0e-3 * omega * vder[EB_PAR_R1]*EB_RSUN;
  vder[EB_PAR_VSYNC2] = 1.0e-3 * omega * vder[EB_PAR_R2]*EB_RSUN;

  /* Eq. 6.1 of Zahn 1977 for stars with convective envelopes except in Gyr */
  vder[EB_PAR_TSYNC] = 1.0e-5 * v[EB_PAR_P]*v[EB_PAR_P]*v[EB_PAR_P]*v[EB_PAR_P]
                     * qpo*qpo / (4*v[EB_PAR_Q]*v[EB_PAR_Q]);

  /* Eq. 6.2 of Zahn 1977 for stars with convective envelopes except in Gyr */
  vder[EB_PAR_TCIRC] = 1.0e-3 * pow((1.0+v[EB_PAR_Q]) / 2, 5.0/3.0)
                     * pow(v[EB_PAR_P], 16.0/3.0) / v[EB_PAR_Q];

  vder[EB_PAR_TSEC] = v[EB_PAR_T0]
                    + v[EB_PAR_P]*eb_phisec(v[EB_PAR_ESINW], v[EB_PAR_ECOSW]);
  
  eb_phicont(v[EB_PAR_ESINW], v[EB_PAR_ECOSW],
             v[EB_PAR_COSI], v[EB_PAR_RASUM], phi);

  dpp = phi[1] - phi[0];
  if(dpp < 0)
    dpp += 1.0;

  dps = phi[3] - phi[2];
  if(dps < 0)
    dps += 1.0;

  vder[EB_PAR_DURPRI] = dpp * v[EB_PAR_P] * 24;
  vder[EB_PAR_DURSEC] = dps * v[EB_PAR_P] * 24;
}
  
/* Fortran wrappers.  For simplicity we define both the single and
   double underscore variants, so we don't need to know how the
   Fortran compiler is set up. */

void eb_phiperi_ (double *esinw, double *ecosw, double *phiperi) {
  *phiperi = eb_phiperi(*esinw, *ecosw);
}

void eb_phiperi__ (double *esinw, double *ecosw, double *phiperi) {
  *phiperi = eb_phiperi(*esinw, *ecosw);
}

void eb_phisec_ (double *esinw, double *ecosw, double *phisec) {
  *phisec = eb_phisec(*esinw, *ecosw);
}

void eb_phisec__ (double *esinw, double *ecosw, double *phisec) {
  *phisec = eb_phisec(*esinw, *ecosw);
}

void eb_phicont_ (double *esinw, double *ecosw, double *cosi,
                  double *d, double *phi) {
  eb_phicont(*esinw, *ecosw, *cosi, *d, phi);
}

void eb_phicont__ (double *esinw, double *ecosw, double *cosi,
                  double *d, double *phi) {
  eb_phicont(*esinw, *ecosw, *cosi, *d, phi);
}

void eb_getvder_ (double *parm, double *gamma, double *ktot, double *vder) {
  eb_getvder(parm, *gamma, *ktot, vder);
}

void eb_getvder__ (double *parm, double *gamma, double *ktot, double *vder) {
  eb_getvder(parm, *gamma, *ktot, vder);
}

