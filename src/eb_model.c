#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "eb.h"
#include "machconst.h"
#include "sincos.h"

/* Temporary */
extern void ELLKE (DATATYPE y, DATATYPE *ke);
extern DATATYPE ELLPI (DATATYPE p, DATATYPE y);

/* These define how far the Newton-Raphson iterations to solve Kepler's
   equation go.  The precision here is probably overkill for most
   applications, so tune according to taste if it's too slow. */
#define KEPLER_PREC    1.0e-13
#define KEPLER_MAXITER 100

/* This structure contains the per-star parameters.  Having them all
   in once place makes it easier to extract the ones for the occulter
   when computing the eclipse. */
struct star {
  DATATYPE rr;      /* radius ratio ("that"/"this") */
  DATATYPE rrsq;

  DATATYPE q;       /* mass ratio ("that"/"this") */
  DATATYPE u1;      /* limb darkening */
  DATATYPE u2;
  DATATYPE gd;
  DATATYPE refl;

  DATATYPE uss;     /* u1 + 2*u2 */
  DATATYPE ldint;   /* LD integral */
  DATATYPE ldnorm;  /* 1/(LD integral) */

  DATATYPE rsq;     /* (r/a)**2 */
  DATATYPE rcb;     /* (r/a)**3 */
  DATATYPE aorsq;   /* (a/r)**2 */

  DATATYPE o;
  DATATYPE delt;
  DATATYPE l;
  DATATYPE rl;
  DATATYPE ol;
  DATATYPE ltot;

  double rot;       /* rotation parameter */
  double qfltt;     /* light travel time mass ratio factor */

  DATATYPE fecs;    /* fraction of spots eclipsed */
  DATATYPE ob;      /* base spottedness out of eclipse */
  DATATYPE oa1;     /* sin(...) coeff */
  DATATYPE ob1;
  DATATYPE oa2;     /* sin(2*...) coeff */
  DATATYPE ob2;
  DATATYPE oav;     /* average out of eclipse level */
};

static inline void bright (struct star *s) {
  DATATYPE eps, rb, fmax, fmin;

  /* Biaxial ellipsoid dimensions and oblateness.
     See Chandrasekhar (1933) and Binnendijk (1974, VA, 16, 61). */
  eps = 1.5 * s->q * s->rcb / (1.0 + s->rcb * (1.0 + 7*s->q) / 6);
  rb = TCBRT(1.0-eps);

  /* Min and max brightness, from Gimenez. */
  fmax = 1.0 - s->u1*(1.0 - 0.4*eps)/3.0 - s->u2*(1.0 - 0.6*eps)/6.0
    + 2.0*s->gd*eps*(1.5 - 13.0*s->u1/30.0 - 0.2*s->u2);
  fmax /= (1.0 - eps);
  fmin = 1.0 - s->u1*(1.0 + 0.8*eps)/3.0 - s->u2*(1.0 + 1.2*eps)/6.0
    + 2.0*s->gd*eps*(1.0 - 7.0*s->u1/15.0 - 4.0*s->u2/15.0);

  s->delt = (fmax-fmin)/fmax;

  /* Unnormalized brightness. */
  s->o = s->rsq * rb*rb * fmax;
}

static inline void ltt (double cltt,
                        double ecc, double roe, double sinw, double cosw,
                        double ma, double ea,
                        double se, double ce, double f,
                        double *rv, double *svw, double *cvw) {
  double dsvw, df, delta;
  int i;

  /* Solve modified Kepler's equation with light travel time.
     We need to add a correction to the mean anomaly of
     dma = rv*sin(v+w) * a omega sin i / c
     cltt = a omega sin i / c = ktot sqrt(1-e^2) / c
     and is precomputed.
     our "svw" variable is actually rv*sin(v+w).  */
  
  /* First iteration using previous results */
  dsvw = ce*cosw*roe - se*sinw;
  
  f -= (*svw)*cltt;
  df = (*rv) - dsvw*cltt;
  
  ea -= f / df;
  
  /* Loop */
  for(i = 0; i < KEPLER_MAXITER; i++) {
    inline_bare_sincos(ea, se, ce);
    
    *rv  = 1.0 - ecc * ce;
    *svw = se*cosw*roe + (ce-ecc)*sinw;
    dsvw = ce*cosw*roe - se*sinw;
    
    f = ea - ecc * se - ma - (*svw)*cltt;
    df = (*rv) - dsvw*cltt;
    
    delta = f / df;
    
    if(fabs(delta) < KEPLER_PREC)
      break;
    
    ea -= delta;
  }

  if(i >= KEPLER_MAXITER)
    fprintf(stderr, "ltt: iteration limit reached\n");

  *cvw = (ce-ecc)*cosw - se*sinw*roe;
}

void FUNC (double *parm, double *t, DATATYPE *ol1, DATATYPE *ol2,
           unsigned char *typ, DATATYPE *out, unsigned char *ieclout,
           int flags, int npt) {
  struct star sp, ss, *s;

  double ecosw0, esinw0, tconj, period, dphi;
  double omega, esq, ecc, roe;
  double maconj;
  double dt, ma, ea, alpha, beta, z, ste;
  double se, ce, f, df, delta, vnorm, sv, cv;
  double dwdt, dw, sdw, cdw, ecosw, esinw, cosw, sinw;
  double svw, cvw, phio, stid, ctid, so, co;
  double cltt, rv1, rv2, svw1, cvw1, svw2, cvw2, dys, dzs;

  DATATYPE jsb, rsum, cosi;
  DATATYPE atid, lth, zp;

  DATATYPE csqi, ssqi, sini;
  DATATYPE rp, norm;

  int p, i;

  DATATYPE svwt, cpsi, csqpsi, tmp;
  DATATYPE d, dsq;
  DATATYPE rr, rrsq, a, b, asq, bsq;
  DATATYPE eta, lam, kap0, kap1, ksrt;
  DATATYPE y, q, ke[2], tpnk, tt, ab, ts, drr;
  DATATYPE area, fecl;
  DATATYPE ltot;

  unsigned char iecl;

  /* Unpack parameter vector into a more convenient form. */
  jsb     = parm[EB_PAR_J];
  rsum    = parm[EB_PAR_RASUM];
  sp.rr   = parm[EB_PAR_RR];
  cosi    = parm[EB_PAR_COSI];
  ecosw0  = parm[EB_PAR_ECOSW];
  esinw0  = parm[EB_PAR_ESINW];
  sp.u1   = parm[EB_PAR_LDLIN1];
  sp.u2   = parm[EB_PAR_LDNON1];
  ss.u1   = parm[EB_PAR_LDLIN2];
  ss.u2   = parm[EB_PAR_LDNON2];
  sp.gd   = parm[EB_PAR_GD1];
  ss.gd   = parm[EB_PAR_GD2];
  sp.refl = parm[EB_PAR_REFL1];
  ss.refl = parm[EB_PAR_REFL2];
  sp.q    = parm[EB_PAR_Q];
  atid    = parm[EB_PAR_TIDANG];
  lth     = parm[EB_PAR_L3];
  dphi    = parm[EB_PAR_PHI0];
  zp      = parm[EB_PAR_M0];
  tconj   = parm[EB_PAR_T0];
  period  = parm[EB_PAR_P];
  cltt    = parm[EB_PAR_CLTT];

  sp.rot  = parm[EB_PAR_ROT1];
  ss.rot  = parm[EB_PAR_ROT2];
  sp.fecs = parm[EB_PAR_FSPOT1];
  ss.fecs = parm[EB_PAR_FSPOT2];

  sp.ob   = parm[EB_PAR_OOE1O];
  sp.oa1  = parm[EB_PAR_OOE11A];
  sp.ob1  = parm[EB_PAR_OOE11B];
  sp.oa2  = parm[EB_PAR_OOE12A];
  sp.ob2  = parm[EB_PAR_OOE12B];

  ss.ob   = parm[EB_PAR_OOE2O];
  ss.oa1  = parm[EB_PAR_OOE21A];
  ss.ob1  = parm[EB_PAR_OOE21B];
  ss.oa2  = parm[EB_PAR_OOE22A];
  ss.ob2  = parm[EB_PAR_OOE22B];

  dwdt    = parm[EB_PAR_DWDT];

  if(flags & EB_FLAG_PHI) {
    /* [0,1] phase, not time */
    tconj = 0.0;
    omega = TWOPI;
  }
  else {
    /* Angular frequency of orbit */
    omega = TWOPI / period;
  }

  /* Convert phase offset to radians */
  dphi *= TWOPI;

  /* Eccentricity */
  esq = ecosw0*ecosw0 + esinw0*esinw0;
  ecc = sqrt(esq);
  roe = sqrt(1.0 - esq);

  /* Convert cltt from ktot/c to a omega sin i / c */
  cltt *= roe;

  /* Mean anomaly at inferior conjunction */
  maconj = atan2(roe * ecosw0, esq + esinw0) - roe * ecosw0 / (1 + esinw0);

  /* Precompute these */
  csqi = cosi*cosi;
  ssqi = 1.0 - csqi;
  sini = TSQRT(ssqi);

  if(atid != 0.0) {
    inline_sincos(atid*M_PI/180.0, stid, ctid);
  }
  else {
    stid = 0;
    ctid = 1;
  }

  /* Radius ratio, both ways up */
  ss.rr = 1.0 / sp.rr;

  /* Square of radius ratio */
  sp.rrsq = sp.rr*sp.rr;
  ss.rrsq = ss.rr*ss.rr;

  /* Limb darkening integral */
  sp.uss = sp.u1 + 2*sp.u2;
  sp.ldint = 1.0 - sp.u1 / 3.0 - sp.u2 / 6.0;
  sp.ldnorm = 1.0 / sp.ldint;

  ss.uss = ss.u1 + 2*ss.u2;
  ss.ldint = 1.0 - ss.u1 / 3.0 - ss.u2 / 6.0;
  ss.ldnorm = 1.0 / ss.ldint;

  /* Radii: this is a bit subtle.  The eclipse calculation needs the
     values to be consistent in the sense of having the correct ratios,
     so the calculation is arranged to compute the quantity for the
     primary star and then multiply by the appropriate power of the
     radius ratio. */

  /* Radius of primary */
  rp = rsum / (1.0+sp.rr);

  /* Compute (r/a)**2 */
  sp.rsq = rp*rp;
  ss.rsq = sp.rsq * sp.rrsq;

  /* Compute (a/r)**2 */
  sp.aorsq = 1.0 / sp.rsq;
  ss.aorsq = sp.aorsq * ss.rrsq;

  /* 1 / (1+q) and -q / (1+q) scaling factors
     used in light travel correction */
  ss.qfltt = 1.0 / (1.0 + sp.q);
  sp.qfltt = -sp.q * ss.qfltt;

  /* Normalized brightness of stars at quadrature */
  if(sp.q > 0) {
    ss.q = 1.0 / sp.q;

    /* (r/a)**3 */
    sp.rcb = sp.rsq * rp;
    ss.rcb = sp.rcb * sp.rrsq*sp.rr;

    bright(&sp);
    bright(&ss);
  }
  else {
    ss.q = 0;  /* plugs up compiler warning */

    sp.o = sp.rsq * sp.ldint;
    ss.o = ss.rsq * ss.ldint;

    sp.delt = 0;
    ss.delt = 0;
  }

  ss.o *= jsb;  /* apply surface brightness ratio */

  norm = 1.0 / (sp.o + ss.o);

  sp.o *= norm;
  ss.o *= norm;

  if(!(flags & EB_FLAG_REFL)) {
    /* "refl" is the albedo, so multiply by l(other) (r/a)^2 */
    csqpsi = ctid*ctid;

    sp.refl *= ss.o * (1.0 - ss.delt * csqpsi) * sp.rsq;
    ss.refl *= sp.o * (1.0 - sp.delt * csqpsi) * ss.rsq;
  }
  /* otherwise, the original EBOP-alike model is used, this is
     usually better when fitting for reflection. */

  /* Average OOE light */
  sp.oav = sp.ob;
  ss.oav = ss.ob;

  if(sp.rot)
    sp.oav += TSQRT(sp.oa1*sp.oa1 + sp.ob1*sp.ob1 +
                    sp.oa2*sp.oa2 + sp.ob2*sp.ob2);
  if(ss.rot)
    ss.oav += TSQRT(ss.oa1*ss.oa1 + ss.ob1*ss.ob1 +
                    ss.oa2*ss.oa2 + ss.ob2*ss.ob2);

  /* Loop over light curve points */
  for(p = 0; p < npt; p++) {
    if(typ[p] == EB_OBS_AVLR) {
      /* Special call to return avg. light ratio (CHECKME!) */
      tmp = 0.5*(1.0 + 0.5*ssqi);

      *out = (ss.o * (1.0 - 0.5*ss.delt*ssqi) * (1.0-ss.oav) + ss.refl*tmp)
           / (sp.o * (1.0 - 0.5*sp.delt*ssqi) * (1.0-sp.oav) + sp.refl*tmp);

      continue;
    }

    /* Time relative to T0 */
    dt = t[p] - tconj;

    /* Apsidal motion, if needed */
    if(dwdt == 0) {
      dw = 0;
      ecosw = ecosw0;
      esinw = esinw0;
    }
    else {
      if(flags & EB_FLAG_PHI)
        dw = dwdt * dt * period;  /* dt is phase */
      else
        dw = dwdt * dt;  /* dt is time */

      inline_sincos(dw, sdw, cdw);

      ecosw = ecosw0 * cdw - esinw0 * sdw;
      esinw = esinw0 * cdw + ecosw0 * sdw;
    }

    /* Mean anomaly, reduced to [-pi, pi] */
    ma = remainder(omega * dt + maconj - dphi - dw, TWOPI);

    /* Initial eccentric anomaly */
    if(ecc > 0) {
      /* For eccentric orbits, use cubic approximation from Mikkola (1987)
         to provide initial guess of eccentric anomaly. */

      /* Eq. 9a */
      tmp = 1.0 / (4 * ecc + 0.5);

      alpha = (1.0 - ecc) * tmp;
      beta = 0.5 * ma * tmp;
      
      /* Eq. 9b */
      z = cbrt(beta + copysign(sqrt(beta*beta + alpha*alpha*alpha), beta));
      
      /* Eq. 9c: initial value of sin(E/3) */
      ste = z - alpha / z;
      
      /* Eq. 7: 5th order correction term */
      ste -= 0.078 * ste*ste*ste*ste*ste / (1.0 + ecc);
      
      /* Eq. 8: eccentric anomaly */
      ea = ma + ecc * ste * (3.0 - 4.0 * ste*ste);
    }
    else
      ea = ma;

    /* Refine solution of Kepler's equation using Newton's method */
    for(i = 0; i < KEPLER_MAXITER; i++) {
      inline_bare_sincos(ea, se, ce);

      f = ea - ecc * se - ma;
      df = 1.0 - ecc * ce;
      delta = f / df;
      
      if(fabs(delta) < KEPLER_PREC) {
	/* I think that's enough... */
	break;
      }

      ea -= delta;
    }

    if(i >= KEPLER_MAXITER)
      fprintf(stderr, "kepler: iteration limit reached\n");

    /* rv * sin, cos of true anomaly */
    sv = se * roe;
    cv = ce - ecc;

    /* Need to substitute these when e=0. */
    if(ecc > 0) {
      cosw = ecosw / ecc;
      sinw = esinw / ecc;
    }
    else {
      cosw = 0.0;
      sinw = 1.0;
    }

    /* rv * sin(v+w) and rv * cos(v+w) */
    svw = sv*cosw + cv*sinw;
    cvw = cv*cosw - sv*sinw;

    if(typ[p] == EB_OBS_LSS) {
      out[p] = svw * sini;
      continue;
    }

    /* Values for separate stars at time light was emitted */
    rv1 = df;
    rv2 = df;
    svw1 = svw;
    svw2 = svw;
    cvw1 = cvw;
    cvw2 = cvw;

    if(cltt) {
      /* Light travel corrected rv, rv*sin(v+w) and rv*cos(v+w) */
      ltt(cltt * sp.qfltt,
          ecc, roe, sinw, cosw, ma, ea, se, ce, f,
          &rv1, &svw1, &cvw1);
      ltt(cltt * ss.qfltt,
          ecc, roe, sinw, cosw, ma, ea, se, ce, f,
          &rv2, &svw2, &cvw2);

      /* Plane of sky coords star 1 - star 2 in CMS */
      dys = cvw1 * sp.qfltt - cvw2 * ss.qfltt;
      dzs = svw1 * sp.qfltt - svw2 * ss.qfltt;
    }
    else {
      dys = cvw;
      dzs = svw;
    }

    /* Handle radial velocity, if requested */
    if(typ[p] == EB_OBS_VRAD1) {
      /* vrad = gamma + K*out[p]
         out[p] = cos(v+w) + e cos w */
      out[p] = cvw1 / rv1 + ecosw;
      continue;  /* skip rest */
    }
    else if(typ[p] == EB_OBS_VRAD2) {
      out[p] = cvw2 / rv2 + ecosw;
      continue;
    }
    /* else: light needed */

    /* 1/rv without light time */
    vnorm = 1.0/df;

    /* sin(v+w+atid) */
    svwt = (ctid*svw + stid*cvw) * vnorm;

    /* Base light from each component including shape (ellipsoidal).
       Following Binnendijk (1974, VA, 16, 61), Eq. 17, approximate:

       L = 1 - 0.5 * e^2 * cos^2 (psi)
         = 1 - delt * cos^2 (psi)

       where psi is the angle between the line of sight and the
       binary major axis, and delt is 1 - Lmin/Lmax as calculated
       above.

       The angle psi is given by Eq. 12:

       cos psi = sin i sin(v+w+atid)

       where we have included the tidal lead/lag term above. */

    /* Eq. 12 squared */
    csqpsi = ssqi * svwt*svwt;

    /* Eq. 17 */
    sp.l = sp.o * (1.0 - sp.delt * csqpsi);
    ss.l = ss.o * (1.0 - ss.delt * csqpsi);

    /* Out of eclipse variations due to spots */
    phio = omega * dt - dphi;

    if(sp.rot) {
      inline_sincos(phio * sp.rot, so, co);
      sp.ol = (sp.oa1*so + sp.ob1*co +
               2*sp.oa2*so*co + sp.ob2*(co+so)*(co-so) - sp.oav);
    }
    else
      sp.ol = -sp.oav;

    if(ol1)
      sp.ol += ol1[p];

    sp.ol *= sp.l;

    if(ss.rot) {
      inline_sincos(phio * ss.rot, so, co);
      ss.ol = (ss.oa1*so + ss.ob1*co +
               2*ss.oa2*so*co + ss.ob2*(co+so)*(co-so) - ss.oav);
    }
    else
      ss.ol = -ss.oav;

    if(ol2)
      ss.ol += ol2[p];

    ss.ol *= ss.l;

    /* Reflection effect using the simple formula from Milne
       (1926, MNRAS, 87, 43) and Russell (1939, ApJ, 90, 641).
       The reflected light is:

       L(refl)_1 = (1/8) L_2 (R_1/a)^2 (1 + 2 cos psi)^2
       L(refl)_2 = (1/8) L_1 (R_2/a)^2 (1 - 2 cos psi)^2

       where L_1 and L_2 are the luminosities of the two stars.
       The direct light to the observer is L/4, so:

       L_1(refl) = (1/2) L_2(direct) (R_1/a)^2 (1 + 2 cos psi)^2

       and similarly for star 2.

       In the default formulation here, the entries in the parameter
       vector are constants multiplying this expression, i.e. the
       albedo, where an albedo of unity corresponds to the Milne and
       Russell formulation.  In EBOP, and here when EB_FLAG_REFL is
       set, the reflection coefficient also includes the
       L_1(direct) (R_1/a)^2 and L_2(direct) (R_2/a)^2 factors. */

    /* Cosine of phase angle for reflection = cos psi, this time
       without tidal lead/lag. */
    cpsi = sini * svw * vnorm;

    /* (1/2) (1 + cos^2 psi) */
    tmp = 0.5*(1.0 + cpsi*cpsi);

    /* (1/2) (1 + cos^2 psi +/- 2 cos psi) = (1/2) (1 + cos psi)^2 */
    sp.rl = sp.refl * (tmp + cpsi);
    ss.rl = ss.refl * (tmp - cpsi);

    /* Total uneclipsed light on each star */
    sp.ltot = sp.l + sp.ol + sp.rl;
    ss.ltot = ss.l + ss.ol + ss.rl;

    /* Which eclipse? */
    if(svw > 0)  /* primary */
      s = &sp;
    else  /* secondary */
      s = &ss;

    /* Plane of sky distance.  Without light travel, this is
       (d/r)^2 = (d/a)^2 (a/r)^2
       (d/a)^2 = cos^2(v+w) + sin^2(v+w) * cos^2 i
       with light travel, we use dys and dzs computed above. */
    dsq = TABS(s->aorsq * (dys*dys + dzs*dzs*csqi));
    d = TSQRT(dsq);

    if(typ[p] == EB_OBS_PSS) {
      out[p] = d;
      continue;  /* skip rest */
    }

    rrsq = s->rrsq;
    rr = s->rr;

    a = d-rr;
    b = d+rr;

    if(typ[p] == EB_OBS_A) {
      out[p] = a;

      if(ieclout) {
        if(a >= 1 || rr == 0)
          ieclout[p] = 0;
        else
          ieclout[p] = 1;
      }

      continue;  /* skip rest */
    }

    if(a >= 1 || rr == 0) {  /* no overlap, integrals trivially zero */
      area = 0;
      fecl = 0;
      iecl = 0;
    }
    else if(a <= -1) {  /* completely eclipsed */
      area = 1;
      fecl = 1;
      iecl = 1;
    }
    else {  /* gotta do some integrations then */
      /* Method of Mandel & Agol (2002).  The notation has been modified
         to suit the author's preferences.  There are also some errors
         in the equations of Sect. 4 in the paper, see the errata at:
         http://www.astro.washington.edu/users/agol/mandel_agol_errata.pdf

         The area integral is lambda^e in their notation, and is given by
         Eq. 1.  The limb darkened integral comes in two parts in the
         paper, lambda^d and eta^d, which are called simply lam and eta
         here, and is given by Eq. 7.  (4*Omega)^-1 is ldnorm in our
         notation and has already been computed.

         The calculation is split depending on whether the occulter is
         entirely within the disk, which corresponds to the condition
         d < 1-rr, equivalently b < 1. */
      asq = a*a;
      bsq = b*b;

      if(TABS(bsq-1.0) < TEPS) {  /* touches limb (case 4, errata 3, 5) */
        area = rrsq;
        eta = rrsq*(0.5*rrsq + dsq);

        lam = TACOS(1 - 2*rr)*2/(3*M_PI) - 
          (3+2*rr-8*rrsq)*TSQRT(rr*(1-rr))*4/(9*M_PI);
      }
      else if(b < 1) {  /* occulter inside disk */
        /* Cases 3, 5, 6, 9, 10 from the paper */
        area = rrsq;
        eta = rrsq*(0.5*rrsq + dsq);

        if(asq <= TSQRTMIN*bsq) {  /* cases 5 and 6: touches center */
          if(d == 0.5) {  /* case 6 */
            lam = 1.0/3 - 4.0/(9*M_PI);
          }
          else {  /* case 5 and erratum 1 */
            ELLKE(1 - 4*rrsq, ke);
            lam = 1.0/3 + ((1-4*rrsq)*ke[0] + 4*(2*rrsq-1)*ke[1])*2/(9*M_PI);
          }
        }
        else if(d == 0) {  /* case 10: dead center */
          tt = 1.0 - rrsq;
          lam = (1.0-tt*TSQRT(tt))*2.0/3.0;
        }
        else {  /* cases 3 and 9 */
          ab = a*b;

          /* Elliptic integrals 1st and 2nd */
          y = (1-bsq)/(1-asq);
          ELLKE(y, ke);

          /* Elliptic integral 3rd */
          q = b / a;
          tpnk = 3*q * ELLPI(TABS(q), y);

          /* lambda^d */
          lam = 2*((1 - 5*dsq + rrsq + ab*ab)*ke[0] +
                   (1-asq)*(dsq + 7*rrsq - 4)*ke[1] + tpnk)
            / (9*M_PI*TSQRT(1-asq));

          if(a < 0)
            lam += 2.0/3.0;  /* heaviside(rr-d) term */
        }
      }
      else {  /* occulter overlaps disk: |1-rr| < d <= 1+rr or b > 1 */
        /* Cases 2, 7 and 8 from the paper */
        tt = 1 - rrsq + dsq;
        ts = rrsq + dsq;
        drr = d*rr;

        kap0 = rrsq*TACOS((ts - 1) / (2*drr));
        kap1 = TACOS(tt / (2*d));
        ksrt = TSQRT(dsq - 0.25*tt*tt);

        area = (kap0 + kap1 - ksrt) / M_PI;  /* Eq. 1 */
        eta  = (kap1 + (ts + dsq)*kap0 -
                0.5*(1 + 5*rrsq + dsq)*ksrt) / (2*M_PI);

        if(asq <= TSQRTMIN) {  /* case 7 and erratum 2: touches center */
          ELLKE(1 - 1.0/(4*rrsq), ke);
          lam = 1/3.0 - ((1-4*rrsq)*(3-8*rrsq)*ke[0]/rr -
                         16*(2*rrsq-1)*ke[1]*rr) / (9*M_PI);
        }
        else {  /* cases 2 and 8 */
          ab = a*b;

          /* Elliptic integrals 1st and 2nd */
          y = (bsq-1)/(bsq-asq);
          ELLKE(y, ke);

          /* Elliptic integral 3rd */
          q = 1.0 / a;
          tpnk = 3*q*b * ELLPI(TABS(q), y);

          /* lambda^d */
          lam = (((1-bsq)*(2*bsq+asq-3) + 3*ab*(bsq-2))*ke[0] +
                 4*drr*(dsq+7*rrsq-4)*ke[1] + tpnk)
            / (9*M_PI*TSQRT(drr));

          if(a < 0)
            lam += 2.0/3.0;  /* heaviside(rr-d) term */
        }
      }

      fecl = (area + s->uss*(lam-area) + s->u2*eta) * s->ldnorm;
      iecl = 1;
    }

    /* Take off the eclipse.  Assumes spots have the same limb
       darkening as the photosphere, seems reasonable. */
    s->ltot -= fecl*s->l + fecl*s->ol*s->fecs + area*s->rl;

    /* Handle light ratio, if requested */
    if(typ[p] == EB_OBS_LRAT)
      out[p] = ss.ltot / sp.ltot;
    else {
      /* Normalized final light */
      ltot = (sp.ltot + ss.ltot) * (1.0 - lth) + lth;
      
      if(typ[p] == EB_OBS_MAG) {
        if(ltot > 0) {
          out[p] = zp - 2.5 * TLOG10(ltot);
        }
        else
          out[p] = zp;
      }
      else
        out[p] = ltot;
    }

    if(ieclout)
      ieclout[p] = iecl;
  }
}

/* Fortran wrappers.  For simplicity we define both the single and
   double underscore variants, so we don't need to know how the
   Fortran compiler is set up. */

void FSUBA (double *parm, double *t, DATATYPE *ol1, DATATYPE *ol2,
            unsigned char *typ, DATATYPE *out, unsigned char *ieclout,
            int *flags, int *npt) {
  FUNC(parm, t, ol1, ol2,
       typ, out, ieclout,
       *flags, *npt);
}

void FSUBB (double *parm, double *t, DATATYPE *ol1, DATATYPE *ol2,
            unsigned char *typ, DATATYPE *out, unsigned char *ieclout,
            int *flags, int *npt) {
  FUNC(parm, t, ol1, ol2,
       typ, out, ieclout,
       *flags, *npt);
}

