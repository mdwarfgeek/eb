#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "ebmc.h"

#define NSAMP 29

void model_resamp (double *parm, double *t, double *ol1, double *ol2,
                   unsigned char *typ, double *y, unsigned char *iecl,
                   int flags, int nmeas,
                   float *texp) {
  double ttmp[nmeas], ytmp[nmeas];
  unsigned char iecltmp[nmeas];

  double dtwin, dt1, dtsamp, dt;
  int imeas, isamp;

  /* Size of window (relative to texp) in same units as "t" */
  if(flags & EB_FLAG_PHI)
    dtwin = 1.0 / (parm[EB_PAR_P] * EB_DAY);
  else
    dtwin = 1.0 / EB_DAY;

  /* Sampling */
  dt1 = -0.5*dtwin;
  dtsamp = dtwin / NSAMP;

  /* Initialize accumulators */
  for(imeas = 0; imeas < nmeas; imeas++) {
    y[imeas] = 0;

    if(iecl)
      iecl[imeas] = 0;
  }

  /* Process samples in parallel */
  for(isamp = 0; isamp < NSAMP; isamp++) {
    /* Sample point */
    dt = dt1 + dtsamp * (isamp+0.5);

    /* Compute times */
    for(imeas = 0; imeas < nmeas; imeas++)
      ttmp[imeas] = t[imeas] + dt * texp[imeas];

    /* Calculate model */
    eb_model_dbl(parm, ttmp, NULL, NULL, typ, ytmp, iecltmp, flags, nmeas);

    /* Accmuluate */
    for(imeas = 0; imeas < nmeas; imeas++) {
      y[imeas] += ytmp[imeas];

      if(iecl)
        iecl[imeas] |= iecltmp[imeas];
    }
  }

  /* Convert to average */
  for(imeas = 0; imeas < nmeas; imeas++)
    y[imeas] /= NSAMP;
}

void fit_func (struct fit_parms *par, int id,
               double *a,
               double *t, double *y, unsigned char *iecl, int nmeas,
               int flags, int ifull, int icor) {

  double *v;
  int iparm, oparm, meas, ot, inormoff;
  double scl;

  unsigned char typ[nmeas];

  ot = par->dlist[id].obstype;

  v = par->vtmp;

  if(a) {
    /* Form full parameter vector */
    for(iparm = 0, oparm = 0; oparm < par->nparm; oparm++) {
      if(par->vary[oparm] == 1 || (ifull && par->vary[oparm] == 2))
        v[oparm] = a[iparm++] * par->vscl[oparm] + par->voff[oparm];
      else
        v[oparm] = par->vinit[oparm];
    }
  }
  else {
    for(oparm = 0; oparm < par->nparm; oparm++)
      v[oparm] = par->vinit[oparm];
  }

  /* Sort out surface brightness ratio and LD for each filter */
  if(par->dlist[id].iband > 0) {
    v[EB_PAR_J] = v[par->pj[par->dlist[id].iband-1]];
    v[EB_PAR_L3] = v[par->pl3[par->dlist[id].iband-1]];

    v[EB_PAR_LDLIN1] = v[par->pldlin1[par->dlist[id].iband-1]];
    v[EB_PAR_LDLIN2] = v[par->pldlin2[par->dlist[id].iband-1]];

    v[EB_PAR_LDNON1] = v[par->pldnon1[par->dlist[id].iband-1]];
    v[EB_PAR_LDNON2] = v[par->pldnon2[par->dlist[id].iband-1]];
  }

  /* Check for LD2 = LD1 case; copy if necessary */
  if(par->ldtype[1] == LD_SAME) {
    v[EB_PAR_LDLIN2] = v[EB_PAR_LDLIN1];
    v[EB_PAR_LDNON2] = v[EB_PAR_LDNON1];
  }

  /* Light travel time parameter */
  if(par->vary[EB_PAR_KTOTC] < 0)
    v[EB_PAR_KTOTC] = v[PAR_KTOT]*1000 / EB_LIGHT;

  /* Priors are easy, no calculations needed */
  if(ot == OBS_PRIOR) {
    for(meas = 0; meas < nmeas; meas++)
      y[meas] = v[par->dlist[id].ipar[meas]];

    return;
  }
  /* else */

  if(ot == OBS_LC) {
    for(meas = 0; meas < nmeas; meas++)
      typ[meas] = EB_OBS_MAG;
  }
  else if(ot == OBS_RV) {
    if(par->dlist[id].component == 2) {
      for(meas = 0; meas < nmeas; meas++)
        typ[meas] = EB_OBS_VRAD2;
    }
    else {
      for(meas = 0; meas < nmeas; meas++)
        typ[meas] = EB_OBS_VRAD1;
    }
  }
  else if(ot == OBS_LRAT) {
    for(meas = 0; meas < nmeas; meas++)
      typ[meas] = t[meas] < 0 ? EB_OBS_AVLR : EB_OBS_LRAT;
  }

  /* Compute model */
  if(par->dlist[id].texp && icor)
    model_resamp(v, t, NULL, NULL, typ, y, iecl,
                 flags, nmeas,
                 par->dlist[id].texp);
  else
    eb_model_dbl(v, t, NULL, NULL, typ, y, iecl,
                 flags, nmeas);

  /* Compute outputs */
  if(ot == OBS_LC) {
    for(meas = 0; meas < nmeas; meas++) {
      if(icor && par->dlist[id].nseg > 0) {  /* segment no. to parameter */
        assert(par->dlist[id].iseg[meas] >= 0);
        inormoff = par->dlist[id].segtbl[par->dlist[id].iseg[meas]];
        assert(inormoff >= 0);
      }
      else
        inormoff = 0;

      y[meas] += v[par->dlist[id].pnorm+inormoff] + par->dlist[id].ymed;

      if(icor && par->dlist[id].fitairm)
        y[meas] += v[par->dlist[id].pairm]*(par->dlist[id].airmass[meas]-1.0);
      if(icor && par->dlist[id].fitcm)
        y[meas] += v[par->dlist[id].pcm]*par->dlist[id].cm[meas];
    }
  }
  else if(ot == OBS_RV) {
    if(par->dlist[id].component == 2)
      scl = -v[PAR_KTOT] / (1.0 + v[EB_PAR_Q]);
    else
      scl = v[PAR_KTOT] * v[EB_PAR_Q] / (1.0 + v[EB_PAR_Q]);

    for(meas = 0; meas < nmeas; meas++)
      y[meas] = v[PAR_GAMMA] + y[meas] * scl;
  }
}

