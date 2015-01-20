#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ebmc.h"

void fit_func (struct fit_parms *par, int id,
               double *a,
               double *t, double *y, unsigned char *iecl, int nmeas,
               int iphi, int ifull, int icor) {

  double *v;
  int iparm, oparm, meas, ot;
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
  if(par->vary[EB_PAR_CLTT] < 0)
    v[EB_PAR_CLTT] = v[PAR_KTOT]*1000 / LIGHT;

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
  eb_model_dbl(v, t, NULL, NULL, typ, y, iecl,
               iphi ? EB_FLAG_PHI : 0, nmeas);

  /* Compute outputs */
  if(ot == OBS_LC) {
    for(meas = 0; meas < nmeas; meas++) {
      y[meas] += v[par->dlist[id].pnorm] + par->dlist[id].ymed;

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

