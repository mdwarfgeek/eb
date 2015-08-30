#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef USE_LEVMAR
#include <levmar.h>
#else
#include <mpfit.h>
#endif

#include "ebmc.h"

#ifdef USE_LEVMAR
static void lm_func (double *a, double *ans, int nparfit, int nmeas,
                     void *data) {
#else
static int lm_func (int nmeas, int nparfit, double *a, double *ans,
                    double **deriv, void *data) {
#endif
  struct fit_parms *par;
  int idat, meas;
  int dptr, mptr;
  double wt;

#if 0
  int ipar;

  for(ipar = 0; ipar < nparfit; ipar++)
    printf(" %f", a[ipar]);

  printf("\n");
#endif

  par = (struct fit_parms *) data;

  /* Compute model */
  for(idat = 0; idat < par->ndata; idat++)
    fit_func(par, idat, a,
             par->dlist[idat].hjd, par->dlist[idat].m, NULL,
             par->dlist[idat].nmeas,
             0, 0, 1);

  /* Extract useful parts */
  for(meas = 0; meas < nmeas; meas++) {
    dptr = par->fit_dptr[meas];
    mptr = par->fit_mptr[meas];

    if(par->dlist[dptr].iflag[mptr]) {
      wt = par->fit_wt[meas];
#ifdef USE_LEVMAR
      ans[meas] = par->dlist[dptr].m[mptr] * wt;
#else
      ans[meas] = (par->fit_y[meas] - par->dlist[dptr].m[mptr]) * wt;
#endif      
    }
    else
      ans[meas] = 0.0;
  }

#ifndef USE_LEVMAR
  return 0;
#endif
}

int do_fit (struct fit_parms *par, int nsigma, FILE *ofp, char *errstr) {
  struct fit_data *dlist;
  int ndata;

  double *a = (double *) NULL, *acov = (double *) NULL;

  double *v;
  double vder[EB_NDER];

#ifdef USE_LEVMAR
  double info[LM_INFO_SZ];
#else
  mp_par *mpp = (mp_par *) NULL;
  mp_config mpc;
  mp_result mpr;
#endif

  int rv;

  int iaddband;

  double *fit_y = (double *) NULL, *fit_wt;
  int *fit_dptr = (int *) NULL, *fit_mptr;
  int idat, nmeasfit;

  int ivparm, iaparm, jvparm, japarm;

  double chisq;
  int nchisq;

  double scl;

  double phisec;
  int meas;

  int iter;

  dlist = par->dlist;
  ndata = par->ndata;
  v = par->vinit;

  /* Allocate parameter vectors for fits */
  a = (double *) malloc(par->nvaryf * sizeof(double));
  acov = (double *) malloc(par->nvaryf*par->nvaryf * sizeof(double));
  if(!a || !acov) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Scaling for period only for the moment (most badly scaled) */
  par->vscl[EB_PAR_P] = v[EB_PAR_P] / DAY;  /* rescales to DAY = initial */
  par->voff[EB_PAR_P] = (DAY - 1.0) * par->vscl[EB_PAR_P];  /* offset to 1 */

  /* Pack parameter vectors for fits */
  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1) {
      a[iaparm] = (v[ivparm] - par->voff[ivparm]) / par->vscl[ivparm];
      iaparm++;
    }

#ifndef USE_LEVMAR
  /* MPFIT parameter init */
  mpp = (mp_par *) calloc(par->nvaryf, sizeof(mp_par));
  if(!mpp) {
    report_syserr(errstr, "calloc");
    goto error;
  }

  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1) {
      /* Manually set some step sizes for parameters that can cause
         trouble for the autoscaler. */
      if(ivparm == EB_PAR_COSI)
        mpp[iaparm].step = 1.0e-3;
      else if(ivparm == EB_PAR_ROT1 ||
              ivparm == EB_PAR_ROT2)
        mpp[iaparm].step = 1.0e-6;

      iaparm++;
    }

  memset(&mpc, 0, sizeof(mpc));
  mpc.maxiter = 1000;
#endif

  /* Allocate temporary packed data vectors for fit */
  fit_y = (double *) malloc(par->nmeastot * 2 * sizeof(double));
  fit_dptr = (int *) malloc(par->nmeastot * 2 * sizeof(int));
  if(!fit_y || !fit_dptr) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  fit_wt = fit_y + par->nmeastot;
  fit_mptr = fit_dptr + par->nmeastot;

  /* Fit */
  par->fit_y = fit_y;
  par->fit_wt = fit_wt;
  par->fit_dptr = fit_dptr;
  par->fit_mptr = fit_mptr;

  for(iter = 0; iter < 5; iter++) {
    nmeasfit = 0;

    for(idat = 0; idat < ndata; idat++)
      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
	if((dlist[idat].yerr[meas] > 0 ||
	    (dlist[idat].obstype == OBS_RV && dlist[idat].fiterr)) &&
	   (iter == 0 || dlist[idat].sigresid == 0 || nsigma <= 0 ||
	    dlist[idat].obstype != OBS_LC ||
	    fabs(dlist[idat].resid[meas]-dlist[idat].medresid) < nsigma*dlist[idat].sigresid)) {
          /* errscale not fit so we can precompute weights */
          fit_wt[nmeasfit] = 1.0 / errscale(dlist+idat, v, meas);
#ifdef USE_LEVMAR
	  fit_y[nmeasfit] = dlist[idat].y[meas] * fit_wt[nmeasfit];
#else
	  fit_y[nmeasfit] = dlist[idat].y[meas];
#endif
	  fit_dptr[nmeasfit] = idat;
	  fit_mptr[nmeasfit] = meas;
	  nmeasfit++;

	  dlist[idat].iflag[meas] = 1;
	}
	else
	  dlist[idat].iflag[meas] = 0;
      }

#ifdef USE_LEVMAR
    rv = dlevmar_dif(lm_func, a, fit_y, par->nvaryf, nmeasfit, 1000,
		     (double *) NULL, info, (double *) NULL, acov,
		     (void *) par);

    tprintf(ofp, "l-m termination %.0f\n", info[6]);
    
    if(rv < 0)
      fprintf(stderr, "WARNING: l-m error %d\n", rv);
#else
    memset(&mpr, 0, sizeof(mpr));
    mpr.covar = acov;

    rv = mpfit(lm_func, nmeasfit, par->nvaryf, a,
               mpp, &mpc, (void *) par, &mpr);

    tprintf(ofp, "mpfit termination %d\n", rv);

    if(rv <= 0)
      fprintf(stderr, "WARNING: mpfit error %d\n", rv);
#endif

    for(iaparm = 0, ivparm = 0; ivparm < par->nparm; ivparm++)
      if(par->vary[ivparm] == 1)
	v[ivparm] = a[iaparm++] * par->vscl[ivparm] + par->voff[ivparm];

    /* Check for LD2 = LD1 case; copy if necessary */
    if(par->ldtype[1] == LD_SAME) {
      v[EB_PAR_LDLIN2] = v[EB_PAR_LDLIN1];
      v[EB_PAR_LDNON2] = v[EB_PAR_LDNON1];

      for(iaddband = 0; iaddband < par->naddband; iaddband++) {
	v[par->pldlin2[iaddband]] = v[par->pldlin1[iaddband]];
	v[par->pldnon2[iaddband]] = v[par->pldnon1[iaddband]];
      }
    }

    /* Light travel time parameter */
    if(par->vary[EB_PAR_CLTT] < 0)
      v[EB_PAR_CLTT] = v[PAR_KTOT]*1000 / EB_LIGHT;

    for(idat = 0; idat < ndata; idat++) {
      chisq = 0.0;
      nchisq = 0;

      fit_func(par, idat, a,
               dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
               0, 0, 1);
    
      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
	dlist[idat].resid[meas] = dlist[idat].y[meas] - dlist[idat].m[meas];
	dlist[idat].m[meas] = dlist[idat].resid[meas];
	
	if(dlist[idat].iflag[meas]) {
	  if(dlist[idat].obstype == OBS_RV)
	    chisq += dlist[idat].resid[meas]*dlist[idat].resid[meas] / (dlist[idat].yerr[meas]*dlist[idat].yerr[meas]+dlist[idat].errguess*dlist[idat].errguess);
	  else
	    chisq += dlist[idat].resid[meas]*dlist[idat].resid[meas] / (dlist[idat].yerr[meas]*dlist[idat].yerr[meas]);
	  nchisq++;
	}
      }

      dlist[idat].chisq = chisq;
      dlist[idat].nchisq = nchisq;

      if(dlist[idat].nmeas > 1)
	dmedsig(dlist[idat].m, dlist[idat].nmeas,
                &(dlist[idat].medresid), &(dlist[idat].sigresid));
      else {
	dlist[idat].medresid = 0;
	dlist[idat].sigresid = 0;
      }
    }

    /* Compute scale factors - LC only, RV will be done in MC */
    for(idat = 0; idat < ndata; idat++) {
      if(dlist[idat].fiterr && dlist[idat].obstype == OBS_LC) {
	if(dlist[idat].nchisq > 10*par->nvaryflc)  /* x10 is a bit of a kludge */
	  v[dlist[idat].perrsc]
            = sqrt(dlist[idat].chisq / (dlist[idat].nchisq-par->nvaryflc));
	else
	  v[dlist[idat].perrsc] = 1.0;
      }
    }
  }

  for(idat = 0; idat < ndata; idat++) {
    if(dlist[idat].obstype == OBS_LC && dlist[idat].nchisq > par->nvaryflc)
      tprintf(ofp,
	      "chisq = %.3f ndof = %.1f reduced chisq = %.3f\n"
	      "multiply errors by %.3f\n"
	      "scatter of residuals = %.4f\n",
	      dlist[idat].chisq, dlist[idat].nchisq-par->nvaryflc,
	      dlist[idat].chisq/(dlist[idat].nchisq-par->nvaryflc),
	      sqrt(dlist[idat].chisq/(dlist[idat].nchisq-par->nvaryflc)),
	      dlist[idat].sigresid);
    else if(dlist[idat].obstype == OBS_RV && dlist[idat].nchisq > par->nvaryfrv)
      tprintf(ofp,
	      "chisq = %.3f ndof = %.1f reduced chisq = %.3f\n"
	      "multiply errors by %.3f\n"
	      "scatter of residuals = %.4f\n",
	      dlist[idat].chisq, dlist[idat].nchisq-par->nvaryfrv,
	      dlist[idat].chisq/(dlist[idat].nchisq-par->nvaryfrv),
	      sqrt(dlist[idat].chisq/(dlist[idat].nchisq-par->nvaryfrv)),
	      dlist[idat].sigresid);
  }

  fflush(ofp);

  tprintf(ofp, "Initial parameters:\n");

  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++) {
    tprintf(ofp, "%-10s %14.6f",
            par->vnames[ivparm],
            (ivparm == EB_PAR_T0 ? v[ivparm]+par->hjdoff : v[ivparm]));
    if(par->vary[ivparm] == 1) {
      tprintf(ofp, " %9.6f",
              sqrtf(acov[iaparm*par->nvaryf+iaparm]) * par->vscl[ivparm]);
      iaparm++;
    }
    else if(par->vary[ivparm] == -1)
      tprintf(ofp, " (computed)");
    else if(par->vary[ivparm] == 0) {
      if(par->ldtype[1] == LD_SAME &&
         (ivparm == EB_PAR_LDLIN2 || ivparm == EB_PAR_LDNON2))
	tprintf(ofp, " (same)");
      else
	tprintf(ofp, " (fixed)");
    }
    else if(par->vary[ivparm] == 2)
      tprintf(ofp, " (fixed initially)");

    tprintf(ofp, " %s\n", par->vunits[ivparm]);
  }

  tprintf(ofp, "Derived parameters:\n");

  eb_getvder(v, v[PAR_GAMMA], v[PAR_KTOT], vder);

  for(ivparm = 0; ivparm < EB_NDER; ivparm++)
    tprintf(ofp, "%-10s %14.6f %s\n",
            eb_dernames[ivparm],
            (ivparm == EB_PAR_TSEC ? vder[ivparm]+par->hjdoff : vder[ivparm]),
            eb_derunits[ivparm]);

  fflush(stdout);

  phisec = eb_phisec(v[EB_PAR_ESINW], v[EB_PAR_ECOSW]);

  tprintf(ofp, "Secondary minimum at %.4f = HJD %.6f\n",
	  phisec, par->hjdoff+v[EB_PAR_T0]+v[EB_PAR_P]*phisec);

  /* Update master parameter vector and covariance matrix */
  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++) {
    if(par->vary[ivparm] == 1) {
      scl = par->vscl[ivparm];

      v[ivparm] = a[iaparm]*scl + par->voff[ivparm];

      for(jvparm = 0, japarm = 0; jvparm < par->nparm; jvparm++) {
        if(par->vary[jvparm] == 1) {
          par->vcov[ivparm*par->nparm+jvparm]
            = acov[iaparm*par->nvaryf+japarm]*scl*par->vscl[jvparm];
          japarm++;
        }
      }

      iaparm++;
    }
  }
  
  /* Reset scaling */
  par->vscl[EB_PAR_P] = 1.0;
  par->voff[EB_PAR_P] = 0.0;

  free((void *) a);
  a = (double *) NULL;
  free((void *) acov);
  acov = (double *) NULL;

#ifndef USE_LEVMAR
  free((void *) mpp);
  mpp = (mp_par *) NULL;
#endif

  free((void *) fit_y);
  fit_y = (double *) NULL;
  free((void *) fit_dptr);
  fit_dptr = (int *) NULL;

  return(0);

 error:
  if(a)
    free((void *) a);
  if(acov)
    free((void *) acov);

#ifndef USE_LEVMAR
  if(mpp)
    free((void *) mpp);
#endif

  if(fit_y)
    free((void *) fit_y);
  if(fit_dptr)
    free((void *) fit_dptr);

  return(-1);
}

