#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ebmc.h"

int make_parms (struct fit_parms *par, FILE *ofp,
                struct fit_data *dlist, int ndata,
                int ldtype[2],
                double *vfix, int *varyfix, double *vsgfix,
                struct fit_pset_entry *pextra, int npextra,
                char *errstr) {

  double *v = (double *) NULL, *voff, *vscl, *vsg, *vtmp;
  double *vcov = (double *) NULL;
  int *vary = (int *) NULL;
  char **vnames = (char **) NULL, **vunits;
  char *vnamesbuf = (char *) NULL, *vunitsbuf;
  int nparm;

  double hjdoff;

  int iaddband, naddband, nband;
  int *pj = (int *) NULL, *pldlin1, *pldlin2, *pldnon1, *pldnon2;
  int *pgenairm = (int *) NULL, *pgencm;

  int ipextra, ifound;

  int idat, meas, nmeastot, iseg;

  int iparm, jparm, nvaryf, nvarym;
  double nvaryflc, nvaryfrv;

  int ntmp;

  /* Init flags and compute median of input data */
  nmeastot = 0;

  for(idat = 0; idat < ndata; idat++) {
    /* Allocate buffer */
    dlist[idat].phitmp = (float *) malloc(2*dlist[idat].nmeas * sizeof(float));
    dlist[idat].resid = (double *) malloc(3*dlist[idat].nmeas * sizeof(double));
    dlist[idat].iflag
      = (unsigned char *) malloc(2*dlist[idat].nmeas * sizeof(unsigned char));

    if(!dlist[idat].phitmp || !dlist[idat].resid || !dlist[idat].iflag) {
      report_syserr(errstr, "malloc");
      goto error;
    }

    dlist[idat].ytmp = dlist[idat].phitmp + dlist[idat].nmeas;

    dlist[idat].m = dlist[idat].resid + dlist[idat].nmeas;
    dlist[idat].corr = dlist[idat].resid + 2*dlist[idat].nmeas;

    dlist[idat].iecl = dlist[idat].iflag + dlist[idat].nmeas;

    ntmp = 0;
    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      if(dlist[idat].yerr[meas] <= 0 &&
	 !(dlist[idat].obstype == OBS_RV && dlist[idat].fiterr))
	dlist[idat].iflag[meas] = 0;
      else
	dlist[idat].iflag[meas] = 1;
      
      dlist[idat].m[ntmp] = dlist[idat].y[meas];
      ntmp++;
    }

    dmedsig(dlist[idat].m, ntmp, &(dlist[idat].ymed), &(dlist[idat].ysig));

    tprintf(ofp, "Median %d: %.3f +/- %.3f\n",
            idat+1, dlist[idat].ymed, dlist[idat].ysig);

    nmeastot += dlist[idat].nmeas;
  }

  /* Subtract off T_0 initial guess to improve stability */
  hjdoff = vfix[EB_PAR_T0];  /* ref = first observation */

  vfix[EB_PAR_T0] -= hjdoff;
  for(idat = 0; idat < ndata; idat++)
    for(meas = 0; meas < dlist[idat].nmeas; meas++)
      if(dlist[idat].obstype != OBS_PRIOR)
	dlist[idat].hjd[meas] -= hjdoff;

  tprintf(ofp, "Subtracted HJD=%f\n", hjdoff);

  /* Count bands */
  naddband = 0;

  for(idat = 0; idat < ndata; idat++)
    if(dlist[idat].obstype == OBS_LC)
      naddband = MAX(dlist[idat].iband, naddband);

  nband = naddband+1;

  /* Figure out how many parameters there are */
  nparm = NPARFIX + 5*naddband;

  for(idat = 0; idat < ndata; idat++) {
    if(dlist[idat].obstype == OBS_LC) {
      if(dlist[idat].nseg > 0)
        nparm += dlist[idat].nseg;  /* zpt per segment */
      else
        nparm++;  /* for zpt */

      if(dlist[idat].fiterr)
	nparm++;

      if(dlist[idat].fitairm > 0)
	nparm++;

      if(dlist[idat].fitcm > 0)
	nparm++;
    }
    else if(dlist[idat].obstype == OBS_RV) {
      if(dlist[idat].fiterr)
	nparm++;
    }
  }

  /* Allocate master parameter vectors */
  v = (double *) malloc(5 * nparm * sizeof(double));
  vcov = (double *) malloc(nparm*nparm * sizeof(double));
  vary = (int *) malloc(nparm * sizeof(int));
  vnames = (char **) malloc(2 * nparm * sizeof(char *));
  vnamesbuf = (char *) malloc(2 * nparm * PARNAME_MAX * sizeof(char));
  pj = (int *) malloc(5 * naddband * sizeof(int));
  pgenairm = (int *) malloc(2 * nband * sizeof(int));
  if(!v || !vcov || !vary || !vnames || !vnamesbuf || !pj) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  voff = v + nparm;
  vscl = v + 2*nparm;
  vsg  = v + 3*nparm;
  vtmp = v + 4*nparm;

  vunits = vnames + nparm;
  vunitsbuf = vnamesbuf + nparm * PARNAME_MAX;

  pldlin1 = pj + naddband;
  pldlin2 = pj + 2*naddband;
  pldnon1 = pj + 3*naddband;
  pldnon2 = pj + 4*naddband;

  pgencm = pgenairm + nband;

  for(iparm = 0; iparm < nparm; iparm++) {
    vnames[iparm] = vnamesbuf + iparm * PARNAME_MAX;
    vunits[iparm] = vunitsbuf + iparm * PARNAME_MAX;  
  }

  /* Pack master parameter vectors */
  for(iparm = 0; iparm < NPARFIX; iparm++) {
    v[iparm] = vfix[iparm];
    vary[iparm] = varyfix[iparm];
    strncpy(vnames[iparm], parnames[iparm], PARNAME_MAX);
    strncpy(vunits[iparm], parunits[iparm], PARNAME_MAX);

    vscl[iparm] = 1.0;
    voff[iparm] = 0;
    vsg[iparm] = vsgfix[iparm];
  }

  for(iaddband = 0; iaddband < naddband; iaddband++) {
    v[iparm] = vfix[EB_PAR_J];
    vary[iparm] = varyfix[EB_PAR_J];
    snprintf(vnames[iparm], PARNAME_MAX, "J_%d", iaddband+1);
    strncpy(vunits[iparm], parunits[EB_PAR_J], PARNAME_MAX);

    vscl[iparm] = 1.0;
    voff[iparm] = 0;
    vsg[iparm] = vsgfix[EB_PAR_J];

    pj[iaddband] = iparm;
    iparm++;

    v[iparm] = vfix[EB_PAR_LDLIN1];
    vary[iparm] = varyfix[EB_PAR_LDLIN1];
    snprintf(vnames[iparm], PARNAME_MAX, "mu_1_f%d", iaddband+1);
    strncpy(vunits[iparm], parunits[EB_PAR_LDLIN1], PARNAME_MAX);

    vscl[iparm] = 1.0;
    voff[iparm] = 0;
    vsg[iparm] = vsgfix[EB_PAR_LDLIN1];

    pldlin1[iaddband] = iparm;
    iparm++;

    v[iparm] = vfix[EB_PAR_LDLIN2];
    vary[iparm] = varyfix[EB_PAR_LDLIN2];
    snprintf(vnames[iparm], PARNAME_MAX, "mu_2_f%d", iaddband+1);
    strncpy(vunits[iparm], parunits[EB_PAR_LDLIN2], PARNAME_MAX);

    vscl[iparm] = 1.0;
    voff[iparm] = 0;
    vsg[iparm] = vsgfix[EB_PAR_LDLIN2];

    pldlin2[iaddband] = iparm;
    iparm++;

    v[iparm] = vfix[EB_PAR_LDNON1];
    vary[iparm] = varyfix[EB_PAR_LDNON1];
    snprintf(vnames[iparm], PARNAME_MAX, "mup_1_f%d", iaddband+1);
    strncpy(vunits[iparm], parunits[EB_PAR_LDNON1], PARNAME_MAX);

    vscl[iparm] = 1.0;
    voff[iparm] = 0;
    vsg[iparm] = vsgfix[EB_PAR_LDNON1];

    pldnon1[iaddband] = iparm;
    iparm++;

    v[iparm] = vfix[EB_PAR_LDNON2];
    vary[iparm] = varyfix[EB_PAR_LDNON2];
    snprintf(vnames[iparm], PARNAME_MAX, "mup_2_f%d", iaddband+1);
    strncpy(vunits[iparm], parunits[EB_PAR_LDNON2], PARNAME_MAX);

    vscl[iparm] = 1.0;
    voff[iparm] = 0;
    vsg[iparm] = vsgfix[EB_PAR_LDNON2];

    pldnon2[iaddband] = iparm;
    iparm++;
  }

  for(idat = 0; idat < ndata; idat++) {
    if(dlist[idat].obstype == OBS_LC) {
      dlist[idat].pnorm = iparm;  /* record this now */

      if(dlist[idat].nseg > 0) {
        for(iseg = 0; iseg < dlist[idat].nseg; iseg++) {
          v[iparm] = 1.0e-3;
          vary[iparm] = 1;
          snprintf(vnames[iparm], PARNAME_MAX, "z_%d_%d", idat+1, iseg+1);
          strncpy(vunits[iparm], "mag", PARNAME_MAX);
          
          vscl[iparm] = 1.0;
          voff[iparm] = 0;
          vsg[iparm] = 1.0e-3;
          
          iparm++;
        }
      }
      else {
        v[iparm] = 1.0e-3;
        vary[iparm] = 1;
        snprintf(vnames[iparm], PARNAME_MAX, "z_%d", idat+1);
        strncpy(vunits[iparm], "mag", PARNAME_MAX);
        
        vscl[iparm] = 1.0;
        voff[iparm] = 0;
        vsg[iparm] = 1.0e-3;
        
        iparm++;
      }      

      if(dlist[idat].fiterr) {
	v[iparm] = 1.0;  /* scale factor for LCs */
	vary[iparm] = 2;
	snprintf(vnames[iparm], PARNAME_MAX, "s_%d", idat+1);
	strncpy(vunits[iparm], "", PARNAME_MAX);

	vscl[iparm] = 1.0;
	voff[iparm] = 0;
	vsg[iparm] = 1.0e-2;

	dlist[idat].perrsc = iparm;
	iparm++;
      }
      if(dlist[idat].fitairm > 0) {
	v[iparm] = 0.001;
	vary[iparm] = 1;
	snprintf(vnames[iparm], PARNAME_MAX, "k_%d", idat+1);
	strncpy(vunits[iparm], "mag/airmass", PARNAME_MAX);

	vscl[iparm] = 1.0;
	voff[iparm] = 0;
	vsg[iparm] = 1.0e-3;

	dlist[idat].pairm = iparm;
	pgenairm[dlist[idat].iband] = iparm;

	iparm++;
      }
      else if(dlist[idat].fitairm < 0) {
	/* Use generic one in same filter */
	dlist[idat].pairm = pgenairm[dlist[idat].iband];
      }
      if(dlist[idat].fitcm > 0) {
	v[iparm] = 1.0;
	vary[iparm] = 1;
	snprintf(vnames[iparm], PARNAME_MAX, "c_%d", idat+1);
	strncpy(vunits[iparm], "", PARNAME_MAX);

	vscl[iparm] = 1.0;
	voff[iparm] = 0;
	vsg[iparm] = 1.0e-2;

	dlist[idat].pcm = iparm;
	pgencm[dlist[idat].iband] = iparm;

	iparm++;
      }
      else if(dlist[idat].fitcm < 0) {
	/* Use generic one in same filter */
	dlist[idat].pcm = pgencm[dlist[idat].iband];
      }
    }
    else if(dlist[idat].obstype == OBS_RV) {
      if(dlist[idat].fiterr) {
	v[iparm] = dlist[idat].errguess;  /* quadrature for RVs */
	vary[iparm] = 2;
	snprintf(vnames[iparm], PARNAME_MAX, "s_%d", idat+1);
	strncpy(vunits[iparm], "km/s", PARNAME_MAX);

	vscl[iparm] = 1.0;
	voff[iparm] = 0;
	vsg[iparm] = 1.0e-1;

	dlist[idat].perrsc = iparm;
	iparm++;
      }
    }
  }

  /* Update initial values etc for non-fixed parameter set
     (e.g. extra filters) */
  for(ipextra = 0; ipextra < npextra; ipextra++) {
    /* Lookup name */
    ifound = -1;
    for(iparm = 0; iparm < nparm; iparm++)
      if(!strcmp(pextra[ipextra].name, vnames[iparm])) {
	ifound = iparm;
	break;
      }

    if(ifound >= 0) {
      if(ifound == EB_PAR_COSI)
	v[ifound] = cos(pextra[ipextra].v * M_PI / 180.0);
      else
	v[ifound] = pextra[ipextra].v;
 
      vary[ifound] = pextra[ipextra].vary;
	  
      if(pextra[ipextra].havesig)
	vsgfix[ifound] = pextra[ipextra].sig;
    }
  }

  /* Check for LD2 = LD1 case; disable varying LD2 if it's on */
  if(ldtype[1] == LD_SAME) {
    v[EB_PAR_LDLIN2] = v[EB_PAR_LDLIN1];
    v[EB_PAR_LDNON2] = v[EB_PAR_LDNON1];
    vary[EB_PAR_LDLIN2] = 0;
    vary[EB_PAR_LDNON2] = 0;
    
    for(iaddband = 0; iaddband < naddband; iaddband++) {
      v[pldlin2[iaddband]] = v[pldlin1[iaddband]];
      v[pldnon2[iaddband]] = v[pldnon1[iaddband]];
      vary[pldlin2[iaddband]] = 0;
      vary[pldnon2[iaddband]] = 0;
    }
  }

  /* Figure out how many parameters we are varying */
  for(iparm = 0, nvaryf = 0, nvaryflc = 0, nvarym = 0; iparm < nparm; iparm++)
    switch(vary[iparm]) {
    case 1:
      nvaryf++;
      nvaryflc++;
    case 2:
      nvarym++;
      break;
    }

  nvaryflc++;  /* for normalization etc. */

  nvaryfrv = 0;
  if(vary[EB_PAR_Q] == 1)
    nvaryfrv++;
  if(vary[EB_PAR_ESINW] == 1)
    nvaryfrv += 0.5;        /* well, sort of */
  if(vary[PAR_KTOT] == 1)
    nvaryfrv++;
  if(vary[PAR_GAMMA] == 1)
    nvaryfrv++;

  nvaryflc -= nvaryfrv;

  /* Initial guess at covariance matrix: diagonal */
  for(iparm = 0; iparm < nparm; iparm++) {
    for(jparm = 0; jparm < iparm; jparm++)
      vcov[iparm*nparm+jparm] = 0;

    vcov[iparm*nparm+iparm] = vsg[iparm]*vsg[iparm];

    for(jparm = iparm+1; jparm < nparm; jparm++)
      vcov[iparm*nparm+jparm] = 0;
  }

  /* Fill in structure */
  par->dlist = dlist;
  par->ndata = ndata;
  par->nmeastot = nmeastot;

  par->vinit = v;
  par->voff = voff;
  par->vscl = vscl;
  par->vsg = vsg;
  par->vtmp = vtmp;
  par->vcov = vcov;
  par->nparm = nparm;

  par->vary = vary;
  par->nvaryflc = nvaryflc;
  par->nvaryfrv = nvaryfrv;
  par->nvaryf = nvaryf;
  par->nvarym = nvarym;

  par->vnames = vnames;
  par->vunits = vunits;

  par->ldtype = ldtype;

  par->hjdoff = hjdoff;

  par->nband = nband;
  par->naddband = naddband;

  par->pj = pj;
  par->pldlin1 = pldlin1;
  par->pldlin2 = pldlin2;
  par->pldnon1 = pldnon1;
  par->pldnon2 = pldnon2;

  return(0);

 error:

  return(-1);
}
