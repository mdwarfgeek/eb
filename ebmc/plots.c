#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "ebmc.h"

#ifdef HAVE_PGPLOT

#include <cpgplot.h>

#define RVEXPAND 0.75

void init_plots (char *pgdev) {
  cpgopen(pgdev);
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
}

void close_plots (void) {
  cpgclos();
}

int do_plots (struct fit_parms *par,
	      FILE *ofp,
              char **filtnamelist, int nfiltname,
              int novertsep,
              char *errstr) {
  struct fit_data *dlist;
  int ndata;
  double *v;

  double phicont[4];

  double tmp;
  float x, y, y1, y2;
  double phi, phifold;
  float tmin, tmax, xmin = 0, xmax = 0, xpri, xsec, ymin, ymax, yrange;
  float residmin, residmax, residrange;
  float xrmag;
  float xminmag[2], xmaxmag[2], xminecl[2], xmaxecl[2];
  float ymagmin[2], ymagmax[2];
  float residmagmin[2], residmagmax[2];
  int ynotdone, ymagnotdone[2];

  float rvxmin, rvxmax, rvmin, rvmax, rvresidmin, rvresidmax;
  int irv, nrv;

  int meas, i;

  float lx[3001], ly[3001];
  double samp, dlx[3001], dly[3001];
  unsigned char iecl[3001];

  float vx1, vx2, vy1, vy2, vh, vw, vxch, vych, vpad;
  int iband, iplot, cyc, cycmin, cycmax, ncyc;
  int *icyclist = (int *) NULL, nicyc, nicycmax;

  float yzp;
  float ch;
  float durpri, dursec;

  int ieclall;
  int idat, ntmp;

  float xch, ych;
  char lab[64], htlab[64];
  int ipanel, npanel;

  double dlyn[3001];
  unsigned char typ[3001];
  float ypad;
  int ispot;

  float tx, txlist[10], ty, tylist[10], yy;
  char ttlist[10][64];
  int ii;

  float lcexpand    = 0.20;
  float residexpand = 0.75;
  int phasefold;

  dlist = par->dlist;
  ndata = par->ndata;
  v = par->vinit;

  /* Stack vertically or simply combine? */
  if(novertsep) {
    lcexpand = 0;
    residexpand = 0;
  }

  /* Decide if we should plot phase-folded by default based on
     rotation parameters */
  if(v[EB_PAR_ROT1] == 1.0 &&
     v[EB_PAR_ROT2] == 1.0)
    phasefold = 1;
  else
    phasefold = 0;
  
  /* Contact points and eclipse durations */
  eb_phicont(v[EB_PAR_ESINW], v[EB_PAR_ECOSW],
             v[EB_PAR_COSI], v[EB_PAR_RASUM],
             phicont);

  if(phicont[0] > phicont[1])
    phicont[0] -= 1.0;

  if(phicont[2] > phicont[3])
    phicont[2] -= 1.0;

  durpri = phicont[1]-phicont[0];
  dursec = phicont[3]-phicont[2];

  xpri = 0.5*(phicont[0]+phicont[1]);
  xsec = 0.5*(phicont[2]+phicont[3]);

  tprintf(ofp,
	  "Primary duration:   %.6f = %.2fh\n"
	  "Secondary duration: %.6f = %.2fh\n",
	  durpri, 24*durpri*v[EB_PAR_P],
	  dursec, 24*dursec*v[EB_PAR_P]);

  snprintf(htlab, sizeof(htlab),
           "HJD-T\\d0\\u (HJD-%.5lf)", par->hjdoff+v[EB_PAR_T0]);

  /* Range for magnified light curve plots */
  xrmag = MAX(durpri, dursec);

  xminmag[0] = xpri-xrmag;
  xmaxmag[0] = xpri+xrmag;
  xminmag[1] = xsec-xrmag;
  xmaxmag[1] = xsec+xrmag;

  if(durpri > 0) {
    xminecl[0] = xpri-0.5*durpri;
    xmaxecl[0] = xpri+0.5*durpri;
  }
  else {
    xminecl[0] = xpri;
    xmaxecl[0] = xpri;
  }

  if(dursec > 0) {
    xminecl[1] = xsec-0.5*dursec;
    xmaxecl[1] = xsec+0.5*dursec;
  }
  else {
    xminecl[1] = xsec;
    xmaxecl[1] = xsec;
  }

  ymagmin[0] = FLT_MAX;
  ymagmax[0] = -FLT_MAX;
  ymagmin[1] = FLT_MAX;
  ymagmax[1] = -FLT_MAX;
  residmagmin[0] = FLT_MAX;
  residmagmax[0] = -FLT_MAX;
  residmagmin[1] = FLT_MAX;
  residmagmax[1] = -FLT_MAX;
  cycmin = 1;
  cycmax = -1;
  ymagnotdone[0] = 1;
  ymagnotdone[1] = 1;

  rvxmin = FLT_MAX;
  rvxmax = -FLT_MAX;
  rvmin = FLT_MAX;
  rvmax = -FLT_MAX;
  rvresidmin = FLT_MAX;
  rvresidmax = -FLT_MAX;
  nrv = 0;

  /* Plot diagonstic and accumulate ranges for magnified */
  cpgsubp(1, 1);
  cpgsch(powf(3, -1.0/3.0));

  for(idat = 0; idat < ndata; idat++) {
    if(dlist[idat].obstype == OBS_LC)
      yzp = v[dlist[idat].pnorm] + dlist[idat].ymed;  /* subtract off zpt */
    else if(dlist[idat].obstype == OBS_RV) {
      yzp = 0;
      nrv++;
    }
    else
      continue;

    cpgpage();
    cpgvstd();
    cpgqvp(0, &vx1, &vx2, &vy1, &vy2);
    vh = (vy2 - vy1) / 5;
    
    cpgsvp(vx1, vx2, vy1+3.1*vh, vy1+5*vh);
    
    tmin = FLT_MAX;
    tmax = -FLT_MAX;
    xmin = FLT_MAX;
    xmax = -FLT_MAX;
    ymin = FLT_MAX;
    ymax = -FLT_MAX;
    ynotdone = 1;
    
    ieclall = 0;

    fit_func(par, idat, NULL,
             dlist[idat].hjd, dlist[idat].m, dlist[idat].iecl,
             dlist[idat].nmeas,
             0, 1, 1);
    
    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      /* Compute un-folded phase */
      tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
      phi = tmp / v[EB_PAR_P];
      phifold = phi - floor(phi);

      cyc = rint(phi-0.5*xsec);
      
      if(meas == 0 || tmp < tmin)
	tmin = tmp;
      if(meas == 0 || tmp > tmax)
	tmax = tmp;
      
      if(meas == 0 || phi < xmin)
	xmin = phi;
      if(meas == 0 || phi > xmax)
	xmax = phi;
      
      /* If data point included... */
      if(dlist[idat].iflag[meas]) {
	/* Accumulate y-range */
	y1 = dlist[idat].y[meas]-errscale(dlist+idat, v, meas)-yzp;
	y2 = dlist[idat].y[meas]+errscale(dlist+idat, v, meas)-yzp;
	
	if(ynotdone || y1 < ymin)
	  ymin = y1;
	if(ynotdone || y2 > ymax)
	  ymax = y2;

	ynotdone = 0;

	/* Check if anything in eclipse */
	if(dlist[idat].iecl[meas])
	  ieclall = 1;

	/* Check if in appropriate folded phase range */
	for(iplot = 0; iplot < 2; iplot++) {
	  if(dlist[idat].obstype == OBS_LC &&
	     ((phifold >= xminmag[iplot] && phifold <= xmaxmag[iplot]) ||
	      (phifold-1 >= xminmag[iplot] && phifold-1 <= xmaxmag[iplot]) ||
	      (phifold+1 >= xminmag[iplot] && phifold+1 <= xmaxmag[iplot]))) {
	    /* Accumulate cycle nos */
	    if((ymagnotdone[0] && ymagnotdone[1]) || cyc < cycmin)
	      cycmin = cyc;
	    if((ymagnotdone[0] && ymagnotdone[1]) || cyc > cycmax)
	      cycmax = cyc;
	    
	    /* Accumulate y-range for magnified */
	    if(ymagnotdone[iplot] || y1 < ymagmin[iplot])
	      ymagmin[iplot] = y1;
	    if(ymagnotdone[iplot] || y2 > ymagmax[iplot])
	      ymagmax[iplot] = y2;
	    
	    ymagnotdone[iplot] = 0;
	  }
	}

	if(dlist[idat].obstype == OBS_RV) {
	  rvxmin = MIN(phi, rvxmin);
	  rvxmax = MAX(phi, rvxmax);
	  rvmin = MIN(y1, rvmin);
	  rvmax = MAX(y2, rvmax);
	}
      }
    }
    
    /* Compute model and expand y-range if needed */
    samp = 1.0 / 1000;
    for(i = 0; i <= 1000; i++)
      dlx[i] = samp*i;

    fit_func(par, idat, NULL,
             dlx, dly, iecl, 1001,
             EB_FLAG_PHI, 1, 0);

    for(i = 0; i <= 1000; i++) {
      x = dlx[i];
      y = dly[i] - yzp;

      phifold = x - floor(x);

      if(dlist[idat].obstype == OBS_RV) {
        if(y < rvmin)
          rvmin = y;
        if(y > rvmax)
          rvmax = y;
      }

      /* Skip eclipses if data have none */
      if(dlist[idat].obstype != OBS_LC || ieclall || !iecl[i]) {
	if(ynotdone || y < ymin)
	  ymin = y;
	if(ynotdone || y > ymax)
	  ymax = y;

	ynotdone = 0;

	/* Check if in appropriate folded phase range */
	for(iplot = 0; iplot < 2; iplot++) {
	  if(dlist[idat].obstype == OBS_LC &&
	     ((phifold >= xminmag[iplot] && phifold <= xmaxmag[iplot]) ||
	      (phifold-1 >= xminmag[iplot] && phifold-1 <= xmaxmag[iplot]) ||
	      (phifold+1 >= xminmag[iplot] && phifold+1 <= xmaxmag[iplot]))) {
	    /* Accumulate y-range for magnified */
	    if(ymagnotdone[iplot] || y1 < ymagmin[iplot])
	      ymagmin[iplot] = y1;
	    if(ymagnotdone[iplot] || y2 > ymagmax[iplot])
	      ymagmax[iplot] = y2;
	    
	    ymagnotdone[iplot] = 0;
	  }
	}
      }
    }
    
    /* If plots are phase-folded, range is always the same */
    if(phasefold) {
      if(dlist[idat].obstype == OBS_LC) {
        xmin = -0.25;
        xmax = 0.75;
      }
      else {
        xmin = 0.0;
        xmax = 1.0;
      }

      rvxmin = 0.0;
      rvxmax = 1.0;
    }
    else {
      xmin -= 0.05*(xmax-xmin);
      xmax += 0.05*(xmax-xmin);
    }

    tmin -= 0.05*(tmax-tmin);
    tmax += 0.05*(tmax-tmin);

    yrange = ymax - ymin;
    ymin -= 0.05*yrange;
    ymax += 0.05*yrange;
    
    if(dlist[idat].obstype == OBS_LC) {
      cpgswin(tmin, tmax, ymax, ymin);
      cpgbox("BCMST", 0.0, 0, "BCNST", 0.0, 0);

      if(dlist[idat].obstype == OBS_LC)
	snprintf(lab, sizeof(lab),
                 "\\gD%s",
                 dlist[idat].iband >= nfiltname ?
                 "Mag" :
                 filtnamelist[dlist[idat].iband]);
      else
	snprintf(lab, sizeof(lab), "RV (km/s)");

      cpglab("", lab, htlab);
      
      /* Full model */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
               0, 1, 1);

      /* Model without systematics corrections for each data point */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].corr, NULL, dlist[idat].nmeas,
               0, 1, 0);
      
      /* Correction */
      for(meas = 0; meas < dlist[idat].nmeas; meas++)
        dlist[idat].corr[meas] = dlist[idat].m[meas] - dlist[idat].corr[meas];

      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
	/* Compute delta time */
	tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
	
	x = tmp;
	y = dlist[idat].y[meas]-yzp;
	y1 = y-errscale(dlist+idat, v, meas);
	y2 = y+errscale(dlist+idat, v, meas);
	
	if(!dlist[idat].iflag[meas])
	  cpgsci(2);
	
	//cpgpt(1, &x, &y, 2);
	cpgerry(1, &x, &y1, &y2, 1);
	
	cpgsci(1);
	
	dlist[idat].phitmp[meas] = tmp;
	dlist[idat].ytmp[meas] = dlist[idat].m[meas] - yzp;
      }
      
      cpgsci(2);
      cpgslw(4);
      cpgqch(&ch);
      cpgsch(0.75*ch);
      //cpgline(dlist[idat].nmeas, dlist[idat].phitmp, dlist[idat].ytmp);
      cpgpt(dlist[idat].nmeas, dlist[idat].phitmp, dlist[idat].ytmp, 17);
      cpgsch(ch);
      cpgslw(1);
      cpgsci(1);
    }
    else {
      /* Full model */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
               0, 1, 1);

      /* Model without systematics corrections for each data point */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].corr, NULL, dlist[idat].nmeas,
               0, 1, 0);
      
      /* Correction */
      for(meas = 0; meas < dlist[idat].nmeas; meas++)
        dlist[idat].corr[meas] = dlist[idat].m[meas] - dlist[idat].corr[meas];
    }

    cpgsvp(vx1, vx2, vy1+vh, vy1+2.9*vh);
    
    if(dlist[idat].obstype == OBS_LC)
      cpgswin(xmin, xmax, ymax, ymin);
    else
      cpgswin(xmin, xmax, ymin, ymax);
    cpgbox("BCST", 0.0, 0, "BCNST", 0.00, 0);
    cpglab("", dlist[idat].obstype == OBS_LC ? "Corrected" : "RV (km/s)", "");
    
    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      /* Compute phase */
      tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
      phi = tmp / v[EB_PAR_P];
      
      if(phasefold)
        x = phi - floor(phi-xmin);
      else
        x = phi;

      y = dlist[idat].y[meas] - yzp;
      
      /* Apply systematics correction */
      y -= dlist[idat].corr[meas];
      
      y1 = y-errscale(dlist+idat, v, meas);
      y2 = y+errscale(dlist+idat, v, meas);
      
      if(!dlist[idat].iflag[meas])
	cpgsci(2);
      
      //cpgpt(1, &x, &y, 2);
      cpgerry(1, &x, &y1, &y2, 1);
      
      cpgsci(1);
    }
    
    samp = (xmax-xmin) / 1000;
    ntmp = 0;
    for(i = 0; i <= 1000; i++)
      dlx[i] = xmin + samp*i;

    fit_func(par, idat, NULL,
             dlx, dly, iecl, 1001,
             EB_FLAG_PHI, 1, 0);

    for(i = 0; i <= 1000; i++) {
      lx[ntmp] = dlx[i];
      ly[ntmp] = dly[i] - yzp;
	
      /* Skip eclipses if data have none */
      if(dlist[idat].obstype != OBS_LC || ieclall || !iecl[i])
	ntmp++;
    }
    
    cpgsci(2);
    cpgslw(4);
    cpgline(ntmp, lx, ly);
    cpgslw(1);
    cpgsci(1);
    
    cpgsvp(vx1, vx2, vy1, vy1+vh);
    
    residmin = FLT_MAX;
    residmax = -FLT_MAX;
    
    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      /* Compute un-folded phase */
      tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
      phi = tmp / v[EB_PAR_P];
      phifold = phi - floor(phi);

      if(dlist[idat].iflag[meas]) {
	y1 = dlist[idat].resid[meas]-errscale(dlist+idat, v, meas);
	y2 = dlist[idat].resid[meas]+errscale(dlist+idat, v, meas);
	
	residmin = MIN(y1, residmin);
	residmax = MAX(y2, residmax);

	/* Check if in appropriate folded phase range */
	for(iplot = 0; iplot < 2; iplot++) {
	  if(dlist[idat].obstype == OBS_LC &&
	     ((phifold >= xminmag[iplot] && phifold <= xmaxmag[iplot]) ||
	      (phifold-1 >= xminmag[iplot] && phifold-1 <= xmaxmag[iplot]) ||
	      (phifold+1 >= xminmag[iplot] && phifold+1 <= xmaxmag[iplot]))) {
	    /* Accumulate range for magnified */
	    residmagmin[iplot] = MIN(y1, residmagmin[iplot]);
	    residmagmax[iplot] = MAX(y2, residmagmax[iplot]);
	  }
	}

	if(dlist[idat].obstype == OBS_RV) {
	  rvresidmin = MIN(y1, rvresidmin);
	  rvresidmax = MAX(y2, rvresidmax);
	}
      }
    }
    
    residrange = residmax - residmin;
    residmin -= 0.05*residrange;
    residmax += 0.05*residrange;
    
    if(dlist[idat].obstype == OBS_LC)
      cpgswin(xmin, xmax, residmax, residmin);
    else
      cpgswin(xmin, xmax, residmin, residmax);
    cpgbox("1BCNST", 0.0, 0, "BCNST", 0.0, 0);
    cpglab("Phase", "Residual", "");
    
    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      /* Compute phase */
      tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
      phi = tmp / v[EB_PAR_P];
      
      if(phasefold)
        x = phi - floor(phi-xmin);
      else
        x = phi;

      y = dlist[idat].resid[meas];
      y1 = y-errscale(dlist+idat, v, meas);
      y2 = y+errscale(dlist+idat, v, meas);
      
      if(!dlist[idat].iflag[meas])
	cpgsci(2);
      
      //cpgpt(1, &x, &y, 2);
      cpgerry(1, &x, &y1, &y2, 1);
      
      cpgsci(1);
    }
    
    lx[0] = xmin;
    lx[1] = xmax;
    ly[0] = 0.0;
    ly[1] = ly[0];
      
    cpgsci(2);
    cpgslw(4);
    cpgline(2, lx, ly);
    cpgslw(1);
    cpgsci(1);
  }

  /* Allocate list */
  ncyc = cycmax-cycmin+1;

  if(ncyc > 0) {
    icyclist = (int *) malloc(ncyc*sizeof(int));
    if(!icyclist) {
      report_syserr(errstr, "malloc");
      goto error;
    }

    /* Figure out max number of curves per panel */
    nicycmax = 0;

    ymin = MIN(ymagmin[0], ymagmin[1]);
    ymax = MAX(ymagmax[0], ymagmax[1]);
    residmin = MIN(residmagmin[0], residmagmin[1]);
    residmax = MAX(residmagmax[0], residmagmax[1]);
  
    yrange = ymax-ymin;
    residrange = residmax-residmin;

    npanel = 0;

    for(iband = 0; iband < par->nband; iband++) 
      for(iplot = 0; iplot < 2; iplot++) {
        /* Clear list */
        for(cyc = 0; cyc < ncyc; cyc++)
          icyclist[cyc] = -1;

        /* Loop through datasets and see who has data here */
        for(idat = 0; idat < ndata; idat++)
          if(dlist[idat].obstype == OBS_LC &&
             dlist[idat].iband == iband)
            for(meas = 0; meas < dlist[idat].nmeas; meas++) {
              /* Compute phase */
              tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
              phi = tmp / v[EB_PAR_P] - floor(tmp / v[EB_PAR_P]);
              cyc = rint(tmp / v[EB_PAR_P] - 0.5*xsec) - cycmin;

              /* Anything in plot range? */
              if((phi >= xminecl[iplot] && phi <= xmaxecl[iplot]) ||
                 (phi-1 >= xminecl[iplot] && phi-1 <= xmaxecl[iplot]) ||
                 (phi+1 >= xminecl[iplot] && phi+1 <= xmaxecl[iplot]))
                icyclist[cyc] = 1;
            }

        /* How many cycles have data? */
        nicyc = 0;

        for(cyc = 0; cyc < ncyc; cyc++)
          if(icyclist[cyc] >= 0)
            nicyc++;

        if(nicyc > nicycmax)
          nicycmax = nicyc;

        if(nicyc > 0)
          npanel++;

        /* Desired ymax and residmax */
        ymax = MAX(ymax, ymagmax[iplot]
             + lcexpand*(nicyc-1)*yrange);
        residmax = MAX(residmax, residmagmax[iplot]
                 + residexpand*(nicyc-1)*residrange);
      }
  
    ymin -= 0.05*yrange;
    ymax += 0.05*yrange;

    residmin -= 0.05*residrange;
    residmax += 0.05*residrange;

    if(npanel < 2)
      npanel = 2;  /* plot smaller */

    /* Plot magnified LC, each band separately */
    cpgpage();
    cpgvstd();
    cpgqvp(0, &vx1, &vx2, &vy1, &vy2);
    cpgsch(powf(npanel, -1.0/3.0));
    cpgqcs(0, &vxch, &vych);
    vpad = 3*vxch;

    vx1 -= vpad;
    vx2 += vpad;

    vh = (vy2 - vy1) / 5;
    vw = (vx2 - vx1) / npanel;

    ipanel = 0;

    for(iplot = 0; iplot < 2; iplot++) {
      for(iband = 0; iband < par->nband; iband++) {
        /* Clear list */
        for(cyc = 0; cyc < ncyc; cyc++)
          icyclist[cyc] = -1;

        /* Loop through datasets and see who has data here */
        for(idat = 0; idat < ndata; idat++)
          if(dlist[idat].obstype == OBS_LC &&
             dlist[idat].iband == iband)
            for(meas = 0; meas < dlist[idat].nmeas; meas++) {
              /* Compute phase */
              tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
              phi = tmp / v[EB_PAR_P] - floor(tmp / v[EB_PAR_P]);
              cyc = rint(tmp / v[EB_PAR_P] - 0.5*xsec) - cycmin;

              /* Anything in plot range? */
              if((phi >= xminecl[iplot] && phi <= xmaxecl[iplot]) ||
                 (phi-1 >= xminecl[iplot] && phi-1 <= xmaxecl[iplot]) ||
                 (phi+1 >= xminecl[iplot] && phi+1 <= xmaxecl[iplot]))
                icyclist[cyc] = 1;
            }

        /* Assign numbers */
        nicyc = 0;

        for(cyc = 0; cyc < ncyc; cyc++)
          if(icyclist[cyc] >= 0) {
            icyclist[cyc] = nicyc;
            nicyc++;
          }

        /* Check if there are any! */
        if(nicyc < 1)
          continue;

        xmin = xminmag[iplot];
        xmax = xmaxmag[iplot];
      
        cpgsvp(vx1+ipanel*vw+vpad, vx1+(ipanel+1)*vw-vpad, vy1+3*vh, vy1+5*vh);
      
        cpgswin(xmin, xmax, ymax, ymin);
        cpgbox("BCST", 0.0, 0, "BCNST", 0.0, 0);

        snprintf(lab, sizeof(lab),
                 "\\gD%s",
                 iband >= nfiltname ? "Mag" : filtnamelist[iband]);
        cpglab("", lab, "");
      
        cpgqcs(4, &xch, &ych);

        for(idat = 0; idat < ndata; idat++) {
          if(dlist[idat].obstype == OBS_LC && dlist[idat].iband == iband)
            yzp = v[dlist[idat].pnorm] + dlist[idat].ymed;  /* subtract zp */
          else
            continue;
	
          fit_func(par, idat, NULL,
                   dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
                   0, 1, 1);

          for(meas = 0; meas < dlist[idat].nmeas; meas++) {
            /* Compute phase */
            tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
            phi = tmp / v[EB_PAR_P] - floor(tmp / v[EB_PAR_P]);
            cyc = rint(tmp / v[EB_PAR_P] - 0.5*xsec) - cycmin;
	  
            x = phi;
            y = dlist[idat].y[meas] + lcexpand*yrange*icyclist[cyc] - yzp;
            y1 = y-errscale(dlist+idat, v, meas);
            y2 = y+errscale(dlist+idat, v, meas);
	  
            if(!dlist[idat].iflag[meas])
              cpgsci(2);
	
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            x = phi-1.0;
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            x = phi+1.0;
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            cpgsci(1);
	  
            if(phi >= xminmag[iplot] && phi <= xmaxmag[iplot])
              dlist[idat].phitmp[meas] = phi;
            else if(phi-1 >= xminmag[iplot] && phi-1 <= xmaxmag[iplot])
              dlist[idat].phitmp[meas] = phi-1;
            else
              dlist[idat].phitmp[meas] = phi+1;

            dlist[idat].ytmp[meas] = dlist[idat].m[meas]
                                   + lcexpand*yrange*icyclist[cyc] - yzp;
          }
	
          cpgsci(2);
          cpgslw(4);
          cpgqch(&ch);
          cpgsch(0.75*ch);
          //cpgline(dlist[idat].nmeas, dlist[idat].phitmp, dlist[idat].ytmp);
          cpgpt(dlist[idat].nmeas, dlist[idat].phitmp, dlist[idat].ytmp, 17);
          cpgsch(ch);
          cpgslw(1);
          cpgsci(1);
        }

        if(!novertsep) {
          for(cyc = cycmin; cyc <= cycmax; cyc++) {
            if(icyclist[cyc-cycmin] >= 0) {
              snprintf(lab, sizeof(lab), "%d", cyc);
              cpgptxt(xmax+0.5*xch,
                      lcexpand*yrange*icyclist[cyc-cycmin]-0.5*ych,
                      0.0, 0.0, lab);
            }
          }
        }

        cpgsvp(vx1+ipanel*vw+vpad, vx1+(ipanel+1)*vw-vpad, vy1+vh, vy1+3*vh);
      
        cpgswin(xmin, xmax, ymax, ymin);
        cpgbox("BCST", 0.0, 0, "BCNST", 0.00, 0);
        cpglab("", "Corrected", "");
      
        for(idat = 0; idat < ndata; idat++) {
          if(dlist[idat].obstype == OBS_LC && dlist[idat].iband == iband)
            yzp = v[dlist[idat].pnorm] + dlist[idat].ymed;  /* subtract zp */
          else
            continue;
	
          for(meas = 0; meas < dlist[idat].nmeas; meas++) {
            /* Compute phase */
            tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
            phi = tmp / v[EB_PAR_P] - floor(tmp / v[EB_PAR_P]);
            cyc = rint(tmp / v[EB_PAR_P] - 0.5*xsec) - cycmin;
	  
            x = phi;
            y = dlist[idat].y[meas] + lcexpand*yrange*icyclist[cyc] - yzp;
	  
            /* Apply systematics correction */
            y -= dlist[idat].corr[meas];
	  
            y1 = y-errscale(dlist+idat, v, meas);
            y2 = y+errscale(dlist+idat, v, meas);
	  
            if(!dlist[idat].iflag[meas])
              cpgsci(2);
	  
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            x = phi-1.0;
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            x = phi+1.0;
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            cpgsci(1);
          }
	
          cpgqcs(4, &xch, &ych);

          for(cyc = cycmin; cyc <= cycmax; cyc++) {
            if(icyclist[cyc-cycmin] >= 0) {
              samp = (xmax-xmin) / 1000;
              ntmp = 0;
              for(i = 0; i <= 1000; i++) {
                x = xmin+samp*i;
                dlx[i] = x+cyc;
                lx[i] = x;
              }
 
              fit_func(par, idat, NULL,
                       dlx, dly, iecl, 1001,
                       EB_FLAG_PHI, 1, 0);

              for(i = 0; i <= 1000; i++) {
                ly[ntmp] = dly[i]
                         + lcexpand*yrange*icyclist[cyc-cycmin] - yzp;
                ntmp++;
              }
	    
              cpgsci(2);
              cpgslw(4);
              cpgline(ntmp, lx, ly);
              cpgslw(1);
              cpgsci(1);

              if(!novertsep) {
                snprintf(lab, sizeof(lab), "%d", cyc);
                cpgptxt(xmax+0.5*xch,
                        lcexpand*yrange*icyclist[cyc-cycmin]-0.5*ych,
                        0.0, 0.0, lab);
              }
            }
          }
	
          if(iplot == 0) {
            lx[0] = phicont[0];
            lx[1] = phicont[0];
            lx[2] = phicont[1];
            lx[3] = phicont[1];
          }
          else if(iplot == 1) {
            lx[0] = phicont[2];
            lx[1] = phicont[2];
            lx[2] = phicont[3];
            lx[3] = phicont[3];
          }
          ly[0] = ymin;
          ly[1] = ymax;
          ly[2] = ymin;
          ly[3] = ly[1];
	
          if(lx[0] != lx[2]) {
            cpgline(2, &(lx[0]), &(ly[0]));
            cpgline(2, &(lx[2]), &(ly[2]));
          }
        }
      
        cpgsvp(vx1+ipanel*vw+vpad, vx1+(ipanel+1)*vw-vpad, vy1, vy1+vh);
      
        cpgswin(xmin, xmax, residmax, residmin);
        cpgbox("1BCNST", 0.0, 0, "BCNST", 0.0, 0);
        cpglab("Phase", "Residual", "");
      
        for(idat = 0; idat < ndata; idat++) {
          if(dlist[idat].obstype == OBS_LC && dlist[idat].iband == iband)
            yzp = v[dlist[idat].pnorm] + dlist[idat].ymed;  /* subtract zp */
          else
            continue;
	
          for(meas = 0; meas < dlist[idat].nmeas; meas++) {
            /* Compute phase */
            tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
            phi = tmp / v[EB_PAR_P] - floor(tmp / v[EB_PAR_P]);
            cyc = rint(tmp / v[EB_PAR_P] - 0.5*xsec) - cycmin;
	  
            x = phi;
            y = dlist[idat].resid[meas] + residexpand*residrange*icyclist[cyc];
            y1 = y-errscale(dlist+idat, v, meas);
            y2 = y+errscale(dlist+idat, v, meas);
	  
            if(!dlist[idat].iflag[meas])
              cpgsci(2);
	  
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            x = phi-1.0;
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            x = phi+1.0;
            //cpgpt(1, &x, &y, 2);
            cpgerry(1, &x, &y1, &y2, 1);
	  
            cpgsci(1);
          }

          cpgqcs(4, &xch, &ych);
	
          for(cyc = cycmin; cyc <= cycmax; cyc++) {
            if(icyclist[cyc-cycmin] >= 0) {
              lx[0] = xmin;
              lx[1] = xmax;
              ly[0] = residexpand*residrange*icyclist[cyc-cycmin];
              ly[1] = ly[0];
	    
              cpgsci(2);
              cpgslw(4);
              cpgline(2, lx, ly);
              cpgslw(1);
              cpgsci(1);

              if(!novertsep) {
                snprintf(lab, sizeof(lab), "%d", cyc);
                cpgptxt(xmax+0.5*xch,
                        residexpand*residrange*icyclist[cyc-cycmin]-0.5*ych,
                        0.0, 0.0, lab);
              }
            }
          }
        }

        ipanel++;
      }
    }

    /* First figure out where we have eclipses */
    for(cyc = 0; cyc < ncyc; cyc++)
      icyclist[cyc] = 0;

    for(idat = 0; idat < ndata; idat++)
      if(dlist[idat].obstype == OBS_LC)
        for(iplot = 0; iplot < 2; iplot++)
          for(meas = 0; meas < dlist[idat].nmeas; meas++) {
            /* Compute phase */
            tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
            phi = tmp / v[EB_PAR_P] - floor(tmp / v[EB_PAR_P]);
            cyc = rint(tmp / v[EB_PAR_P] - 0.5*xsec) - cycmin;
	  
            /* Anything in plot range? */
            if((phi >= xminecl[iplot] && phi <= xmaxecl[iplot]) ||
               (phi-1 >= xminecl[iplot] && phi-1 <= xmaxecl[iplot]) ||
               (phi+1 >= xminecl[iplot] && phi+1 <= xmaxecl[iplot]))
              icyclist[cyc] |= (1<<iplot);
          }

    /* Plot OOE */
    for(idat = 0; idat < ndata; idat++) {
      if(dlist[idat].obstype == OBS_LC)
        yzp = v[dlist[idat].pnorm] + dlist[idat].ymed;  /* subtract off zpt */
      else
        continue;

      tmin = FLT_MAX;
      tmax = -FLT_MAX;
      xmin = FLT_MAX;
      xmax = -FLT_MAX;
      ymin = FLT_MAX;
      ymax = -FLT_MAX;
      ynotdone = 1;
    
      ieclall = 0;
    
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].m,
               dlist[idat].iecl, dlist[idat].nmeas,
               0, 1, 1);

      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute un-folded phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
        phi = tmp / v[EB_PAR_P];
        phifold = phi - floor(phi);

        cyc = rint(phi-0.5*xsec);
      
        if(meas == 0 || tmp < tmin)
          tmin = tmp;
        if(meas == 0 || tmp > tmax)
          tmax = tmp;
      
        if(meas == 0 || phi < xmin)
          xmin = phi;
        if(meas == 0 || phi > xmax)
          xmax = phi;
      
        /* If data point included... */
        if(dlist[idat].iflag[meas]) {
          /* Accumulate y-range */
          y1 = dlist[idat].y[meas]-errscale(dlist+idat, v, meas)-yzp;
          y2 = dlist[idat].y[meas]+errscale(dlist+idat, v, meas)-yzp;
	
          if(ynotdone || y1 < ymin)
            ymin = y1;
          if(ynotdone || y2 > ymax)
            ymax = y2;

          ynotdone = 0;

          /* Check if anything in eclipse */
          if(dlist[idat].iecl[meas])
            ieclall = 1;
        }
      }

      /* Skip if has eclipses */
      if(ieclall)
        continue;

      /* Compute model and expand y-range if needed */
      samp = 1.0 / 1000;
      for(i = 0; i <= 1000; i++)
        dlx[i] = samp*i;

      fit_func(par, idat, NULL,
               dlx, dly, iecl, 1001,
               EB_FLAG_PHI, 1, 0);

      for(i = 0; i <= 1000; i++) {
        x = dlx[i];
        phifold = x - floor(x);

        y = dly[i] - yzp;
      
        /* Skip eclipses  */
        if(!iecl[i]) {
          if(ynotdone || y < ymin)
            ymin = y;
          if(ynotdone || y > ymax)
            ymax = y;

          ynotdone = 0;
        }
      }

      cpgpage();
      cpgvstd();
      cpgqvp(0, &vx1, &vx2, &vy1, &vy2);
      vh = (vy2 - vy1) / 5;
    
      cpgsvp(vx1, vx2, vy1+3.1*vh, vy1+5*vh);
    
      xmin -= 0.05*(xmax-xmin);
      xmax += 0.05*(xmax-xmin);

      tmin -= 0.05*(tmax-tmin);
      tmax += 0.05*(tmax-tmin);

      yrange = ymax - ymin;
      ymin -= 0.05*yrange;
      ymax += 0.05*yrange;
    
      cpgswin(tmin, tmax, ymax, ymin);
      cpgbox("BCMST", 0.0, 0, "BCNST", 0.0, 0);
    
      snprintf(lab, sizeof(lab),
               "\\gD%s",
               dlist[idat].iband >= nfiltname ?
               "Mag" :
               filtnamelist[dlist[idat].iband]);
      cpglab("", lab, htlab);
    
      /* Mark where we got eclipses */
      cpgslw(4);
      for(cyc = 0; cyc < ncyc; cyc++) {
        if(icyclist[cyc] & 0x01) {
          /* Primary */
          lx[0] = v[EB_PAR_P]*(cycmin+cyc);
          lx[1] = lx[0];
          ly[0] = ymin;
          ly[1] = ymin+0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
          ly[0] = ymax;
          ly[1] = ymax-0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
        }
        if(icyclist[cyc] & 0x02) {
          /* Secondary */
          lx[0] = v[EB_PAR_P]*(cycmin+cyc+xsec);
          lx[1] = lx[0];
          ly[0] = ymin;
          ly[1] = ymin+0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
          ly[0] = ymax;
          ly[1] = ymax-0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
        }
      }
      cpgslw(1);

      /* Full model */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
               0, 1, 1);

      /* Model without systematics corrections for each data point */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].corr, NULL, dlist[idat].nmeas,
               0, 1, 0);

      /* Correction */
      for(meas = 0; meas < dlist[idat].nmeas; meas++)
        dlist[idat].corr[meas] = dlist[idat].m[meas] - dlist[idat].corr[meas];

      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
      
        x = tmp;
        y = dlist[idat].y[meas]-yzp;
        y1 = y-errscale(dlist+idat, v, meas);
        y2 = y+errscale(dlist+idat, v, meas);
      
        if(!dlist[idat].iflag[meas])
          cpgsci(2);
      
        //cpgpt(1, &x, &y, 2);
        cpgerry(1, &x, &y1, &y2, 1);
      
        cpgsci(1);
      
        dlist[idat].phitmp[meas] = tmp;
        dlist[idat].ytmp[meas] = dlist[idat].m[meas] - yzp;
      }
    
      cpgsci(2);
      cpgslw(4);
      cpgqch(&ch);
      cpgsch(0.75*ch);
      //cpgline(dlist[idat].nmeas, dlist[idat].phitmp, dlist[idat].ytmp);
      cpgpt(dlist[idat].nmeas, dlist[idat].phitmp, dlist[idat].ytmp, 17);
      cpgsch(ch);
      cpgslw(1);
      cpgsci(1);

      cpgsvp(vx1, vx2, vy1+vh, vy1+2.9*vh);
    
      cpgswin(xmin, xmax, ymax, ymin);
      cpgbox("BCST", 0.0, 0, "BCNST", 0.00, 0);
      cpglab("",
             dlist[idat].obstype == OBS_LC ? "Corrected" : "RV (km/s)",
             "");
    
      /* Mark where we got eclipses */
      cpgslw(4);
      for(cyc = 0; cyc < ncyc; cyc++) {
        if(icyclist[cyc] & 0x01) {
          /* Primary */
          lx[0] = cycmin+cyc;
          lx[1] = lx[0];
          ly[0] = ymin;
          ly[1] = ymin+0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
          ly[0] = ymax;
          ly[1] = ymax-0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
        }
        if(icyclist[cyc] & 0x02) {
          /* Secondary */
          lx[0] = cycmin+cyc+xsec;
          lx[1] = lx[0];
          ly[0] = ymin;
          ly[1] = ymin+0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
          ly[0] = ymax;
          ly[1] = ymax-0.1*(ymax-ymin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
        }
      }
      cpgslw(1);

      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
        phi = tmp / v[EB_PAR_P];
      
        x = phi;
        y = dlist[idat].y[meas] - yzp;
      
        /* Apply systematics correction */
        y -= dlist[idat].corr[meas];

        y1 = y-errscale(dlist+idat, v, meas);
        y2 = y+errscale(dlist+idat, v, meas);
      
        if(!dlist[idat].iflag[meas])
          cpgsci(2);
      
        //cpgpt(1, &x, &y, 2);
        cpgerry(1, &x, &y1, &y2, 1);
      
        cpgsci(1);
      }
    
      samp = (xmax-xmin) / 1000;
      ntmp = 0;
      for(i = 0; i <= 1000; i++)
        dlx[i] = xmin + samp*i;

      fit_func(par, idat, NULL,
               dlx, dly, iecl, 1001,
               EB_FLAG_PHI | EB_FLAG_NOEC, 1, 0);

      for(i = 0; i <= 1000; i++) {
        lx[ntmp] = dlx[i];
        ly[ntmp] = dly[i] - yzp;
        ntmp++;
      }
    
      cpgsci(2);
      cpgslw(4);
      cpgline(ntmp, lx, ly);
      cpgslw(1);
      cpgsci(1);
    
      cpgsvp(vx1, vx2, vy1, vy1+vh);
    
      residmin = FLT_MAX;
      residmax = -FLT_MAX;
    
      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute un-folded phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
        phi = tmp / v[EB_PAR_P];
        phifold = phi - floor(phi);

        if(dlist[idat].iflag[meas]) {
          y1 = dlist[idat].resid[meas]-errscale(dlist+idat, v, meas);
          y2 = dlist[idat].resid[meas]+errscale(dlist+idat, v, meas);
	
          residmin = MIN(y1, residmin);
          residmax = MAX(y2, residmax);
        }
      }
    
      residrange = residmax - residmin;
      residmin -= 0.05*residrange;
      residmax += 0.05*residrange;
    
      cpgswin(xmin, xmax, residmax, residmin);
      cpgbox("1BCNST", 0.0, 0, "BCNSTX", 0.0, 0);
      cpglab("Phase", "Residual", "");
    
      /* Mark where we got eclipses */
      cpgslw(4);
      for(cyc = 0; cyc < ncyc; cyc++) {
        if(icyclist[cyc] & 0x01) {
          /* Primary */
          lx[0] = cycmin+cyc;
          lx[1] = lx[0];
          ly[0] = residmin;
          ly[1] = residmin+0.1*(residmax-residmin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
          ly[0] = residmax;
          ly[1] = residmax-0.1*(residmax-residmin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
        }
        if(icyclist[cyc] & 0x02) {
          /* Secondary */
          lx[0] = cycmin+cyc+xsec;
          lx[1] = lx[0];
          ly[0] = residmin;
          ly[1] = residmin+0.1*(residmax-residmin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
          ly[0] = residmax;
          ly[1] = residmax-0.1*(residmax-residmin);
          cpgarro(lx[0], ly[0], lx[1], ly[1]);
        }
      }
      cpgslw(1);

      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
        phi = tmp / v[EB_PAR_P];
      
        x = phi;
        y = dlist[idat].resid[meas];
        y1 = y-errscale(dlist+idat, v, meas);
        y2 = y+errscale(dlist+idat, v, meas);
      
        if(!dlist[idat].iflag[meas])
          cpgsci(2);
      
        //cpgpt(1, &x, &y, 2);
        cpgerry(1, &x, &y1, &y2, 1);
      
        cpgsci(1);
      }
    
      lx[0] = xmin;
      lx[1] = xmax;
      ly[0] = 0.0;
      ly[1] = ly[0];
      
      cpgsci(2);
      cpgslw(4);
      cpgline(2, lx, ly);
      cpgslw(1);
      cpgsci(1);
    }

    free((void *) icyclist);
    icyclist = (int *) NULL;
  }
  else {
    /* Set some reasonable sized defaults */
    npanel = 2;
  }

  if(nrv > 0) {
    /* Plot combined RV */
    ymin = rvmin;
    ymax = rvmax;
    yrange = ymax-ymin;
    
    ymin -= 0.05*yrange;
    ymax += 0.05*yrange;
    
    residmin = rvresidmin;
    residmax = rvresidmax;
    residrange = residmax-residmin;
    
    residmin -= 0.05*residrange;
    residmax += 0.05*residrange;
    
    residmin -= (nrv-1)*residrange;
    
    cpgpage();
    cpgvstd();
    cpgqvp(0, &vx1, &vx2, &vy1, &vy2);
    cpgsch(powf(npanel, -1.0/3.0));  /* same size as before */
    cpgqcs(0, &vxch, &vych);
    vpad = 3*vxch;
    
    vh = (vy2 - vy1) / 5;
    vw = (vx2 - vx1) / 1;
    
    xmin = rvxmin;
    xmax = rvxmax;
    
    if(!phasefold) {
      xmin -= 0.05*(rvxmax-rvxmin);
      xmax += 0.05*(rvxmax-rvxmin);
    }    

    cpgsvp(vx1+2.0*vpad/3.0, vx1+vw-2.0*vpad/3.0, vy1+3*vh, vy1+5*vh);
    
    cpgswin(xmin, xmax, ymin, ymax);
    cpgbox("BCST", 0.0, 0, "BCNST", 0.0, 0);
    cpglab("", "RV (km/s)", "");
    
    irv = 0;
    
    cpgqch(&ch);

    for(idat = 0; idat < ndata; idat++) {
      if(dlist[idat].obstype == OBS_RV)
        yzp = 0;
      else
        continue;
      
      /* Full model */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
               0, 1, 1);

      /* Model without systematics corrections for each data point */
      fit_func(par, idat, NULL,
               dlist[idat].hjd, dlist[idat].corr, NULL, dlist[idat].nmeas,
               0, 1, 0);

      /* Correction */
      for(meas = 0; meas < dlist[idat].nmeas; meas++)
        dlist[idat].corr[meas] = dlist[idat].m[meas] - dlist[idat].corr[meas];

      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
        phi = tmp / v[EB_PAR_P];

        if(phasefold)
          x = phi - floor(phi);
        else
          x = phi;

        y = dlist[idat].y[meas] - yzp;

        /* Apply exposure time correction */
        y -= dlist[idat].corr[meas];
        
        y1 = y-errscale(dlist+idat, v, meas);
        y2 = y+errscale(dlist+idat, v, meas);
        
        if(!dlist[idat].iflag[meas])
          cpgsci(2);
        
        cpgsch(2.5*ch);
        cpgpt(1, &x, &y, 16+irv);
        cpgsch(ch);
        cpgerry(1, &x, &y1, &y2, 1);
        
        cpgsci(1);
      }
      
      samp = (xmax-xmin) / 3000;
      for(i = 0; i <= 3000; i++)
        dlx[i] = xmin + samp*i;
      
      fit_func(par, idat, NULL,
               dlx, dly, NULL, 3001,
               EB_FLAG_PHI, 1, 0);
      
      for(i = 0; i <= 3000; i++) {
        lx[i] = dlx[i];
        ly[i] = dly[i];
      }
      
      cpgsci(2);
      cpgslw(4);
      cpgline(3001, lx, ly);
      cpgslw(1);
      cpgsci(1);
      
      irv++;
    }

    cpgsvp(vx1+2.0*vpad/3.0, vx1+vw-2.0*vpad/3.0, vy1+2*vh, vy1+3*vh);
    
    cpgswin(xmin, xmax, residmin, residmax);
    cpgbox("1BCNST", 0.0, 0, "BCNST", 0.0, 0);
    cpglab("Phase", "Residual", "");
    
    irv = 0;
    
    cpgqch(&ch);
    
    for(idat = 0; idat < ndata; idat++) {
      if(dlist[idat].obstype == OBS_RV)
        yzp = 0;
      else
        continue;
      
      for(meas = 0; meas < dlist[idat].nmeas; meas++) {
        /* Compute phase */
        tmp = dlist[idat].hjd[meas] - v[EB_PAR_T0];
        phi = tmp / v[EB_PAR_P];

        if(phasefold)
          x = phi - floor(phi);
        else
          x = phi;

        y = dlist[idat].resid[meas]-RVEXPAND*residrange*irv;
        y1 = y-errscale(dlist+idat, v, meas);
        y2 = y+errscale(dlist+idat, v, meas);
        
        if(!dlist[idat].iflag[meas])
          cpgsci(2);
        
        cpgsch(1.4*ch);
        cpgpt(1, &x, &y, 16+irv);
        cpgsch(ch);
        cpgerry(1, &x, &y1, &y2, 1);
        
        cpgsci(1);
      }
    
      lx[0] = xmin;
      lx[1] = xmax;
      ly[0] = -RVEXPAND*residrange*irv;
      ly[1] = ly[0];

      cpgsci(2);
      cpgslw(4);
      cpgline(2, lx, ly);
      cpgslw(1);
      cpgsci(1);
      
      cpgqcs(4, &xch, &ych);
      
      snprintf(lab, sizeof(lab), "%d", irv+1);
      cpgptxt(xmax+0.5*xch,
              -RVEXPAND*residrange*irv-0.5*ych,
              0.0, 0.0, lab);
      
      irv++;
    }
  }

  cpgpage();
  cpgvstd();
  cpgqvp(0, &vx1, &vx2, &vy1, &vy2);
  cpgsch(powf(npanel, -1.0/3.0));  /* same size as before */
  cpgqcs(0, &vxch, &vych);
  vpad = 3*vxch;
  ypad = 3*vych;

  vh = (vy2 - vy1) / 2;
  vw = (vx2 - vx1) / 2;

  for(iplot = 1; iplot <= 2; iplot++) {
    if(iplot == 1) {
      xmin = xpri-xrmag;
      xmax = xpri+xrmag;
    }
    else if(iplot == 2) {
      xmin = xsec-xrmag;
      xmax = xsec+xrmag;
    }

    ymin = -0.005;
    ymax = 0.025;

    cpgsvp(vx1+2.0*vpad/3.0, vx1+vw-2.0*vpad/3.0,
           vy2-iplot*vh+ypad, vy2-(iplot-1)*vh-ypad);
    
    cpgswin(xmin, xmax, ymax, ymin);
    cpgbox("1BCNST", 0.0, 0, "BCNST", 0.0, 0);
    cpglab("Phase", "Delta mag", "");

    cpgqcs(4, &xch, &ych);

    for(cyc = 4; cyc < 10; cyc++) {
      for(ispot = 1; ispot >= 0; ispot--) {
	samp = (xmax-xmin) / 1000;
	for(i = 0; i <= 1000; i++) {
	  x = xmin + samp*i;

          dlx[i] = x+cyc;
          lx[i] = x;
          typ[i] = EB_OBS_MAG;
        }

        eb_model_dbl(v, dlx, NULL, NULL, typ, dly, NULL, EB_FLAG_PHI, 1001);
        
        memcpy(par->vtmp, v, sizeof(par->vtmp));
        par->vtmp[EB_PAR_OOE1O] = 0;
        par->vtmp[EB_PAR_OOE11A] = 0;
        par->vtmp[EB_PAR_OOE11B] = 0;
        par->vtmp[EB_PAR_OOE12A] = 0;
        par->vtmp[EB_PAR_OOE12B] = 0;
        par->vtmp[EB_PAR_OOE2O] = 0;
        par->vtmp[EB_PAR_OOE21A] = 0;
        par->vtmp[EB_PAR_OOE21B] = 0;
        par->vtmp[EB_PAR_OOE22A] = 0;
        par->vtmp[EB_PAR_OOE22B] = 0;

        eb_model_dbl(par->vtmp, dlx, NULL, NULL, typ, dlyn, NULL, EB_FLAG_PHI, 1001);

        for(i = 0; i <= 1000; i++)
	  ly[i] = dly[i]-dlyn[i]; // - 0.02*ispot;

	//cpgsci(1+ispot);
	cpgsls(1+ispot);
	cpgline(1001, lx, ly);
	cpgsls(1);
	//cpgsci(1);
      }

      snprintf(lab, sizeof(lab), "%d", cyc);
      tx = xmax+0.5*xch;
      ty = ly[1000]-0.25*ych;

      for(ii = 4; ii < cyc; ii++) {
	if(ttlist[ii][0] != '\0') {
	  if(fabsf(ty - tylist[ii]) < 0.1*fabsf(ych)) {
	    /* Concatenate */
	    if(strlen(ttlist[ii])+strlen(lab)+2 < sizeof(ttlist[ii])) {
	      strcat(ttlist[ii], ",");
	      strcat(ttlist[ii], lab);
	    }
	    else
	      warning("buffer not large enough");

	    lab[0] = '\0';

	    break;
	  }
	  else if(fabsf(ty - tylist[ii]) < fabsf(ych)) {
	    /* Offset */
	    yy = 0.5*(ty+tylist[ii]);
	    if(ty < tylist[ii]) {
	      ty = yy - fabsf(ych)/2;
	      tylist[ii] = yy + fabsf(ych)/2;
	    }
	    else {
	      ty = yy + fabsf(ych)/2;
	      tylist[ii] = yy - fabsf(ych)/2;
	    }
	  }
	}
      }

      txlist[cyc] = tx;
      tylist[cyc] = ty;
      strncpy(ttlist[cyc], lab, sizeof(ttlist[cyc]));
      ttlist[cyc][sizeof(ttlist[cyc])-1] = '\0';
    }

    for(ii = 4; ii < 10; ii++)
      if(ttlist[ii][0] != '\0' && tylist[ii] >= ymin && tylist[ii] <= ymax)
	cpgptxt(txlist[ii], tylist[ii], 0.0, 0.0, ttlist[ii]);

    if(iplot == 1) {
      lx[0] = phicont[0];
      lx[1] = phicont[0];
      lx[2] = phicont[1];
      lx[3] = phicont[1];
    }
    else if(iplot == 2) {
      lx[0] = phicont[2];
      lx[1] = phicont[2];
      lx[2] = phicont[3];
      lx[3] = phicont[3];
    }
    ly[0] = ymin;
    ly[1] = ymax;
    ly[2] = ymin;
    ly[3] = ly[1];
    
    if(lx[0] != lx[2]) {
      cpgline(2, &(lx[0]), &(ly[0]));
      cpgline(2, &(lx[2]), &(ly[2]));
    }
  }

  return(0);

 error:
  if(icyclist)
    free((void *) icyclist);

  return(-1);
}

#define NHIST 20

void plothist (double *a, int nsamp,
	       char *xtitle, int tx, int ty, int doinit,
	       double offset, double scfac,
               double *perc,
               double xbest, int havebest) {
  int isamp, bin;
  float xlist[NHIST], hist[NHIST], binsize, xl, xh, yh;
  float xtmp, xls, xmed, xhs;
  float xx[2], yy[2];

  /* Decide plotting range */
  xl = 0.0;
  xh = 0.0;

  for(isamp = 0; isamp < nsamp; isamp++) {
    xtmp = offset+a[isamp]*scfac;

    if(isamp == 0 || xtmp < xl)
      xl = xtmp;
    if(isamp == 0 || xtmp > xh)
      xh = xtmp;
  }

  xls = offset + perc[0]*scfac;
  xmed = offset + perc[1]*scfac;
  xhs = offset + perc[2]*scfac;

  xbest = offset + xbest*scfac;

  binsize = (xh-xl) / (NHIST-1);

  /* Initialise */
  for(bin = 0; bin < NHIST; bin++)
    hist[bin] = 0;

  /* Bin up */
  for(isamp = 0; isamp < nsamp; isamp++) {
    bin = floor((offset + a[isamp] * scfac - xl)/binsize);
    if(bin >= 0 && bin < NHIST)
      hist[bin]++;
  }

  /* Renormalise */
  yh = 0;

  for(bin = 0; bin < NHIST; bin++) {
    xlist[bin] = xl+(bin+0.5)*binsize;

    /* Unit integral */
    hist[bin] /= (nsamp*binsize);

    if(hist[bin] > yh)
      yh = hist[bin];
  }

  /* Plot */
  if(doinit) {
    cpgpage();
    cpgvstd();
  }

  cpgswin(xl, xh, 0, yh*1.05);
  if(doinit) {
    cpgbox(tx ? "BCNST" : "BCST", 0.0, 0,
           (doinit && ty) ? "BCNST" : "BCST", 0.0, 0);
    cpglab(tx ? xtitle : "", doinit ? "Probability density" : "", "");
  }

  cpgbin(NHIST, xlist, hist, 1);

  yy[0] = 0;
  yy[1] = yh*1.05;
  if(havebest) {
    cpgsci(3);
    xx[0] = xbest;
    xx[1] = xbest;
    cpgline(2, xx, yy);
  }
  cpgsci(2);
  xx[0] = xmed;
  xx[1] = xmed;
  cpgline(2, xx, yy);
  cpgsls(2);
  xx[0] = xls;
  xx[1] = xls;
  cpgline(2, xx, yy);
  xx[0] = xhs;
  xx[1] = xhs;
  cpgline(2, xx, yy);
  cpgsls(1);
  cpgsci(1);
}

#define NHISTX 20
#define NHISTY 20

static void plot2d (double *x, double *y, int nsamp,
                    char *xtitle, char *ytitle,
                    int tx, int ty, int doinit, int boxonly,
                    double offx, double scfacx,
                    double offy, double scfacy) {
  int isamp;
  float xl, xh, yl, yh, xtmp, ytmp;
  float hist2d[NHISTX*NHISTY];
  float binsizex = 0.0, binsizey = 0.0, bmax = 0.0;
  float tr[6];
  int bx, by;

  /* Decide plotting range */
  xl = 0.0;
  xh = 0.0;
  yl = 0.0;
  yh = 0.0;

  for(isamp = 0; isamp < nsamp; isamp++) {
    xtmp = offx+x[isamp]*scfacx;
    ytmp = offy+y[isamp]*scfacy;

    if(isamp == 0 || xtmp < xl)
      xl = xtmp;
    if(isamp == 0 || xtmp > xh)
      xh = xtmp;
    if(isamp == 0 || ytmp < yl)
      yl = ytmp;
    if(isamp == 0 || ytmp > yh)
      yh = ytmp;
  }

  if(!boxonly) {
    /* Compute bin sizes */
    binsizex = (xh-xl)/NHISTX;
    binsizey = (yh-yl)/NHISTY;
    
    /* Initialize map */
    for(bx = 0; bx < NHISTX; bx++)
      for(by = 0; by < NHISTY; by++)
	hist2d[by*NHISTX+bx] = 0.0;
    
    /* Bin up */
    for(isamp = 0; isamp < nsamp; isamp++) {
      xtmp = offx+x[isamp]*scfacx;
      ytmp = offy+y[isamp]*scfacy;
      
      bx = floor((xtmp - xl)/binsizex);
      by = floor((ytmp - yl)/binsizey);
      if(bx >= 0 && bx < NHISTX &&
	 by >= 0 && by < NHISTY)
	hist2d[by*NHISTX+bx]++;
    }
    
    /* Compute range */
    bmax = 0;
    
    for(bx = 0; bx < NHISTX; bx++)
      for(by = 0; by < NHISTY; by++) {
	hist2d[by*NHISTX+bx] /= nsamp;
	
	if(hist2d[by*NHISTX+bx] > bmax)
	  bmax = hist2d[by*NHISTX+bx];
      }
  }

  /* Plot */
  if(doinit) {
    cpgpage();
    cpgvstd();
  }
  cpgswin(xl, xh, yl, yh);

  if(!boxonly) {
    cpgsitf(2);
    
    tr[0] = xl-0.5*binsizex;
    tr[1] = binsizex;
    tr[2] = 0.0;
    tr[3] = yl-0.5*binsizey;
    tr[4] = 0.0;
    tr[5] = binsizey;
    cpggray(hist2d, NHISTX, NHISTY, 1, NHISTX, 1, NHISTY, bmax, 0.0, tr);
    
#if 0
    cpgsch(2.0);
    cpgwedg("RG", 1.0, 2.0, bmax, 0.0, "");
    cpgsch(1.4);
#endif
  }

  cpgbox(tx ? "BCNST" : "BCST", 0.0, 0, ty ? "BCNST" : "BCST", 0.0, 0);
  cpglab(tx ? xtitle : "", ty ? ytitle : "", "");
  

  if(!boxonly) {
#if 0
    for(isamp = 0; isamp < nsamp; isamp++) {
      xtmp = offx+x[isamp]*scfacx;
      ytmp = offy+y[isamp]*scfacy;
      
      cpgpt(1, &xtmp, &ytmp, 1);
    }
#endif
  }
}

void plot_fried_eggs (double *ainit, char **anames, int nvary,
                      double *mc_res, int nalloc, int nsimd,
                      double *aperc, double *abest) {
  float vx1, vx2, vy1, vy2, vw, vh;
  double *mc_xptr, *mc_yptr;
  int xparm, yparm;

  /* Fried egg matrix */
  cpgsubp(1, 1);
  cpgsch(powf(nvary, -1.0/3.0));

  cpgpage();
  cpgvstd();
  cpgwnad(0.0, 1.0, 0.0, 1.0);  /* set to 1:1 aspect ratio */

  cpgqvp(0, &vx1, &vx2, &vy1, &vy2);
  vw = (vx2 - vx1) / nvary;
  vh = (vy2 - vy1) / nvary;

  for(yparm = 0; yparm < nvary; yparm++) {
    mc_yptr = mc_res+yparm*nalloc;

    for(xparm = 0; xparm <= yparm; xparm++) {
      mc_xptr = mc_res+xparm*nalloc;

      cpgsvp(vx1+xparm*vw, vx1+(xparm+1)*vw,
             vy1+(nvary-1-yparm)*vh, vy1+(nvary-yparm)*vh);

      plot2d(mc_xptr, mc_yptr, nsimd,
             anames[xparm], anames[yparm],
             yparm == nvary-1, xparm == 0, 0, xparm == yparm,
             ainit[xparm], 1.0, ainit[yparm], 1.0);      
    }

    /* xparm = yparm */
    plothist(mc_yptr, nsimd,
             anames[yparm], yparm == 0, yparm == 0, 0,
             ainit[yparm], 1.0,
             aperc + 3*yparm,
             abest ? abest[yparm] : 0.0, abest ? 1 : 0);
  }
}

void plot_derived_hist (double *mc_der, int nalloc, int nsimd,
                        double *derperc, double *vderbest) {
  int iparm;
  double *mc_ptr;

  cpgsubp(3, 3);
  cpgsch(1.4);

  for(iparm = 0; iparm < EB_NDER; iparm++) {
    mc_ptr = mc_der+iparm*nalloc;

    plothist(mc_ptr, nsimd,
             eb_dernames[iparm], 1, 1, 1,
             0.0, 1.0,
             derperc + 3*iparm,
	     vderbest ? vderbest[iparm] : 0.0, vderbest ? 1 : 0);
  }
}

#else  /* !HAVE_PGPLOT */

/* PGPLOT not available so all plotting functions are no-ops */

void init_plots (char *pgdev) {

}

void close_plots (void) {

}

int do_plots (struct fit_parms *par,
	      FILE *ofp,
              char **filtnamelist, int nfiltname,
              char *errstr) {
  return(0);
}

void plot_fried_eggs (double *ainit, char **anames, int nvary,
                      double *mc_res, int nalloc, int nsimd,
                      double *aperc, double *abest) {

}

void plot_derived_hist (double *mc_der, int nalloc, int nsimd,
                        double *derperc, double *vderbest) {

}

#endif  /* HAVE_PGPLOT */
