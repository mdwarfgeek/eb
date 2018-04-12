#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#include "ebmc.h"

static inline int ldok (double *v, int ldtype[2], int sec);
static inline int ooeok (double *v, int sec);

/* This is shared and does the output (plots, printed, vector update) */
static int mc_out (struct fit_parms *par, double *ainit,
                   double *mc_res, double *mc_der, double *mc_lrat,
                   int nalloc, int nsimd,
                   double *abest, double *vderbest,
                   FILE *ofp, FILE *tfp, char *errstr) {
  double *vinit;

  double *mc_ptr;
  double *mc_tmp = (double *) NULL;

  size_t ind[3];
  double *aperc = (double *) NULL, *perc;
  double derperc[3*EB_NDER];

  char **anames = (char **) NULL;

  int ivparm, iaparm, ndp;
  double med, err;

  double theomega, domega;
  int isimd;

  char parstr[PARNAME_MAX];
  char *unit;

  double *lratperc = (double *) NULL;
  char **lratnames = (char **) NULL;
  char *lratnamebuf = (char *) NULL;
  int iband;

  vinit = par->vinit;

  /* Allocate workspace */
  mc_tmp = (double *) malloc(nsimd * sizeof(double));
  aperc = (double *) malloc(3*par->nvarym * sizeof(double));
  anames = (char **) malloc(par->nvarym * sizeof(char *));
  lratperc = (double *) malloc(3*par->nband * sizeof(double));
  lratnames = (char **) malloc(par->nband * sizeof(char *));
  lratnamebuf = (char *) malloc(par->nband * PARNAME_MAX);
  if(!mc_tmp || !aperc || !anames ||
     !lratperc || !lratnames || !lratnamebuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  for(iband = 0; iband < par->nband; iband++)
    lratnames[iband] = lratnamebuf + iband * PARNAME_MAX;

#ifdef ASYMM
  /* Set up desired percentiles */
  ind[0] = 0.1585 * nsimd;
  ind[1] = 0.5    * nsimd;
  ind[2] = 0.8415 * nsimd;
#else
  /* "New" way of doing it.  Analysis assumes nsimd is very large
     so the fractional samples we neglect don't matter. */
  ind[0] = 0.5    * nsimd;
  ind[1] = 0.683  * nsimd;
#endif

  /* Begin output */
  tprintf(ofp, "MC parameters:\n");
  fprintf(tfp,
          "\\hline\n"
          "MC parameters\\\\\n"
          "\\hline\n");

  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2) {
      mc_ptr = mc_res+iaparm*nalloc;
      memcpy(mc_tmp, mc_ptr, nsimd * sizeof(double));

      perc = &(aperc[3*iaparm]);

#ifdef ASYMM
      dmultquickselect(mc_tmp, nsimd,
                       ind, 3,
                       perc);

      tprintf(ofp,
	      "%-10s %.6f %.6f +%.6f\n", par->vnames[ivparm],
	      perc[1]+ainit[iaparm]+(ivparm == EB_PAR_T0 ? par->hjdoff : 0.0),
	      perc[0]-perc[1],
	      perc[2]-perc[1]);
#else
      med = dquickselect(mc_tmp, ind[0], nsimd);

      for(isimd = 0; isimd < nsimd; isimd++)
        mc_tmp[isimd] = fabs(mc_tmp[isimd] - med);

      err = dquickselect(mc_tmp, ind[1], nsimd);

      perc[0] = med - err;
      perc[1] = med;
      perc[2] = med + err;

      if(err > 0) {
        ndp = 1 - floor(log10(err));
        if(ndp < 0)
          ndp = 0;
      }
      else
        ndp = 0;

      tprintf(ofp,
	      "%-10s %.*f +/- %.*f\n", par->vnames[ivparm],
              ndp < 6 ? 6 : ndp,
	      med+ainit[iaparm]+(ivparm == EB_PAR_T0 ? par->hjdoff : 0.0),
              ndp < 6 ? 6 : ndp,
	      err);

      unit = par->vunits[ivparm];

      if(*unit)
        snprintf(parstr, sizeof(parstr),
                 "$%s$ (%s)", par->vtexsym[ivparm], unit);
      else
        snprintf(parstr, sizeof(parstr),
                 "$%s$", par->vtexsym[ivparm]);
      
      fprintf(tfp,
              "%-36s & $%.*f \\pm %.*f$ \\\\\n",
              parstr,
              ndp,
              med+ainit[iaparm]+(ivparm == EB_PAR_T0 ? par->hjdoff : 0.0),
              ndp,
              err);
#endif

      anames[iaparm] = par->vnames[ivparm];
      iaparm++;
    }

  /* Update full parameter vector only so we can compute omega */
  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2) {
      vinit[ivparm] = ainit[iaparm] + aperc[3*iaparm+1];  /* median */
      iaparm++;
    }

  /* Check for LD2 = LD1 case; copy if necessary */
  if(par->ldtype[1] == LD_SAME) {
    vinit[EB_PAR_LDLIN2] = vinit[EB_PAR_LDLIN1];
    vinit[EB_PAR_LDNON2] = vinit[EB_PAR_LDNON1];
  }

  /* Light travel time parameter */
  if(par->vary[EB_PAR_CLTT] < 0)
    vinit[EB_PAR_CLTT] = vinit[PAR_KTOT]*1000 / EB_LIGHT;

  /* Argument of periastron needs special treatment */
  theomega = atan2(vinit[EB_PAR_ESINW], vinit[EB_PAR_ECOSW]) * 180.0/M_PI;
  if(theomega < 0.0)  /* wrap to conventional range [0,360) */
    theomega += 360.0;

  fprintf(tfp,
          "\\hline\n"
          "Derived parameters\\\\\n"
          "\\hline\n");

  for(ivparm = 0; ivparm < EB_NDER; ivparm++) {
    mc_ptr = mc_der+ivparm*nalloc;
    memcpy(mc_tmp, mc_ptr, nsimd * sizeof(double));

    perc = &(derperc[3*ivparm]);

    /* Force correct revolution for omega by wrapping delta into (-180,180] */
    if(ivparm == EB_PAR_OMEGA) {
      for(isimd = 0; isimd < nsimd; isimd++) {
        domega = mc_tmp[isimd] - theomega;
        mc_tmp[isimd] = theomega + remainder(domega, 360.0);
      }
    }

#ifdef ASYMM
    dmultquickselect(mc_tmp, nsimd,
                     ind, 3,
                     perc);

    tprintf(ofp,
	    "%-10s %.6f %.6f +%.6f\n", eb_dernames[ivparm],
	    perc[1]+(ivparm == EB_PAR_TSEC ? par->hjdoff : 0.0),
	    perc[0]-perc[1],
	    perc[2]-perc[1]);
#else
    med = dquickselect(mc_tmp, ind[0], nsimd);

    for(isimd = 0; isimd < nsimd; isimd++)
      mc_tmp[isimd] = fabs(mc_tmp[isimd] - med);
    
    err = dquickselect(mc_tmp, ind[1], nsimd);
    
    perc[0] = med - err;
    perc[1] = med;
    perc[2] = med + err;

    if(err > 0) {
      ndp = 1 - floor(log10(err));
      if(ndp < 0)
        ndp = 0;
    }
    else
      ndp = 0;

    tprintf(ofp,
	    "%-10s %.*f +/- %.*f\n", eb_dernames[ivparm],
            ndp < 6 ? 6 : ndp,
	    med+(ivparm == EB_PAR_TSEC ? par->hjdoff : 0.0),
            ndp < 6 ? 6 : ndp,
	    err);

    unit = eb_derunits[ivparm];

    if(*unit) {
      if(!strcmp(unit, "Msol"))
        unit = "$\\msol$";
      else if(!strcmp(unit, "Rsol"))
        unit = "$\\rsol$";

      snprintf(parstr, sizeof(parstr),
               "$%s$ (%s)", eb_dertexsym[ivparm], unit);
    }
    else
      snprintf(parstr, sizeof(parstr),
               "$%s$", eb_dertexsym[ivparm]);

    fprintf(tfp,
	    "%-36s & $%.*f \\pm %.*f$ \\\\\n",
            parstr,
            ndp,
	    med+(ivparm == EB_PAR_TSEC ? par->hjdoff : 0.0),
            ndp,
	    err);
#endif
  }

  /* Light ratios */
  for(iband = 0; iband < par->nband; iband++) {
    mc_ptr = mc_lrat+iband*nalloc;
    memcpy(mc_tmp, mc_ptr, nsimd * sizeof(double));

    if(iband)
      snprintf(lratnames[iband], PARNAME_MAX,
               "L_2/L_1(%d)", iband);
    else
      snprintf(lratnames[iband], PARNAME_MAX,
               "L_2/L_1");

    perc = lratperc + 3*iband;

#ifdef ASYMM
    dmultquickselect(mc_tmp, nsimd,
                     ind, 3,
                     perc);

    tprintf(ofp,
	    "%-10s %.6f %.6f +%.6f\n", lratnames[iband],
	    perc[1]+(ivparm == EB_PAR_TSEC ? par->hjdoff : 0.0),
	    perc[0]-perc[1],
	    perc[2]-perc[1]);
#else
    med = dquickselect(mc_tmp, ind[0], nsimd);

    for(isimd = 0; isimd < nsimd; isimd++)
      mc_tmp[isimd] = fabs(mc_tmp[isimd] - med);
    
    err = dquickselect(mc_tmp, ind[1], nsimd);
    
    perc[0] = med - err;
    perc[1] = med;
    perc[2] = med + err;

    if(err > 0) {
      ndp = 1 - floor(log10(err));
      if(ndp < 0)
        ndp = 0;
    }
    else
      ndp = 0;

    tprintf(ofp,
	    "%-10s %.*f +/- %.*f\n", lratnames[iband],
            ndp < 6 ? 6 : ndp,
	    med+(ivparm == EB_PAR_TSEC ? par->hjdoff : 0.0),
            ndp < 6 ? 6 : ndp,
	    err);
#endif

  }

  /* Plots */
  plot_fried_eggs(ainit, anames, par->nvarym,
                  mc_res, nalloc, nsimd,
                  aperc, abest);

  plot_derived_hist(mc_der, nalloc, nsimd,
                    derperc, vderbest);

  plot_lrat_hist(mc_lrat, par->nband, nalloc, nsimd,
                 lratperc, lratnames);

  /* Update fit parameter vector now */
  for(iaparm = 0; iaparm < par->nvarym; iaparm++)
    ainit[iaparm] += aperc[3*iaparm+1];  /* median */

  free((void *) mc_tmp);
  mc_tmp = (double *) NULL;
  free((void *) aperc);
  aperc = (double *) NULL;
  free((void *) anames);
  anames = (char **) NULL;
  free((void *) lratperc);
  lratperc = (double *) NULL;
  free((void *) lratnames);
  lratnames = (char **) NULL;
  free((void *) lratnamebuf);
  lratnamebuf = (char *) NULL;

  return(0);

 error:
  if(mc_tmp)
    free((void *) mc_tmp);
  if(aperc)
    free((void *) aperc);
  if(anames)
    free((void *) anames);
  if(lratperc)
    free((void *) lratperc);
  if(lratnames)
    free((void *) lratnames);
  if(lratnamebuf)
    free((void *) lratnamebuf);

  return(-1);
}

/* This first due to emacs issues with indentation in "mc" */
int read_mc (struct fit_parms *par, FILE *ofp, FILE *tfp,
             char **filelist, int nfiles,
             char *errstr) {
  double *vinit;

  double *ainit = (double *) NULL, *a, *abest;

  FILE *mfp;
  int ifile, iwantv, iwanti;
  char line[16384], *p, *ep;
  double dtmp;
  int vtmp, rv;

  int ninit, nburn, nsim, nsimd, nalloc, isimd, isimthis, isim;
  int *nsimdfile = (int *) NULL, *nskipfile;
  int nvaryfile;

  double *mc_res = (double *) NULL, *mc_der;
  double *mc_lrat = (double *) NULL;
  int ivparm, iaparm;

  double *v = (double *) NULL;
  double *vltmp = (double *) NULL;
  double *voff = (double *) NULL;
  double vder[EB_NDER], vderbest[EB_NDER];

  double bestnlap = 0, nlap;

  double phitmp;
  unsigned char typtmp;
  int iband;

  vinit = par->vinit;

  /* Allocate parameter vector */
  v = (double *) malloc(par->nparm * sizeof(double));
  vltmp = (double *) malloc(par->nparm * sizeof(double));
  voff = (double *) malloc(par->nparm*nfiles * sizeof(double));
  if(!v || !vltmp || !voff) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Allocate lists to store file props */
  nsimdfile = (int *) malloc(2*nfiles * sizeof(int));
  if(!nsimdfile) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  nskipfile = nsimdfile + nfiles;

  /* Loop through file list reading headers */
  nsimd = 0;

  for(ifile = 0; ifile < nfiles; ifile++) {
    /* Open file */
    mfp = fopen(filelist[ifile], "r");
    if(!mfp) {
      report_syserr(errstr, "open: %s", filelist[ifile]);
      goto error;
    }

    /* Read interesting part of header */
    iwantv = 0;
    iwanti = 0;
    while(fgets(line, sizeof(line), mfp)) {
      /* We want comments */
      p = strchr(line, '#');
      if(!p)
	break;

      /* Now skip past it */
      p++;

      /* Clean up the mess */
      p = sstrip(p);

      if(*p == '\0')
	continue;

      if(iwantv == 1) {  /* Is this the line we want for the vary? */
	for(ivparm = 0, nvaryfile = 0; ivparm < par->nparm; ivparm++) {
	  vtmp = strtol(p, &ep, 0);
	  if(p == ep || !(isspace(*ep) || *ep == '\0')) {
	    report_err(errstr, "could not understand: %s", p);
	    goto error;
	  }
	  
	  if(vtmp == 1 || vtmp == 2)
	    nvaryfile++;

          if(ifile) {
            if(par->vary[ivparm] != vtmp) {
              report_err(errstr, "vary mismatch: %d != %d",
                         par->vary[ivparm], vtmp);
              goto error;
            }
          }
          else
            par->vary[ivparm] = vtmp;

	  p = ep;
	}

        par->nvarym = nvaryfile;

	iwantv = 2;
      }
      else if(iwanti == 1) {  /* Is this the line we want for the names? */
	for(ivparm = 0; ivparm < par->nparm; ivparm++) {
          /* Skip whitespace */
          while(*p != '\0' && isspace(*p))
            p++;
          
          /* Figure out how many non-whitespace */
          ep = p;
          while(*ep != '\0' && !isspace(*ep))
            ep++;
          
          /* Check name agrees */
          if(strncmp(par->vnames[ivparm], p, ep-p)) {
            report_err(errstr, "parameter name mismatch: %s != %s",
                       par->vnames[ivparm], p);
            goto error;
          }
          
          p = ep;
        }

	iwanti = 2;
      }
      else if(iwanti == 2) {  /* Is this the line we want for the vinit? */
	for(ivparm = 0; ivparm < par->nparm; ivparm++) {
	  dtmp = strtod(p, &ep);
	  if(p == ep || !isspace(*ep)) {
	    report_err(errstr, "could not understand: %s", p);
	    goto error;
	  }
	  
          voff[ifile*par->nparm+ivparm] = dtmp;

	  p = ep;
	}

        if(ifile == 0) {
          /* Init "best" */
          for(ivparm = 0; ivparm < EB_NDER; ivparm++) {
            dtmp = strtod(p, &ep);
            if(p == ep || !isspace(*ep)) {
              report_err(errstr, "could not understand: %s", p);
              goto error;
            }
            
            vderbest[ivparm] = dtmp;
            
            p = ep;
          }
          
          dtmp = strtod(p, &ep);
          if(p == ep || !isspace(*ep)) {  /* isspace ok */
            report_err(errstr, "could not understand: %s", p);
            goto error;
          }
          
          p = ep;
          
          bestnlap = strtod(p, &ep);
          if(p == ep || !(isspace(*ep) || *ep == '\0')) {
            report_err(errstr, "could not understand: %s", p);
            goto error;
          }
        }

	iwanti = 3;
      }
      else {
	/* Look for keywords */
	if(!strncmp(p, "Simulations:", 12)) {
	  p += 12;
	  rv = sscanf(p, "%d covinit %d burn %d total", &ninit, &nburn, &nsim);
	  if(rv != 3) {
	    report_err(errstr, "could not understand: %s", p);
	    goto error;
	  }
	}
        if(!strncmp(p, "hjdoff:", 7)) {
          p += 7;
          par->hjdoff = strtod(p, &ep);
          if(p == ep || !(isspace(*ep) || *ep == '\0')) {
            report_err(errstr, "could not understand: %s", p);
	    goto error;
	  }
        }
	else if(!strncmp(p, "Vary:", 5))
	  iwantv = 1;
        else if(!strncmp(p, "Initial", 7))
          iwanti = 1;
      }
    }

    nskipfile[ifile] = MAX(ninit, nburn);
    nsimdfile[ifile] = nsim - nskipfile[ifile];

    nsimd += nsimdfile[ifile];

    fclose(mfp);
  }

  /* Allocate memory */
  ainit = (double *) malloc(3*par->nvarym * sizeof(double));
  if(!ainit) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  a = ainit + par->nvarym;
  abest = ainit + 2*par->nvarym;

  nalloc = nsimd;

  mc_res = (double *) malloc(nalloc * (par->nvarym+EB_NDER) * sizeof(double));
  mc_lrat = (double *) malloc(nalloc * par->nband * sizeof(double));
  if(!mc_res || !mc_lrat) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  mc_der = mc_res + nalloc*par->nvarym;

  /* Pack parameter vectors for fits */
  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++) {
    vinit[ivparm] = voff[ivparm];
    
    if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2) {
      ainit[iaparm] = vinit[ivparm];
      abest[iaparm] = 0;
      iaparm++;
    }
  }

  /* Correction to system of first file */
  for(ifile = 0; ifile < nfiles; ifile++)
    for(ivparm = 0; ivparm < par->nparm; ivparm++)
      voff[ifile*par->nparm+ivparm] -= vinit[ivparm];

  /* Loop through file list reading data */
  isimd = 0;

  for(ifile = 0; ifile < nfiles; ifile++) {
    /* Open file */
    mfp = fopen(filelist[ifile], "r");
    if(!mfp) {
      report_syserr(errstr, "open: %s", filelist[ifile]);
      goto error;
    }

    /* Read away */
    isimthis = 0;

    while(fgets(line, sizeof(line), mfp)) {
      p = strchr(line, '\n');
      if(!p)  /* skip unterminated lines.  XXX - buffer undersizing! */
	continue;

      p = strchr(line, '#');
      if(p)
	*p = '\0';

      /* Clean up the whitespace */
      p = sstrip(line);

      if(*p == '\0')
	continue;

      /* Simulation number */
      isim = strtol(p, &ep, 0);
      if(p == ep || !isspace(*ep)) {
	report_err(errstr, "could not understand: %s", p);
	goto error;
      }

      isim--;  /* 0-indexed */

      if(isim >= nskipfile[ifile]) {
	p = ep;
	
	/* First check we aren't going to blow up */
	if(isimthis >= nsimdfile[ifile] || isimd >= nalloc) {
	  /* Something gone wrong */
	  report_err(errstr, "too many simulations at %d", isim+1);
	  goto error;
	}
	/* too few is not a problem */

	/* Data */
	for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++) {
          if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2) {
            dtmp = strtod(p, &ep);
            if(p == ep || !isspace(*ep)) {  /* isspace ok */
              report_err(errstr, "could not understand: %s", p);
              goto error;
            }
	  
            p = ep;

            /* Correct for any differences in initial parameters */
            dtmp += voff[ifile*par->nparm+ivparm];

	    mc_res[iaparm*nalloc+isimd] = dtmp;
            a[iaparm] = dtmp;
            v[ivparm] = vinit[ivparm] + dtmp;
	    iaparm++;
	  }
          else
            v[ivparm] = vinit[ivparm];
	}

        dtmp = strtod(p, &ep);
        if(p == ep || !isspace(*ep)) {  /* isspace ok */
          report_err(errstr, "could not understand: %s", p);
          goto error;
        }
	
        p = ep;

        nlap = strtod(p, &ep);
        if(p == ep || !(isspace(*ep) || *ep == '\0')) {
          report_err(errstr, "could not understand: %s", p);
          goto error;
        }
	
        /* Check for LD2 = LD1 case; copy if necessary */
        if(par->ldtype[1] == LD_SAME) {
          v[EB_PAR_LDLIN2] = v[EB_PAR_LDLIN1];
          v[EB_PAR_LDNON2] = v[EB_PAR_LDNON1];
        }
        
        /* Light travel time parameter */
        if(par->vary[EB_PAR_CLTT] < 0)
          v[EB_PAR_CLTT] = v[PAR_KTOT]*1000 / EB_LIGHT;

        /* Derived parameters */
        eb_getvder(v, v[PAR_GAMMA], v[PAR_KTOT], vder);

	for(ivparm = 0; ivparm < EB_NDER; ivparm++)
	  mc_der[ivparm*nalloc+isimd] = vder[ivparm];

        /* Light ratios */
        for(iband = 0; iband < par->nband; iband++) {
          memcpy(vltmp, v, par->nparm*sizeof(double));
        
          if(iband > 0) {
            vltmp[EB_PAR_J] = v[par->pj[iband-1]];
            vltmp[EB_PAR_LDLIN1] = v[par->pldlin1[iband-1]];
            vltmp[EB_PAR_LDLIN2] = v[par->pldlin2[iband-1]];
            vltmp[EB_PAR_LDNON1] = v[par->pldnon1[iband-1]];
            vltmp[EB_PAR_LDNON2] = v[par->pldnon2[iband-1]];
          }
          
          phitmp = 0;
          typtmp = EB_OBS_AVLR;
          eb_model_dbl(vltmp, &phitmp, NULL, NULL, &typtmp,
                       mc_lrat + iband*nalloc+isimd,
                       NULL, 0, 1);
        }

        if(nlap < bestnlap) {
          bestnlap = nlap;

          for(iaparm = 0; iaparm < par->nvarym; iaparm++)
            abest[iaparm] = a[iaparm];
          
          for(ivparm = 0; ivparm < EB_NDER; ivparm++)
            vderbest[ivparm] = vder[ivparm];
        }

	isimthis++;
	isimd++;
      }
    }

    fclose(mfp);
  }

  /* Update number of simulations in case files had less (which is fine) */
  nsimd = isimd;

  /* Now compute and populate parameter vectors */
  if(mc_out(par, ainit,
            mc_res, mc_der, mc_lrat,
            nalloc, nsimd,
            abest, vderbest,
            ofp, tfp, errstr))
     goto error;

  free((void *) ainit);
  ainit = (double *) NULL;
  free((void *) v);
  v = (double *) NULL;
  free((void *) vltmp);
  vltmp = (double *) NULL;
  free((void *) voff);
  voff = (double *) NULL;
  free((void *) mc_res);
  mc_res = (double *) NULL;
  free((void *) mc_lrat);
  mc_lrat = (double *) NULL;
  free((void *) nsimdfile);
  nsimdfile = (int *) NULL;

  return(0);

 error:
  if(ainit)
    free((void *) ainit);
  if(v)
    free((void *) v);
  if(vltmp)
    free((void *) vltmp);
  if(voff)
    free((void *) voff);
  if(mc_res)
    free((void *) mc_res);
  if(mc_lrat)
    free((void *) mc_lrat);
  if(nsimdfile)
    free((void *) nsimdfile);

  return(-1);
}

/* 
We briefly summarize the Metropolis-Hastings algorithm here.
Starting from an initial point in parameter space, the algorithm takes
the most recent set of parameters and perturbs one or more parameters
by a random Gaussian deviate.  If the perturbed parameter set has a
lower $\chi^2$ than its progenitor, it is accepted as a new point in
the chain.  If it has a larger $\chi^2$, it is accepted with a
probability $\exp(-\Delta \chi^2/2)$.  If it is not accepted, the
original point is repeated in the chain.  The size of the
perturbations were adjusted so that $20-30\%$ of the proposed points
were accepted.
*/

#define EPS 1.0e-100

int do_mc (struct fit_parms *par, FILE *ofp, FILE *tfp,
           char *filename, int nsim, int iseed, char *errstr) {
  struct fit_data *dlist;
  int ndata;
  double *vinit;

  FILE *mfp;
  struct rng_state rngs;

  double *ainit = (double *) NULL, *a, *atrial, *aatrial;
  double *amean, *ameanprev, *ameantmp, *abest;
  double *acov = (double *) NULL, *awork;

  double *v = (double *) NULL, *vltmp, *vtrial;

  int ivparm, iaparm, jvparm, japarm;

  double *mc_res = (double *) NULL, *mc_der;
  double *mc_lrat = (double *) NULL;
  double *mc_init = (double *) NULL, sum;

  int meas, isim, nburn, ninit, nskip, nsimd, nacc, itmp;
  double chisq, newchisq, nlap, bestnlap, esq, newesq, sd;
  double nll, nllnew, baye, bayenew, bf;
  double vder[EB_NDER], vderbest[EB_NDER];
  int idat, nchisq;
  int rvok;

  double rvvar, rvavsig;
  int nrvvar;
  double swt;

  int haveerr;

  double phitmp;
  unsigned char typtmp;
  int iband;

  dlist = par->dlist;
  ndata = par->ndata;
  vinit = par->vinit;

  nburn = nsim/2;
  ninit = MAX(1000, nsim/10);

  nskip = MAX(ninit, nburn);

  sd = 2.4*2.4 / par->nvarym;

  /* Allocate parameter vectors for fits */
  ainit = (double *) malloc(8*par->nvarym * sizeof(double));
  acov = (double *) malloc(2*par->nvarym*par->nvarym * sizeof(double));
  v = (double *) malloc(3*par->nparm * sizeof(double));
  if(!ainit || !acov || !v) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  a = ainit + par->nvarym;
  atrial = ainit + 2*par->nvarym;
  aatrial = ainit + 3*par->nvarym;
  amean = ainit + 4*par->nvarym;
  ameanprev = ainit + 5*par->nvarym;
  ameantmp = ainit + 6*par->nvarym;
  abest = ainit + 7*par->nvarym;

  awork = acov + par->nvarym*par->nvarym;

  vtrial = v + par->nparm;
  vltmp = v + 2*par->nparm;

  /* Pack parameter vectors for fits */
  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2) {
      ainit[iaparm] = vinit[ivparm];

      a[iaparm] = 0.0;  /* zero-mean system for better numerics */
      
      ameantmp[iaparm] = 0;
      
      amean[iaparm] = 0;
      ameanprev[iaparm] = 0;
      
      /* Initial covariance from Levenberg-Marquardt */
      for(jvparm = 0, japarm = 0; jvparm < par->nparm; jvparm++) {
        if(par->vary[jvparm] == 1 || par->vary[jvparm] == 2) {
          acov[iaparm*par->nvarym+japarm]
            = sd*par->vcov[ivparm*par->nparm+jvparm];
          japarm++;
        }
      }

      iaparm++;
    }

  /* Init output file */
  mfp = fopen(filename, "w");
  if(!mfp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  /* Compute avg RV variance for use in priors on semiamplitude.
   * This used to be purely empirical (scatter about fit), but
   * it was rather susceptible to the usual L-M convergence
   * problems, so changed to something purely a function of the
   * operator's inputs.
   */
  rvvar = 0;
  nrvvar = 0;

  for(idat = 0; idat < ndata; idat++)
    if(dlist[idat].obstype == OBS_RV) {
      /* Errors on individual points */
      for(meas = 0; meas < dlist[idat].nmeas; meas++)
	rvvar += dlist[idat].yerr[meas]*dlist[idat].yerr[meas];
      
      /* Initial guess of quadrature cpt */
      if(dlist[idat].fiterr)
	rvvar += dlist[idat].nmeas*dlist[idat].errguess*dlist[idat].errguess;
      
      nrvvar += dlist[idat].nmeas;
    }

  if(nrvvar > 0)
    rvvar /= nrvvar;

  rvavsig = sqrt(rvvar);  /* rough estimate */

  tprintf(ofp, "rvavsig = %.3f\n", rvavsig);

  /* Arrays to store results */
  mc_res = (double *) malloc(nsim * (par->nvarym+EB_NDER) * sizeof(double));
  mc_lrat = (double *) malloc(nsim * par->nband * sizeof(double));
  mc_init = (double *) malloc(ninit * par->nvarym * sizeof(double));
  if(!mc_res || !mc_lrat || !mc_init) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  mc_der = mc_res + nsim*par->nvarym;

  /* Initial residuals and chi^2 */
  chisq = 0;
  nchisq = 0;

  for(idat = 0; idat < ndata; idat++) {
    fit_func(par, idat, ainit,
             dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
             0, 1, 1);

    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      dlist[idat].resid[meas] = dlist[idat].y[meas] - dlist[idat].m[meas];
      if(dlist[idat].iflag[meas]) {
        chisq += dlist[idat].resid[meas]*dlist[idat].resid[meas] / varscale(dlist+idat, vinit, meas);
	nchisq++;
      }
    }
  }

  esq = vinit[EB_PAR_ECOSW]*vinit[EB_PAR_ECOSW] + 
	vinit[EB_PAR_ESINW]*vinit[EB_PAR_ESINW];  

  /* Negative log likelihood ("unlikelihood") */
  nll = chisq / 2;
  
  for(idat = 0; idat < ndata; idat++)
    if(dlist[idat].fiterr) {
      if(dlist[idat].obstype == OBS_LC) {
	/* Factor of 1/S^N in likelihood */
	nll += dlist[idat].nchisq*log(vinit[dlist[idat].perrsc]);
      }
      else if(dlist[idat].obstype == OBS_RV) {
	/* Factor of product of 1/sqrt(sigma_i^2+s^2) over all i */
	for(meas = 0; meas < dlist[idat].nmeas; meas++)
	  if(dlist[idat].iflag[meas])
	    nll += 0.5*log(dlist[idat].yerr[meas]*dlist[idat].yerr[meas]+vinit[dlist[idat].perrsc]*vinit[dlist[idat].perrsc] / dlist[idat].yinfwt[meas]);
      }
    }
  
  /* Now compute log of denominator of Bayes factor */
  baye = 0.0;
  
  /* Eccentricity - convert to uniform prior in e.
   * the "if" statement prevents the obvious KABOOM if circular.
   */
  if(esq != 0)
    baye += 0.5*log(esq);
  
  /* Modified Jeffreys prior in semiamplitude, transition at 10% of
     empirical scatter from above */
  baye += log(vinit[PAR_KTOT]+0.1*rvavsig);

  /* Jeffreys priors in error scaling parameters
     Modified for RV with transition at 10% of empirical scatter. */
  for(idat = 0; idat < ndata; idat++)
    if(dlist[idat].fiterr) {
      if(dlist[idat].obstype == OBS_LC)
	baye += log(vinit[dlist[idat].perrsc]);
      else if(dlist[idat].obstype == OBS_RV) {
	/* Sum of squared reciprocals */
	swt = 0;

	haveerr = 0;

	for(meas = 0; meas < dlist[idat].nmeas; meas++)
	  if(dlist[idat].iflag[meas]) {
	    if(dlist[idat].yerr[meas] > 0)
	      haveerr = 1;

	    swt += 1.0/(dlist[idat].yerr[meas]*dlist[idat].yerr[meas]+vinit[dlist[idat].perrsc]*vinit[dlist[idat].perrsc] / dlist[idat].yinfwt[meas]);
	  }

	/* Bayes factor - reduces to |s| when the sigmas are zero */
	if(haveerr)  /* modified Jeffreys prior */
	  baye -= log((vinit[dlist[idat].perrsc] +
                       0.1*dlist[idat].sigresid) * swt);
	else  /* Jeffreys prior */
	  baye -= log(vinit[dlist[idat].perrsc] * swt);
      }
    }

  for(ivparm = 0; ivparm < par->nparm; ivparm++)
    v[ivparm] = vinit[ivparm];
  
  eb_getvder(v, v[PAR_GAMMA], v[PAR_KTOT], vder);

  /* Init best pset */
  for(iaparm = 0; iaparm < par->nvarym; iaparm++)
    abest[iaparm] = a[iaparm];

  for(ivparm = 0; ivparm < EB_NDER; ivparm++)
    vderbest[ivparm] = vder[ivparm];

  bestnlap = nll + baye;

  /* Write out header to file */
  fprintf(mfp,
	  "# Initial solution:\n"
	  "#");

  for(ivparm = 0; ivparm < par->nparm; ivparm++)
    fprintf(mfp, " %s", par->vnames[ivparm]);

  for(ivparm = 0; ivparm < EB_NDER; ivparm++)
    fprintf(mfp, " %s", eb_dernames[ivparm]);

  fprintf(mfp,
          " chisq nlap\n"
          "#");

  for(ivparm = 0; ivparm < par->nparm; ivparm++)
    fprintf(mfp, " %.*e", DBL_DIG, v[ivparm]);

  for(ivparm = 0; ivparm < EB_NDER; ivparm++)
    fprintf(mfp, " %.*e", DBL_DIG, vder[ivparm]);
  
  fprintf(mfp, " %.*e %.*e\n", DBL_DIG, chisq, DBL_DIG, bestnlap);

  fprintf(mfp,
	  "# Vary:\n"
	  "#");

  for(ivparm = 0; ivparm < par->nparm; ivparm++)
    fprintf(mfp, " %d", par->vary[ivparm]);

  fprintf(mfp, "\n");

  fprintf(mfp,
	  "# Simulations: %d covinit %d burn %d total\n",
	  ninit, nburn, nsim);

  fprintf(mfp,
          "# hjdoff: %.*e\n", DBL_DIG, par->hjdoff);

  fprintf(mfp,
	  "# Columns:\n"
	  "#");

  fprintf(mfp, " isim");

  for(ivparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2)
      fprintf(mfp, " %s", par->vnames[ivparm]);

  fprintf(mfp, " chisq nlap\n");

  /* RNG */
  rng_init(&rngs, iseed);

  /* MC */
  nsimd = 0;
  nacc = 0;

  for(isim = 0; isim < nsim; isim++) {
    if((isim+1) % 100 == 0)
      fprintf(stderr, "\rSimulation %d of %d", isim+1, nsim);

#if 0 /* periodic rather than continuous updates */
    if(isim > 0 && (isim % ninit) == 0 && isim <= nburn) {
#else
    if(isim == ninit) {
#endif
      /* Means */
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	ameanprev[iaparm] = (ameantmp[iaparm] - a[iaparm]) / (ninit-1);

      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	amean[iaparm] = ameantmp[iaparm] / ninit;

      /* Covariance */
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	for(japarm = 0; japarm <= iaparm; japarm++) {
	  sum = 0;

	  for(itmp = 0; itmp < ninit; itmp++)
	    sum += (mc_init[iaparm*ninit+itmp]-amean[iaparm])
                 * (mc_init[japarm*ninit+itmp]-amean[japarm]);

	  acov[iaparm*par->nvarym+japarm] = sd*(sum / (ninit-1) + EPS);
	}

#if 0  /* periodic rather than continuous updates */
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	ameantmp[iaparm] = 0;
#endif
    }
#if 1  /* continuous rather than periodic updates */
    else if(isim > ninit) {
      /* Save this for later */
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	ameanprev[iaparm] = amean[iaparm];

      /* Update mean - easiest the old-fashioned way */
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	amean[iaparm] = ameantmp[iaparm] / isim;

      /* Update covariance using recursion relation */
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	for(japarm = 0; japarm <= iaparm; japarm++)
	  acov[iaparm*par->nvarym+japarm]
            = (acov[iaparm*par->nvarym+japarm] * (isim-2) +
               sd * ((isim-1) * ameanprev[iaparm]*ameanprev[japarm] -
                     isim * amean[iaparm]*amean[japarm] +
                     a[iaparm]*a[japarm] +
                     EPS)) / (isim-1);
    }
#endif

#if 1
    /* Compute proposal - multivariate Gaussian deviate */
    if(rng_fetch_mvgauss(&rngs, a, acov, awork, atrial, par->nvarym)) {
      report_err(errstr, "mvgauss failure");
      goto error;
    }
#else
    /* Generate trial parameters */
    for(iaparm = 0; iaparm < par->nvarym; iaparm++)
      atrial[iaparm] = a[iaparm];

    /* Choose one and perturb it */
    iaparmadj = isim % par->nvarym;  /* one at a time round robin */
    atrial[iaparmadj] += sqrt(acov[iaparmadj*par->nvarym+iaparmadj])
                      * rng_fetch_one_gauss(&rngs);
#endif

    /* Non-offset vector for chi^2 */
    for(iaparm = 0; iaparm < par->nvarym; iaparm++)
      aatrial[iaparm] = ainit[iaparm]+atrial[iaparm];

    /* Full vector for bounds checks */
    for(iaparm = 0, ivparm = 0; ivparm < par->nparm; ivparm++) {
      if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2)
	vtrial[ivparm] = aatrial[iaparm++];
      else
	vtrial[ivparm] = vinit[ivparm];
    }

    /* Check for LD2 = LD1 case; copy if necessary */
    if(par->ldtype[1] == LD_SAME) {
      vtrial[EB_PAR_LDLIN2] = vtrial[EB_PAR_LDLIN1];
      vtrial[EB_PAR_LDNON2] = vtrial[EB_PAR_LDNON1];
    }

    /* Light travel time parameter */
    if(par->vary[EB_PAR_CLTT] < 0)
      vtrial[EB_PAR_CLTT] = vtrial[PAR_KTOT]*1000 / EB_LIGHT;

    /* Check RV scale factors are non-negative */
    rvok = 1;

    for(idat = 0; idat < ndata; idat++)
      if(dlist[idat].fiterr &&
	 dlist[idat].obstype == OBS_RV &&
	 vtrial[dlist[idat].perrsc] < 0)
	rvok = 0;

    /* Check boundaries on parameters */
    if(vtrial[EB_PAR_RASUM] >= 0 && vtrial[EB_PAR_RR] >= 0 &&
       vtrial[EB_PAR_COSI] >= 0 && vtrial[EB_PAR_L3] >= 0 &&
       rvok &&
       ldok(vtrial, par->ldtype, 0) &&
       ldok(vtrial, par->ldtype, 1) &&
       vtrial[EB_PAR_FSPOT1] >= 0 && vtrial[EB_PAR_FSPOT1] <= 1 &&
       vtrial[EB_PAR_FSPOT2] >= 0 && vtrial[EB_PAR_FSPOT2] <= 1 &&
       ooeok(vtrial, 0) &&
       ooeok(vtrial, 1)) {

      /* Evaluate chi squared */
      newchisq = 0.0;

      for(idat = 0; idat < ndata; idat++) {
        fit_func(par, idat, aatrial,
                 dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
                 0, 1, 1);

	for(meas = 0; meas < dlist[idat].nmeas; meas++) {
	  dlist[idat].resid[meas] = dlist[idat].y[meas] - dlist[idat].m[meas];
          if(dlist[idat].iflag[meas])
            newchisq += dlist[idat].resid[meas]*dlist[idat].resid[meas] / varscale(dlist+idat, vtrial, meas);
        }
      }
      
      newesq = vtrial[EB_PAR_ECOSW]*vtrial[EB_PAR_ECOSW] + 
	       vtrial[EB_PAR_ESINW]*vtrial[EB_PAR_ESINW];

      /* Negative log likelihood ("unlikelihood") */
      nllnew = newchisq / 2;

      for(idat = 0; idat < ndata; idat++)
	if(dlist[idat].fiterr) {
	  if(dlist[idat].obstype == OBS_LC) {
	    /* Factor of 1/S^N in likelihood */
	    nllnew += dlist[idat].nchisq*log(vtrial[dlist[idat].perrsc]);
	  }
	  else if(dlist[idat].obstype == OBS_RV) {
	    /* Factor of product of 1/sqrt(sigma_i^2+s^2) over all i */
	    for(meas = 0; meas < dlist[idat].nmeas; meas++)
	      if(dlist[idat].iflag[meas])
		nllnew += 0.5*log(dlist[idat].yerr[meas]*dlist[idat].yerr[meas]+vtrial[dlist[idat].perrsc]*vtrial[dlist[idat].perrsc] / dlist[idat].yinfwt[meas]);
	  }
	}

      /* Now compute log of denominator of Bayes factor */
      bayenew = 0.0;

      /* Eccentricity - convert to uniform prior in e.
       * the "if" statement prevents the obvious KABOOM if circular.
       */
      if(newesq != 0)
	bayenew += 0.5*log(newesq);

      /* Modified Jeffreys prior in semiamplitude, transition at 10% of empirical scatter */
      bayenew += log(vtrial[PAR_KTOT]+0.1*rvavsig);

      /* Jeffreys priors in error scaling parameters
       * Modified for RV with transition at 10% of empirical scatter.
       */
      for(idat = 0; idat < ndata; idat++)
	if(dlist[idat].fiterr) {
	  if(dlist[idat].obstype == OBS_LC)
	    bayenew += log(vtrial[dlist[idat].perrsc]);
	  else if(dlist[idat].obstype == OBS_RV) {
	    /* Sum of squared reciprocals */
	    swt = 0;

	    haveerr = 0;

	    for(meas = 0; meas < dlist[idat].nmeas; meas++)
	      if(dlist[idat].iflag[meas]) {
		if(dlist[idat].yerr[meas] > 0)
		  haveerr = 1;

		swt += 1.0/(dlist[idat].yerr[meas]*dlist[idat].yerr[meas]+vtrial[dlist[idat].perrsc]*vtrial[dlist[idat].perrsc] / dlist[idat].yinfwt[meas]);
	      }

	    if(haveerr)  /* modified Jeffreys prior */
	      bayenew -= log((vtrial[dlist[idat].perrsc] +
                              0.1*dlist[idat].sigresid) * swt);
	    else  /* Jeffreys prior */
	      bayenew -= log(vtrial[dlist[idat].perrsc] * swt);
	  }
	}

      bf = exp(nll+baye - (nllnew+bayenew));

      /* Metropolis-Hastings acceptance criterion */
      if(bf > 1 || rng_fetch_one_uniform(&rngs) < bf) {
	/* Accept new point */
	for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	  a[iaparm] = atrial[iaparm];

	for(ivparm = 0; ivparm < par->nparm; ivparm++)
	  v[ivparm] = vtrial[ivparm];

	chisq = newchisq;  /* ?! */
	nll = nllnew;
	baye = bayenew;
	esq = newesq;

	nacc++;
      }
      /* else repeat old point */
    }
    /* else repeat old point */

    /* Accumulate for empirical mean and covariance */
    for(iaparm = 0; iaparm < par->nvarym; iaparm++) {
      ameantmp[iaparm] += a[iaparm];
      
      if(isim < nburn)
	mc_init[iaparm*ninit+(isim % ninit)] = a[iaparm];
    }

    /* Derived parameters */
    eb_getvder(v, v[PAR_GAMMA], v[PAR_KTOT], vder);
    
    nlap = nll+baye;

    /* Store */
    if(isim >= nskip) {
      for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	mc_res[iaparm*nsim+nsimd] = a[iaparm];

      for(ivparm = 0; ivparm < EB_NDER; ivparm++)
	mc_der[ivparm*nsim+nsimd] = vder[ivparm];

      /* Light ratios */
      for(iband = 0; iband < par->nband; iband++) {
        memcpy(vltmp, v, par->nparm*sizeof(double));
        
        if(iband > 0) {
          vltmp[EB_PAR_J] = v[par->pj[iband-1]];
          vltmp[EB_PAR_LDLIN1] = v[par->pldlin1[iband-1]];
          vltmp[EB_PAR_LDLIN2] = v[par->pldlin2[iband-1]];
          vltmp[EB_PAR_LDNON1] = v[par->pldnon1[iband-1]];
          vltmp[EB_PAR_LDNON2] = v[par->pldnon2[iband-1]];
        }
        
        phitmp = 0;
        typtmp = EB_OBS_AVLR;
        eb_model_dbl(vltmp, &phitmp, NULL, NULL, &typtmp,
                     mc_lrat + iband*nsim+nsimd,
                     NULL, 0, 1);
      }

      /* Accumulate maximum a-posteriori */
      if(nlap < bestnlap) {
	bestnlap = nlap;

	for(iaparm = 0; iaparm < par->nvarym; iaparm++)
	  abest[iaparm] = a[iaparm];

	for(ivparm = 0; ivparm < EB_NDER; ivparm++)
	  vderbest[ivparm] = vder[ivparm];
      }

      nsimd++;
    }

    fprintf(mfp, "%d", isim+1);

    for(iaparm = 0; iaparm < par->nvarym; iaparm++)
      fprintf(mfp, " %.*e", DBL_DIG, a[iaparm]);

    fprintf(mfp, " %.*e %.*e\n", DBL_DIG, chisq, DBL_DIG, nlap);
  }

  fprintf(stderr, "\n");

  fclose(mfp);

  tprintf(ofp,
	  "Acceptance rate: %d of %d = %.1f%%\n",
	  nacc, nsim, (100.0*nacc)/nsim);

  /* ACF */
  for(ivparm = 0, iaparm = 0; ivparm < par->nparm; ivparm++)
    if(par->vary[ivparm] == 1 || par->vary[ivparm] == 2) {
      acf(mc_res+iaparm*nsim, nsimd, par->vnames[ivparm], 1.0, ofp);
      iaparm++;
    }

  /* Call output subroutine */
  if(mc_out(par, ainit,
            mc_res, mc_der, mc_lrat,
            nsim, nsimd,
            abest, vderbest,
            ofp, tfp, errstr))
    goto error;

  free((void *) ainit);
  ainit = (double *) NULL;
  free((void *) acov);
  acov = (double *) NULL;
  free((void *) v);
  v = (double *) NULL;

  free((void *) mc_res);
  mc_res = (double *) NULL;
  free((void *) mc_lrat);
  mc_lrat = (double *) NULL;

  return(0);

 error:
  if(ainit)
    free((void *) ainit);
  if(acov)
    free((void *) acov);
  if(v)
    free((void *) v);
  if(mc_res)
    free((void *) mc_res);
  if(mc_lrat)
    free((void *) mc_lrat);

  return(-1);
}

static inline int ldok (double *v, int ldtype[2], int sec) {
  int ildlin, ildnon, ildtype, ldok;

  /* Figure out which coefficients */
  if(sec && ldtype[sec] != LD_SAME) {
    ildlin = EB_PAR_LDLIN2;
    ildnon = EB_PAR_LDNON2;
    ildtype = ldtype[sec];
  }
  else {
    ildlin = EB_PAR_LDLIN1;
    ildnon = EB_PAR_LDNON1;
    ildtype = ldtype[0];
  }

  /* Require specific intensity to be non-negative and monotonically
   * decreasing from centre to limb.  Note that these assumptions
   * are NOT guaranteed to be valid - for example, negative specific
   * intensities are produced in U-band with a linear LD law.
   * In reality, this should be treated by clipping them at zero
   * while integrating to generate the light curve, but the model
   * doesn't do this and will instead dutifully subtract flux - so
   * the non-negative constraint is justifiable in this case.  The
   * monotonic decrease may not be so reasonable.  According to
   * Claret (2000), some of the lowest temperature PHOENIX model
   * atmospheres (in the brown dwarf domain) showed limb brightening,
   * so these assumptions may need to be revised if fitting objects
   * below the hydrogen burning limit.
   */
  ldok = 0;

  switch(ildtype) {
  case LD_LIN:
    if(v[ildlin] >= 0 && v[ildlin] <= 1)
      ldok = 1;
    break;
  case LD_QUAD:
    if((v[ildlin]+v[ildnon]) >= 0 &&
       (v[ildlin]+v[ildnon]) <= 1 &&
       v[ildlin] >= 0 &&
       (v[ildlin]+2*v[ildnon]) >= 0)  /* this one often neglected?! */
      ldok = 1;
  case LD_SQRT:
    if((v[ildlin]+v[ildnon]) >= 0 &&
       (v[ildlin]+v[ildnon]) <= 1 &&
       (2*v[ildlin]+v[ildnon]) >= 0 &&
       v[ildnon] >= 0)  /* only way to prevent -ve divergence at mu=0 */
      ldok = 1;
    break;
    /* LOG not normalized in light so no point - do not use! */
  }

  return ldok;
}

static inline int ooeok (double *v, int sec) {
  int io, i1a, i1b, i2a, i2b;

  if(sec) {
    io = EB_PAR_OOE2O;
    i1a = EB_PAR_OOE21A;
    i1b = EB_PAR_OOE21B;
    i2a = EB_PAR_OOE22A;
    i2b = EB_PAR_OOE22B;
  }
  else {
    io = EB_PAR_OOE1O;
    i1a = EB_PAR_OOE11A;
    i1b = EB_PAR_OOE11B;
    i2a = EB_PAR_OOE12A;
    i2b = EB_PAR_OOE12B;
  }

  /* Make sure out of eclipse light never goes negative */
  return(v[io] +
         2*sqrt(v[i1a]*v[i1a]+
                v[i1b]*v[i1b]+
                v[i2a]*v[i2a]+
                v[i2b]*v[i2b]) < 1);
}
