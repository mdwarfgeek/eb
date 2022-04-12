#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "ebmc.h"

static char **read_list (char *str, int *n_r, char *errstr);
static int read_lc (char *filename, struct fit_data *d, int iband, int oldlc,
                    char *errstr);
static int read_rv (char *filename, struct fit_data *d,
		    int component, int obstype, int iband, char *errstr);
static int read_prior (char *filename, struct fit_data *d, char *errstr);

static int read_pset (char *filename, int ldtype[2],
		      double *vfix, int *varyfix, double *vsgfix,
		      struct fit_pset_entry **pextra_r, int *npextra_r,
		      char *errstr);
static int write_pset (char *filename, 
		       double *vfix, int *varyfix, double *vsgfix,
		       char *errstr);

/* Getopt stuff */
extern char *optarg;
extern int optind;

static void usage (char *av) {
  fprintf(stderr,
	  "Usage:\t%s [-f names] [-k nsig] [-m mcfile] [-l mcfile] [-n trials] [-o file] [-O] [-p pgdev] pset lcfile[,texp,fitairm,fitcm] [lc?=] [rv?=] [lrat?=] [prior=] [...]\n\n"
	  "General options:\n"
	  "\t-f names   Set filter names.  Comma separated list.\n"
	  "\t-k nsig    Outlier clip data at nsig*1.48*MAD.  (default: 10)\n"
          "\t-o file    Write final parameters to file.\n"
          "\t-O         Use old light curve format.  (default: new format)\n\n"
	  "Monte Carlo options:\n"
	  "\t-l mcfile  Load samples from 'mcfile' and do output.\n"
	  "\t-m mcfile  Run Monte Carlo, direct samples to 'mcfile'.\n"
	  "\t           Also directs screen output to 'mcfile'.out, and\n"
          "\t           plots to 'mcfile'.ps.\n"
	  "\t-n ntrials Use 'ntrials' Monte Carlo trials.\n"
          "\t           (default: 10,000 for testing - you will want a\n"
          "\t            lot more for a final solution!)\n"
          "\t-s seed    Set RNG seed\n\n"
          "Plotting options:\n"
	  "\t-p pgdev   Plot to PGPLOT device 'pgdev'.  (default: asks)\n"
          "\t-R         Omit raw data panels from light curve plots\n"
          "\t-S         Don't separate light curves vertically on plots\n\n",
	  av);
  exit(1);
}  

int main (int argc, char *argv[]) {
  char *pn = (char *) NULL, *avzero;
  int c;

  struct fit_data *dlist = (struct fit_data *) NULL;
  char *p, *cp, *ep;
  int iband, iobj;

  char errstr[ERRSTR_LEN];
  int f;

  int dofit = 1;

  char mcfile[1024], outfile[1024], texfile[1024];
  char **mcfilelist = (char **) NULL;
  int domc = 0;
  int nmc = 10000;
  int nmcfile = 0;

  char outparfile[1024];
  int dooutpar = 0;

  char pgdev[1024] = "?";
  int pgdev_set = 0;
  int noraw = 0;
  int novertsep = 0;

  char **filtnamelist = (char **) NULL, filtnamebuf[1024];
  int nfiltname = 0;

  FILE *ofp = (FILE *) NULL;
  FILE *tfp = (FILE *) NULL;

  int nsigma = 10;
  int iseed = 42;

  int ldtype[2] = { LD_QUAD, LD_QUAD };
  double vfix[NPARFIX], vsgfix[NPARFIX];
  int varyfix[NPARFIX];

  struct fit_pset_entry *pextra = (struct fit_pset_entry *) NULL;
  int npextra = 0;

  int oldlc = 0;

  int iparm;

  struct fit_parms par;

  int idat, meas;

  double chisq;
  int nchisq;

  double phi, lrat;
  unsigned char typ;

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "ebmc";

  setprogname(pn);

  avzero = argv[0];

  /* Init constants */
  init_const();

  /* Extract command-line arguments */
  while((c = getopt(argc, argv, "f:k:l:m:n:o:Op:s:RS")) != -1)
    switch(c) {
    case 'f':
      strncpy(filtnamebuf, optarg, sizeof(filtnamebuf)-1);
      filtnamebuf[sizeof(filtnamebuf)-1] = '\0';

      filtnamelist = read_list(filtnamebuf, &nfiltname, errstr);
      if(!filtnamelist)
	fatal(1, "read_list: %s", errstr);

      break;
    case 'k':
      nsigma = (int) strtol(optarg, &ep, 0);
      if(*ep != '\0')
	fatal(1, "invalid n-sigma: %s", optarg);

      break;
    case 'l':
      if(domc)
	fatal(1, "can't use -l and -m together");

      strncpy(mcfile, optarg, sizeof(mcfile)-1);
      mcfile[sizeof(mcfile)-1] = '\0';

      mcfilelist = read_list(mcfile, &nmcfile, errstr);
      if(!mcfilelist)
	fatal(1, "read_list: %s", errstr);

      domc = -1;
      dofit = 0;

      break;
    case 'm':
      if(domc)
	fatal(1, "can't use -l and -m together");

      strncpy(mcfile, optarg, sizeof(mcfile)-1);
      mcfile[sizeof(mcfile)-1] = '\0';
      if(!pgdev_set)
	snprintf(pgdev, sizeof(pgdev), "%s.ps/cps", optarg);
      snprintf(outfile, sizeof(outfile), "%s.out", optarg);
      snprintf(texfile, sizeof(texfile), "%s.tex", optarg);
      domc = 1;

      break;
    case 'n':
      nmc = (int) strtol(optarg, &ep, 0);
      if(*ep != '\0')
	fatal(1, "invalid ntrials: %s", optarg);

      break;
    case 'o':
      strncpy(outparfile, optarg, sizeof(outparfile)-1);
      outparfile[sizeof(outparfile)-1] = '\0';
      dooutpar = 1;

      break;
    case 'O':
      oldlc = 1;

      break;
    case 'p':
      strncpy(pgdev, optarg, sizeof(pgdev)-1);
      pgdev[sizeof(pgdev)-1] = '\0';
      pgdev_set = 1;

      break;
    case 's':
      iseed = (int) strtol(optarg, &ep, 0);
      if(*ep != '\0' || iseed < 0)
        fatal(1, "invalid seed: %s", optarg);

      break;
    case 'R':
      noraw = 1;
      break;
    case 'S':
      novertsep = 1;
      break;
    case '?':
    default:
      usage(avzero);
    }

  argc -= optind;
  argv += optind;

  /* Get remaining arguments */
  if(argc < 2)
    usage(avzero);

  /* Init pset */
  for(iparm = 0; iparm < NPARFIX; iparm++) {
    vfix[iparm] = 0;
    varyfix[iparm] = 0;
    vsgfix[iparm] = default_sigguess[iparm];
  }

  vfix[EB_PAR_GD1] = 0.32;   /* conv. envelopes from Lucy 1967 */
  vfix[EB_PAR_GD2] = 0.32;
  vfix[EB_PAR_REFL1] = 0.4;  /* standard value for albedo */
  vfix[EB_PAR_REFL2] = 0.4;
  vfix[EB_PAR_BEAM1] = 1.0;  /* standard value for beaming */
  vfix[EB_PAR_BEAM2] = 1.0;
  vfix[EB_PAR_TIDANG] = 0.0;
  vfix[EB_PAR_INTEG] = 1.0;  /* for compatibility only, not used */

  varyfix[EB_PAR_KTOTC] = -1;  /* default is to compute */

  if(read_pset(argv[0],
               ldtype, vfix, varyfix, vsgfix,
               &pextra, &npextra, errstr))
    fatal(1, "read_pset: %s", errstr);

  argv++;
  argc--;

  /* Allocate curves array */
  dlist = (struct fit_data *) malloc(argc * sizeof(struct fit_data));
  if(!dlist)
    error(1, "malloc");

  /* Ensure initialized to zero */
  memset(dlist, 0, argc * sizeof(struct fit_data));

  for(f = 0; f < argc; f++) {
    p = strchr(argv[f], '=');
    if(p)
      *p = '\0';

    if(!p || !strncasecmp(argv[f], "lc", 2)) {
      if(p) {
	cp = argv[f]+2;
	iband = (int) strtol(cp, &ep, 0);
	if(*ep != '\0')
	  fatal(1, "could not understand arg: %s", argv[f]);
	else if(cp == ep)
	  iband = 0;
      }
      else
	iband = 0;

      if(read_lc(p ? p+1 : argv[f], dlist+f, iband, oldlc, errstr))
	fatal(1, "read_lc: %s", errstr);
    }
    else if(!strncasecmp(argv[f], "rv", 2)) {
      cp = argv[f]+2;
      iobj = (int) strtol(cp, &ep, 0);
      if(*ep != '\0' || cp == ep)
	fatal(1, "could not understand arg: %s", argv[f]);

      if(read_rv(p+1, dlist+f, iobj, OBS_RV, 0, errstr))
	fatal(1, "read_rv: %s", errstr);
    }
    else if(!strncasecmp(argv[f], "lrat", 4)) {
      cp = argv[f]+4;
      iband = (int) strtol(cp, &ep, 0);
      if(*ep != '\0')
	fatal(1, "could not understand arg: %s", argv[f]);
      else if(cp == ep)
	iband = 0;

      if(read_rv(p+1, dlist+f, 0, OBS_LRAT, iband, errstr))
	fatal(1, "read_rv: %s", errstr);
    }
    else if(!strcasecmp(argv[f], "prior")) {
      if(read_prior(p+1, dlist+f, errstr))
	fatal(1, "read_prior: %s", errstr);
    }
    else
      fatal(1, "unknown data type: %*s", p-argv[f], argv[f]);
  }

  init_plots(pgdev);

  /* Open output file if there is one */
  if(domc == 1) {
    ofp = fopen(outfile, "w");
    if(!ofp)
      error(1, "open: %s", outfile);
  }

  /* Set up parameter vectors */
  if(make_parms(&par, ofp,
                dlist, argc,
                ldtype,
                vfix, varyfix, vsgfix,
                pextra, npextra,
                errstr))
    fatal(1, "make_parms: %s", errstr);

  if(dofit) {
    /* Fit */
    if(do_fit(&par, nsigma, ofp, errstr))
      fatal(1, "do_fit: %s", errstr);
  }

  if(domc == 1) {
    tfp = fopen(texfile, "w");
    if(!tfp)
      error(1, "open: %s", texfile);

    /* Now MC */
    if(do_mc(&par, ofp, tfp,
             mcfile, nmc, iseed, errstr))
      fatal(1, "do_mc: %s", errstr);

    fclose(tfp);
  }
  else if(domc == -1) {
    /* Read MC */
    if(read_mc(&par, ofp, (FILE *) NULL,
               mcfilelist, nmcfile, errstr))
      fatal(1, "read_mc: %s", errstr);
  }

  /* Recompute residual and chi squared */
  chisq = 0.0;
  nchisq = 0;

  for(idat = 0; idat < par.ndata; idat++) {
    fit_func(&par, idat, NULL,
             dlist[idat].hjd, dlist[idat].m, NULL, dlist[idat].nmeas,
             0, 1, 1);

    for(meas = 0; meas < dlist[idat].nmeas; meas++) {
      dlist[idat].resid[meas] = dlist[idat].y[meas] - dlist[idat].m[meas];
      if(dlist[idat].iflag[meas]) {
        chisq += dlist[idat].resid[meas]*dlist[idat].resid[meas] / varscale(dlist+idat, par.vinit, meas);
        nchisq++;
      }
    }
  }

  tprintf(ofp, "Final chisq = %.3f ndof = %d reduced chisq = %.3f\n",
	  chisq, nchisq-par.nvarym, chisq/(nchisq-par.nvarym));

  /* Work out light ratio for various bands */
  for(iband = 0; iband < par.nband; iband++) {
    memcpy(par.vtmp, par.vinit, par.nparm*sizeof(double));

    if(iband > 0) {
      par.vtmp[EB_PAR_J] = par.vinit[par.pj[iband-1]];
      par.vtmp[EB_PAR_L3] = par.vinit[par.pl3[iband-1]];
      par.vtmp[EB_PAR_LDLIN1] = par.vinit[par.pldlin1[iband-1]];
      par.vtmp[EB_PAR_LDLIN2] = par.vinit[par.pldlin2[iband-1]];
      par.vtmp[EB_PAR_LDNON1] = par.vinit[par.pldnon1[iband-1]];
      par.vtmp[EB_PAR_LDNON2] = par.vinit[par.pldnon2[iband-1]];
    }

    phi = 0;
    typ = EB_OBS_AVLR;
    eb_model_dbl(par.vtmp, &phi, NULL, NULL, &typ, &lrat, NULL, 0, 1);

    if(nfiltname > iband)
      printf("L_2/L_1 band %d (%s): %.4f\n",
             iband, filtnamelist[iband], lrat);
    else
      printf("L_2/L_1 band %d: %.4f\n", iband, lrat);
  }

  /* Plot */
  if(do_plots(&par, ofp,
              filtnamelist, nfiltname,
              noraw, novertsep,
              errstr))
    fatal(1, "do_plots: %s", errstr);

  if(ofp)
    fclose(ofp);

  close_plots();

  if(dooutpar) {
    /* Put this back */
    par.vinit[EB_PAR_T0] += par.hjdoff;

    for(iparm = 0; iparm < NPARFIX; iparm++) {
      vfix[iparm] = par.vinit[iparm];
      varyfix[iparm] = par.vary[iparm];
    }

    if(write_pset(outparfile, vfix, varyfix, vsgfix, errstr))
      fatal(1, "write_pset: %s", errstr);
  }

  return(0);
}

static char **read_list (char *str, int *n_r, char *errstr) {
  char **list = (char **) NULL, *p, *cp;
  int n, iend;

  p = str;

  n = 0;
  iend = 0;
  do {
    cp = strchr(p, ',');
    if(cp)
      *cp = '\0';
    else
      iend = 1;

    list = (char **) realloc(list, (n+1)*sizeof(char *));
    if(!list) {
      report_syserr(errstr, "realloc");
      goto error;
    }

    list[n] = p;
    n++;
    
    /* Update */
    if(!iend)
      p = cp+1;
    
  } while(!iend);

  *n_r = n;

  return(list);

 error:
  if(list)
    free((void *) list);

  return((char **) NULL);
}

static int read_lc (char *filename, struct fit_data *d, int iband, int oldlc,
                    char *errstr) {
  FILE *fp;
  char line[1024], *p, *cp;
  int rv, iend;

  double *hjd = (double *) NULL;
  float *mag = (float *) NULL, *magerr = (float *) NULL;
  float *airmass = (float *) NULL, *cm = (float *) NULL;
  int nmeas = 0;

  int *iseg = (int *) NULL, thisiseg;
  int *haveseg = (int *) NULL, nhaveseg = 0, iclr, inew;

  float thisextinc, thisfwhm, thisell, thisx, thisy, thistheta;
  float thissky, thispeak;
  int thisv, thisr, thismerid, thisflag;

  int dotexp = 0;
  float *texp = (float *) NULL;

  /* Defaults for options */
  d->noacf = 0;
  d->fiterr = 0;
  d->fitairm = 0;
  d->fitcm = 0;

  /* Pull out optional bits */
  p = strchr(filename, ',');
  if(p) {
    *p = '\0';
    p++;

    iend = 0;
    do {
      cp = strchr(p, ',');
      if(cp)
	*cp = '\0';
      else
	iend = 1;

      if(!strncasecmp(p, "texp", 4))
	dotexp = 1;
      else if(!strncasecmp(p, "noacf", 5))
	d->noacf = 1;
      else if(!strncasecmp(p, "fitair", 6) || !strncasecmp(p, "air", 3))
	d->fitairm = 1;
      else if(!strncasecmp(p, "fitc", 4) || !strncasecmp(p, "c", 1))
	d->fitcm = 1;
      else if(!strncasecmp(p, "fiterr", 6) || !strncasecmp(p, "err", 3))
	d->fiterr = 1;
      else if(!strncasecmp(p, "useair", 6))
	d->fitairm = -1;
      else if(!strncasecmp(p, "usec", 4))
	d->fitcm = -1;
      else {
	report_err(errstr, "unknown option: %s", p);
	goto error;
      }

      /* Update */
      if(!iend)
	p = cp+1;

    } while(!iend);
  }

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  while(fgets(line, sizeof(line), fp)) {
    p = strchr(line, '#');
    if(p)
      *p = '\0';

    p = sstrip(line);

    if(*p != '\0') {
      /* Allocate array element */
      hjd = (double *) realloc(hjd, (nmeas+1) * sizeof(double));
      mag = (float *) realloc(mag, (nmeas+1) * sizeof(float));
      magerr = (float *) realloc(magerr, (nmeas+1) * sizeof(float));
      airmass = (float *) realloc(airmass, (nmeas+1) * sizeof(float));
      cm = (float *) realloc(cm, (nmeas+1) * sizeof(float));
      iseg = (int *) realloc(iseg, (nmeas+1) * sizeof(int));
      texp = (float *) realloc(texp, (nmeas+1) * sizeof(float));

      /* Supply defaults for these */
      iseg[nmeas] = -1;
      texp[nmeas] = 0.0;

      if(oldlc)
        rv = sscanf(p,
                    "%lf %f %f %f %f %f %d %f %f %d %f",
                    hjd+nmeas, mag+nmeas, magerr+nmeas,
                    &thisextinc, &thisfwhm, airmass+nmeas,
                    &thismerid, &thisx, &thisy, &thisflag, cm+nmeas);
      else
        rv = sscanf(p,
                    "%lf %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f",
                    hjd+nmeas, mag+nmeas, magerr+nmeas,
                    texp+nmeas, &thisextinc, &thisfwhm, &thisell,
                    airmass+nmeas, &thisx, &thisy, &thistheta, &thissky,
                    &thispeak, iseg+nmeas, &thisv, &thisr, &thisflag,
                    cm+nmeas);

      if(rv < 3) {
	report_err(errstr, "could not understand: %s", p);
	goto error;
      }

      if(hjd[nmeas] < 2400000.5)
	hjd[nmeas] += 2400000.5;

      /* Make a table to keep track of which segment numbers have data */
      if(iseg[nmeas] > 0) {
        thisiseg = iseg[nmeas];

        if(thisiseg+1 > nhaveseg) {
          /* Allocate extra */
          haveseg = (int *) realloc(haveseg, (thisiseg+1) * sizeof(int));
          if(!haveseg)
            goto error;

          /* Clear new elements */
          for(iclr = nhaveseg; iclr <= thisiseg; iclr++)
            haveseg[iclr] = 0;

          nhaveseg = thisiseg+1;
        }

        haveseg[thisiseg] = 1;
      }

      if(magerr[nmeas] > 0)
	nmeas++;
    }
  }

  if(ferror(fp)) {
    report_syserr(errstr, "%s: read", filename);
    goto error;
  }

  fclose(fp);

  d->hjd = hjd;
  d->y = mag;
  d->yerr = magerr;
  d->nmeas = nmeas;
  d->component = 0;
  d->obstype = OBS_LC;
  d->airmass = airmass;
  d->cm = cm;
  d->iband = iband;

  if(dotexp)
    d->texp = texp;
  else
    d->texp = (float *) NULL;

  if(haveseg) {
    /* Renumber segments to ensure gapless when assigning
       new parameters for individual zero points. */
    inew = 0;

    for(thisiseg = 0; thisiseg < nhaveseg; thisiseg++) {
      if(haveseg[thisiseg]) {
        haveseg[thisiseg] = inew;
        inew++;
      }
      else
        haveseg[thisiseg] = -1;
    }

    d->iseg = iseg;
    d->segtbl = haveseg;
    d->nseg = inew;
  }
  else {
    d->iseg = (int *) NULL;
    d->segtbl = (int *) NULL;
    d->nseg = 0;
  }

  return(0);

 error:
  if(hjd)
    free((void *) hjd);
  if(mag)
    free((void *) mag);
  if(magerr)
    free((void *) magerr);
  if(airmass)
    free((void *) airmass);
  if(iseg)
    free((void *) iseg);
  if(haveseg)
    free((void *) haveseg);
  if(cm)
    free((void *) cm);
  if(texp)
    free((void *) texp);

  return(-1);
}

static int read_rv (char *filename, struct fit_data *d,
		    int component, int obstype, int iband, char *errstr) {
  FILE *fp;
  char line[1024], *p, *cp, *np, *ep;
  int rv, iend;

  double *hjd = (double *) NULL;
  float *y = (float *) NULL, *yerr = (float *) NULL, *yinfwt = (float *) NULL;
  int nmeas = 0;

  int dotexp = 0;
  float *texp = (float *) NULL;

  /* Defaults for options */
  d->fiterr = 0;
  d->errguess = 0;
  d->havewt = 0;  /* error column should be treated as weights multiplying
                   * a computed error instead */

  /* Pull out optional bits */
  p = strchr(filename, ',');
  if(p) {
    *p = '\0';
    p++;

    iend = 0;
    do {
      cp = strchr(p, ',');
      if(cp)
	*cp = '\0';
      else
	iend = 1;

      if(!strncasecmp(p, "texp", 4))
	dotexp = 1;
      else if(!strncasecmp(p, "fiterr", 6) || !strncasecmp(p, "err", 3)) {
	d->fiterr = 1;

	np = strchr(p, '=');
	if(np) {
	  np++;  /* skip the '=' */

	  d->errguess = (float) strtod(np, &ep);
	  if(*ep != '\0') {
	    /* Not understood */
	    report_err(errstr, "problems understanding: %s", p);
	    goto error;
	  }
	}
      }
      else if(!strncasecmp(p, "wt", 2) || !strncasecmp(p, "weight", 6))
	d->havewt = 1;
      else {
	report_err(errstr, "unknown option: %s", p);
	goto error;
      }

      /* Update */
      if(!iend)
	p = cp+1;

    } while(!iend);
  }

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  while(fgets(line, sizeof(line), fp)) {
    p = strchr(line, '#');
    if(p)
      *p = '\0';

    p = sstrip(line);

    if(*p != '\0') {
      /* Allocate array element */
      hjd = (double *) realloc(hjd, (nmeas+1) * sizeof(double));
      y = (float *) realloc(y, (nmeas+1) * sizeof(float));
      yerr = (float *) realloc(yerr, (nmeas+1) * sizeof(float));
      yinfwt = (float *) realloc(yinfwt, (nmeas+1) * sizeof(float));
      texp = (float *) realloc(texp, (nmeas+1) * sizeof(float));

      /* Init errors to zero - they are optional and will be
       * computed if not given.
       */
      yerr[nmeas] = 0.0;
      yinfwt[nmeas] = 1.0;
      texp[nmeas] = 0.0;

      /* Read */
      rv = sscanf(p, "%lf %f %f %f",
                  hjd+nmeas,
                  y+nmeas,
                  (d->havewt ? yinfwt : yerr)+nmeas,
                  texp+nmeas);
      if(rv < 2) {
        report_err(errstr, "could not understand: %s", p);
        goto error;
      }

      if(hjd[nmeas] > 0 && hjd[nmeas] < 2400000.5)
	hjd[nmeas] += 2400000.5;

      nmeas++;
    }
  }

  if(ferror(fp)) {
    report_syserr(errstr, "%s: read", filename);
    goto error;
  }

  fclose(fp);

  d->hjd = hjd;
  d->y = y;
  d->yerr = yerr;
  d->yinfwt = yinfwt;
  d->nmeas = nmeas;
  d->component = component;
  d->obstype = obstype;
  d->iband = iband;

  if(dotexp)
    d->texp = texp;
  else
    d->texp = (float *) NULL;

  return(0);

 error:
  if(hjd)
    free((void *) hjd);
  if(y)
    free((void *) y);
  if(yerr)
    free((void *) yerr);
  if(texp)
    free((void *) texp);
  
  return(-1);
}

static int read_prior (char *filename, struct fit_data *d, char *errstr) {
  FILE *fp;
  char line[1024], name[1024], *p;
  int rv, i, ifound;

  int *ipar = (int *) NULL;
  double *hjd = (double *) NULL;
  float *y = (float *) NULL, *yerr = (float *) NULL;
  int nmeas = 0;

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  while(fgets(line, sizeof(line), fp)) {
    p = strchr(line, '#');
    if(p)
      *p = '\0';

    p = sstrip(line);

    if(*p != '\0') {
      /* Allocate array element */
      ipar = (int *) realloc(ipar, (nmeas+1) * sizeof(int));
      hjd = (double *) realloc(hjd, (nmeas+1) * sizeof(double));
      y = (float *) realloc(y, (nmeas+1) * sizeof(float));
      yerr = (float *) realloc(yerr, (nmeas+1) * sizeof(float));

      rv = sscanf(p, "%1023s %f %f",
		  name, y+nmeas, yerr+nmeas);
      if(rv != 3) {
	report_err(errstr, "could not understand: %s", p);
	goto error;
      }

      /* Lookup name */
      ifound = -1;
      for(i = 0; i < NPARFIX; i++)
	if(!strcmp(name, parnames[i])) {
	  ifound = i;
	  break;
	}
    
      if(ifound >= 0) {
	ipar[nmeas] = ifound;
	hjd[nmeas] = 0;  /* dummy, wasteful but I'm too lazy to fix the rest! */
      }
      else {
	report_err(errstr, "unknown parameter: %s", name);
	goto error;
      }

      nmeas++;
    }
  }

  if(ferror(fp)) {
    report_syserr(errstr, "%s: read", filename);
    goto error;
  }

  fclose(fp);

  d->ipar = ipar;
  d->hjd = hjd;
  d->y = y;
  d->yerr = yerr;
  d->nmeas = nmeas;
  d->obstype = OBS_PRIOR;
  d->texp = (float *) NULL;

  return(0);

 error:
  if(ipar)
    free((void *) ipar);
  if(hjd)
    free((void *) hjd);
  if(y)
    free((void *) y);
  if(yerr)
    free((void *) yerr);
  
  return(-1);
}

static int read_pset (char *filename, int ldtype[2],
		      double *vfix, int *varyfix, double *vsgfix,
		      struct fit_pset_entry **pextra_r, int *npextra_r,
		      char *errstr) {
  FILE *fp;
  char line[1024], name[1024], *p, *ep;
  int rv, i, ifound, istar;

  double yinit, ysigguess;
  int ivary;

  struct fit_pset_entry *pextra = (struct fit_pset_entry *) NULL;
  int npextra = 0;

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  while(fgets(line, sizeof(line), fp)) {
    p = strchr(line, '#');
    if(p)
      *p = '\0';

    p = sstrip(line);

    if(*p != '\0') {
      /* First check for LD type - special case (ugly, I know... sorry folks) */
      if(!strncmp(p, "ldtype", 6)) {
	p += 6;

	istar = strtol(p, &ep, 0);
	if(!isspace((unsigned char) *ep)) {
	  report_err(errstr, "could not understand: %s", p);
	  goto error;
	}

	if(istar < 1 || istar > 2) {
	  report_err(errstr, "invalid star number: %d", istar);
	  goto error;
	}

	p = ep;

	while(*p != '\0' && isspace((unsigned char) *p))
	  p++;

	if(!strcmp(p, "lin"))
	  ldtype[istar-1] = LD_LIN;
	else if(!strcmp(p, "quad"))
	  ldtype[istar-1] = LD_QUAD;
	else if(!strcmp(p, "same")) {
	  if(istar != 2) {
	    report_err(errstr, "LD_SAME only works on star 2");
	    goto error;
	  }

	  ldtype[istar-1] = LD_SAME;
	}
	else {
	  report_err(errstr, "unknown LD type: %s", p);
	  goto error;
	}
      }
      else {
	/* Normal arg */
	rv = sscanf(p, "%1023s %lf %d %lf",
		    name, &yinit, &ivary, &ysigguess);
	if(rv < 3) {
	  report_err(errstr, "could not understand: %s", p);
	  goto error;
	}
	
	/* Lookup name */
	ifound = -1;
	for(i = 0; i < NPARFIX; i++)
	  if(!strcmp(name, parnames[i])) {
	    ifound = i;
	    break;
	  }
	
	if(ifound >= 0) {
	  if(ifound == EB_PAR_COSI)
	    vfix[ifound] = cos(yinit * M_PI / 180.0);
	  else
	    vfix[ifound] = yinit;
	  
	  varyfix[ifound] = ivary;
	  
	  if(rv > 3)
	    vsgfix[ifound] = ysigguess;
	}
	else {
	  /* Assume it's an extra one - stash */
	  pextra =
            (struct fit_pset_entry *) realloc(pextra,
                                              (npextra+1)*
                                              sizeof(struct fit_pset_entry));
	  if(!pextra) {
	    report_syserr(errstr, "malloc");
	    goto error;
	  }

	  memcpy(pextra[npextra].name, name, sizeof(pextra[npextra].name));
	  pextra[npextra].v = yinit;
	  pextra[npextra].vary = ivary;

	  if(rv > 3) {
	    pextra[npextra].sig = ysigguess;
	    pextra[npextra].havesig = 1;
	  }
	  else
	    pextra[npextra].havesig = 0;

	  npextra++;
	}
      }
    }
  }

  if(ferror(fp)) {
    report_syserr(errstr, "%s: read", filename);
    goto error;
  }

  fclose(fp);

  *pextra_r = pextra;
  *npextra_r = npextra;

  return(0);

 error:
  return(-1);
}

static int write_pset (char *filename,
		       double *vfix, int *varyfix, double *vsgfix,
		       char *errstr) {
  FILE *fp;
  int i, len, lmax;
  double y;

  /* Open file */
  fp = fopen(filename, "w");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  lmax = 0;
  for(i = 0; i < NPARFIX; i++) {
    len = strlen(parnames[i]);
    if(len > lmax)
      lmax = len;
  }

  for(i = 0; i < NPARFIX; i++) {
    if(i == EB_PAR_COSI)
      y = acos(vfix[i]) * 180.0 / M_PI;
    else
      y = vfix[i];

    fprintf(fp, "%-*s %16.8lf %2d %11.8lf\n",
	    lmax, parnames[i], y, varyfix[i], vsgfix[i]);
  }

  if(ferror(fp)) {
    report_syserr(errstr, "%s: write", filename);
    goto error;
  }

  if(fclose(fp)) {
    report_syserr(errstr, "%s: close", filename);
    goto error;
  }

  return(0);

 error:
  return(-1);
}
