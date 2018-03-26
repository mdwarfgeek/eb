#ifndef EBMC_H
#define EBMC_H

#include <stdio.h>

#include <lfa.h>
#include <eb.h>
#include <util.h>

/* Add radial velocity parameters to all internal vectors */
#define PAR_KTOT  (EB_NPAR)
#define PAR_GAMMA (EB_NPAR+1)
#define NPARFIX   (EB_NPAR+2)

/* LD types */
#define LD_SAME      0
#define LD_LIN       1
#define LD_LOG       2  /* not supported */
#define LD_SQRT      3  /* not supported */
#define LD_QUAD      4

/* Observation types */
#define OBS_LC       0
#define OBS_RV       1
#define OBS_LRAT     2
#define OBS_TMIN     3
#define OBS_PRIOR    4  /* prior on any parameter in 'v' */

/* Max parameter name length */
#define PARNAME_MAX 64

/* Structure for each input dataset */
struct fit_data {
  double *hjd;
  float *y;
  float *yerr;
  float *yinfwt;   /* weights from input file */
  float *airmass;  /* airmass for LCs */
  float *cm;       /* common mode for LCs */
  int *ipar;       /* parameter no. for priors */
  int nmeas;

  int component;   /* component - 0: both, 1: pri, 2: sec */
  int obstype;     /* see eb.h */
  int iband;       /* passband number for LC / light ratio */

  float *texp;     /* integration time */
  float errguess;  /* initial guess at error inflation parameter */

  double *resid;
  double *m;
  double *corr;

  double ymed;
  double ysig;
  double medresid;
  double sigresid;

  double chisq;
  int nchisq;

  /* Temp arrays for plotting */
  float *phitmp;
  float *ytmp;

  unsigned char *iflag;
  unsigned char *iecl;

  int noacf;       /* LC only: no ACF for covariance */
  int fiterr;      /* LC/RV only: adjust errors? */
  int fitairm;     /* LC only: fit for airmass coeff? */
  int fitcm;       /* LC only: fit for common mode coeff? */
  int havewt;      /* RV only: weight the error inflation parameter? */

  int pnorm;       /* pointer to normalisation in parameter vector */
  int perrsc;      /* pointer to error adjustment parameter */
  int pairm;       /* pointer to airmass coefficient */
  int pcm;         /* pointer to common-mode coefficient */
};

/* Structure for passing around fit data and parameters */
struct fit_parms {
  struct fit_data *dlist;
  int ndata;
  int nmeastot;

  double *vinit;
  double *voff;
  double *vscl;
  double *vsg;
  double *vtmp;
  double *vcov;
  int nparm;
  
  int *vary;
  double nvaryflc;
  double nvaryfrv;
  int nvaryf;
  int nvarym;

  char **vnames;
  char **vunits;

  int *ldtype;

  double hjdoff;

  int nband;
  int naddband;

  int *pj;
  int *pldlin1;
  int *pldlin2;
  int *pldnon1;
  int *pldnon2;

  double *fit_y;
  double *fit_wt;
  int *fit_dptr;
  int *fit_mptr;
};

/* Structure for passing around "extra" parameters */
struct fit_pset_entry {
  char name[PARNAME_MAX];
  double v;
  double sig;
  int vary;
  int havesig;
};

/* Global lists of parameter names, units and default rms */
extern char *parnames[NPARFIX];
extern char *parunits[NPARFIX];
extern double default_sigguess[NPARFIX];

/* Inline functions for error scaling */

static inline double errscale (struct fit_data *d, double *v, int meas) {
  double rv = d->yerr[meas];

  if(d->fiterr) {
    if(d->obstype == OBS_LC)
      rv *= v[d->perrsc];
    else if(d->obstype == OBS_RV)
      rv = sqrt(rv*rv + v[d->perrsc]*v[d->perrsc] / d->yinfwt[meas]);
  }

  return rv;
}

static inline double varscale (struct fit_data *d, double *v, int meas) {
  double rv = d->yerr[meas]*d->yerr[meas];

  if(d->fiterr) {
    if(d->obstype == OBS_LC)
      rv *= v[d->perrsc]*v[d->perrsc];
    else if(d->obstype == OBS_RV)
      rv += v[d->perrsc]*v[d->perrsc] / d->yinfwt[meas];
  }

  return rv;
}

/* -- common.c -- */

void acf (double *a, int nsamp, char *xtitle, double scfac, FILE *ofp);
void tprintf (FILE *fp, const char *fmt, ...);

/* -- constants.c -- */

void init_const (void);

/* -- fit.c -- */

int do_fit (struct fit_parms *par, int nsigma, FILE *ofp, char *errstr);

/* -- func.c -- */

void fit_func (struct fit_parms *par, int id,
               double *a,
               double *t, double *y, unsigned char *iecl, int nmeas,
               int iphi, int ifull, int icor);

/* -- mc.c -- */

int do_mc (struct fit_parms *par, FILE *ofp,
           char *filename, int nsim, int iseed, char *errstr);

int read_mc (struct fit_parms *par, FILE *ofp,
             char **mcfilelist, int nmcfile,
             char *errstr);

/* -- parms.c -- */

int make_parms (struct fit_parms *par, FILE *ofp,
                struct fit_data *dlist, int ndata,
                int ldtype[2],
                double *vfix, int *varyfix, double *vsgfix,
                struct fit_pset_entry *pextra, int npextra,
                char *errstr);

/* -- plots.c -- */

void init_plots (char *pgdev);
void close_plots (void);

int do_plots (struct fit_parms *par,
	      FILE *ofp,
              char **filtnamelist, int nfiltname,
              char *errstr);

void plot_fried_eggs (double *ainit, char **anames, int nvary,
                      double *mc_res, int nalloc, int nsimd,
                      double *aperc, double *abest);

void plot_derived_hist (double *mc_der, int nalloc, int nsimd,
                        double *derperc, double *vderbest);

#endif  /* EBMC_H */
