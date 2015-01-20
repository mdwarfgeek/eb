#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>

#include "ebmc.h"

void acf (double *a, int nsamp, char *xtitle, double scfac, FILE *ofp) {
  int isamp, jsamp, istart, iend, icl;
  double mean, var, sum, ac;

  /* Compute mean and variance */
  mean = 0.0;

  for(isamp = 0; isamp < nsamp; isamp++)
    mean += a[isamp];

  mean /= nsamp;

  var = 0.0;

  for(isamp = 0; isamp < nsamp; isamp++)
    var += (a[isamp]-mean)*(a[isamp]-mean);

  var /= nsamp;

  /* Step out */
  istart = 0;
  iend = 1;

  while(iend < nsamp) {
    /* Compute ACF */
    sum = 0;

    for(jsamp = 0; jsamp < nsamp-iend; jsamp++)
      sum += a[jsamp]*a[jsamp+iend];

    sum /= (nsamp-iend);

    ac = (sum - mean*mean) / var;

    if(ac > 0.5) {
      istart = iend;
      iend *= 2;
    }
    else
      break;
  }

  /* Binary chop to find correlation length */
  icl = (istart+iend)/2;
  ac = 0;

  while((iend-istart) >= 2) {
    /* Compute ACF */
    sum = 0;

    for(jsamp = 0; jsamp < nsamp-icl; jsamp++)
      sum += a[jsamp]*a[jsamp+icl];

    sum /= (nsamp-icl);

    ac = (sum - mean*mean) / var;

    if(ac > 0.5) {
      istart = icl;
      icl = (istart + iend)/2;
    }
    else {
      iend = icl;
      icl = (istart + iend)/2;
    }
  }

  tprintf(ofp, "%s correlation length: %d ac=%.3f\n", xtitle, icl, ac);

  if(icl > nsamp/100)
    tprintf(ofp, "WARNING: correlation length exceeds 1%% of simulation\n");
}

void tprintf (FILE *fp, const char *fmt, ...) {
  va_list ap;
  char buf[1024];

  va_start(ap, fmt);
  vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);

  fputs(buf, stdout);
  if(fp)
    fputs(buf, fp);
}

