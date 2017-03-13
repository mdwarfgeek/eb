#ifndef EB_H
#define EB_H

#include <math.h>

/* Fundamental (defining) constants, mostly IAU and IERS */
#define EB_GMSUN  1.3271244e20      /* m^3 / s^2, IAU 2015 Resol B3 */
#define EB_AU     149597870700.0    /* m, IAU 2009 system */
#define EB_LIGHT  2.99792458e8      /* m/s, definition */

/* Other astrophysical quantities */
#define EB_RSUN   6.957e8        /* m, IAU 2015 Resol B3 */

/* This is pretty useful too */
#ifndef TWOPI
#define TWOPI     (2.0*M_PI)
#endif

/* Parameters in the "parm" array.  Order is for compatibility with jktebop. */
#define EB_PAR_J        0
#define EB_PAR_RASUM    1
#define EB_PAR_RR       2
#define EB_PAR_LDLIN1   3
#define EB_PAR_LDLIN2   4
#define EB_PAR_COSI     5
#define EB_PAR_ECOSW    6 
#define EB_PAR_ESINW    7
#define EB_PAR_GD1      8
#define EB_PAR_GD2      9
#define EB_PAR_REFL1   10  /* albedo (default) or reflection */
#define EB_PAR_REFL2   11
#define EB_PAR_Q       12
#define EB_PAR_TIDANG  13  /* deg */
#define EB_PAR_L3      14
#define EB_PAR_PHI0    15  /* [0,1] */
#define EB_PAR_M0      16
#define EB_PAR_INTEG   17  /* not used */
#define EB_PAR_P       18
#define EB_PAR_T0      19  /* epoch of inferior conjunction if phi0=0 */
#define EB_PAR_LDNON1  20
#define EB_PAR_LDNON2  21
#define EB_PAR_CLTT    22  /* ktot / c if want light travel time corr */
#define EB_PAR_ROT1    23  /* rotation parameter */
#define EB_PAR_ROT2    24
#define EB_PAR_FSPOT1  25  /* fraction of spots eclipsed */
#define EB_PAR_FSPOT2  26
#define EB_PAR_OOE1O   27  /* base spottedness out of eclipse star 1 */
#define EB_PAR_OOE11A  28  /* sin(...) coeff star 1 */
#define EB_PAR_OOE11B  29  /* cos(...) coeff star 1 */
#define EB_PAR_OOE12A  30  /* sin(2*...) coeff star 1 */
#define EB_PAR_OOE12B  31  /* cos(2*...) coeff star 1 */
#define EB_PAR_OOE2O   32  /* base spottedness out of eclipse star 2 */
#define EB_PAR_OOE21A  33  /* sin(...) coeff star 2 */
#define EB_PAR_OOE21B  34  /* cos(...) coeff star 2 */
#define EB_PAR_OOE22A  35  /* sin(2*...) coeff star 2 */
#define EB_PAR_OOE22B  36  /* cos(2*...) coeff star 2 */
#define EB_PAR_DWDT    37  /* apsidal precession rate (rad/day) */
#define EB_NPAR        38

/* Flags */
#define EB_FLAG_REFL 0x01  /* reflection rather than albedo */
#define EB_FLAG_PHI  0x02  /* phase rather than time */

/* Observation types we can request from the light curve generator. */
#define EB_OBS_MAG      0  /* magnitude */
#define EB_OBS_LIGHT    1  /* total normalized light (as Agol) */
#define EB_OBS_LRAT     2  /* L_2 / L_1 */
#define EB_OBS_AVLR     3  /* orbit averaged L_2 / L_1 */
#define EB_OBS_VRAD1    4  /* cos(v+w) + e cos w for star 1 */
#define EB_OBS_VRAD2    5  /* cos(v+w) + e cos w for star 2 */
#define EB_OBS_PSS      6  /* plane of sky separation ("d" in our notation) */
#define EB_OBS_A        7  /* d - rr, used to test if in eclipse (< 1) */
#define EB_OBS_LSS      8  /* line of sight separation */

/* Function prototype for the generator.  Fills "out" with model
   observables.  The array "t" is a time array in the same units as P
   and T0, and "typ" is an array of the OBS_* parameters specifying
   the observables to be computed.  These may appear in any combination,
   but note that a lot of unnecessary computations are done at each
   function call if all you want are radial velocities. */
void eb_model_dbl (double *parm, double *t, double *ol1, double *ol2,
                   unsigned char *typ, double *out, unsigned char *iecl,
                   int flags, int npt);
void eb_model_flt (double *parm, double *t, float *ol1, float *ol2,
                   unsigned char *typ, float *out, unsigned char *iecl,
                   int flags, int npt);

/* Utility subroutines from ebutil.c */
double eb_phisec (double esinw, double ecosw);
void eb_phicont (double esinw, double ecosw, double cosi,
                 double d, double *phi);

/* Parameters in the "vder" (derived) array */
#define EB_PAR_I        0
#define EB_PAR_R1A      1
#define EB_PAR_R2A      2
#define EB_PAR_E        3
#define EB_PAR_OMEGA    4
#define EB_PAR_A        5
#define EB_PAR_MTOT     6
#define EB_PAR_M1       7
#define EB_PAR_M2       8
#define EB_PAR_RTOT     9
#define EB_PAR_R1      10
#define EB_PAR_R2      11
#define EB_PAR_LOGG1   12
#define EB_PAR_LOGG2   13
#define EB_PAR_VSYNC1  14
#define EB_PAR_VSYNC2  15
#define EB_PAR_TSYNC   16
#define EB_PAR_TCIRC   17
#define EB_PAR_TSEC    18
#define EB_PAR_DURPRI  19
#define EB_PAR_DURSEC  20

#define EB_NDER        21

void eb_getvder (double *v, double gamma, double ktot, double *vder);

/* Arrays of names and units for parameters */
extern char *eb_parnames[EB_NPAR];
extern char *eb_parunits[EB_NPAR];
extern char *eb_dernames[EB_NDER];
extern char *eb_derunits[EB_NDER];

#endif  /* EB_H */
