#include "ebmc.h"

/* Names and units, extended vectors including radial velocity */
char *parnames[NPARFIX];
char *parunits[NPARFIX];

/* Initial guesses at standard deviation of parameters.  This is used
 * to initialize the covariance matrix in fitting mode 2, when the
 * initial guesses are too close to zero.  In mode 1, the covariance
 * from the L-M is used.
 */
double default_sigguess[NPARFIX] = {
  0.01,    /* J */
  0.01,    /* r_1+r_2 */
  0.01,    /* R_2/R_1 */
  0.1,     /* mu_1 */
  0.1,     /* mu_2 */
  0.01,    /* cosi */
  0.001,   /* ecosw */
  0.01,    /* esinw */
  0.01,    /* grav_1 */
  0.01,    /* grav_2 */
  0.001,   /* refl_1 */
  0.001,   /* refl_2 */
  0.01,    /* q */
  0.01,    /* tidang */
  0.01,    /* L_3 */
  0.001,   /* Phi_0 */
  0.001,   /* m_0 */
  0,       /* integ (cannot be varied) */
  0.00001, /* P */
  0.001,   /* T_0 */
  0.1,     /* mup_1 */
  0.1,     /* mup_2 */
  0,       /* cltt (not varied) */
  0.001,   /* Rot_1 */
  0.001,   /* Rot_2 */
  0,       /* Fspot_1 (cannot be varied) */
  0,       /* Fspot_2 (cannot be varied) */
  0.1,     /* o_1 */
  0.001,   /* a_11 */
  0.001,   /* b_11 */
  0.001,   /* a_12 */
  0.001,   /* b_12 */
  0.1,     /* o_2 */
  0.001,   /* a_21 */
  0.001,   /* b_21 */
  0.001,   /* a_22 */
  0.001,   /* b_22 */
  1.0,     /* K_1+K_2 */
  1.0      /* gamma */
};

/* These strings need to be global to ensure they don't disappear */
static char *name_ktot = "K_1+K_2";
static char *name_gamma = "gamma";
static char *km_per_sec = "km/s";

void init_const (void) {
  int i;

  for(i = 0; i < EB_NPAR; i++) {
    parnames[i] = eb_parnames[i];
    parunits[i] = eb_parunits[i];
  }

  parnames[PAR_KTOT] = name_ktot;
  parnames[PAR_GAMMA] = name_gamma;

  parunits[PAR_KTOT] = km_per_sec;
  parunits[PAR_GAMMA] = km_per_sec;
}

