#include <stdlib.h>
#include <math.h>

/* Coefficients for the n=4 case from Cody (1965) "Chebyshev Approximations
   for the Complete Elliptic Integrals K and E".  According to this source
   the error in the approximation itself should be below 1.57e-8 in the
   worst case (E).  This is close to machine epsilon.  */

static const float fellke_coeff[20] = {
  /*     K              E      */
  1.451338556e-2, 1.736314854e-2,  /* a,c_4 */
  4.418398230e-3, 5.263789328e-3,  /* b,d_4 */
  3.742539571e-2, 4.757404429e-2,  /* a,c_3 */
  3.328521016e-2, 4.069468414e-2,  /* b,d_3 */
  3.589980090e-2, 6.260761942e-2,  /* a,c_2 */
  6.880295505e-2, 9.200109374e-2,  /* b,d_2 */
  9.666338350e-2, 4.432515145e-1,  /* a,c_1 */
  1.249859468e-1, 2.499836641e-1,  /* b,d_1 */
  2*M_LN2,        1.0           ,  /* a,c_0 */
  0.5,            0.0              /* b,d_0 */
};

void fellke_gen (float y, float v[2]) {
  float sa, sb, sc, sd, ly;
  int i = 0;

  ly = logf(y);

  /* Sum polynomials */
  sa = fellke_coeff[i++] * y;
  sc = fellke_coeff[i++] * y;
  sb = fellke_coeff[i++] * y;
  sd = fellke_coeff[i++] * y;

  while(i < 16) {
    sa = (sa + fellke_coeff[i++]) * y;
    sc = (sc + fellke_coeff[i++]) * y;
    sb = (sb + fellke_coeff[i++]) * y;
    sd = (sd + fellke_coeff[i++]) * y;
  }

  sa += fellke_coeff[i++];
  sc += fellke_coeff[i++];
  sb += fellke_coeff[i++];
  /* d0 = 0, so omitted */

  /* Result */
  v[0] = sa - sb*ly;
  v[1] = sc - sd*ly;
}
