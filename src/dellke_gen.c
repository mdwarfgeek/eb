#include <stdlib.h>
#include <math.h>

/* Coefficients for the n=9 case from Cody (1965) "Chebyshev Approximations
   for the Complete Elliptic Integrals K and E".  According to this source
   the error in the approximation itself should be below 1.45e-16 in the
   worst case (E).  This is close to machine epsilon.

   The peculiar arrangement of this array is optimized for summing
   both polynomials simultaneously using bunched, back-to-back
   additions and multiplications with Horner's method.  By doing it
   in this order, we take advantage of pipelining by avoiding
   dependencies between the consecutive multiply instructions.
   This is more for the assembly language implementations than the
   C implementation.  On x87 this causes the compiler to emit a 
   sequence of fmul and fxch instructions, which is okay on the
   Pentium (where the fxch gets combined with the fmul) but may
   be slower than summing the polynomials one at a time on earlier
   CPUs or other architectures. */

static const double dellke_coeff[40] = {
  /*         K                       E         */
  3.00725199036864838e-4, 3.25192015506390418e-4,  /* a,c_9 */
  6.66317524646073151e-5, 7.20316963457154599e-5,  /* b,d_9 */
  3.96847090209897819e-3, 4.30253777479311659e-3,  /* a,c_8 */
  1.72161470979865212e-3, 1.86453791840633632e-3,  /* b,d_8 */
  1.07959904905916349e-2, 1.17858410087339355e-2,  /* a,c_7 */
  9.28116038296860419e-3, 1.00879584943751004e-2,  /* b,d_7 */
  1.05899536209893585e-2, 1.18419259955012494e-2,  /* a,c_6 */
  2.06902400051008404e-2, 2.26603098916041221e-2,  /* b,d_6 */
  7.51938672180838102e-3, 9.03552773754088184e-3,  /* a,c_5 */
  2.95037293486887130e-2, 3.28110691727210618e-2,  /* b,d_5 */
  8.92664629455646620e-3, 1.17167669446577228e-2,  /* a,c_4 */
  3.73355466822860296e-2, 4.26725101265917523e-2,  /* b,d_4 */
  1.49420291422820783e-2, 2.18361314054868967e-2,  /* a,c_3 */
  4.88271550481180099e-2, 5.85927071842652739e-2,  /* b,d_3 */
  3.08851730018997099e-2, 5.68052233293082895e-2,  /* a,c_2 */
  7.03124954595466082e-2, 9.37499951163670673e-2,  /* b,d_2 */
  9.65735903017425285e-2, 4.43147180583368137e-1,  /* a,c_1 */
  1.24999999997640658e-1, 2.49999999997461423e-1,  /* b,d_1 */
  2*M_LN2,                1.0                   ,  /* a,c_0 */
  0.5,                    0.0                      /* b,d_0 */
};

/* Computes K(k) and E(e) given y=1-k^2.  v[0] receives K and v[1] receives E.
   The output vector format is to maintain compatibility with the assembly
   language routine for SSE2. */

void dellke_gen (double y, double v[2]) {
  double sa, sb, sc, sd, ly;
  int i = 0;

  ly = log(y);

  /* Sum polynomials */
  sa = dellke_coeff[i++] * y;
  sc = dellke_coeff[i++] * y;
  sb = dellke_coeff[i++] * y;
  sd = dellke_coeff[i++] * y;

  while(i < 36) {
    sa = (sa + dellke_coeff[i++]) * y;
    sc = (sc + dellke_coeff[i++]) * y;
    sb = (sb + dellke_coeff[i++]) * y;
    sd = (sd + dellke_coeff[i++]) * y;
  }

  sa += dellke_coeff[i++];
  sc += dellke_coeff[i++];
  sb += dellke_coeff[i++];
  /* d0 = 0, so omitted */

  /* Result */
  v[0] = sa - sb*ly;
  v[1] = sc - sd*ly;
}

