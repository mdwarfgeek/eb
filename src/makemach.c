#include <stdio.h>
#include <math.h>
#include <float.h>

/* Write out a header file containing some precomputed machine-specific
   constants used by the routines to save on runtime computation. */

int main (int argc, char *argv[]) {
  printf("#ifndef MACHCONST_H\n"
         "#define MACHCONST_H\n"
         "\n"
         "#define SQRT_FLT_MIN     %.*e\n"
         "#define SQRT_DBL_MIN     %.*le\n"
         "\n"
         "#define SQRT_FLT_EPSILON %.*e\n"
         "#define SQRT_DBL_EPSILON %.*le\n"
         "\n"
         "#endif  /* MACHCONST_H */\n",
         FLT_DIG+4,  /* on the safe side! */
         sqrtf(FLT_MIN),
         DBL_DIG+4,
         sqrt(DBL_MIN),
         FLT_DIG+4,
         sqrtf(FLT_EPSILON),
         DBL_DIG+4,
         sqrt(DBL_EPSILON));

  return(0);
}
