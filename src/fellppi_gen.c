/* These tolerances should be enough to reach machine precision. */
#define CTOL 0.04
#define JTOL 0.05

#define DATATYPE float
#define FUNC     fellppi_gen
#define TSQRT    sqrtf
#define TABS     fabsf

#include "ellppi_gen.c"
