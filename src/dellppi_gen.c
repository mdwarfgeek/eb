/* These tolerances should be enough to reach machine precision. */
#define CTOL 1.0e-3
#define JTOL 1.0e-3

#define DATATYPE double
#define FUNC     dellppi_gen
#define TSQRT    sqrt
#define TABS     fabs

#include "ellppi_gen.c"
