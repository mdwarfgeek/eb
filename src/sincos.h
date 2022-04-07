#ifndef SINCOS_H
#define SINCOS_H

#include <math.h>

#if defined(__GLIBC__)
#include <features.h>
#endif

#if defined(__GLIBC__) && defined(_GNU_SOURCE)

/* Try to use glibc provided routines in math.h */
#define inline_sincos(a, s, c) sincos((a), &(s), &(c))
#define inline_bare_sincos(a, s, c) sincos((a), &(s), &(c))

#else

/* Generic implementation using standard math.h functions */
#define inline_sincos(a, s, c) (s) = sin(a); (c) = cos(a)
#define inline_bare_sincos(a, s, c) (s) = sin(a); (c) = cos(a)

#endif

#endif  /* SINCOS_H */
