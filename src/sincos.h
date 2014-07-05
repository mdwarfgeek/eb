#ifndef SINCOS_H
#define SINCOS_H

#include <math.h>

#if defined(__GLIBC__)
#include <features.h>
#endif

#if defined(__GNUC__) && (defined(__i386) || defined(__amd64))

/* Inline x86 assembler routines, GNU C syntax */
#define inline_sincos(a, s, c)   \
  asm("fsincos"           "\n\t" \
      "fnstsw  %%ax"      "\n\t" \
      "btw     $10, %%ax" "\n\t" \
      "jnc     1f"        "\n\t" \
      "fldpi"             "\n\t" \
      "fadd    %%st(0)"   "\n\t" \
      "fxch    %%st(1)"   "\n\t" \
"2: " "fprem1"            "\n\t" \
      "fnstsw  %%ax"      "\n\t" \
      "btw     $10, %%ax" "\n\t" \
      "jc      2b"        "\n\t" \
      "fstp    %%st(1)"   "\n\t" \
      "fsincos"           "\n\t" \
"1: "                            \
      : "=t" (c), "=u" (s)       \
      : "0" (a)                  \
      : "ax", "cc")

#define inline_bare_sincos(a, s, c) \
  asm("fsincos" : "=t" (c), "=u" (s) : "0" (a))

#elif defined(__GLIBC__) && defined(_GNU_SOURCE)

/* Try to use glibc provided routines in math.h */
#define inline_sincos(a, s, c) sincos((a), &(s), &(c))
#define inline_bare_sincos(a, s, c) sincos((a), &(s), &(c))

#else

/* Generic implementation using standard math.h functions */
#define inline_sincos(a, s, c) (s) = sin(a); (c) = cos(a)
#define inline_bare_sincos(a, s, c) (s) = sin(a); (c) = cos(a)

#endif

#endif  /* SINCOS_H */
