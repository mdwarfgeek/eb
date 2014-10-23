eb
==

A bare bones implementation of a model for well-detached eclipsing
binaries.  The model is very similar to the widely-used EBOP and
JKTEBOP models, and is essentially a reimplementation of the light
curve generator used in Irwin et al. (2011) on LSPM J1112+7626, which
was based on JKTEBOP version 14 (since that time, many of the features
I added myself have been implemented independently in JKTEBOP).

The light curve generator has been rewritten in C, using the analytic
method of Mandel & Agol (2002) for the quadratic limb darkening law.
The implementation aims for nearly full machine precision, although
does not quite achieve this at the present time.  The treatments of
ellipsoidal variation and the reflection effect are intentionally very
similar to JKTEBOP for compatibility, as is the composition of the
parameter vector.  Some features I didn't use were not reimplemented,
including the other non-linear limb darkening types, which cannot be
done using elliptic integrals.

Most of the additions made in the 2011 paper have been carried over,
including the spot model and the corrections for the "classical" light
travel time across the system (in the solar system, the effect is
called the Romer delay).  Reflection was modified to be able to use
albedo rather than the original EBOP-like reflection parameter if
desired.  The model can also compute light ratios and radial
velocities.

A Python wrapper is provided, and is intended to be the primary way to
use the model.  The original Monte Carlo program (in C) from the paper
will also be provided as an example, once the source has been cleaned
up (the original is messy and has many dependencies).

Dependencies
============

Light curve generator (directory src)
-------------------------------------

Only a C99 compiler should be needed.  Some x86-specific assembly
optimizations are included, which need a compiler supporting GNU
syntax for inline assembler.  If another compiler or architecture is
used, these are currently disabled.  More assembly language
optimizations (particularly of the elliptic integrals) are in the
works.

Python module (directory python)
--------------------------------

Python >= 2.6  (for test program, module may work with earlier versions)
Numpy  >= 1.4  (advised; earlier versions may work)

Python 3 is not yet supported.

Building
========

The standard Makefile builds the light curve generator and Python
module by default when you type "make".

Automatic dependency generation for header files is supported, and
recommended if you plan to update using git or make modifications.
Run "make depend" before "make" to get this.

Acknowledgments
===============

The author gratefully acknowledges the original JKTEBOP and EBOP
authors John Southworth and Paul B. Etzel.  I have tried to avoid
reusing any code from other authors, but the full equations for
ellipsoidal variation and reflection are not available in any
published source to my knowledge, so it was necessary to get them from
the JKTEBOP Fortran code.  The series expansion at the end of the
calculation for R_J is based on the one in the public domain SLATEC
implementation, file drj.f by Bille Carlson.

