eb
==

A bare bones implementation of a model for well-detached eclipsing
binaries.  The model is essentially a reimplementation of the light
curve generator used in Irwin et al. (2011) on LSPM J1112+7626, which
was based on software by Southworth et al. (2004 - 2009) and Etzel
(1981).

The light curve generator has been rewritten in C, using the analytic
method of Mandel & Agol (2002) for the quadratic limb darkening law.
The implementation aims for nearly full machine precision, although
does not quite achieve this at the present time.  The treatment of
ellipsoidal variation is based on Binnendijk (1974) and Etzel (1981),
and the reflection effect on Milne (1926), Russell (1939), and Etzel
(1981).

Most of the additions made in the 2011 paper have been carried over,
including the spot model and the corrections for the "classical" light
travel time across the system (in the solar system, the effect is
called the Romer delay).  The model can also compute light ratios and
radial velocities.

A Python wrapper is provided, and is intended to be the primary way to
use the model.  The original Monte Carlo program (in C) from the paper
will also be provided as an example, once the source has been cleaned
up (the original is messy and has many dependencies).

The file API.txt describes how to use the light curve generator, and
gives the meanings of the parameters.

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

Docstrings are embedded in the module, and should (hopefully) be
sufficient to figure out how to use the Python API in conjunction with
the examples (*.py) and the C API documentation.  The "eb_" or "EB_"
prefix used to avoid namespace pollution in C is not needed in the
Python bindings and is therefore omitted.

Building
========

The standard Makefile builds the light curve generator and Python
module by default when you type "make".

Automatic dependency generation for header files is supported, and
recommended if you plan to update using git or make modifications.
Run "make depend" before "make" to get this.

Acknowledgments
===============

The author gratefully acknowledges the work of Paul B. Etzel, Daniel
Popper, John Southworth, and Alvaro Gimenez.  These authors are
thanked for their contributions to this research area, and for making
their software available, upon which this model is based.

