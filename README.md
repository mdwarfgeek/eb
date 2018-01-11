eb
==

A bare bones implementation of a model for well-detached eclipsing
binaries.  The model is essentially a reimplementation of the light
curve generator used in Irwin et al. (2011) on LSPM J1112+7626, which
was a modified version of the Fortran subroutine "light" by Southworth
et al. (2004 - 2009) and Etzel (1981).

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
use the model, but it can also be called directly from C and via a
thin wrapper from Fortran.  An experimental IDL wrapper is provided in
the subdirectory "idl", but I don't use IDL so it has had very little
testing.

An updated version of the Monte Carlo program (in C) from the 2011
paper is also provided in "ebmc".  It has a number of dependencies
that may make it tricky to get working and is therefore not included
in the default set of Makefile targets.  This is very bare bones, not
general, messy, and has no user interface to speak of.  It is not
intended to be used as-is, but may serve as a place to borrow ideas
from.  It has accumulated a large number of changes since the version
used for the paper, and therefore does not exactly reproduce the
published results.

The file API.txt describes how to use the light curve generator, and
gives the meanings of the parameters.  Some practical notes on fitting
real data intended to assist less experienced researchers are also
included in the file HINTS.txt.  Note that these documents are not an
introductory text on eclipsing binary models, and it is assumed the
reader is already familiar with this class of model.

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

C fitter / Monte Carlo program (directory ebmc)
-----------------------------------------------

My C subroutine library from "lib" on github.  The Makefiles currently
assume it's located in a directory "lib" at the same level as the "eb"
directory (i.e. ../lib if this README is in the current directory).

MPFIT, C version
http://www.physics.wisc.edu/~craigm/idl/cmpfit.html

The original program used levmar instead of MPFIT to provide the
Levenberg-Marquardt minimizer.  This is still supported (for now) as
an alternative to MPFIT by defining USE_LEVMAR (see make.inc).
http://users.ics.forth.gr/~lourakis/levmar/

PGPLOT >= 5.1.1 (optional, plots disabled if unavailable)
http://www.astro.caltech.edu/~tjp/pgplot/

The colour Postscript ("cps") driver is assumed to be available.
Everything else is optional.  Successfully linking to PGPLOT can be
tricky, particularly on Mac OS X.  Some suggestions are included as
comments in make.inc.  To disable the plots, override PGPLOT_INC to
remove the -DHAVE_PGPLOT.

Building
========

The standard Makefile builds the light curve generator and Python
module by default when you type "make".  To build the IDL wrapper,
type "make idl" (there is no need to type "make" first).

To build the fitter / Monte Carlo program, change to the ebmc
directory and run "make".  You may need to alter the top level
make.inc file or override some of the variables so it can find the
dependencies or to disable some features.

Automatic dependency generation for header files is supported, and
recommended if you plan to update using git or make modifications.
Run "make depend" before "make" to get this.

Acknowledgments
===============

The author gratefully acknowledges the work of Paul B. Etzel, Daniel
Popper, John Southworth, and Alvaro Gimenez.  These authors are
thanked for their contributions to this research area, and for making
their software available, upon which this model is based.  Guillermo
Torres is thanked for his extensive contributions to the development,
testing and documentation of the present software.

