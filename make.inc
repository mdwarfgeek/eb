### Constants: edit to suit your system ###

# C compiler
#CC=gcc

# Optimization flags

# DEBUG:
#OPT=-g

# OPT:
OPT=-g -O3 -ffast-math

# Compiler flags
CFLAGS=-std=gnu99 -Wall $(OPT) -fPIC

### Dependencies for "ebmc".  Not needed for the model itself ###
### (directory "src") or for the Python module.               ###

# Fortran compiler.  Used to link "ebmc" due to the Fortran dependencies
# it has.  Should be the same compiler used to build PGPLOT.  "make"
# usually defaults to FC=f77, whereas you probably want to use gfortran
# on recent installs.
#FC=gfortran

# For L-M, we need MPFIT (C version) or Levmar.
# I prefer MPFIT.
LM_INC?=
LM_LIBS?=-lmpfit
# Levmar (used by older versions)
#LM_INC?=-DUSE_LEVMAR
#LM_LIBS?=-llevmar -llapack -lf77blas -lcblas -latlas

# PGPLOT.  These should be sufficient to link on standard UNIX platforms.
# If you have trouble, see LIBS in $PGPLOT_DIR/makefile.  Running ldd
# (or on Mac OS X, otool -L) on one of the PGPLOT demo programs usually
# found under $PGPLOT_DIR can also be a useful trick.
#
# Mac OS X installs seem to be particularly troublesome.  I'm using
# Fink 0.36.4 on mine, and there the following additions work for me:
# -L/sw/lib -laquaterm -Wl,-framework -Wl,Foundation -lpng
# 
PGPLOT_INC?=-I$(PGPLOT_DIR) -DHAVE_PGPLOT
PGPLOT_LIBS?=-L$(PGPLOT_DIR) -lpgplot -lcpgplot -L/usr/X11R6/lib -lX11

### End constants section ###

# Rules for building

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

