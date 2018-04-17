#!/usr/bin/env python
#
# A noddy example to exercise most of the features in the "eb" module.
# Demonstrates how I recommend filling in the parameter vector - this
# way internal rearrangements of the vector as the model evolves won't
# break all of your scripts.
#

import numpy
import eb
import matplotlib.pyplot as plt

# Allocate main parameter vector, init to zero.
parm = numpy.zeros(eb.NPAR, dtype=numpy.double)

# These are the basic parameters of the model.
parm[eb.PAR_J]      =  0.799853  # J surface brightness ratio
parm[eb.PAR_RASUM]  =  0.015548  # (R_1+R_2)/a
parm[eb.PAR_RR]     =  0.785989  # R_2/R_1
parm[eb.PAR_COSI]   =  0.004645  # cos i

# Mass ratio is used only for computing ellipsoidal variation and
# light travel time.  Set to zero to disable ellipsoidal.
parm[eb.PAR_Q]      =  0.6904

# Light travel time coefficient.
ktot = 55.602793  # K_1+K_2 in km/s
cltt = 1000*ktot / eb.LIGHT

# Set to zero if you don't need light travel correction (it's fairly slow
# and can often be neglected).
parm[eb.PAR_CLTT]   =  cltt      # ktot / c

# Radiative properties of star 1.
parm[eb.PAR_LDLIN1] =  0.2094    # u1 star 1
parm[eb.PAR_LDNON1] =  0.6043    # u2 star 1
parm[eb.PAR_GD1]    =  0.32      # gravity darkening, std. value
parm[eb.PAR_REFL1]  =  0.4       # albedo, std. value

# Spot model.  Assumes spots on star 1 and not eclipsed.
parm[eb.PAR_ROT1]   =  0.636539  # rotation parameter (1 = sync.)
parm[eb.PAR_FSPOT1] =  0.0       # fraction of spots eclipsed
parm[eb.PAR_OOE1O]  =  0.0       # base spottedness out of eclipse
parm[eb.PAR_OOE11A] =  0.006928  # *sin
parm[eb.PAR_OOE11B] =  0.005088  # *cos

# PAR_OOE12* are sin(2*rot*omega) on star 1,
# PAR_OOE2* are for spots on star 2.

# Assume star 2 is the same as star 1 but without spots.
parm[eb.PAR_LDLIN2] = parm[eb.PAR_LDLIN1]
parm[eb.PAR_LDNON2] = parm[eb.PAR_LDNON1]
parm[eb.PAR_GD2]    = parm[eb.PAR_GD1]
parm[eb.PAR_REFL2]  = parm[eb.PAR_REFL1]

# Orbital parameters.
parm[eb.PAR_ECOSW]  =  0.152408  # ecosw
parm[eb.PAR_ESINW]  =  0.182317  # esinw
parm[eb.PAR_P]      = 41.032363  # period
parm[eb.PAR_T0]     = 2455290.046183  # T0 (epoch of primary eclipse)

# OTHER NOTES:
#
# To do standard transit models (a'la Mandel & Agol),
# set J=0, q=0, cltt=0, albedo=0.
# This makes the secondary dark, and disables ellipsoidal and reflection.
#
# The strange parameterization of radial velocity is to retain the
# flexibility to be able to model just light curves, SB1s, or SB2s.
#
# For improved precision, it's best to subtract most of the "DC offset"
# from T0 and the time array (e.g. take off the nominal value of T0 or
# the midtime of the data array) and add it back on at the end when
# printing parm[eb.PAR_T0] and vder[eb.PAR_TSEC].  Likewise the period
# can cause scaling problems in minimization routines (because it has
# to be so much more precise than the other parameters), and may need
# similar treatment.

# Simple (but not astronomer friendly) dump of model parameters.
print("Model parameters:")

for name, value, unit in zip(eb.parnames, parm, eb.parunits):
  print("{0:<10} {1:14.6f} {2}".format(name, value, unit))

# Derived parameters.
vder = eb.getvder(parm, -61.070553, ktot)

print("Derived parameters:")

for name, value, unit in zip(eb.dernames, vder, eb.derunits):
  print("{0:<10} {1:14.6f} {2}".format(name, value, unit))

# Phases of contact points.
(ps, pe, ss, se) = eb.phicont(parm)
if ps > 0.5:
  ps -= 1.0

# Use max(duration) for sampling range.
pdur=pe-ps
sdur=se-ss
if pdur > sdur:
  mdur = pdur
else:
  mdur = sdur

# ...centered on the average of start and end points.
pa=0.5*(ps+pe)
sa=0.5*(ss+se)

# Generate phase array: primary, secondary, and out of eclipse
# in leading dimension.
phi = numpy.empty([3, 1000], dtype=numpy.double)
phi[0] = numpy.linspace(pa-mdur, pa+mdur, phi.shape[1])
phi[1] = numpy.linspace(sa-mdur, sa+mdur, phi.shape[1])
phi[2] = numpy.linspace(-0.25, 1.25, phi.shape[1])

# All magnitudes.
typ = numpy.empty_like(phi, dtype=numpy.uint8)
typ.fill(eb.OBS_MAG)

# These calls both do the same thing.  First, phase.
y = eb.model(parm, phi, typ, eb.FLAG_PHI)
# Alternative using time.
#t = parm[eb.PAR_T0] + phi*parm[eb.PAR_P]
#y = eb.model(parm, t, typ)

# Plot eclipses in top 2 panels.  Manual y range, forced same on
# both plots.  x range is already guaranteed to be the same above.
ymin = y.min()
ymax = y.max()
yrange = ymax - ymin

plt.subplot(2, 2, 1)
plt.ylim(ymax+0.05*yrange, ymin-0.05*yrange)
plt.plot(phi[0], y[0])

plt.plot([ps, ps], [ymax, ymin], linestyle='--')
plt.plot([pe, pe], [ymax, ymin], linestyle='--')

plt.subplot(2, 2, 2)
plt.ylim(ymax+0.05*yrange, ymin-0.05*yrange)
plt.plot(phi[1], y[1])

plt.plot([ss, ss], [ymax, ymin], linestyle='--')
plt.plot([se, se], [ymax, ymin], linestyle='--')

# Out of eclipse plot across the bottom with 5*sigma clipping, robust
# MAD estimator.
median = numpy.median(y[2])
adiff = numpy.absolute(y[2]-median)
sigma = 1.48*numpy.median(adiff)
tmp = numpy.compress(adiff < 5*sigma, y[2])

ymin = tmp.min()
ymax = tmp.max()
yrange = ymax - ymin

plt.subplot(2, 1, 2)
plt.ylim(ymax+0.05*yrange, ymin-0.05*yrange)
plt.plot(phi[2], y[2])

plt.show()
