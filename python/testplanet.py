#!/usr/bin/env python
#
# Simple planet and fitting example, works in normalized phase,
# (t - T0) / period, where T0 is the transit epoch.  The phases
# don't need to be wrapped to [0,1] range, and if using the out
# of transit modulation feature with a non-synchronized spin and
# orbit, they shouldn't be wrapped.
#
# The parameters are for Kepler-37b.
#
# This example uses the "phase offset" (really a mean anomaly
# adjustment) parameter to adjust the time of the eclipse for
# simplicity.  In practice we would presumably want to determine
# the ephemeris (P and T0) and would fit for those instead of
# adjusting PHI0.
#
# For non-linear least-squares, I have used "lmmin" in pwkit
# by Peter Williams in this example.  It can be found at:
#   https://github.com/pkgw/pwkit

import numpy
import eb
import matplotlib.pyplot as plt

try:
  from pwkit.lmmin import ResidualProblem
  do_fit=1
except ImportError:
  do_fit=0

# Allocate main parameter vector, init to zero.
parm = numpy.zeros(eb.NPAR, dtype=numpy.double)

# These are the basic parameters of the model.
parm[eb.PAR_RASUM]  = 0.03518      # (R_1+R_2)/a
parm[eb.PAR_RR]     = 0.00360      # R_2/R_1
parm[eb.PAR_COSI]   = 0.0239       # cos i

# Radiative properties of star.
parm[eb.PAR_LDLIN1] = 0.4323       # u1
parm[eb.PAR_LDNON1] = 0.2777       # u2

# Orbital parameters.
parm[eb.PAR_ECOSW]  = -0.47        # ecosw
parm[eb.PAR_ESINW]  = -0.54        # esinw
parm[eb.PAR_P]      = 13.367308    # period

# Simple (but not astronomer friendly) dump of model parameters.
print("Model parameters:")

for name, value, unit in zip(eb.parnames, parm, eb.parunits):
  print("{0:<10} {1:14.6f} {2}".format(name, value, unit))

# Phases of contact points.
(ps, pe, ss, se) = eb.phicont(parm)
if ps > 0.5:
  ps -= 1.0

# Transit duration.
pdur=pe-ps

# Centre for plots.
pa=0.5*(ps+pe)

# Generate phase array, making sure there is a reasonable amount
# of baseline available for fitting.
phi = numpy.linspace(pa-2*pdur, pa+2*pdur, 5000)

# Tell the model to compute magnitudes (author's preference when
# fitting, can be changed to OBS_LIGHT to compute normalized flux).
typ = numpy.empty_like(phi, dtype=numpy.uint8)
typ.fill(eb.OBS_MAG)

# Compute model.
y = eb.model(parm, phi, typ, eb.FLAG_PHI)

# Are we fitting?
if do_fit:
  # Make sure we always get the same random deviates.
  numpy.random.seed(42)

  # Add noise.
  nois = 30.0e-6  # 30 ppm
  
  yobs = y + numpy.random.normal(loc=0.0, scale=nois, size=y.shape)

  # Observational uncertainties for fitter.  We'll be nice here and
  # use the correct answer.
  e_yobs = numpy.empty_like(y, dtype=numpy.double)
  e_yobs.fill(nois)

  # A simple way to interface this to "lmmin".
  def fit_func (trial_parm, ymod):
    eb.model(trial_parm, phi, typ, eb.FLAG_PHI, out=ymod)

  # New problem.
  p = ResidualProblem(eb.NPAR, yobs, 1.0/e_yobs, fit_func, None)

  # Start with all parameters fixed.
  for iparm, vparm in enumerate(parm):
    p.p_value(iparm, vparm, fixed=True)

  # Fit these ones, perturbed a bit.
  tofit = [eb.PAR_M0, eb.PAR_PHI0, eb.PAR_RASUM, eb.PAR_RR, eb.PAR_COSI]

  print("Initial guesses:")

  for iparm in tofit:
    # Decide perturbation.  Perturbing zero won't work too well in
    # relative so use absolute for those parameters.
    if iparm == eb.PAR_M0:
      pertur=0.1  # absolute
    elif iparm == eb.PAR_PHI0:
      pertur=0.01  # absolute
    else:
      pertur=0.5*parm[iparm]
      
    # Generate our initial guesses (intentionally bad).
    initval = numpy.random.normal(parm[iparm], scale=pertur)

    # Write out what we're going to use.
    print("{0:<10} {1:14.6f} {2}".format(eb.parnames[iparm],
                                         initval,
                                         eb.parunits[iparm]))

    # Set in the problem.
    p.p_value(iparm, parm[iparm], fixed=False)

  # Phase seems to confuse the automatic step size choice.
  # I found this to work a bit better.
  p.p_step(eb.PAR_PHI0, pdur*0.1)

  # Run fit.
  sol = p.solve()

  # Final parameters.
  print("Final parameters:")

  for iparm in tofit:
    print("{0:<10} {1:14.6f} +/- {2:8.6f} {3}".format(eb.parnames[iparm],
                                                      sol.params[iparm],
                                                      sol.perror[iparm],
                                                      eb.parunits[iparm]))
    
  print("chi squared = {0:f} ndof = {1:d}".format(sol.fnorm, sol.ndof))

  # Model (without weights).
  ymod = eb.model(sol.params, phi, typ, eb.FLAG_PHI)

  # Plot data and fit.
  plt.errorbar(phi, yobs, e_yobs)
  plt.plot(phi, ymod)

# Plot original model.
plt.plot(phi, y)

# Bright up!
plt.gca().invert_yaxis()

plt.show()
