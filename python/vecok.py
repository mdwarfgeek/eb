# This example implements some of the common boundary conditions needed
# on the "eb" parameter vector in an MCMC or similar solution.

import math
import numpy

import eb

def vecok (parm):
  ecc = math.hypot(parm[eb.PAR_ECOSW], parm[eb.PAR_ESINW])

  ol1 = parm[eb.PAR_OOE1O] + 2.0*math.sqrt(parm[eb.PAR_OOE11A]**2+
                                           parm[eb.PAR_OOE11B]**2+
                                           parm[eb.PAR_OOE12A]**2+
                                           parm[eb.PAR_OOE12B]**2)

  ol2 = parm[eb.PAR_OOE2O] + 2.0*math.sqrt(parm[eb.PAR_OOE21A]**2+
                                           parm[eb.PAR_OOE21B]**2+
                                           parm[eb.PAR_OOE22A]**2+
                                           parm[eb.PAR_OOE22B]**2)

  # Maximum cosine of inclination where there is (barely) a
  # primary eclipse at inferior conjunction.
  # Note that it is possible to have an eclipse around but not
  # at inferior conjunction for some peculiar configurations
  # of eccentric orbits, so we should really figure out how to
  # do a better job of this (root finding on plane of sky sep?).
  pmcosi = parm[eb.PAR_RASUM] * (1.0 + parm[eb.PAR_ESINW]) / (1.0-ecc*ecc)

  # Similar for secondary eclipse.
  smcosi = parm[eb.PAR_RASUM] * (1.0 - parm[eb.PAR_ESINW]) / (1.0-ecc*ecc)

  # Maximum cosine of inclination, here we require there is
  # an eclipse (primary or secondary).  This prevents the
  # solution getting stuck when it strays into regions of
  # parameter space with no eclipse, at which point the
  # derivative vanishes.
  max_cosi = max(pmcosi, smcosi)

  return(parm[eb.PAR_RASUM] >= 0 and
         parm[eb.PAR_RR] >= 0 and
         parm[eb.PAR_COSI] >= 0 and
         parm[eb.PAR_COSI] < max_cosi and  # has an eclipse
         parm[eb.PAR_LDLIN1]+parm[eb.PAR_LDNON1] >= 0 and  # LD triangle
         parm[eb.PAR_LDLIN1]+parm[eb.PAR_LDNON1] <= 1 and  # star 1
         parm[eb.PAR_LDLIN1] >= 0 and
         parm[eb.PAR_LDLIN1]+2*parm[eb.PAR_LDNON1] >= 0 and
         parm[eb.PAR_FSPOT1] >= 0 and parm[eb.PAR_FSPOT1] <= 1 and
         parm[eb.PAR_LDLIN2]+parm[eb.PAR_LDNON2] >= 0 and  # LD triangle
         parm[eb.PAR_LDLIN2]+parm[eb.PAR_LDNON2] <= 1 and  # star 2
         parm[eb.PAR_LDLIN2] >= 0 and
         parm[eb.PAR_LDLIN2]+2*parm[eb.PAR_LDNON2] >= 0 and
         parm[eb.PAR_FSPOT2] >= 0 and parm[eb.PAR_FSPOT2] <= 1 and
         ecc < 1.0 and
         ol1 < 1.0 and
         ol2 < 1.0)
