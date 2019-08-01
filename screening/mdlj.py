# routine accociated with non-bonded forces


import numpy
import math

# ----------------------------------------------------------
def zero_momentum(masses,vel):  # zero the liniar momentum
    mom = masses*vel # get momentum
    tmom = numpy.sum(mom,axis=0)/numpy.sum(masses,axis=0) # total mom/mass
    vel -= tmom # zero out

# ----------------------------------------------------

