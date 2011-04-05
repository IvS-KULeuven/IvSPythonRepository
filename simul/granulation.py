import sys
from numpy import *
from numpy.random import *



#------------------------------------------
def granulation(time, timescale, varscale):
#------------------------------------------

  #-- Input: time[0..Ntime-1]: time points
  #--        timescale[0..Ncomp-1]: time scale tau of each component
  #--                               of the granulation/magnetic activity
  #--        varscale[0..Ncomp-1]: variation scale of each component
  #--                              of the granulation/magnetic activity
  #--                              in the appropriate passband

  Ntime = len(time)
  Ncomp = len(timescale)

  print >> sys.stderr, "Simulating", Ncomp, "component(s)."

  # Set the kick (= reexcitation) timestep to be one 100th of the
  # shortest granulation time scale (i.e. kick often enough).

  kicktimestep = min(timescale) / 100.0
  print >> sys.stderr, "Kicktimestep = ", kicktimestep

  # Predefine some arrays

  signal = zeros(Ntime)
  granul = zeros(Ncomp)
  mu = zeros(Ncomp)
  sigma = sqrt(kicktimestep/timescale)*varscale

  # Warm up the first-order autoregressive process

  print >> sys.stderr, "Warming up..."

  for i in range(2000):
    granul = granul * (1.0 - kicktimestep / timescale) + normal(mu, sigma)

  # Start simulating the granulation time series

  print >> sys.stderr, "Simulating granulation signal"

  delta = 0.0
  currenttime = time[0] - kicktimestep

  for i in range(Ntime):

    # Compute the contribution of each component separately.
    # First advance the time series right *before* the time point i,

    while((currenttime + kicktimestep) < time[i]):
      granul = granul * (1.0 - kicktimestep / timescale) + normal(mu,sigma)
      currenttime = currenttime + kicktimestep

    # Then advance the time series with a small time step right *on* time[i]

    delta = time[i] - currenttime
    granul = granul * (1.0 - delta / timescale)    \
               + normal(mu, sqrt(delta/timescale)*varscale)
    currenttime = time[i]

    # Add the different components to the signal. 

    signal[i] = sum(granul)


  # That's it!

  return(signal)

