import sys
from numpy import *
from numpy.random import *





def solarosc(time, freq, ampl, eta):

  """
  # PURPOSE: compute time series of stochastically excited damped modes
  #
  # INPUT: . time[0..Ntime-1]: time points    (unit: e.g. Ms)
  #        . freq[0..Nmodes-1]: frequencies   (unit: e.g. microHz)
  #        . ampl[0..Nmodes-1]: amplitude of each of the oscillation modes
  #                             rms amplitude = ampl / sqrt(2)
  #        . eta[0..Nmodes-1]: damping rates  (unit: e.g. (Ms)^{-1})
  # 
  # OUTPUT: . signal[0..Ntime-1]
  """
  
  Ntime = len(time)
  Nmode = len(freq)

  print >> sys.stderr, "Simulating", Nmode, "mode(s)"

  # Set the kick (= reexcitation) timestep to be one 100th of the
  # shortest damping time. (i.e. kick often enough).

  kicktimestep = (1.0 / max(eta)) / 100.0
  print >> sys.stderr, "Kicktimestep = ", kicktimestep

  #-- Init start values of amplitudes, and the kicking amplitude
  #-- so that the amplitude of the oscillator will be on average be
  #-- constant and equal to the user given amplitude

  amplcos = 0.0
  amplsin = 0.0
  kick_amplitude = ampl * sqrt(kicktimestep * eta)

  #-- Warm up the stochastic excitation simulator to forget the
  #-- initial conditions. Do this during the longest damping time.
  #-- But put a maximum on the number of kicks, as there might
  #-- be almost-stable modes with damping time = infinity
  
  damp = exp(-eta * kicktimestep)
  Nwarmup = min(20000, int(floor(1.0 / min(eta) / kicktimestep)))

  print >> sys.stderr, Nwarmup, "kicks for warm up..."

  for i in range(Nwarmup):
    amplsin = damp * amplsin + normal(zeros(Nmode), kick_amplitude)
    amplcos = damp * amplcos + normal(zeros(Nmode), kick_amplitude)


  #-- Initialize the last kick times for each mode to be randomly chosen
  #-- a little before the first user time point. This is to avoid that
  #-- the kicking time is always exactly the same for all of the modes.

  last_kicktime = uniform(time[0] - kicktimestep, time[0], Nmode)
  next_kicktime = last_kicktime + kicktimestep


  #--Start simulating the time series.

  print >> sys.stderr, "Simulating stochastic oscillations"

  signal = zeros(Ntime)

  for j in range(Ntime):

    # Compute the contribution of each mode separatly

    for i in range(Nmode):

      # Let the oscillator evolve until right before 'time[j]'

      while (next_kicktime[i] <= time[j]):

        deltatime = next_kicktime[i] - last_kicktime[i]
        damp = exp(-eta[i] * deltatime)
        amplsin[i] = damp * amplsin[i] + kick_amplitude[i] * normal(0.,1.)
        amplcos[i] = damp * amplcos[i] + kick_amplitude[i] * normal(0.,1.)
        last_kicktime[i] = next_kicktime[i]
        next_kicktime[i] = next_kicktime[i] + kicktimestep

      # Now make the last small step until 'time[j]'

      deltatime = time[j] - last_kicktime[i]
      damp = exp(-eta[i] * deltatime)
      signal[j] = signal[j] + damp * (amplsin[i] * sin(2*pi*freq[i]*time[j])    \
                                      + amplcos[i] * cos(2*pi*freq[i]*time[j]))

  #-- Return the resulting signal

  return(signal)

