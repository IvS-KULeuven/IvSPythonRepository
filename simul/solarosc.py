"""
Simulation of a solarlike oscillations

Author: Joris De Ridder

"""

import numpy as np
from numpy.random import uniform, normal
from math import sin,cos, floor, pi


def solarosc(time, freq, ampl, eta, logger=None):
    
    """
    Compute time series of stochastically excited damped modes
    
    See also De Ridder et al., 2006, MNRAS 365, pp. 595-605.
    
    @param time: time points [0..Ntime-1] (unit: e.g. Ms)
    @type time: ndarray
    @param freq: oscillation freqs [0..Nmodes-1] (unit: e.g. microHz)
    @type freq: ndarray
    @param ampl: amplitude of each oscillation mode
                 rms amplitude = ampl / sqrt(2.)
    @type ampl: ndarray
    @param eta: damping rates (unit: e.g. (Ms)^{-1})
    @type eta: ndarray
    @param logger: standard python logging facility
    @type logger: logging instance
    @return: ndarray, signal[0..Ntime-1]

    Example:
    
    >>> time = np.linspace(0, 40, 100)      # in Ms
    >>> freq = np.array([23.0, 23.5])       # in microHz
    >>> ampl = np.array([100.0, 110.0])     # in ppm
    >>> eta = np.array([1.e-6, 3.e-6])      # in 1/Ms
    >>> oscsignal = solarosc(time, freq, ampl, eta)
    >>> flux = 1000000.0                      # average flux level
    >>> signal = flux * (1.0 + oscsignal)
    >>> # The same with a logger
    >>> import sys, logging, logging.handlers
    >>> myLogger = logging.getLogger('myLogger')
    >>> myLogger.addHandler(logging.StreamHandler(sys.stdout))
    >>> myLogger.setLevel(logging.INFO)
    >>> oscsignal = solarosc(time, freq, ampl, eta, myLogger)
    Simulating 2 modes
    Oscillation kicktimestep: 3333.333333
    300 kicks for warm up for oscillation signal
    Simulating stochastic oscillations
    
    """
  
    Ntime = len(time)
    Nmode = len(freq)

    if logger is not None:
        logger.info("Simulating %d modes" % Nmode)

    # Set the kick (= reexcitation) timestep to be one 100th of the
    # shortest damping time. (i.e. kick often enough).

    kicktimestep = (1.0 / max(eta)) / 100.0
    
    if logger is not None:
        logger.info("Oscillation kicktimestep: %f" % kicktimestep)

    # Init start values of amplitudes, and the kicking amplitude
    # so that the amplitude of the oscillator will be on average be
    # constant and equal to the user given amplitude

    amplcos = 0.0
    amplsin = 0.0
    kick_amplitude = ampl * np.sqrt(kicktimestep * eta)

    # Warm up the stochastic excitation simulator to forget the
    # initial conditions. Do this during the longest damping time.
    # But put a maximum on the number of kicks, as there might
    # be almost-stable modes with damping time = infinity
  
    damp = np.exp(-eta * kicktimestep)
    Nwarmup = min(20000, int(floor(1.0 / min(eta) / kicktimestep)))

    if logger is not None:
        logger.info("%d kicks for warm up for oscillation signal" % Nwarmup)
        
    for i in range(Nwarmup):
        amplsin = damp * amplsin + normal(np.zeros(Nmode), kick_amplitude)
        amplcos = damp * amplcos + normal(np.zeros(Nmode), kick_amplitude)


    # Initialize the last kick times for each mode to be randomly chosen
    # a little before the first user time point. This is to avoid that
    # the kicking time is always exactly the same for all of the modes.

    last_kicktime = uniform(time[0] - kicktimestep, time[0], Nmode)
    next_kicktime = last_kicktime + kicktimestep


    # Start simulating the time series.

    if logger is not None:
        logger.info("Simulating stochastic oscillations")
        
    signal = np.zeros(Ntime)

    for j in range(Ntime):

        # Compute the contribution of each mode separatly

        for i in range(Nmode):

            # Let the oscillator evolve until right before 'time[j]'

            while (next_kicktime[i] <= time[j]):

                deltatime = next_kicktime[i] - last_kicktime[i]
                damp = np.exp(-eta[i] * deltatime)
                amplsin[i] = damp * amplsin[i] + kick_amplitude[i] * normal(0.,1.)
                amplcos[i] = damp * amplcos[i] + kick_amplitude[i] * normal(0.,1.)
                last_kicktime[i] = next_kicktime[i]
                next_kicktime[i] = next_kicktime[i] + kicktimestep

            # Now make the last small step until 'time[j]'

            deltatime = time[j] - last_kicktime[i]
            damp = np.exp(-eta[i] * deltatime)
            signal[j] = signal[j] + damp * (amplsin[i] * sin(2*pi*freq[i]*time[j])    \
                                  + amplcos[i] * cos(2*pi*freq[i]*time[j]))

    # Return the resulting signal

    return(signal)

