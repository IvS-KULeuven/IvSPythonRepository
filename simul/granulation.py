"""
Simulation of a granulation signal

Author: Joris De Ridder

"""

import numpy as np
from numpy.random import normal



def granulation(time, timescale, varscale, logger=None):

    """
    Simulates a time series showing granulation variations
    
    A first-order autoregressive process is used, as this gives a Harvey
    model in the frequency domain. See also:
    De Ridder et al., 2006, MNRAS 365, pp. 595-605.
    
    @param time: time points
    @type time: ndarray
    @param timescale: array of time scale "tau_i" of each granulation component
                      of the granulation/magnetic activity. Same units as 'time'.
    @type timescale: ndarray
    @param varscale: array of variation scale "sigma_i" of each component of the 
                     granulation/magnetic activity in the appropriate passband.
                     Same size as the timescale array. Unit: ppm
    @type varscale: ndarray
    @param logger: logger (optional). logger.info() will be called to write info strings.
    @type logger: logging instance 
    @return: ndarray containing the granulation signal
    
    Example:
    
    >>> time = np.linspace(0,100,200)             # E.g. in days
    >>> timescale = np.array([5.0, 20.])          # time scales in days
    >>> varscale = np.array([10.0, 50.0])         # variation scale in ppm
    >>> gransignal = granulation(time, timescale, varscale)
    >>> flux = 100000.0                           # mean flux level
    >>> signal = flux * (1.0 + gransignal)        # signal in flux
    
    """
    
    
    Ntime = len(time)
    Ncomp = len(timescale)

    if logger is not None:
        logger.info("Simulating %d granulation components\n" % Ncomp)
        
    # Set the kick (= reexcitation) timestep to be one 100th of the
    # shortest granulation time scale (i.e. kick often enough).

    kicktimestep = min(timescale) / 100.0
    
    if logger is not None:
        logger.info("Kicktimestep = %f\n" % kicktimestep)

    # Predefine some arrays

    signal = np.zeros(Ntime)
    granul = np.zeros(Ncomp)
    mu = np.zeros(Ncomp)
    sigma = np.sqrt(kicktimestep/timescale)*varscale

    # Warm up the first-order autoregressive process

    if logger is not None:
        logger.info("Granulation process warming up...\n")
        
    for i in range(2000):
        granul = granul * (1.0 - kicktimestep / timescale) + normal(mu, sigma)

    # Start simulating the granulation time series

    if logger is not None:
        logger.info("Simulating granulation signal.\n")

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
                 + normal(mu, np.sqrt(delta/timescale)*varscale)
        currenttime = time[i]

        # Add the different components to the signal. 

        signal[i] = sum(granul)


    # That's it!

    return(signal)

