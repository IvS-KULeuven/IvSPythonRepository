# -*- coding: utf-8 -*-
"""
Frequency Analysis Routines.

Author: Joris De Ridder
"""
import numpy as np
from ivs.timeseries import fit
from ivs.timeseries import evaluate


#{ Convenience functions

def find_frequency(times,signal,
            f0=None,fn=None,df=None,threads=1,
            max_loops=20, scale_region=0.1, scale_df=0.20, nyq_stat=np.min,
            **kwargs):
    """
    Find one frequency, automatically going to maximum precision and parameters
    & error estimates.
    
    This routine will make the frequency grid finer until it is well below the
    estimated error on the frequency. After that, it will compute harmonic
    parameters and estimated errors.
    
    There is a possibility to escape this optimization by setting dfscale=0 or
    freqregscale=0.
    
    Possible Extra keywords: see definition of the used periodogram function,
    evaluating and fitting functions, etc...
    
    Example keywords:
        - 'correlation_correction', default=True
        - 'freqregscale', default=0.5: factor for zooming in on frequency
        - 'dfscale', default = 0.25: factor for optimizing frequency resolution
    
    There is a possibility to enter frequency search range in units of the
    Nyquist frequency, by setting the keyword 'units' to 'relative'.
    
    There is a possibiliity to include a nonlinear least square update of the
    parameters, by setting the keyword 'nllsq' to True.
    
    Example usage:
    
    Import necessary modules:
        >>> from pylab import plot,figure,title
        >>> from sigproc.base import stats
    
    Generate test data:
        >>> times = linspace(0,150,10000)
        >>> pars = array([[0,1,1/20.,0],[0,0.3,1/6.,0]])
        >>> signal = teval.sine_eval(times,pars) + random.normal(size=len(times),scale=3.5)
    
    Search for the first frequency (linear and nonlinear):
        >>> o1 = find_frequency(times,signal,norm='amplitude',fn=1)
        >>> o2 = find_frequency(times,signal,norm='amplitude',fn=1,nllsq=True)
    
    Output the results:
        >>> p=figure();p=title("test:FIND_FREQUENCY: fit result")
        >>> p=plot(times,signal,'ko')
        >>> p=plot(times,teval.sine_eval(times,o1['pars']),'r-',linewidth=2)
        >>> p=plot(times,teval.sine_eval(times,o2['pars']),'b--',linewidth=2)
        >>> p=figure();p=title("test:FIND_FREQUENCY: periodogram result")
        >>> p=plot(o1['pergram'][0],o1['pergram'][1],'k-')
        >>> p=plot(o1['pergram_fine'][0],o1['pergram_fine'][1],'r-')
        >>> varred_linear = stats.varred(signal,teval.sine_eval(times,o1['pars']))
        >>> varred_nonlinear = stats.varred(signal,teval.sine_eval(times,o2['pars']))
        >>> print varred_linear < varred_nonlinear
        True
    
    @rtype: dict
    @return: [frequency, peak], full_periodogram, fine_periodogram, stats,
             parameters, errors, rho
    """
    #-- make sure the series is sorted
    gaps = np.diff(times)

    if np.any(gaps<=0):
        logger.warning("Time series not sorted or equal timepoints: sorting now")
        sa = np.argsort(times)
        times = times[sa]
        signal = signal[sa]
        
    T = times[-1]-times[0]
    nyquist = getNyquist(times,nyq_stat=nyq_stat)
    N = len(times)
    
    #-- basic input: start/stop frequency, and frequency step
    if f0 is None: f0 = 0.1/T
    if fn is None: fn = nyquist
    if df is None: df = 0.1/T
    
    #-- initial values
    e_f = 0
    freq_diff = inf
    prev_freq = -inf
    counter = 0
    
    #-- calculate periodogram until frequency precision is
    #   under 1/10th of correlation corrected version
    while freq_diff>e_f/10.:
        #-- calculate periodogram
        freqs,ampls = scargle(times,signal,f0=f0,fn=fn,df=df,**kwargs)
        frequency = freqs[np.argmax(ampls)]
        #-- estimate parameters and calculate a fit, errors and residuals
        params = fit.sine(times,signal,frequency,**kwargs)
        errors = fit.e_sine(times,signal,params,**kwargs)
        #-- improve precision
        freq_diff = abs(frequency-prev_freq)
        prev_freq = frequency
        freq_region = fn-f0
        f0 = max(0      ,frequency-freq_region*freqregscale/2.)
        fn = min(nyquist,frequency+freq_region*freqregscale/2.)
        df *= dfscale
        #-- possibility to escape optimization
        if freqregscale==0 or dfscale==0: break
        if counter >= max_loops:
            logger.error("Frequency precision not reached in %d steps, breaking loop"%(max_loops))
        counter += 1

    return frequency, parameters


















def Zwavelet(time, signal, freq, position, sigma=10.0):

    """
    Computes "Weighted Wavelet Z-Transform"
    which is a type of time-frequency diagram more suitable
    for non-equidistant time series.
    See: G. Foster, 1996, Astronomical Journal, 112, 1709 
    
    The result can be plotted as a colorimage (imshow in pylab).
    
    Mind the max, min order for position. It's often useful to try log(Z)
    and/or different sigmas.
    
    @param time: time points [0..Ntime-1]
    @type time: ndarray
    @param signal: observed data points [0..Ntime-1]
    @type signal: ndarray
    @param freq: frequencies (omega/2pi in Foster's paper) [0..Nfreq-1]
                 the array should not contain 0.0
    @type freq: ndarray
    @param position: time parameter: tau in Foster's paper [0..Npos-1]
    @type position: ndarray
    @param sigma: smoothing parameter in time domain: sigma in Foster's paper
    @type sigma: float
    @return: Z[0..Npos-1, 0..Nfreq-1]: the Z-transform: time-freq diagram
    @rtype: ndarray
    
    """
  
    y = np.zeros(3)
    S = np.zeros([3,3])
    Z = np.zeros([len(position),len(freq)])
  
    for i in range(len(position)):
        
        tau = position[i]
        arg = 2.0*pi*(time-tau)
 
        for j in range(len(freq)):
        
            nu = freq[j]
      
            # Compute statistical weights akin the Morlet wavelet
      
            weight = np.exp(-(time-tau)**2 * (nu / 2./sigma)**2)
            W = np.sum(weight)
      
            # Compute the base functions. A 3rd base function is the constant 1.
      
            cosine = np.cos(arg*nu)
            sine = np.sin(arg*nu)
          
            print arg, nu
            print sine, cosine
          
            # Compute the innerproduct of the base functions
            # phi_0 = 1 (constant), phi_1 = cosine, phi_2 = sine
      
            S[0,0] = 1.0
            S[0,1] = S[1,0] = np.sum(weight * 1.0 * cosine) / W
            S[0,2] = S[2,0] = np.sum(weight * 1.0 * sine) / W
            S[1,1] = np.sum(weight * cosine * cosine) / W
            S[1,2] = S[2,1] = np.sum(weight * cosine * sine) / W
            S[2,2] = np.sum(weight * sine * sine) / W
            print S
            invS = np.linalg.inv(S)
      
            # Determine the best-fit coefficients y_k of the base functions
      
            for k in range(3):
                y[k] =   np.sum(weight * 1.0 * signal) / W * invS[k,0]     \
                       + np.sum(weight * cosine * signal) / W * invS[k,1]  \
                       + np.sum(weight * sine * signal) / W * invS[k,2]      
      
            # Compute the best-fit model
      
            model = y[0] + y[1] * cosine + y[2] * sine
      
            # Compute the weighted variation of the signal and the model functions
      
            Vsignal = np.sum(weight * signal**2) / W -  (np.sum(weight * signal) / W)**2
            Vmodel = np.sum(weight * model**2) / W -  (np.sum(weight * model) / W)**2
      
            # Calculate the weighted Wavelet Z-Transform
      
            Neff = W**2 / np.sum(weight**2)
            Z[i,j] = (Neff - 3) * Vmodel / 2. / (Vsignal - Vmodel)
      
    # That's it!
  
    return Z
    
  
  



  



def getNyquist(times,nyq_stat=np.inf):
    """
    Calculate Nyquist frequency.
    
    Typical use is minimum or median of time points differences.
    
    If C{nyq_stat} is not callable, it is assumed to be a number and that number
    will just be returned: this you can do to search for frequencies above the
    nyquist frequency
    
    @param times: sorted array containing time points
    @type times: numpy array
    @param nyq_stat: statistic to use or absolute value of the Nyquist frequency
    @type nyq_stat: callable or float
    @return: Nyquist frequency
    @rtype: float
    """
    if not hasattr(nyq_stat,'__call__'):
        nyquist = nyq_stat
    else:
        nyquist = 1/(2.*nyq_stat(np.diff(times)))
    return nyquist
    