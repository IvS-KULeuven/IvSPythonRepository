# -*- coding: utf-8 -*-
"""
Frequency Analysis Routines.

Author: Joris De Ridder
"""


import numpy as np
from math import ceil, fmod, pi
  



def DFTpower(time, signal, startfreq, stopfreq, stepfreq):

    """
    Computes the modulus square of the fourier transform. 
    Unit: square of the unit of signal. Time points need not be equidistant.
    The normalisation is such that a signal A*sin(2*pi*nu_0*t)
    gives power A^2 at nu=nu_0
    
    @param time: time points [0..Ntime-1] 
    @type time: ndarray
    @param signal: signal [0..Ntime-1]
    @type signal: ndarray
    @param startfreq: the power is computed for the frequencies
                      freq = arange(startfreq,stopfreq,stepfreq)
    @type startfreq: float
    @param stopfreq: see startfreq
    @type stopfreq: float
    @param stepfreq: see startfreq
    @type stepfreq: float
    @return: power spectrum of the signal
    @rtype: ndarray 
    
    """
  
    Ntime = len(time)
    Nfreq = int(ceil((stopfreq-startfreq)/stepfreq))
  
    A = np.exp(1j*2.*pi*startfreq*time) * signal
    B = np.exp(1j*2.*pi*stepfreq*time)
    ft = np.zeros(Nfreq, complex) 
    ft[0] = A.sum()

    for k in range(1,Nfreq):
        A *= B
        ft[k] = np.sum(A)
    
    return (ft.real**2 + ft.imag**2) * 4.0 / Ntime**2
  
  
  
  
  
  
  
def FFTpower(signal, timestep):

    """
    Computes power spectrum of an equidistant time series 'signal'
    using the FFT algorithm. The length of the time series need not
    be a power of 2 (zero padding is done automatically). 
    Normalisation is such that a signal A*sin(2*pi*nu_0*t)
    gives power A^2 at nu=nu_0  (IF nu_0 is in the 'freq' array)
    
    @param signal: the time series [0..Ntime-1]
    @type signal: ndarray
    @param timestep: time step fo the equidistant time series
    @type timestep: float
    @return: frequencies and the power spectrum
    @rtype: (ndarray,ndarray)
    
    """
  
    # Compute the FFT of a real-valued signal. If N is the number 
    # of points of the original signal, 'Nfreq' is (N/2+1).
  
    fourier = np.fft.rfft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)
  
    # Compute the power
  
    power = np.abs(fourier)**2 * 4.0 / Ntime**2
  
    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep
  
    freq = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep
  
    # That's it!
  
    return (freq, power)





  
  
  
def FFTpowerdensity(signal, timestep):
  
    """
    Computes the power density of an equidistant time series 'signal',
    using the FFT algorithm. The length of the time series need not
    be a power of 2 (zero padding is done automatically). 

    @param signal: the time series [0..Ntime-1]
    @type signal: ndarray
    @param timestep: time step fo the equidistant time series
    @type timestep: float
    @return: frequencies and the power density spectrum
    @rtype: (ndarray,ndarray)

    """
  
    # Compute the FFT of a real-valued signal. If N is the number 
    # of points of the original signal, 'Nfreq' is (N/2+1).
  
    fourier = np.fft.rfft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)
  
    # Compute the power density
  
    powerdensity = np.abs(fourier)**2 / Ntime * timestep
  
    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep
  
    freq = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep
  
    # That's it!
  
    return (freq, powerdensity)

  


  
  



def weightedpower(time, signal, weight, freq):

    """
    Compute the weighted power spectrum of a time signal.
    For each given frequency a weighted sine fit is done using
    chi-square minimization.
    
    @param time: time points [0..Ntime-1] 
    @type time: ndarray
    @param signal: observations [0..Ntime-1]
    @type signal: ndarray
    @param weight: 1/sigma_i^2 of observation
    @type weight: ndarray
    @param freq: frequencies [0..Nfreq-1] for which the power 
                 needs to be computed
    @type freq: ndarray
    @return: weighted power [0..Nfreq-1]
    @rtype: ndarray
    
    """

    result = np.zeros(len(freq))

    for i in range(len(freq)):
        if (freq[i] != 0.0):
            sine   = np.sin(2.0*pi*freq[i]*time)
            cosine = np.cos(2.0*pi*freq[i]*time)
            a11= np.sum(weight*sine*sine)
            a12 = np.sum(weight*sine*cosine)
            a21 = a12
            a22 = np.sum(weight*cosine*cosine)
            b1 = np.sum(weight*signal*sine)
            b2 = np.sum(weight*signal*cosine)
            denominator = a11*a22-a12*a21
            A = (b1*a22-b2*a12)/denominator
            B = (b2*a11-b1*a21)/denominator
            result[i] = A*A+B*B
        else:
            result[i] = np.sum(signal)/len(signal)

    return(result)









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
  
  
  
  
           
           
           
           
           
           
def windowfunction(time, freq):

    """
    Computes the modulus square of the window function of a set of 
    time points at the given frequencies. The time point need not be 
    equidistant. The normalisation is such that 1.0 is returned at 
    frequency 0.
    
    @param time: time points  [0..Ntime-1]
    @type time: ndarray       
    @param freq: frequency points. Units: inverse unit of 'time' [0..Nfreq-1]
    @type freq: ndarray       
    @return: |W(freq)|^2      [0..Nfreq-1]
    @rtype: ndarray
    
    """
  
    Ntime = len(time)
    Nfreq = len(freq)
    winkernel = np.empty_like(freq)

    for i in range(Nfreq):
        winkernel[i] = np.sum(np.cos(2.0*pi*freq[i]*time))**2     \
                     + np.sum(np.sin(2.0*pi*freq[i]*time))**2

    # Normalise such that winkernel(nu = 0.0) = 1.0 

    return winkernel/Ntime**2








def pdm(time, signal, freq, Nbin=10, Ncover=5):

    """
    Computes the theta-statistics to do a Phase Dispersion Minimisation.
    See Stellingwerf R.F., 1978, ApJ, 224, 953)
    
    @param time: time points  [0..Ntime-1]
    @type time: ndarray       
    @param signal: observed data points [0..Ntime-1]
    @type signal: ndarray
    @param freq: frequency points. Units: inverse unit of 'time' [0..Nfreq-1]
    @type freq: ndarray
    @param Nbin: the number of phase bins (with length 1/Nbin)
    @type Nbin: integer
    @param Ncover: the number of covers (i.e. bin shifts)
    @type Ncover: integer
    @return: theta-statistic for each given frequency [0..Nfreq-1]
    @rtype: ndarray

    """
  
    Ntime = len(time)
    Nfreq = len(freq)
  
    binsize = 1.0 / Nbin
    covershift = 1.0 / (Nbin * Ncover)
  
    theta = np.zeros(Nfreq)
  
    for i in range(Nfreq):
  
        # Compute the phases in [0,1[ for all time points
    
        phase = np.fmod((time - time[0]) * freq[i], 1.0)
    
        # Reset the number of (shifted) bins without datapoints
    
        Nempty = 0
    
        # Loop over all Nbin * Ncover (shifted) bins
    
        for k in range(Nbin):
            for n in range(Ncover):
        
                # Determine the left and right boundary of one such bin
                # Note that due to the modulo, right may be < left. So instead
                # of 0-----+++++------1, the bin might be 0++-----------+++1 .
        
                left = fmod(k * binsize + n * covershift, 1.0) 
                right = fmod((k+1) * binsize + n * covershift, 1.0) 

                # Select all data points in that bin
        
                if (left < right):
                    bindata = np.compress((left <= phase) & (phase < right), signal)
                else:
                    bindata = np.compress(~((right <= phase) & (phase < left)), signal)

                # Compute the contribution of that bin to the theta-statistics  
          
                if (len(bindata) != 0):
                    theta[i] += (len(bindata) - 1) * bindata.var()
                else:
                    Nempty += 1
  
        # Normalize the theta-statistics        

        theta[i] /= Ncover * Ntime - (Ncover * Nbin - Nempty)     
    
    # Normalize the theta-statistics again
  
    theta /= signal.var()  
    
    # That's it!
 
    return theta
    
  
  



  
            
def phasediagram(time, signal, nu0, D=0, t0=None, return_indices=False,
                 chronological=False):
    """
    Construct a phasediagram, using frequency nu0.
    
    Possibility to include a frequency shift D.
    
    If needed, the indices can be returned. To 'unphase', just do:
    
    Example usage:
    
    >>> original = random.normal(size=10)
    >>> indices = argsort(original)
    >>> inv_indices = argsort(indices)
    >>> all(original == take(take(original,indices),inv_indices))
    True
    
    Optionally, the zero time point can be given, it will be subtracted from all
    time points
    
    Joris De Ridder, Pieter Degroote
    
    @param time: time points [0..Ntime-1]
    @type time: ndarray
    @param signal: observed data points [0..Ntime-1]
    @type signal: ndarray
    @param nu0: frequency to compute the phasediagram with 
                (inverse unit of time)
    @type nu0: float
    @param t0: zero time point, if None: t0 = time[0]
    @type t0: float
    @param D: frequency shift
    @type D: float
    @param chronological: set to True if you want to have a list of every phase
    separately
    @type chronological: boolean
    @param return_indices: flag to return indices
    @type return_indices: boolean
    @return: phase points (sorted), corresponding observations
    @rtype: ndarray,ndarray(, ndarray)
    """
    if (t0 is None): t0 = time[0]
    phase = np.fmod(nu0 * (time-t0) + D/2. * (time-t0)**2,1.0)
    phase = np.where(phase<0,phase+1,phase)
    
    if chronological:
        #-- keep track of the phase number, and prepare lists to add the points
        chainnr = np.floor((time-t0)*nu0)
        myphase = [[]]
        mysignl = [[]]
        currentphase = chainnr[0]
        #-- divide the signal into separate phase bins
        for i in xrange(len(signal)):
            if chainnr[i]==currentphase:
                myphase[-1].append(phase[i])
                mysignl[-1].append(signal[i])
            else:
                myphase[-1] = np.array(myphase[-1])
                mysignl[-1] = np.array(mysignl[-1])
                currentphase = chainnr[i]
                myphase.append([phase[i]])
                mysignl.append([signal[i]])
        myphase[-1] = np.array(myphase[-1])
        mysignl[-1] = np.array(mysignl[-1])
    else:
        indices = np.argsort(phase)
    
    #-- possibly only return phase and signal
    if not return_indices and not chronological:
        return phase[indices], signal[indices]
    #-- maybe return phase,signal and indices
    elif not chronological:
        return phase[indices], signal[indices], indices
    #-- or return seperate phase bins
    else:
        return myphase,mysignl