# -*- coding: utf-8 -*-
"""
Functions to compute the mean large separation using the envelope auto-correlation function.

"""


__all__ = ['eacf', 'meanLargeSeparation']


import numpy as np



def DFTpower(time, signal, freqs):

    """
    Computes the power spectrum of a signal using a discrete Fourier transform.

    Auxiliary function. Not to be imported. It is re-implemented here because
    the standard power spectrum functions in the repository only allow for equidistant frequencies.

    @param time: time points, not necessarily equidistant
    @type time: ndarray
    @param signal: signal corresponding to the given time points
    @type signal: ndarray
    @param freqs: frequencies for which the power spectrum will be computed. Unit: inverse of 'time'.
    @type freqs: ndarray
    @return: power spectrum. Unit: square of unit of 'signal'
    @rtype: ndarray
    """
    
    powerSpectrum = np.zeros(len(freqs))

    for i, freq in enumerate(freqs):
        arg = 2.0 * np.pi * freq * time
        powerSpectrum[i] = np.sum(signal * np.cos(arg))**2 + np.sum(signal * np.sin(arg))**2

    powerSpectrum = powerSpectrum * 4.0 / len(time)**2
    return(powerSpectrum)
    


    
    

def eacf(freqs, spectrum, spacings, kernelWidth, minFreq=None, maxFreq=None):

    """
    Compute the Envelope Auto-Correlation Function (EACF) of a signal.

    The EACF is derived by computing the power spectrum of a smoothed version of the 
    power spectrum of the time series. It allows to identify equispaced patterns in the 
    power spectrum

    Source: 
        - Mosser & Appourchaux, 2009, A&A 508, p. 877
        - Mosser, 2010, Astron. Nachr. 331, p. 944
    
    Example:
    
    >>>
    >>>
    
    @param freqs: frequencies in which the power specturm is given. Assumed to be equidistant.
    @type freqs: ndarray
    @param spectrum: power spectrum corresponding to the given frequency points.
    @type spectrum: ndarray
    @param spacings: array of frequency spacings for which the EACF needs to be computed
    @type spacings: ndarray
    @param kernelWidth: the width (in freq units) of the convolution kernel used to first smooth the power spectrum
    @type kernelWidth: float
    @param minFreq: only the part [minFreq, maxFreq] of the spectrum is taken into account.
                    If it is None, the 1st point of the 'freqs' array is taken.
    @type minFreq: float
    @param maxFreq: only the part [minFreq, maxFreq] of the spectrum is taken into account
                    If it is None, the last point of the 'freqs' array is taken.
    @type maxFreq: float
    @return: autoCorrelation, croppedFreqs, smoothedSpectrum
                - autoCorrelation:  the EACF evaluated in the values of 'spacings'
                - croppedFreqs:     the frequencies of the selected part [minFreq, maxFreq]
                - smoothedSpectrum: the smoothed power spectrum used to compute the EACF. Same length as croppedFreqs.    
    @rtype: (ndarray, ndarray, ndarray)
    """
    
    # Set the default values
    
    if minFreq == None: minFreq = freq[0]
    if maxFreq == None: maxFreq = freq[-1]
    freqStep = freqs[1]-freqs[0]
   
    # Crop the spectrum to the specified range
    
    croppedSpectrum = spectrum[(freqs >= minFreq) & (freqs <= maxFreq)]
    croppedFreqs = freqs[(freqs >= minFreq) & (freqs <= maxFreq)]
    
    # Set up a normalized Hann smoothing kernel.
    # Note: a Hann window of size N has FWHM = (N-1)/2. The FWHM is given by 
    # the user in frequency units, so the corresponding N can be derived.
     
    kernelSize = 1 + 2 * int(round(kernelWidth/freqStep))
    kernel = np.hanning(kernelSize)
    kernel = kernel/kernel.sum()
    
    # Smooth the spectrum using a convolution
    
    smoothedSpectrum = np.convolve(kernel, croppedSpectrum, mode='same')
    smoothedSpectrum -= smoothedSpectrum.mean()
    
    # Compute the power spectrum of the power spectrum
    
    autoCorrelation = DFTpower(croppedFreqs, smoothedSpectrum, 1.0/spacings)
    
    # That's it.
    
    return autoCorrelation, croppedFreqs, smoothedSpectrum 
    










def meanLargeSeparation(freq, spec, nuMax, minFreq=None, maxFreq=None):

    """
    Estimate the mean large separation from oscillation peaks in a power spectrum.

    @param freq:  frequencies in which the power specturm is given. Assumed to be equidistant.
    @type freq: ndarray
    @param spec:  power spectrum corresponding to the given frequency points.
    @type spec: ndarray
    @param nuMax: Estimate of the freq of maximum power of the oscillation power excess. Unit: same as 'freqs'.
    @type nuMax: float
    @param minFreq: only the part [minFreq, maxFreq] of the spectrum is taken into account.
                    If it is None, the 1st point of the 'freqs' array is taken.
    @type minFreq: float
    @param maxFreq: only the part [minFreq, maxFreq] of the spectrum is taken into account
                    If it is None, the last point of the 'freqs' array is taken.
    @type maxFreq: float
    @return: mean large separation. Unit: same as 'freqs'.
    @rtype: float
    """
    
    firstEstimateLargeSep = 0.24 * nuMax**0.78
    kernelWidth = firstEstimateLargeSep / 8.0
    spacings = np.arange(1.0, 15.0, 0.005)
    autoCorrelation = eacf(freq, spec, spacings, kernelWidth, minFreq, maxFreq)
    freqOfHighestPeak = spacings[np.argmax(autoCorrelation)]
    
    return 2*freqOfHighestPeak
    

