
__all__ = ['eacf', 'meanLargeSeparation']


import numpy as np



def DFTpower(time, signal, freqs):

    """
    Auxiliary function. Not to be imported.

    @param:
    @type:
    @param: 
    @type: 
    @return: 
    @rtype: 
    
    """
    
    powerSpectrum = np.zeros(len(freqs))

    for i, freq in enumerate(freqs):
        arg = 2.0 * np.pi * freq * time
        powerSpectrum[i] = np.sum(signal * np.cos(arg))**2 + np.sum(signal * np.sin(arg))**2

    powerSpectrum = powerSpectrum * 4.0 / len(time)**2
    return(powerSpectrum)
    


    
    

def eacf(freqs, spectrum, spacings, kernelWidth, minFreq=None, maxFreq=None):

    """
    Compute the Envelope AutoCorrelation Function
    
    Source: 
    Example:
    
    >>>
    >>>
    
    @param:
    @type:
    @param: 
    @type: 
    @return: 
    @rtype: 
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
    
    return autoCorrelation  #, croppedFreqs, smoothedSpectrum 
    





def meanLargeSeparation(freq, spec, nuMax, minFreq, maxFreq):

    """
    
    """
    
    firstEstimateLargeSep = 0.24 * nuMax**0.78
    kernelWidth = firstEstimateLargeSep / 8.0
    spacings = np.arange(1.0, 15.0, 0.005)
    autoCorrelation = eacf(freq, spec, spacings, kernelWidth, minFreq, maxFreq)
    freqOfHighestPeak = spacings[np.argmax(autoCorrelation)]
    
    return 2*freqOfHighestPeak
    




def doubleEacf(freqs, spectrum, spacings, kernelWidth1, kernelWidth2, minFreq=None, maxFreq=None):

    """
    
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
     
    kernelSize1 = 1 + 2 * int(round(kernelWidth1/freqStep))
    kernelSize2 = 1 + 2 * int(round(kernelWidth2/freqStep))
    if kernelSize1 <= kernelSize2:
        shift = (kernelSize2 - kernelSize1) / 2
        kernel = np.hanning(kernelSize2)
        kernel[shift:shift+kernelSize1] *= np.hanning(kernelSize1)
    else:
        shift = (kernelSize1 - kernelSize2) / 2
        kernel = np.hanning(kernelSize1)
        kernel[shift:shift+kernelSize2] *= np.hanning(kernelSize2)
        
    kernel = kernel/kernel.sum()
    
    # Smooth the spectrum using a convolution
    
    smoothedSpectrum = np.convolve(kernel, croppedSpectrum, mode='same')
    smoothedSpectrum -= smoothedSpectrum.mean()
    
    # Compute the power spectrum of the power spectrum
    
    autoCorrelation = DFTpower(croppedFreqs, smoothedSpectrum, 1.0/spacings)
    
    # That's it.
    
    return autoCorrelation
    
    
    
