# -*- coding: utf-8 -*-
"""
Computes the envelope auto-correlation function. This function is estimated by computing the power spectrum of 
a smoothed power spectrum of a time series. It is ideally suited to discover periodicities in the power spectrum,
as is often the case for solar-like oscillators (cf the asymptotic relation of Tassoul), and can thus be used to 
estimate the mean large separation (mean frequency spacing between two consecutive radial modes).

As an example we take the Kepler red giant KIC3744043, which shows a beautiful spectrum of solar-like oscillations:

]include figure]]powerspec_KIC3744043.png]

To compute the EACF, we first import:

>>> import numpy as np
>>> from ivs.timeseries.eacf import eacf

Red giants have spacings in the power spectrum between, say, 1.0 and 15.0, so this is the interval in which we will
compute the EACF:

>>> spacings = np.arange(1.0, 15.0, 0.005)

The figure above shows that the oscillations are only present between (roughly) 80 microHz and 135 microHz, so we will limit ourselves to this frequency interval:

>>> minFreq, maxFreq = 80.0, 135.0

The EACF is now computed as:

>>> autoCorrelation, croppedFreqs, smoothedSpectrum = eacf(freq, spec, spacings, 2.0, minFreq, maxFreq)

The value '2.0' is related to the FWHM of the (Hann) kernel that is used to smooth the power spectrum of the signal, on which we will come back later.

The EACF looks as follows:

]include figure]]eacf_smooth15.png]

The highest peak is caused by the main periodicity in the power spectrum. This does not, however, correspond directly 
to the mean large separation, but to half the large separation. The reason is that almost in the middle between each 
two radial modes there is a dipole mode. To compute the large separation we therefore have to enter:

>>> meanLargeSeparation = 2.0 * spacings[np.argmax(autoCorrelation)]
>>> meanLargeSeparation
9.8599999999998325

How did we choose the value 2.0 (microHz) for the FWHM of the smoothing kernel? The smoothing is done by convolving
the power spectrum with a Hann kernel with a certain width. If you take the width too large, the spectrum will be too
heavily smoothed so that no periodicities can be detected. If you take the width too small, the EACF will contain too 
many features, which makes it difficult to recognize the right peak. As a rule-of-thumb you can take the kernel width
to be 1/8 of your initial estimate of the mean large separation. If you don't have such an initial estimate, you can
estimate nu_max, which is the frequency of maximum power (= the location of the gaussian envelope of the 
power excess), and derive from that value a first estimate for the large separation:

>>> nuMax = 120.0                                      # for KIC3744043
>>> initialEstimateLargeSep = 0.24 * nuMax**0.78       # general approximation for red giants
>>> kernelWidth = initialEstimateLargeSep / 8.0        # rule of thumb
>>> autoCorrelation, croppedFreqs, smoothedSpectrum = eacf(freq, spec, spacings, kernelWidth, minFreq, maxFreq)

"""


__all__ = ['eacf']


import numpy as np
from ivs.timeseries.pergrams import DFTpower2 as DFTpower

def eacf(freqs, spectrum, spacings, kernelWidth, minFreq=None, maxFreq=None, doSanityCheck=False):

    """
    Compute the Envelope Auto-Correlation Function (EACF) of a signal.

    The EACF is derived by computing the power spectrum of a smoothed version of the 
    power spectrum of the time series. It allows to identify equispaced patterns in the 
    power spectrum

    Source: 
        - Mosser & Appourchaux, 2009, A&A 508, p. 877
        - Mosser, 2010, Astron. Nachr. 331, p. 944
    
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
    @param doSanityCheck: if True: make a few sanity checks of the arguments (e.g. equidistancy of freqs).
                          If False, do nothing.
    @type doSanityCheck: boolean
    @return: autoCorrelation, croppedFreqs, smoothedSpectrum
                - autoCorrelation:  the EACF evaluated in the values of 'spacings'
                - croppedFreqs:     the frequencies of the selected part [minFreq, maxFreq]
                - smoothedSpectrum: the smoothed power spectrum used to compute the EACF. Same length as croppedFreqs.    
    @rtype: (ndarray, ndarray, ndarray)
    """
    
    # If requested, perform sanity checks
    # The 1000. in the check on equidistancy was needed to let pass the arrays created by 
    # linspace() and arange().

    if doSanityCheck:
        if len(freqs) != len(spectrum):
            raise ValueError("freqs and spectrum don't have the same length")
        if np.fabs(np.diff(np.diff(freqs)).max()) > 1000. * np.finfo(freqs.dtype).eps:
            raise ValueError("freqs array is not equidistant")
        if np.sometrue(spacings <= 0.0):
            raise ValueError("spacings are not all strictly positive")
        if kernelWidth <= 0.0:
            raise ValueError("kernel width is not > 0.0")
        if minFreq > freqs[-1]:
            raise ValueError("minimum frequency 'minFreq' is larger than largest frequency in freqs array")
        if maxFreq < freqs[0]:
            raise ValueError("maximum frequency 'maxFreq' is smaller than smallest frequency in freqs array")
        if minFreq >= maxFreq:
            raise ValueError("minFreq >= maxFreq")


    # Set the default values
    
    if minFreq == None: minFreq = freqs[0]
    if maxFreq == None: maxFreq = freqs[-1]
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
    

