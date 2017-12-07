# -*- coding: utf-8 -*-
"""
Window functions.

A window with a high dynamic range is a window that can distinguish peaks in
a broad range of amplitudes. Usually, this goes together with loss in resolution
and detection power (higher noise).

Example usage: calculate the window function of timeseries, using different
windows.
    >>> import tfreq
    >>> from pylab import figure,show,legend,plot,xlim
    >>> times = linspace(0,150,10000)
    >>> signal = ones(len(times))
    >>> windows = ['rectangular','hann','hamming','cosine','bartlett']
    >>> windows += ['nuttall','blackman-harris']
    >>> for window in windows:
    ...   pergram = tfreq.scargle(times,signal,fn=1.5,df=0.0005,norm='amplitude',window=window)[1]
    ...   p=figure(1)
    ...   p=plot(pergram[0],pergram[1],'-',label=window)
    ...   p=figure(2)
    ...   p=plot(log10(pergram[0]),log10(pergram[1]),'-',label=window)
    >>> p=figure(1);p=legend();p=xlim(0,1.5)
    >>> p=figure(2);p=legend(loc='lower left');p=xlim(-2.5,-0.3)
"""
import numpy as np
import logging

logger = logging.getLogger("IVS.TS.WINF")

#{ Main wrapper
def getWindowFunction(name,times):
    return globals()[name.lower()](times)

#}

#{ Low dynamic range
def rectangular(times):
    """
    Rectangular (no) window
    """
    logger.debug("Selected rectangular window")
    return np.ones(len(times))
#}

#{ Moderate dynamic range
def hamming(times):
    """
    Hamming window
    """
    a = 0.54836
    window_width = times.ptp()
    window = a - (1-a) * np.cos(2*np.pi*(times-times[0])/window_width)
    logger.debug("Selected Hamming window")
    return window

def hann(times):
    """
    Hann or Hanning window
    """
    n = times-times[0]
    N_1 = times.ptp()
    window = 0.5* (1-np.cos(2*np.pi*n)/N_1)
    logger.debug("Selected Hann window")
    return window

def cosine(times):
    """
    Cosine window
    """
    n = times-times[0]
    N_1 = times.ptp()
    window = np.cos(np.pi*n/N_1 - np.pi/2)
    logger.debug("Selected Cosine window")
    return window

def bartlett(times):
    """
    Bartlett window (zero valued end-points)
    """
    n = times-times[0]
    N_1 = times.ptp()
    window = 2./N_1 * (N_1/2. - np.abs(n-N_1/2.))
    logger.debug("Selected Bartlett window")
    return window

#}
#{ High dynamic range
def nuttall(times):
    """
    Nuttall window
    """
    a0=0.355768
    a1=0.487396
    a2=0.144232
    a3=0.012604
    n = times-times[0]
    N_1 = times.ptp()
    window = a0 - a1*np.cos(2*np.pi*n/N_1) + a2*np.cos(4*np.pi*n/N_1) - a3*np.cos(6*np.pi*n/N_1)
    logger.debug("Selected Nuttall window")
    return window

def blackman_harris(times):
    """
    Blackman Harris window
    """
    a0 = 0.3635819
    a1 = 0.4891775
    a2 = 0.1365995
    a3 = 0.0106411
    n = times-times[0]
    N_1 = times.ptp()
    window = a0 - a1*np.cos(2*np.pi*n/N_1) + a2*np.cos(4*pi*n/N_1) - a3*np.cos(6*np.pi*n/N_1)
    logger.debug("Selected Blackman-Harris window")
    return window
#}

def test():
    """
        >>> from pylab import show
        >>> p=show()
    """
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    test()
