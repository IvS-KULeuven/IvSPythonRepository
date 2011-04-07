# -*- coding: utf-8 -*-
"""
Functions relevant for photometric calibration

Example usage:
    
>>> from pylab import plot,show
>>> for band in ['J','H','KS']:
...    p = plot(*get_response('2MASS.%s'%(band)))
>>> p = show()
"""
import os
import numpy as np

from ivs.misc.decorators import memoized
from ivs.io import ascii

@memoized
def get_response(photband):
    """
    Retrieve the response curve of a photometric system 'SYSTEM.FILTER'
    
    Example usage:
    
    >>> from pylab import plot,show
    >>> for band in ['J','H','KS']:
    ...    p = plot(*get_response('2MASS.%s'%(band)))
    >>> p = show()
    
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @return: (wavelength [A], response)
    @rtype: (array, array)
    """
    photfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),'calibration',photband.upper())
    wave, response = ascii.read2array(photfile).T
    return wave,response

def eff_wave(photband):
    """
    Return the effective wavelength of a photometric passband.
    
    The effective wavelength is defined as the average wavelength weighed wit
    the response curve.
    
    >>> eff_wave('2MASS.J')
    12412.136241640892
    
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @return: effective wavelength [A]
    @rtype: float
    """
    wave,response = get_response(photband)
    return np.average(wave,weights=response)

if __name__=="__main__":
    import doctest
    doctest.testmod()