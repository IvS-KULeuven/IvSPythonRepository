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
import glob
import numpy as np

from ivs.misc.decorators import memoized
from ivs.io import ascii
from ivs import config

basedir = os.path.dirname(os.path.abspath(__file__))

#{ response curves
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
    photfile = os.path.join(basedir,'calibration',photband.upper())
    wave, response = ascii.read2array(photfile).T
    return wave,response





def list_response(name='*',wave_range=(-np.inf,+np.inf)):
    """
    List available response curves.
    
    Specify a glob string C{name} and/or a wavelength range to make a selection
    of all available curves. If nothing is supplied, all curves will be returned.
    
    @param name: list all curves containing this string
    @type name: str 
    @param wave_range: list all curves within this wavelength range (A)
    @type wave_range: (float, float)
    @return: list of curve files
    @rtype: list of str
    """
    #-- collect all curve files
    curve_files = sorted(glob.glob(os.path.join(basedir,glob_string.upper())))
    #-- select in correct wavelength range
    curve_files = [os.path.basename(curve_file) for curve_file in curvefiles \
      if (wave_range[0]<=effwave(curve_files)<=wave_range[1])]
    #-- log to the screen and return
    for curve_file in curve_files: logger.info(curve_file)
    return curve_files







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



#}

if __name__=="__main__":
    import doctest
    doctest.testmod()