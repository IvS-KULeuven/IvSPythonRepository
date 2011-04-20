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
import pyfits
import logging
import numpy as np

from ivs import config
from ivs.misc.decorators import memoized
from ivs.misc import loggers
from ivs.io import ascii

basedir = os.path.dirname(__file__)

logger = logging.getLogger("CAT.VIZIER")
logger.addHandler(loggers.NullHandler)

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
    photfile = os.path.join(basedir,'filters',photband.upper())
    wave, response = ascii.read2array(photfile).T[:2]
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
    @type photband: str ('SYSTEM.FILTER') or array/list of str
    @return: effective wavelength [A]
    @rtype: float or numpy array
    """
    #-- if photband is a string, it's the name of a photband
    if isinstance(photband,str):
        wave,response = get_response(photband)
        my_eff_wave = np.average(wave,weights=response)
    #-- else, it is a container
    else:
        my_eff_wave = []
        for iphotband in photband:
            try:
                wave,response = get_response(iphotband)
                my_eff_wave.append(np.average(wave,weights=response))
            except IOError:
                my_eff_wave.append(np.nan)
            
        my_eff_wave = np.array(my_eff_wave,float)
            
    return my_eff_wave


@memoized
def get_info(photbands=None):
    """
    Return a record array containing all filter information.
    """
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp = ascii.read2recarray(zp_file)
    zp = zp[np.argsort(zp['photband'])]
    
    if photbands is not None:
        keep = np.searchsorted(zp['photband'],photbands)
        zp = zp[keep]
    
    return zp

def update_info(zp):
    """
    Update information in zeropoint file, e.g. after calibration.
    
    @param zp: updated contents from C{zeropoints.dat}
    @type zp: recarray
    """
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp_,comms = ascii.read2recarray(zp_file,return_comments=True)
    ascii.write_array(zp,'zeropoints.dat',header=True,auto_width=True,comments=['#'+line for line in comms[:-2]])
    


if __name__=="__main__":
    import doctest
    doctest.testmod()