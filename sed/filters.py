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
    #print photfile
    wave, response = ascii.read2array(photfile).T[:2]
    sa = np.argsort(wave)
    return wave[sa],response[sa]





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
    #-- collect all curve files but remove human eye responses
    if not '*' in name:
        name = '*' + name + '*'
    curve_files = sorted(glob.glob(os.path.join(basedir,'filters',name.upper())))
    curve_files = [cf for cf in curve_files if not ('HUMAN' in cf or 'EYE' in cf) ]
    #-- select in correct wavelength range
    curve_files = [os.path.basename(curve_file) for curve_file in curve_files if (wave_range[0]<=eff_wave(os.path.basename(curve_file))<=wave_range[1])]
    #-- log to the screen and return
    for curve_file in curve_files: logger.info(curve_file)
    return curve_files



def is_color(photband):
    """
    Return true if the photometric passband is actually a color.
    
    @param photband: name of the photometric passband
    @type photband: string
    """
    if '-' in photband.split('.')[1]:
        return True
    elif photband.split('.')[1].upper() in ['M1','C1']:
        return True
    else:
        return False




def eff_wave(photband,model=None):
    """
    Return the effective wavelength of a photometric passband.
    
    The effective wavelength is defined as the average wavelength weighed with
    the response curve.
    
    >>> eff_wave('2MASS.J')
    12412.136241640892
    
    If you give model fluxes as an extra argument, the wavelengths will take
    these into account to calculate the `true' effective wavelength (e.g., 
    Van Der Bliek, 1996), eq 2.
    
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER') or array/list of str
    @param model: model wavelength and fluxes
    @type model: tuple of 1D arrays (wave,flux)
    @return: effective wavelength [A]
    @rtype: float or numpy array
    """
    
    #-- if photband is a string, it's the name of a photband: put it in a container
    #   but unwrap afterwards
    if isinstance(photband,str):
        single_band = True
        photband = [photband]
    #-- else, it is a container
    else:
        single_band = False
        
    my_eff_wave = []
    for iphotband in photband:
        try:
            wave,response = get_response(iphotband)
            if model is None:
                this_eff_wave = np.average(wave,weights=response)
            else:
                #-- interpolate response curve onto higher resolution model and
                #   take weighted average
                is_response = response>1e-10
                start_response,end_response = wave[is_response].min(),wave[is_response].max()
                fluxm = 10**np.interp(np.log10(wave),np.log10(model[0]),np.log10(model[1]))
                this_eff_wave = np.trapz(wave*fluxm*response,x=wave) / np.trapz(fluxm*response,x=wave)
        #-- if the photband is not defined:
        except IOError:
            this_eff_wave = np.nan
        my_eff_wave.append(this_eff_wave)
    
    if single_band:
        my_eff_wave = my_eff_wave[0]
    else:
        my_eff_wave = np.array(my_eff_wave,float)
    
    return my_eff_wave

@memoized
def get_info(photbands=None):
    """
    Return a record array containing all filter information.
    
    The record arrays contains following columns:
        - photband
        - eff_wave
        - type
        - vegamag, vegamag_lit
        - ABmag, ABmag_lit
        - STmag, STmag_lit
        - Flam0, Flam0_units, Flam0_lit
        - Fnu0, Fnu0_units, Fnu0_lit,
        - source
    
    @param photbands: list of photbands to get the information from. The input
    order is equal to the output order. If C{None}, all filters are returned.
    @type photbands: iterable container (list, tuple, 1Darray)
    @return: record array containing all information on the requested photbands.
    @rtype: record array
    """
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp = ascii.read2recarray(zp_file)
    zp = zp[np.argsort(zp['photband'])]
    
    #-- list photbands in order given, and remove those that do not have
    #   zeropoints etc.
    if photbands is not None:
        order = np.searchsorted(zp['photband'],photbands)
        zp = zp[order]
        keep = (zp['photband']==photbands)
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
    ascii.write_array(zp,'zeropoints.dat',header=True,auto_width=True,comments=['#'+line for line in comms[:-2]],use_float='%g')
    


if __name__=="__main__":
    import doctest
    doctest.testmod()