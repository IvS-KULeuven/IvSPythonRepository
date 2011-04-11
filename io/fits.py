# -*- coding: utf-8 -*-
"""
Read and write FITS files.
"""
import gzip
import logging
import os
import pyfits

import numpy as np
from ivs.misc import loggers

logger = logging.getLogger("IO.FITS")
logger.addHandler(loggers.NullHandler)

#{ Input

def read_spectrum(filename, return_header=False):
    """
    Read a standard 1D spectrum from the primary HDU of a FITS file.
    
    @param filename: FITS filename
    @type filename: str
    @param return_header: return header information as dictionary
    @type return_header: bool
    @return: wavelength, flux(, header)
    @rtype: array, array(, dict)
    """
    flux = pyfits.getdata(filename)
    header = pyfits.getheader(filename)
    
    nu0 = float(header["CRVAL1"])
    dnu = float(header["CDELT1"])
    nun = nu0 + len(flux)*dnu
    wave = np.linspace(nu0,nun,len(flux))
    #-- fix wavelengths for logarithmic sampling
    if 'ctype1' in header and header['CTYPE1']=='log(wavelength)':
        wave = np.exp(wavelength)
    
    logger.debug('Read spectrum %s'%(filename))
    
    if return_header:
        return wave,flux,header
    else:
        return wave,flux
    
#}
