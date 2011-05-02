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
from ivs.units import conversions

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


def read_corot(fits_file,  return_header=False, type_data='hel',
                         channel='sismo', remove_flagged=True):
    """
    Read CoRoT data from a CoRoT FITS file.
    
    Currently only SISMO data is implemented.
    
    type_data is one of:
        - type_data='raw'
        - type_data='hel': heliocentric unequidistant
        - type_data='helreg': heliocentric equidistant
    
    @param fits_file: CoRoT FITS file name
    @type fits_file: string
    @param return_header: return header information as dictionary
    @type return_header: bool
    @param type_data: type of data to return
    @type type_data: string (one of 'raw','hel' or 'helreg')
    @param channel: type of CoRoT data
    @type channel: string (one of 'sismo' or 'exo')
    @param remove_flagged: remove flagged datapoints
    @type remove_flagged: True
    @return: CoRoT data (times, flux, error, flags)
    @type: array, array, array, array(, header)
    """
    #-- read in the FITS file
    # headers: ['DATE', 'DATEJD', 'DATEHEL', 'STATUS', 'WHITEFLUX', 'WHITEFLUXDEV', 'BG', 'CORREC']
    fits_file_    = pyfits.open(fits_file)
    times,flux,error,flags = fits_file_[type_data].data.field(0),\
                             fits_file_[type_data].data.field(1),\
                             fits_file_[type_data].data.field(2),\
                             fits_file_[type_data].data.field(3)
    # extract the header if asked
    if return_header:
        header = fits_file_[0].header
    fits_file_.close()
    
    logger.debug('Read CoRoT file %s'%(fits_file))
    
    # remove flagged datapoints if asked
    if remove_flagged:
        keep = (flags==0) & (error!=-1)
        times,flux,error,flags = times[keep], flux[keep], error[keep], flags[keep]
    
    # convert times to heliocentric JD
    times = conversions.convert('MJD','JD',times,jtype='corot')
    
    if return_header:
        return times, flux, error, flags, header
    else:
        return times, flux, error, flags

#}
