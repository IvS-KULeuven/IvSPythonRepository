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
        wave = np.exp(wave)
    
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
    @rtype: array, array, array, array(, header)
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
        keep1 = (flags==0)
        keep2 = (error!=-1)
        logger.info('Remove: flagged (%d) no valid error (%d) datapoints (%d)'%(len(keep1)-sum(keep1),len(keep2)-sum(keep2),len(keep1)))
        keep = keep1 & keep2
        times,flux,error,flags = times[keep], flux[keep], error[keep], flags[keep]
    
    # convert times to heliocentric JD
    times = conversions.convert('MJD','JD',times,jtype='corot')
    
    if return_header:
        return times, flux, error, flags, header
    else:
        return times, flux, error, flags

#}


#{ Output

def write_recarray(recarr,filename,header_dict={},units={},ext='new'):
    """
    Write or add a record array to a FITS file.
    
    If 'filename' refers to an existing file, the record array will be added
    (ext='new') to the HDUlist or replace an existing HDU (ext=integer). Else,
    a new file will be created.
    
    Units can be given as a dictionary with keys the same names as the column
    names of the record array.
    
    A header_dictionary can be given, it is used to update an existing header
    or create a new one if the extension is new.
    """
    if not os.path.isfile(filename):
        primary = np.array([[0]])
        hdulist = pyfits.HDUList([pyfits.PrimaryHDU(primary)])
        hdulist.writeto(filename)
    
    hdulist = pyfits.open(filename,mode='update')
    
    #-- create the table HDU
    cols = []
    for i,name in enumerate(recarr.dtype.names):
        format = recarr.dtype[i].str.lower().replace('|','').replace('s','a').replace('>','')
        format = format.replace('b1','L').replace('<','')
        unit = name in units and units[name] or 'NA'
        cols.append(pyfits.Column(name=name,format=format,array=recarr[name],unit=unit))
    tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
    
    #   put it in the right place
    if ext=='new':
        hdulist.append(tbhdu)
        ext = -1
    else:
        hdulist[ext] = tbhdu
    
    #-- take care of the header:
    if len(header_dict):
        for key in header_dict:
            hdulist[ext].header.update(key,header_dict[key])
    
    hdulist.close()

def write_array(arr,filename,names=(),units=(),header_dict={},ext='new',close=True):
    """
    Write or add an array to a FITS file.
    
    If 'filename' refers to an existing file, the list of arrays will be added
    (ext='new') to the HDUlist or replace an existing HDU (ext=integer). Else,
    a new file will be created.
    
    Names and units should be given as a list of strings, in the same order as
    the list of arrays.
    
    A header_dictionary can be given, it is used to update an existing header
    or create a new one if the extension is new.
    """
    if close:
        if not os.path.isfile(filename):
            primary = np.array([[0]])
            hdulist = pyfits.HDUList([pyfits.PrimaryHDU(primary)])
            hdulist.writeto(filename)
        hdulist = pyfits.open(filename,mode='update')
    else:
        hdulist = filename
    
    #-- create the table HDU
    cols = []
    for i,name in enumerate(names):
        format = arr[i].dtype.str.lower().replace('|','').replace('s','a').replace('>','')
        format = format.replace('b1','L').replace('<','')
        if format=='f8':
            format = 'D'
        if isinstance(units,dict):
            unit = name in units and units[name] or 'NA'
        elif len(units)>i:
            unit = units[i]
        cols.append(pyfits.Column(name=name,format=format,array=arr[i],unit=unit))
    tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
    
    #   put it in the right place
    if ext=='new':
        hdulist.append(tbhdu)
        ext = -1
    else:
        hdulist[ext] = tbhdu
    
    #-- take care of the header:
    if len(header_dict):
        for key in header_dict:
            hdulist[ext].header.update(key,header_dict[key])
    
    if close:
        hdulist.close()
    else:
        return hdulist

#}