# -*- coding: utf-8 -*-
"""
Read and write FITS files.
"""
import logging
import os
import astropy.io.fits as pf

import numpy as np
from ivs.aux import loggers
from ivs.units import conversions

logger = logging.getLogger("IO.FITS")
logger.addHandler(loggers.NullHandler())

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
    flux = pf.getdata(filename)
    header = pf.getheader(filename)

    #-- Make the equidistant wavelengthgrid using the Fits standard info
    #   in the header
    ref_pix = int(header["CRPIX1"])-1
    dnu = float(header["CDELT1"])
    nu0 = float(header["CRVAL1"]) - ref_pix*dnu
    nun = nu0 + (len(flux)-1)*dnu
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
                         remove_flagged=True):
    """
    Read CoRoT data from a CoRoT FITS file.

    Both SISMO and EXO data are recognised and extracted accordingly.

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
    @param remove_flagged: remove flagged datapoints
    @type remove_flagged: True
    @return: CoRoT data (times, flux, error, flags)
    @rtype: array, array, array, array(, header)
    """
    #-- read in the FITS file
    # headers: ['DATE', 'DATEJD', 'DATEHEL', 'STATUS', 'WHITEFLUX', 'WHITEFLUXDEV', 'BG', 'CORREC']
    fits_file_    = pf.open(fits_file)
    if fits_file_[0].header['hlfccdid'][0]=='A':
        times,flux,error,flags = fits_file_[type_data].data.field(0),\
                                 fits_file_[type_data].data.field(1),\
                                 fits_file_[type_data].data.field(2),\
                                 fits_file_[type_data].data.field(3)
        # extract the header if asked
        if return_header:
            header = fits_file_[0].header
        fits_file_.close()

        logger.debug('Read CoRoT SISMO file %s'%(fits_file))
    elif fits_file_[0].header['hlfccdid'][0]=='E':
        times = fits_file_['bintable'].data.field('datehel')
        if 'blueflux' in fits_file_['bintable'].columns.names:
            blueflux,e_blueflux = fits_file_['bintable'].data.field('blueflux'),fits_file_['bintable'].data.field('bluefluxdev')
            greenflux,e_greenflux = fits_file_['bintable'].data.field('greenflux'),fits_file_['bintable'].data.field('greenfluxdev')
            redflux,e_redflux = fits_file_['bintable'].data.field('redflux'),fits_file_['bintable'].data.field('redfluxdev')
            #-- chromatic light curves
            if type_data=='colors':
                flux = np.column_stack([blueflux,greenflux,redflux])
                error = np.column_stack([e_blueflux,e_greenflux,e_redflux]).min(axis=1)
            #-- white light curves
            else:
                flux = blueflux + greenflux + redflux
                error = np.sqrt(e_blueflux**2 + e_greenflux**2 + e_redflux**2)
        else:
            flux,error = fits_file_['bintable'].data.field('whiteflux'),fits_file_['bintable'].data.field('whitefluxdev')
        flags = fits_file_['bintable'].data.field('status')

    # remove flagged datapoints if asked
    if remove_flagged:
        keep1 = (flags==0)
        keep2 = (error!=-1)
        logger.info('Remove: flagged (%d) no valid error (%d) datapoints (%d)'%(len(keep1)-keep1.sum(),len(keep2)-keep2.sum(),len(keep1)))
        keep = keep1 & keep2
        times,flux,error,flags = times[keep], flux[keep], error[keep], flags[keep]

    # convert times to heliocentric JD
    times = conversions.convert('MJD','JD',times,jtype='corot')

    if return_header:
        return times, flux, error, flags, header
    else:
        return times, flux, error, flags


def read_fuse(ff,combine=True,return_header=False):
    """
    Read FUSE spectrum.

    Modified JD: JD-2400000.5.

    Do 'EXPEND'-'EXPSTART'

    V_GEOCEN,V_HELIO

    ANO: all night only: data obtained during orbital night (highest SNR when airglow is not an issue)
    ALL: all: highest SNR with minimal airglow contamination

    Preference of ANO over ALL for science purpose.

    Use TTAGfcal files.
    """
    ff = pf.open(ff)
    hdr = ff[0].header
    if hdr['SRC_TYPE']=='EE':
        logger.warning("Warning: %s is not thrustworty (see manual)"%(ff))
    waves,fluxs,errrs = [],[],[]
    for seg in range(1,len(ff)):
        if ff[seg].data is None: continue
        waves.append(ff[seg].data.field('WAVE'))
        fluxs.append(ff[seg].data.field('FLUX'))
        errrs.append(ff[seg].data.field('ERROR'))
    ff.close()

    if combine:
        waves = np.hstack(waves)
        fluxs = np.hstack(fluxs)
        errrs = np.hstack(errrs)
        sa = np.argsort(waves)
        waves,fluxs,errrs = waves[sa],fluxs[sa],errrs[sa]
        keep = fluxs>0
        waves,fluxs,errrs = waves[keep],fluxs[keep],errrs[keep]


    if return_header:
        return waves,fluxs,errrs,hdr
    else:
        return waves,fluxs,errrs

def read_iue(filename,return_header=False):
    """
    Read IUE spectrum

    Instrumental profiles: http://starbrite.jpl.nasa.gov/pds/viewInstrumentProfile.jsp?INSTRUMENT_ID=LWR&INSTRUMENT_HOST_ID=IUE

    Better only use .mxlo for reliable absolute calibration!!

    LWP

    - Large-aperture spectral resolution is best between 2700 and 2900 A
    with an average FWHM of 5.2 A and decreases to approximately 8.0 A on
    either side of this range. Small-aperture resolution is optimal
    between 2400 and 3000 A with an average FWHM of 5.5 A and decreases to
    8.1 A at the extreme wavelengths.

    SWP

    - The best resolution occurs around 1200 A, with a FWHM of 4.6 A in the
    large aperture and 3.0 A in the small aperture, and gradually worsens
    towards longer wavelengths: 6.7 A at 1900 A in the large aperture and
    6.3 A in the small. On average, the small-aperture resolution is
    approximately 10% better than the large-aperture resolution.


    """
    ff = pf.open(filename)
    header = ff[0].header
    if os.path.splitext(filename)[1]=='.mxlo':
        try:
            flux = ff[1].data.field('flux')[0]
            error = ff[1].data.field('sigma')[0]
        except:
            flux = ff[1].data.field('abs_cal')[0]
            error = flux
        nu0 = ff[1].data.field('wavelength')[0]
        dnu = ff[1].data.field('deltaw')[0]
        nun = nu0 + len(flux)*dnu
        wavelength = np.linspace(nu0,nun,len(flux))
    elif os.path.splitext(filename)[1]=='.mxhi':
        orders_w = []
        orders_f = []
        orders_e = []
        for i in range(len(ff[1].data.field('abs_cal'))):
            flux = ff[1].data.field('abs_cal')[i]
            error = ff[1].data.field('noise')[i]
            quality = ff[1].data.field('quality')[i]
            npoints = ff[1].data.field('npoints')[i]
            startpix = ff[1].data.field('startpix')[i]
            flux = flux[startpix:startpix+npoints+1]
            error = error[startpix:startpix+npoints+1]
            quality = quality[startpix:startpix+npoints+1]
            flux[quality!=0] = 0
            nu0 = ff[1].data.field('wavelength')[i]
            dnu = ff[1].data.field('deltaw')[i]
            nun = nu0 + len(flux)*dnu
            wavelength = np.linspace(nu0,nun,len(flux))
            if len(orders_f)>0:
                orders_f[-1][wavelength[0]<orders_w[-1]] = 0
            orders_w.append(wavelength)
            orders_f.append(flux)
            orders_e.append(error)
        wavelength = np.hstack(orders_w)
        flux = np.hstack(orders_f)
        error = np.hstack(orders_e)

    wavelength = wavelength[flux!=0]
    error = error[flux!=0]
    flux = flux[flux!=0]

    logger.info("IUE spectrum %s read"%(filename))
    ff.close()
    if not return_header:
        return wavelength,flux,error
    else:
        return wavelength,flux,error,header

#}

#{ Generic reading

def read2recarray(fits_file,ext=1,return_header=False):
    """
    Read the contents of a FITS file to a record array.

    Should add a test that the strings were not chopped of...
    """
    dtype_translator = dict(L=np.bool,D=np.float64,E=np.float32,J=np.int)
    if isinstance(fits_file,str):
        ff = pf.open(fits_file)
    elif isinstance(fits_file,pf.HDUList):
        ff = fits_file
    data = ff[ext].data
    names = ff[ext].columns.names
    formats = ff[ext].columns.formats
    dtypes = []
    for name,dtype in zip(names,formats):
        if 'A' in dtype:
            dtypes.append((name,'S60'))
        else:
            dtypes.append((name,dtype_translator[dtype]))
    dtypes = np.dtype(dtypes)
    data = [np.cast[dtypes[i]](data.field(name)) for i,name in enumerate(names)]
    data = np.rec.array(data,dtype=dtypes)
    header = {}
    for key in ff[ext].header.keys():
        if 'TTYPE' in key: continue
        if 'TUNIT' in key: continue
        if 'TFORM' in key: continue
        header[key] = ff[ext].header[key]
    if isinstance(fits_file,str):
        ff.close()
    if not return_header:
        return data
    else:
        return data,header



#}


#{ Output

def write_primary(filename,data=None,header_dict={}):
    """
    Initiate a FITS file by writing to the primary HDU.

    If data is not given, a 1x1 zero array will be added.
    """
    if data is None:
        data = np.array([[0]])
    hdulist = pf.HDUList([pf.PrimaryHDU(data)])
    for key in header_dict:
        hdulist[0].header.update(key,header_dict[key])
    hdulist.writeto(filename)
    hdulist.close()
    return filename


def write_recarray(recarr,filename,header_dict={},units={},ext='new',close=True):
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
    is_file = isinstance(filename,str) and os.path.isfile(filename)
    if isinstance(filename,str) and not os.path.isfile(filename):
        primary = np.array([[0]])
        hdulist = pf.HDUList([pf.PrimaryHDU(primary)])
        hdulist.writeto(filename)
        hdulist.close()

    if is_file or isinstance(filename,str):
        hdulist = pf.open(filename,mode='update')
    else:
        hdulist = filename


    #-- create the table HDU
    cols = []
    for i,name in enumerate(recarr.dtype.names):
        format = recarr.dtype[i].str.lower().replace('|','').replace('>','')
        format = format.replace('b1','L').replace('<','')
        if 's' in format:                                                                                                              # Changes to be compatible with Pyfits version 3.3
            format = format.replace('s','') + 'A'                                                                                       # Changes to be compatible with Pyfits version 3.3
        unit = name in units and units[name] or 'NA'
        cols.append(pf.Column(name=name,format=format,array=recarr[name],unit=unit))
    tbhdu = pf.BinTableHDU.from_columns(pf.ColDefs(cols))

    #-- take care of the header:
    if len(header_dict):
        for key in header_dict:
            if (len(key)>8) and (not key in tbhdu.header.keys()) and (not key[:9]=='HIERARCH'):
                key_ = 'HIERARCH '+key
            else:
                key_ = key
            tbhdu.header[key_] = header_dict[key]
        if ext!='new':
            tbhdu.header['EXTNAME'] = ext


    #   put it in the right place
    extnames = [iext.header['EXTNAME'] for iext in hdulist if ('extname' in iext.header.keys()) or ('EXTNAME' in iext.header.keys())]
    if ext=='new' or not ext in extnames:
        logger.info('Creating new extension %s'%(ext))
        hdulist.append(tbhdu)
        ext = -1
    else:
        logger.info('Overwriting existing extension %s'%(ext))
        hdulist[ext] = tbhdu

    if close:
        hdulist.close()
        return filename
    else:
        return hdulist

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

    Instead of writing the file, you can give a hdulist and append to it.
    Supply a HDUList for 'filename', and set close=False
    """
    if isinstance(filename,str) and not os.path.isfile(filename):
            primary = np.array([[0]])
            hdulist = pf.HDUList([pf.PrimaryHDU(primary)])
            hdulist.writeto(filename)

    if isinstance(filename,str):
        hdulist = pf.open(filename,mode='update')
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
        else:
            unit = 'NA'
        cols.append(pf.Column(name=name,format=format,array=arr[i],unit=unit))
    tbhdu = pf.BinTableHDU.from_columns(pf.ColDefs(cols))                                                                            # Changes to be compatible with Pyfits version 3.3

    #   put it in the right place
    if ext=='new' or ext==len(hdulist):
        hdulist.append(tbhdu)
        ext = -1
    else:
        hdulist[ext] = tbhdu

    #-- take care of the header:
    if len(header_dict):
        for key in header_dict:
            hdulist[ext].header[key] = header_dict[key]

    if close:
        hdulist.close()
    else:
        return hdulist

#}
