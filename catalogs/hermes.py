# -*- coding: utf-8 -*-
"""
Interface the spectra from the Hermes spectrograph.

The most important function is L{search}. This looks in SIMBAD for the coordinates
of a given object, and finds all spectra matching those within a given radius.
If the object's name is not recognised, it will look for correspondence between
the given name and the contents of the FITS header's keyword C{object}. L{search}
returns a record array, so that you can easily access all columns by their names.

Note that Hermes spectra retrieved in B{logscale} should B{not be corrected} for
the B{barycentric velocity} (the pipeline does it). The spectra retrieved in B{normal
wavelength scale should} be corrected. These two cases are outlined below.

Section 1. Data lookup and reading
==================================

B{Example usage:} retrieve all data on HD170580

>>> mydata = search('HD170580')

Keep only those with a long enough exposure time:

>>> myselection = mydata[mydata['exptime']>500]

Now read in all the data, and plot the spectra. First, we need an extra module
to read the FITS file and the plotting package.

>>> from ivs.inout import fits
>>> import pylab as pl

Then we can easily plot the relevant data:

>>> for fname in myselection['filename']:
...     wave,flux = fits.read_spectrum(fname)
...     p = pl.plot(wave,flux)
>>> p = pl.ylim(0,45000)
>>> p = pl.xlim(4511,4513.5)

]include figure]]ivs_catalogs_hermes_HD170580.png]

Note that you can easily shift them according to some radial velocity: as an
extra example, we retrieve the data in wavelength scale, shift according to the
barycentric velocity, convert to velocity space, and plot them:

First start a new figure and add extra modules:

>>> p = pl.figure()
>>> from ivs.spectra.tools import doppler_shift         # to apply doppler shift
>>> from ivs.units import conversions     # to convert to velocity space

Then get the spectra again, but now not in log scale. Make the same selection
on exposure time.

>>> mydata = search('HD170580',data_type='cosmicsremoved_wavelength')
>>> myselection = mydata[mydata['exptime']>500]

Then read them all in and shift according to the barycentric velocity. Also,
convert the spectra to velocity space afterwards for plotting purposes.

>>> for rv,fname in zip(myselection['bvcor'],myselection['filename']):
...     wave,flux = fits.read_spectrum(fname)
...     wave_shifted = doppler_shift(wave,rv)
...     velo_shifted = conversions.convert('angstrom','km/s',wave_shifted,wave=(4512.3,'angstrom'))
...     p = pl.plot(velo_shifted,flux)
>>> p = pl.ylim(0,45000)
>>> p = pl.xlim(-70,70)

]include figure]]ivs_catalogs_hermes_HD170580_velo.png]

Section 2. Extracting information
=================================

If you want to list the observers of a certain target, and the amount of spectra
they took, you can do the following:

>>> data = search('HD50230')
>>> print ([(i,list(data['observer']).count(i)) for i in set(data['observer'])])
[('Robin Lombaert', 4), ('Steven Bloemen', 6), ('Pieter Degroote', 25), ('Michel Hillen', 3)]

Section 3. Radial velocity computation with the HERMES DRS
==========================================================

Make sure you have sourced the Hermes DRS (C{source ~mercator/HermesDRS.rc}) and
that you make a backup of your config file before proceeding.

In this example, we first make a mask file:

>>> from ivs.spectra import linelists
>>> teff = 12500
>>> logg = 4.0
>>> ll = linelists.get_lines(teff,logg)
>>> mask_file = '%.0f_%.2f.fits'%(teff,logg)
>>> # make_mask_file(ll['wavelength'],ll['depth'],filename=mask_file) The Hermes pipeline is Py2.7

And run the DRS:

>>> #CCFList('HD170200',mask_file=mask_file)


Section 4. Hermes overview file
===============================

To ensure a fast lookup of datafiles, an overview file C{HermesFullDataOverview.tsv}
is created via L{make_data_overview}. The file resides in one of the C{IVSdata}
directories, and there should also exist a copy at C{/STER/mercator/hermes/}.

The best strategy to keep the file up-to-date is by running this module in the
background in a terminal, via (probably on pleiad22)::

    $:> python hermes.py update

This script will look for a C{HermesFullDataOverview.tsv} in one of the data
directories (see L{config} module), check until what date data files are stored,
and add entries for HERMES data files created since the last update. If you
want a complete update of all the HERMES data folders, you will need to remove
the file and run the update. Indeed, if C{HermesFullDataOverview.tsv} does not
exist, it will create a new file in the current working directory from scratch,
and copy it to C{/STER/mercator/hermes/}. It is the user's responsibility to
copy the file also to one of the IVSdata directories where the user has write
permission, so that next time the update is performed, the script can start from
the previous results.

When running the above command in the terminal, you will notice that the script
does not terminate until a manual C{CTRL+C} is performed. Indeed, when left running
the script will perform the update once a day. So once it runs, as long as the
filepaths do not change or computers do not shut down, you don't need to run it
again. If a computer does stop running, just restart again with::

    $:> python hermes.py update

and everything should be fine.

@var hermesDir: Path to the directory in which the hermes folders (hermesAnalysis,
                hermesRun, hermesDebug) are located. By default this is set to
                '/home/user/'.

@var tempDir: Path to temporary directory where files nessessary for hermesVR to run
              can be copied to. You need write permission in this directory. The
              default is set to '/scratch/user/'
"""
import os
import sys
import glob
import time
import shutil
import logging
import getpass
import datetime
import subprocess
import numpy as np
import astropy.io.fits as pf

import copy
# from lxml import etree
# from xml.etree import ElementTree as ET
from collections import defaultdict

from ivs.catalogs import sesame
from ivs.inout import ascii, fits
from ivs.aux import loggers
from ivs.observations import airmass
from ivs.observations.barycentric_correction import helcorr
from ivs.units import conversions
from ivs import config
from ivs.spectra import tools as sptools


hermesDir = os.path.expanduser('~/')
tempDir = '/scratch/%s/'%(getpass.getuser())

logger = logging.getLogger("CAT.HERMES")
logger.addHandler(loggers.NullHandler())

#{ User functions

def search(ID=None,time_range=None,prog_ID=None,data_type='cosmicsremoved_log',
           radius=1.,filename=None):
    """
    Retrieve datafiles from the Hermes catalogue.

    B{If C{ID} is given}: A string search is performed to match the 'object'
    field in the FITS headers. The coordinates are pulled from SIMBAD. If the
    star ID is recognised by SIMBAD, an additional search is done based only on
    the coordinates. The union of both searches is the final result.

    B{If C{time_range} is given}: The search is confined within the defined
    range. If you only give one day, the search is confined to the observations
    made during the night starting at that day. If C{ID} is not given, all
    observations will be returned of the given datatype.

    B{If C{prog_ID} is given}: The search is performed to match the number of
    the program. Individual stars are not queried in SIMBAD, so any information
    that is missing in the header will not be corrected.

    If you don't give either ID or time_range, the info on all data will be
    returned. This is a huge amount of data, so it can take a while before it
    is returned. Remember that the header of each spectrum is read in and checked.

    Data type can be any of:
        1. cosmicsremoved_log: return log merged without cosmics
        2. cosmicsremoved_wavelength: return wavelength merged without cosmics
        3. ext_log: return log merged with cosmics
        4. ext_wavelength: return wavelength merged with cosmics
        5. raw: raw files (also TECH..., i.e. any file in the raw directory)

    This functions needs a C{HermesFullDataOverview.tsv} file located in one
    of the datadirectories from C{config.py}, and subdirectory C{catalogs/hermes}.

    If this file does not exist, you can create it with L{make_data_overview}.

    If you want a summary file with the data you search for, you can give
    C{filename} as an extra keyword argument. The results will be saved to that
    file.

    The columns in the returned record array are listed in L{make_data_overview},
    but are repeated here (capital letters are directly retrieved from the
    fits header, small letters are calculated values. The real header strings
    are all small capitals):

        1.  UNSEQ
        2.  PROG_ID
        3.  OBSMODE
        4.  BVCOR
        5.  OBSERVER
        6.  OBJECT
        7.  RA
        8.  DEC
        9.  BJD
        10. EXPTIME
        11. PMTOTAL
        12. DATE-AVG
        13. OBJECT
        14. airmass
        15. filename

    The column C{filename} contains a string with the absolute location of the
    file. If you need any extra information from the header, you can easily
    retrieve it.

    If BVCOR or BJD are not available from the FITS header, this function will
    attempt to calculate it. It will not succeed if the object's name is not
    recognised by SIMBAD.

    Example usage: retrieve all data on HD50230

    >>> mydata = search('HD50230')

    Keep only those with a long enough exposure time:

    >>> myselection = mydata[mydata['exptime']>500]

    Look up the 'telalt' value in the FITS headers of all these files via a fast
    list comprehension:

    >>> telalts = [pf.getheader(fname)['telalt'] for fname in myselection['filename']]

    Search for all data of HD50230 taken in the night of 22 September 2009:

    >>> data = search('HD50230',time_range='2009-9-22')

    Or within an interval of a few days:

    >>> data = search('HD50230',time_range=('2009-9-22','2009-9-30'))

    Search for all data observed in a given night:

    >>> data = search(time_range='2009-9-22')

    B{Warning:} the heliocentric correction is not calculated when no ID is given,
    so make sure it is present in the header if you need it, or calculate it yourself.

    @param ID: ID of the star, understandable by SIMBAD
    @type ID: str
    @param time_range: range of dates to confine the search to
    @type time_range: tuple strings of type '2009-09-23T04:24:35.712556' or '2009-09-23'
    @param data_type: if None, all data will be returned. Otherwise, subset
    'cosmicsremoved', 'merged' or 'raw'
    @type data_type: str
    @param radius: search radius around the coordinates (arcminutes)
    @type radius: float
    @param filename: write summary to outputfile if not None
    @type filename: str
    @return: record array with summary information on the observations, as well
    as their location (column 'filename')
    @rtype: numpy rec array
    """
    #-- read in the data from the overview file, and get SIMBAD information
    #   of the star
    ctlFile = '/STER/mercator/hermes/HermesFullDataOverview.tsv'
    data = ascii.read2recarray(ctlFile, splitchar='\t')
    #data = ascii.read2recarray(config.get_datafile(os.path.join('catalogs','hermes'),'HermesFullDataOverview.tsv'),splitchar='\t')
    keep = np.array(np.ones(len(data)),bool)
    #-- confined search within given time range
    if time_range is not None:
        if isinstance(time_range,str):
            time_range = _timestamp2datetime(time_range)
            time_range = (time_range,time_range+datetime.timedelta(days=1))
        else:
            time_range = (_timestamp2datetime(time_range[0]),_timestamp2datetime(time_range[1]))
        keep = keep & np.array([(time_range[0]<=_timestamp2datetime(i)<=time_range[1]) for i in data['date-avg']],bool)
        info = None


    #-- search on ID
    if ID is not None:
        info = sesame.search(ID)

        #-- first search on object name only
        ID = ID.replace(' ','').replace('.','').replace('+','').replace('-','').replace('*','')
        match_names = np.array([objectn.replace(' ','').replace('.','').replace('+','').replace('-','').replace('*','') for objectn in data['object']],str)
        keep_id = [((((ID in objectn) or (objectn in ID)) and len(objectn)) and True or False) for objectn in match_names]
        keep_id = np.array(keep_id)
        #   if we found the star on SIMBAD, we use its RA and DEC to match the star
        if info:
            ra,dec = info['jradeg'],info['jdedeg']
            keep_id = keep_id | (np.sqrt((data['ra']-ra)**2 + (data['dec']-dec)**2) < radius/60.)
        keep = keep & keep_id

    if prog_ID is not None:
        keep = keep & (data['prog_id']==prog_ID)

    #-- if some data is found, we check if the C{data_type} string is contained
    #   with the file's name. If not, we remove it.
    if np.any(keep):
        data = data[keep]

        if data_type is not None:
            data_type == data_type.lower()
            #-- now derive the location of the 'data_type' types from the raw
            #   files
            if not data_type=='raw':
                data['filename'] = [_derive_filelocation_from_raw(ff,data_type) for ff in data['filename']]
                existing_files = np.array([ff!='naf' for ff in data['filename']],bool)
                data = data[existing_files]
            seqs = sorted(set(data['unseq']))
        logger.info('ID={}/prog_ID={}: Found {:d} spectra (data type={} with unique unseqs)'.format(ID,prog_ID,len(seqs),data_type))
    else:
        data = data[:0]
        logger.info('%s: Found no spectra'%(ID))

    #-- we now check if the barycentric correction was calculated properly.
    #   If not, we calculate it here, but only if the object was found in
    #   SIMBAD. Else, we have no information on the ra and dec (if bvcorr was
    #   not calculated, ra and dec are not in the header).
    for obs in data:
        if ID is not None and info:
            try:
                jd  = _timestamp2jd(obs['date-avg'])
            except ValueError:
                logger.info('Header probably corrupted for unseq {}: no info on time or barycentric correction'.format(obs['unseq']))
                jd = np.nan
            # the previous line is equivalent to:
            # day = dateutil.parser.parse(header['DATE-AVG'])
            # BJD = ephem.julian_date(day)
            bvcorr, hjd = helcorr(ra/360.*24, dec, jd)
        else:
            break
        if np.isnan(obs['bvcor']):
            logger.info("Corrected 'bvcor' for unseq {} (missing in header)".format(obs['unseq']))
            obs['bvcor'] = float(bvcorr)
        if np.isnan(obs['bjd']):
            logger.info("Corrected 'bjd' for unseq {} (missing in header)".format(obs['unseq']))
            obs['bjd'] = float(hjd)


    #-- do we need the information as a file, or as a numpy array?
    if filename is not None:
        ascii.write_array(data,filename,auto_width=True,header=True)
    else:
        return data

def merge_hermes_spectra(objlist, wscalelist=None, **kwargs):
    """
    Combines HERMES spectra with sigma clipping to remove cosmics. The input spectra
    are given as a list of filenames for the objects and for the wavelength scales.
    The output is the wavelength scale of the first object, the merged and sigma clipped
    flux, and the adapted header of the object. In this adapted header, the exposure
    time is the sum of all individual exposure times, and the observing time is averaged.
    Uses the L{ivs.spectra.tools.merge_cosmic_clipping} method to merge the spectra.

    >>> data = search('KIC9540226')
    >>> objlist = data['filename']
    >>> wave, flux, header = merge_hermes_spectra(objlist, wscalelist=None)

    @param objlist: list of OBJ or order merged filenames
    @param wscalelist: list of wavelengthscale filenames
    @param vrads: list of radial velocities (optional)
    @param vrad_units: units of the radial velocities
    @param sigma: value used for sigma clipping
    @param window: window size used in median filter
    @param runs: number of iterations through the spectra

    @return: The combined spectrum, cosmic clipped. (wave, flux, header)
    @rtype: [array, array, dict]
    """
    kwargs['full_output'] = False

    if wscalelist == None:
        #-- Order Merged spectra are straight forward
        exptime, bjd = 0, np.zeros(len(objlist))
        waves, fluxes = [],[]
        for i, ofile in enumerate(objlist):
            w, f, h = fits.read_spectrum(ofile, return_header=True)
            waves.append(w)
            fluxes.append(f)
            exptime += h['exptime']
            #bjd[i] = h['bjd']

        waves, fluxes = np.array(waves), np.array(fluxes)
        mwave, mflux = sptools.merge_cosmic_clipping(waves, fluxes, **kwargs)

    else:
        #-- The OBJ spectra need to be merged order per order
        # read all spectra
        exptime, bjd = 0, np.zeros(len(objlist))
        waves, fluxes = np.zeros((55,4608, len(objlist))), np.zeros((55,4608, len(wscalelist)))
        for i, (ofile, wfile) in enumerate(zip(objlist, wscalelist)):
            odata = pf.getdata(ofile, 0).T # The flux has format (4608, 55) thus is transposed
            oheader = pf.getheader(ofile, 0)
            wdata = pf.getdata(wfile, 0)
            fluxes[:,:,i] = odata
            waves[:,:,i] = wdata
            exptime += oheader['exptime']
            #bjd[i] = oheader['bjd']

        # merge the spectra
        mwave, mflux = waves[:,:,0], np.zeros((55, 4608))
        for i in range(55):
            wave, flux = sptools.merge_cosmic_clipping(waves[i,:,:].T, fluxes[i,:,:].T, **kwargs)
            mflux[i] = flux

    header = pf.getheader(objlist[0], 0)
    header['exptime'] = exptime
    #header['bjd'] = np.average(bjd)

    return mwave, mflux, header

def make_list_star(ID,direc=None):
    """
    Mimics HermesTool MakeListStar.

    This should work as input for HermesTool CCFList.py

    The result is a file in the current working directory with name C{ID.list}.
    If you have specified the C{direc} keyword, the file will be written inside
    that directory. Make sure you have write permission.

    The contents of the file is:

    unseq, date-avg, ID, bjd, bvcor, prog_id, exptime, airmass, pmtotal

    @param ID: name of the star, understandable by SIMBAD.
    @type ID: string
    @param direc: directory to write the file into (defaults to current working
    directory)
    @type direc: string
    """
    if direc is None:
        direc = os.getcwd()
    data = search(ID)
    fname = os.path.join(direc,'%s.list'%(ID))
    ascii.write_array([data['unseq'],data['date-avg'],[ID for i in data['unseq']],
                       data['bjd'],data['bvcor'],data['prog_id'],data['exptime'],
                       data['airmass'],data['pmtotal']],fname,axis0='cols')

def make_mask_file(wavelength,depth,filename='mymask.fits'):
    """
    Make a mask file for RV calculations with the Hermes pipeline.

    See L{ivs.units.linelists} to select appropriate masks. You readily use the
    columns C{wavelength} and C{depth} as input for this function.

    @param wavelength: wavelength in angstrom
    @type wavelength: 1D numpy array
    @param depth: strenght of the line (1-normalised flux minimum)
    @type depth: 1D numpy array
    """
    data = np.array([wavelength,depth]).T
    hdulist = pf.PrimaryHDU(data)
    hdulist.writeto(filename)


def calibrate(wave,flux,header=None):
    """
    Rough calibration of Hermes spectrum.

    Thanks to Roy Ostensen for computing the response function via detailed
    modeling of Feige 66.
    """
    #-- read instrument response function
    instrument = config.get_datafile('catalogs/hermes','hermes_20110319.resp')
    iwave,iflux = ascii.read2array(instrument).T
    keep = (iwave.min()<=wave) & (wave<=iwave.max())
    #-- interpolate on observed wavelength array
    iflux = np.interp(wave[keep],iwave,iflux)
    return wave[keep],flux[keep]/iflux,iflux
#}

#{ Read / Write

def read_hermesVR_velocities(unique=True, return_latest=True, unseq=None, object=None, **kwargs):
    """
    Read the velocities.data file with the resuls from hermesVR. The results are
    returned as a numpy record array. By setting unique to True, the method will
    only return a unique set of sequence numbers. When return_latest is True it will
    pick the latest added version, otherwise the first added is picked.

    You can supply a list of unseq numbers, and then only those will be returned,
    same goes for object. In both cases the output is ordered in the same way as
    the list of unseq or object that you provided. These lists can be used together
    with the unique option.

    This method will search for 'velocities.data' in 'hermesDir/hermesAnalyses/'.
    You can provide another file location by using the keyword filename.

    The field in the recored array are:

        - unseq: unique number
        - object: object name
        - hjd: HJD
        - exptime: exposure time
        - bvcorr: bary_corr
        - rvdrift: RV drift (2F frames only)
        - telloff: telluric offset
        - tellofferr: error on telluric offset
        - telloffwidth: width telluric fit
        - vrad: Vrad(55-74) (bvcorr applied)
        - vraderr: err on Vrad
        - nlines: number of lines used
        - ccfdepth: depth_CCF
        - ccfdeptherr: error on depth_CCF
        - ccfsigma: sigma_CCF
        - ccfsigmaerr: err on sigma_CCF
        - sn: signal to noise in spectrum (orders 55-74)
        - gaussoc: O-C of Gaussian fit
        - wvlfile: wavelength calibration file
        - maskfile: mask file
        - fffile: Flatfield file

    @param unique: True is only return unique sequence numbers
    @type unique: bool
    @param return_latest: True to pick the last added file
    @type return_latest: bool
    @param unseq: list of sequence numbers to return
    @type unseq: list
    @param object: list of object names to return
    @type object: list
    @param filename: Path to the velocities.data file
    @type filename: str

    @return: array with the requested result
    @rtype: record array
    """

    #-- read data from velocities.data
    velocity_file = kwargs.pop('filename', hermesDir + 'hermesAnalyses/velocities.data')
    if not os.path.isfile(velocity_file):
        logger.warning('velocities.data file not found!')
        return []
    rawdata = ascii.read2list(velocity_file)[0]

    #-- parse the data and store in rec array
    for i, line in enumerate(rawdata):
        # here are some options to catch problems before creating the recarray.
        if len(line) > 21:
            line[20] = " ".join(line[20:])
            rawdata[i] = tuple(line[0:21])

    dtype = [('unseq','int'),('object','U20'),('hjd','f8'),('exptime','f8'),('bvcor','f8'),
             ('rvdrift','f8'),('telloff','f8'),('tellofferr','f8'),('telloffwidth','f8'),
             ('vrad','f8'),('vraderr','f8'),('nlines','int'),('ccfdepth','f8'),
             ('ccfdeptherr','f8'), ('ccfsigma','f8'),('ccfsigmaerr','f8'),('sn','f8'),
             ('gaussoc','f8'),('wvlfile','U80'),('maskfile','U80'),('fffile','U80')]
    data = np.rec.array(rawdata, dtype=dtype)

    #-- select which lines to keep
    if unique:
        unique_sequence = set(data['unseq'])
        order = -1 if return_latest else 0
        select = []
        for seq in unique_sequence:
            select.append(np.where(data['unseq']==seq)[0][order])
        data = data[(select,)]

    if unseq != None:
        # also keep the order in which unseq is provided
        select = np.ravel(np.hstack([np.where(data['unseq']==sq) for sq in unseq]))
        data = data[(select,)]


    if object != None:
        select = np.ravel(np.hstack([np.where(data['object']==ob) for ob in object]))
        data = data[select]

    return data

def read_hermesVR_AllCCF(unseq):
    """
    Read the AllCCF.fits files written by hermesVR with the individual cross
    correlation functions of each order. The cross-correlation functions are stored
    in a HermesCCF object from which they can be requested.

    You can read *_AllCCF.fits file by providing only the unique sequence number as
    an integer. In this case the file will be searched in the hermesDir/hermesAnalyses/
    folder. You can also provide the complete filename instead as a string fx.

    >>> #read_hermesVR_AllCCF(457640)
    >>> #read_hermesVR_AllCCF('/home/jorisv/hermesAnalysis/00457640_AllCCF.fits')

    @param unseq: unique sequence number or filename to read
    @type unseq: int or str

    @return: The cross correlation functions
    @rtype: HermesCCF object
    """

    if type(unseq) == int:
        filename = hermesDir + 'hermesAnalyses/*{:.0f}_AllCCF.fits'.format(unseq)
        files = glob.glob( filename )
        if len(files) > 0:
            filename = files[0]
    else:
        filename = unseq

    if not os.path.isfile(filename):
        logger.warning("Filename ''{:}'' is NOT valid!, returning empty object".format(filename))
        return HermesCCF(filename=None)

    return HermesCCF(filename=filename)

def write_hermesConfig(**kwargs):
    """
    Update the hermesConfig.xml file with the specified keywords. This method will
    return the original content of the file as a dictionary. You can use that dict
    as an argument to this function if you want to restore the hermesConfig file to
    the original state.

    write_hermesConfig will search for hermesConfig.xml in hermesDir/hermesRun/
    unless a path is specified in the kwarg filename.

    example use:

    >>> #old = write_hermesConfig(CurrentNight="/STER/mercator/hermes/20101214")
    >>> #" Do some stuff "
    >>> #write_hermesConfig(**old)

    @attention: This function does not check the validity of the keywords that you
    provide, that is your own responsibility.

    @param filename: complete path to the hermesConfig.xml file
    @type filename: str

    @return: The original content of hermesConfig.xml before modification.
    @rtype: dict
    """

    config_file = kwargs.pop('filename', hermesDir + 'hermesRun/hermesConfig.xml')

    #-- read current config file to dictionary
    old_config = _etree_to_dict( ET.parse(config_file).getroot() )['hermes']

    if len(list(kwargs.keys())) == 0.0:
        logger.debug('Nothing written to hermesConfig.xml, returning current config.')
        return old_config

    #-- update the config dictionary
    new_config = copy.deepcopy(old_config)
    new_config.update(kwargs)

    hermes = etree.Element("hermes")
    for key in list(new_config.keys()):
        child = etree.SubElement(hermes, key)
        child.text = new_config[key]

    #-- write to file
    out = etree.ElementTree(element=hermes)
    out.write(config_file, pretty_print=True)

    logger.debug('Written new hermesConfig to file:\n'+etree.tostring(hermes, pretty_print=True))

    return old_config


#}

#{ Hermes DRS wrappers

def run_hermesVR(filename, mask_file=None, wvl_file=None, cosmic_clipping=True,
                 flatfield=False, vrange='default', vrot=False, version='release',
                 verbose=False, timeout=500, **kwargs):
    """
    Wrapper around hermesVR that will run hermesVR on any given filename. Thus not on
    sequence number. This can be usefull if you first merge together 2 back-to-back
    spectra and then want to calculate the radial velocity of the merged spectrum.

    The original file will not be modified, it is coppied to the tempDir, and the name
    will be changed according to the supplied flags. During the process hermesConfig.xml
    is adapted, but after hermesVR is finished, it is restored to its old state.

    When providing a wavelength scale file, provide the full path, but cut off the
    '_HRF_TH_ext_wavelengthScale.fits' part of the filename. Thus if you want to use:
    '/folder/2451577_HRF_TH_ext_wavelengthScale.fits', use:

    >>> #run_hermesVR(filename, wvl_file = '/folder/2451577')

    The results of hermesVR are saved to the standard paths, as provided in the
    hermesConfig file. This method only returns the unseq number of the filename
    so you can find the corresponding output file, the output of hermesVR as a
    string, and a unix return signal of running hermesVR.
    The returncode is 0 if no errors are found. Any negative number means hermesVR
    failed (codes can be consulted here:  U{unix signals
    <http://people.cs.pitt.edu/~alanjawi/cs449/code/shell/UnixSignals.htm>}).
    If returncode = -35, then hermesVR completed but didn't write any output. In
    this case consult the logs or string output to find out what went wrong.

    @param filename: complete path to the input spectrum for hermesVR
    @type filename: str
    @param mask_file: complete path to the mask file to use or None
    @type mask_file: str
    @param wvl_file: path to wavelengthscale file or None
    @type wvl_file: str
    @param cosmic_clipping: Doesn't do anything for now
    @type cosmic_clipping: bool
    @param flatfield: Use flatfield devision or not
    @type flatfield: bool
    @param vrange: Size of velocity window ('default', 'large' (x2), 'xlarge' (x5))
    @type vrange: str
    @param vrot: Fit CCF with rotationaly broadend profile instead of gaussian
    @type vrot: bool
    @param version: which version of the pipeline to use: 'release' or 'trunk'
    @type version: str
    @param verbose: print output of hermesVR
    @type verbose: bool
    @param timeout: Maximum time hermesVR is allowed to run in seconds
    @type timeout: int

    @return: unseq used, text output of hermesVR, unix return code
    @rtype: (int, str, int)

    """
    #-- Create running directory
    vrDir = tempDir + 'hermesvr/'
    redDir = tempDir + 'hermesvr/reduced/'
    delRedDir, delVrDir = False, False
    if not os.path.isdir(vrDir):
        delVrDir = True
        os.mkdir(vrDir)
    if not os.path.isdir(redDir):
        delRedDir = True
        os.mkdir(redDir)

    #-- Get unseq.
    header = pf.getheader(filename)
    if 'unseq' in kwargs:
        unseq = kwargs.pop('unseq')
        logger.debug('Using provided sequence number: {:}'.format(unseq))
    elif 'unseq' in header:
        unseq = header['unseq']
        logger.debug('Using sequence number from fits header: {:}'.format(unseq))
    else:
        unseq = 1
        logger.debug('Using default sequence number: {:}'.format(unseq))

    #-- Construct filename and copy
    vrFile = redDir+ '{:.0f}_HRF_OBJ_ext.fits'.format(unseq)
    shutil.copyfile(filename, vrFile )
    logger.debug('Coppied input file to: {:}'.format(vrFile))

    #-- Which version to use and update hermesConfig:
    if version == 'release':
        cmd = 'python /STER/mercator/mercator/Hermes/releases/hermes5/pipeline/run/hermesVR.py'
        oldConfig = write_hermesConfig(Nights=tempDir, CurrentNight=vrDir, Reduced=tempDir)
    else:
        cmd = 'python /STER/mercator/mercator/Hermes/trunk/hermes/pipeline/run/hermesVR.py'
        oldConfig = write_hermesConfig(Nights=vrDir, CurrentNight=vrDir, Reduced=tempDir)

    #-- Delete old hermesVR results so it is possible to check if hermesVR completed
    hermesOutput = oldConfig['AnalysesResults'] + "/{:.0f}_AllCCF.data".format(unseq)
    if os.path.isfile(hermesOutput):
        os.remove(hermesOutput)
        logger.debug('Removing existing hermesRV result: {:}'.format(hermesOutput))

    #-- Create command and run hermesVR
    runError = None
    try:
        cmd += " -i {:.0f}".format(unseq)

        if mask_file:
            cmd += " -m {:}".format(mask_file)

        if wvl_file:
            cmd += " -w {:}".format(wvl_file)

        if flatfield and version=='release' or not flatfield and version=='trunk':
            cmd += " -f"

        if vrange == 'large':
            cmd += " -L"
        elif vrange == 'xlarge':
            cmd += " -LL"

        if vrot and version == 'trunk':
            cmd += " -ROT"

        # RUN
        logger.info('running hermesVR {:} version, with command:\n\t{:}'.format(version, cmd))
        out, err, returncode = _subprocess_execute(cmd.split(), timeout)

        # Error Handling
        if returncode == -1:
            logger.warning('Timeout, hermesVR was terminated after' \
                           +' {:.0f}s'.format(timeout))
        elif returncode < 0:
            logger.warning('hermesVR terminated with signal: ' \
                           +' {:.0f}s'.format(abs(returncode)))
        elif not os.path.isfile(hermesOutput):
            logger.warning('hermesVR finished but did not succeed, check log for details')
            returncode = -35
            logger.debug(out)

        if verbose: print(out)

    except Exception as e:
        runError = e

    #-- restore hermesConfig and delete coppied files
    write_hermesConfig(**oldConfig)
    os.remove(vrFile)
    if delRedDir: os.rmdir(redDir)
    if delVrDir: os.rmdir(vrDir)

    #-- Now raise possible errors
    if runError: raise runError

    return unseq, out, returncode

def CCFList(ID,config_dir=None,out_dir=None,mask_file=None,cosmic_clipping=True,flatfield=False):
    """
    Calculate radial velocities for Hermes stars.

    @warning: you have to have the DRS installed locally!

    @warning: this function changes your hermesConfig.xml file, you'd better
    backup it before use
    """
    #-- build the command
    make_list_star(ID)
    cmd = 'python /STER/mercator/mercator/HermesTools/CCFList.py -i %s.list'%(ID)
    if mask_file is not None:
        cmd += ' -m %s'%(os.path.abspath(mask_file))
    if not cosmic_clipping:
        cmd += ' -nc'
    if flatfield:
        cmd += ' -f'

    #-- change configFile
    if config_dir is None:
        config_dir = os.path.expanduser('~/hermesRun')
    if out_dir is None:
        out_dir = os.getcwd()
    out_dir = os.path.abspath(out_dir)

    config_file = r"""<hermes>
<Nights>/STER/mercator/hermes</Nights>
<CurrentNight>/STER/mercator/hermes/20110501</CurrentNight>
<Reduced>/STER/mercator/hermes/20110501/reduced</Reduced>
<AnalysesResults>%s</AnalysesResults>
<DebugPath>%s</DebugPath>
<Raw>raw</Raw>
<ConsoleLogSeverity>debug</ConsoleLogSeverity>
<LogFileName>%s</LogFileName>
<LogFileFormater>%%(asctime)s  ::  %%(levelname)s  ::
  %%(message)s</LogFileFormater>
</hermes>"""%(out_dir,out_dir,os.path.join(out_dir,'hermes.log'))

    #   write to file
    config = open(os.path.join(config_dir,'hermesConfig.xml'),'w')
    config.write(config_file)
    config.close()

    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            logger.error("Child (%s) was terminated by signal %d"%(cmd,-retcode))
        else:
            logger.error("Child (%s) returned %d"%(cmd,retcode))
    except OSError as e:
        logger.error("Execution (%s) failed: %s"%(cmd,e))







#}

class HermesCCF(object):
    """
    Object that holds all information stored in the *_AllCCF.fits files written by hermesVR.
    The original information is stored in the following variables that can be accessed
    directly:

        - filename: The path to the *_AllCCF.fits file
        - header: The header of the original spectra that was the input to hermesVR
        - groups: The summary that hermesVR prints of the 5 preset order groups stored
                  as a rec array
        - ccfs: All cross-correlation functions stored as a array (index 0-54)
        - vr: the velocity scale of the cross correlation functions (array)

    Appart from these attributes, two functions allow you to read the ccf by order or
    a list of orders where the order numbers are the original hermes orders (40-94).
    """

    def __init__(self, filename=None):
        """
        Read the provided file and store all information.

        @param filename: The complete path to the AllCCF.fits file to read.
        @type filename: string
        """
        self.header = None
        self.vr = None
        self.ccfs = None
        self.groups = None

        if filename != None:
            self._read_file(filename)
        else:
            filename = None

    #{ Interaction

    def ccf(self, order):
        """
        Return the ccf function belonging to the requested order. The order number
        is the real order number (40 <= order <= 94).

        @param order: The hermes order of which you want the ccf
        @type order: int

        @return: The cross correlation function of the requested order
        @rtype: array

        @raise IndexError: When the order is out of bounds.
        """
        order = 94 - order

        if order < 0 or order > 54:
            raise IndexError('Order {:.0f} is out of bounds (40 - 94)'.format(94 - order))

        return self.ccfs[order]

    def combine_ccf(self, orders):
        """
        Sum the ccf function belonging to the requested order numbers. The order
        numbers can be given as a list of integers or in a string representation.
        When using the string representation, you can use the ':' sign to give a
        range of order which will include the starting and ending order. Fx the
        following two commands will return the same summed ccf:

        >>> #combine_ccf([54,55,56,57,75,76])
        >>> #combine_ccf('54:57,75,76')

        @param orders: The hermes orders of which you want the summed ccf
        @type orders: list or string

        @return: The summed cross correlation function of the requested orders
        @rtype: array

        @raise IndexError: When an order is out of bounds.
        """

        if type(orders) == str:
            orders = orders.split(',')
            o_list = []
            for o in orders:
                if ':' in o:
                    o = o.split(':')
                    o_list.extend( list(range( int(o[0]), int(o[1])+1)) )
                else:
                    o_list.append( int(o) )
            logger.debug("converted order string: ''{:}'' to list: {:}".format(orders, o_list))
            orders = o_list

        if np.any([(o<40) | (o>94) for o in orders]):
            raise IndexError('Some/All orders {:} are out of bounds (40 - 94)'.format(orders))

        ccf = np.zeros_like(self.vr)
        for o in orders:
            ccf += self.ccf(o)

        return ccf

    #}

    #{ Internal

    def _read_file(self, filename):
        "Read the *_AllCCF.fits file "
        self.filename = filename

        hdu = pf.open(filename)

        #-- header of original observation
        self.header = hdu[0].header

        #-- individual ccfs
        self.ccfs = hdu[0].data.T
        self.vr = hdu[2].data.field('VR')

        #-- data of preselected groups
        self.groups = np.array(hdu[1].data, dtype = hdu[1].data.dtype)

        logger.debug('Read ccf file: {:}'.format(filename))

        return

    def __str__(self):
        "String representation of this object "
        if self.header == None:
            return "< empty Hermes CCF object >"
        else:
            obj, sq = self.header['object'], self.header['unseq']
            return "< CCF of {:} with unseq {:.0f} >".format(obj, sq)

    #}

#{ Administrator functions

def make_data_overview():
    """
    Summarize all Hermes data in a file for easy data retrieval.

    The file is located in one of date data directories (see C{config.py}), in
    subdirectories C{catalogs/hermes/HermesFullDataOverview.tsv}. If it doesn't
    exist, it will be created. It contains the following columns, which are
    extracted from the Hermes FITS headers (except C{filename}:

        1.  UNSEQ
        2.  PROG_ID
        3.  OBSMODE
        4.  BVCOR
        5.  OBSERVER
        6.  OBJECT
        7.  RA
        8.  DEC
        9.  BJD
        10. EXPTIME
        11. PMTOTAL
        12. DATE-AVG
        13. OBJECT
        14. airmass
        15. filename

    This file can most easily be read with the L{ivs.inout.ascii} module and the
    command:

    >>> hermes_file = config.get_datafile(os.path.join('catalogs','hermes'),'HermesFullDataOverview.tsv')
    >>> data = ascii.read2recarray(hermes_file,splitchar='\\t')

    """
    logger.info('Collecting files...')
    #-- all hermes data directories
    dirs = sorted(glob.glob(os.path.join(config.ivs_dirs['hermes'],'20??????')))
    dirs = [idir for idir in dirs if os.path.isdir(idir)]
    obj_files = []
    #-- collect in those directories the raw and relevant reduced files
    for idir in dirs:
        obj_files += sorted(glob.glob(os.path.join(idir,'raw','*.fits')))
        #obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*wavelength_merged.fits')))
        #obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*wavelength_merged_c.fits')))
        #obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*log_merged.fits')))
        #obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*log_merged_c.fits')))

    #-- keep track of what is already in the file, if it exists:
    try:
        overview_file = config.get_datafile('catalogs/hermes','HermesFullDataOverview.tsv')
        #overview_file = config.get_datafile(os.path.join('catalogs','hermes'),'HermesFullDataOverview.tsv')
        overview_data = ascii.read2recarray(overview_file,splitchar='\t')
        outfile = open(overview_file,'a')
        logger.info('Found %d FITS files: appending to overview file %s'%(len(obj_files),overview_file))
    #   if not, begin a new file
    except IOError:
        overview_file = 'HermesFullDataOverview.tsv'
        outfile = open(overview_file,'w')
        outfile.write('#unseq prog_id obsmode bvcor observer object ra dec bjd exptime pmtotal date-avg airmass filename\n')
        outfile.write('#i i a20 >f8 a50 a50 >f8 >f8 >f8 >f8 >f8 a30 >f8 a200\n')
        overview_data = {'filename':[]}
        logger.info('Found %d FITS files: starting new overview file %s'%(len(obj_files),overview_file))

    #-- and summarize the contents in a tab separated file (some columns contain spaces)
    existing_files = np.sort(overview_data['filename'])
    for i,obj_file in enumerate(obj_files):
        sys.stdout.write(chr(27)+'[s') # save cursor
        sys.stdout.write(chr(27)+'[2K') # remove line
        sys.stdout.write('Scanning %5d / %5d FITS files'%(i+1,len(obj_files)))
        sys.stdout.flush() # flush to screen

        #-- maybe this file is already processed: forget about it then
        index = existing_files.searchsorted(obj_file)
        if index<len(existing_files) and existing_files[index]==obj_file:
            sys.stdout.write(chr(27)+'[u') # reset cursor
            continue

        #-- keep track of: UNSEQ, PROG_ID, OBSMODE, BVCOR, OBSERVER,
        #                  OBJECT, RA, DEC, BJD, EXPTIME, DATE-AVG, PMTOTAL,
        #                  airmass and filename (not part of fitsheader)
        contents = dict(unseq=-1,prog_id=-1,obsmode='nan',bvcor=np.nan,observer='nan',
                        object='nan',ra=np.nan,dec=np.nan,
                        bjd=np.nan,exptime=np.nan,pmtotal=np.nan,airmass=np.nan,
                        filename=os.path.realpath(obj_file))
        contents['date-avg'] = 'nan'
        header = pf.getheader(obj_file)
        for key in contents:
            if key in header and key in ['unseq','prog_id']:
                try: contents[key] = int(header[key])
                except: pass
            elif key in header and key in ['obsmode','observer','object','date-avg']:
                contents[key] = str(header[key])
            elif key in header and key in ['ra','dec','exptime','pmtotal','bjd','bvcor']:
                contents[key] = float(header[key])
            elif key=='airmass' and 'telalt' in header:
                if float(header['telalt'])<90:
                    try:
                        contents[key] = airmass.airmass(90-float(header['telalt']))
                    except ValueError:
                        pass

        outfile.write('%(unseq)d\t%(prog_id)d\t%(obsmode)s\t%(bvcor)f\t%(observer)s\t%(object)s\t%(ra)f\t%(dec)f\t%(bjd)f\t%(exptime)f\t%(pmtotal)f\t%(date-avg)s\t%(airmass)f\t%(filename)s\n'%contents)
        outfile.flush()
        sys.stdout.write(chr(27)+'[u') # reset cursor
    outfile.close()
    return overview_file

def _derive_filelocation_from_raw(rawfile,data_type):
    """
    Derive the location of a reduced file from the raw file.
    """
    redfile = rawfile.replace('raw','reduced')
    redfiledir,redfilebase = os.path.dirname(redfile),os.path.basename(redfile)
    base,ext = os.path.splitext(redfilebase)
    if data_type.lower()=='cosmicsremoved_log':
        redfile = os.path.join(redfiledir,'_'.join([base,'ext','CosmicsRemoved','log','merged','c']))+'.fits'
    elif data_type.lower()=='cosmicsremoved_wavelength':
        redfile = os.path.join(redfiledir,'_'.join([base,'ext','CosmicsRemoved','wavelength','merged','c']))+'.fits'
    elif data_type.lower()=='log':
        redfile = os.path.join(redfiledir,'_'.join([base,'ext','log','merged']))+'.fits'
    elif data_type.lower()=='wavelength':
        redfile = os.path.join(redfiledir,'_'.join([base,'ext','wavelength','merged']))+'.fits'
    else:
        redfile = rawfile
    #-- extra check to see if it exists
    if not os.path.isfile(redfile):
        return 'naf'
    return redfile

def _timestamp2jd(timestamp):
    """
    Convert the time stamp from a HERMES FITS 'date-avg' to Julian Date.

    @param timestamp: string from 'date-avg'
    @type timestamp: string
    @return: julian date
    @rtype: float
    """
    date, hour = timestamp.split("T")
    year, month, day = date.split("-")
    hour, minute, second = hour.split(":")
    year   = float(year)
    month  = float(month)
    day    = float(day)
    hour   = float(hour)
    minute = float(minute)
    second = float(second)
    return conversions.convert("CD","JD",(year, month, day, hour, minute, second))

def _timestamp2datetime(timestamp):
    """
    Convert the time stamp from a HERMES FITS 'date-avg' to a datetime object.

    @param timestamp: string from 'date-avg'
    @type timestamp: string
    @return: datetime object
    @rtype: datetime
    """
    if timestamp=='nan':
        timestamp = '1000-01-01T00:00:00'

    timestamp = [int(i) for i in ' '.join(' '.join(' '.join(' '.join(timestamp.split('-')).split('T')).split(':')).split('.')).split()]
    #-- only the day is given, make sure to switch it to mid-day
    if len(timestamp)==3:
        timestamp += [12,0,0]
    return datetime.datetime(*timestamp)

def _etree_to_dict(t):
    """
    Convert a xml tree to a dictionary.

    @param t: the xml tree
    @type t: ElementTree root

    @return: the xml tree as a dictionary
    @rtype: dict
    """
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(_etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
              d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d

def _subprocess_execute(command, time_out=100):
    """executing the command with a watchdog"""

    # launching the command
    c = subprocess.Popen(command, stdout=subprocess.PIPE)

    # now waiting for the command to complete
    t = 0
    while t < time_out and c.poll() is None:
        time.sleep(1)  # (comment 1)
        t += 1

    # there are two possibilities for the while to have stopped:
    if c.poll() is None:
        # in the case the process did not complete, we kill it
        c.terminate()
        # and fill the return code with some error value
        out, err = c.communicate()
        returncode = -1

    else:
        # in the case the process completed normally
        out, err = c.communicate()
        returncode = c.poll()

    return out, err, returncode

#}

if __name__=="__main__":

    #from ivs.aux import loggers
    #logger = loggers.get_basic_logger(clevel='debug')

    #filename = '/home/jorisv/sdB/Uli/hermesvr/reduced/2451576_HRF_OBJ_ext_CosmicsRemoved.fits'
    #wvl_file = '/home/jorisv/sdB/Uli/hermesvr/reduced/2451577'

    #unseq,output, error = run_hermesVR(filename, wvl_file=wvl_file, verbose=False, version='trunk', timeout=500)

    #data = read_hermesVR_velocities(unseq=[unseq])

    #print(data['vrad'], ' +- ', data['vraderr'])

    #import time
    #import sys
    import doctest
    #import shutil
    import pylab as pl

    #if len(sys.argv[1:])==0:
    doctest.testmod()
    pl.show()

    #elif sys.argv[1].lower()=='update':
        #logger = loggers.get_basic_logger()

        #while 1:
            #source = make_data_overview()
            #destination = '/STER/mercator/hermes/HermesFullDataOverview.tsv'
            #if os.path.isfile(destination):
                #original_size = os.path.getsize(destination)
                #logger.info("Original file size: %.6f MB"%(original_size/1.0e6))
            #else:
                #logger.info('New file will be created')
            #new_size = os.path.getsize(source)
            #logger.info("New file size: %.6f MB"%(new_size/1.0e6))
            #os.system('cp %s %s'%(source,destination))
            #logger.info('Copied %s to %s'%(source,destination))

            #logger.info('Going to bed know... see you tomorrow!')
            #time.sleep(24*3600)
            #logger.info('Rise and shine!')


    #elif sys.argv[1].lower()=='copy':
        #while 1:
            #source = '/STER/pieterd/IVSDATA/catalogs/hermes/HermesFullDataOverview.tsv'
            #destination = '/STER/mercator/hermes/HermesFullDataOverview.tsv'
            #if os.path.isfile(destination):
                #original_size = os.path.getsize(destination)
                #logger.info("Original file size: %.5f kB"%(original_size/1000.))
            #else:
                #logger.info('New file will be created')
            #new_size = os.path.getsize(source)
            #logger.info("New file size: %.5f kB"%(new_size/1000.))
            #shutil.copy(source,destination)
            #logger.info('Copied %s to %s'%(source,destination))
            #time.sleep(24*3600)

    #else:
        #logger = loggers.get_basic_logger()
        #for target in sys.argv[1:]:
            #make_list_star(target)
