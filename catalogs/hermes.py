# -*- coding: utf-8 -*-
"""
Interface the spectra from the Hermes spectrograph.

Developer's to-do: remove doubles entries from the catalog, and I{derive}
the path names of the different datafiles (?).

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

>>> from ivs.io import fits
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
>>> from ivs.spectra import model         # to apply doppler shift
>>> from ivs.units import conversions     # to convert to velocity space

Then get the spectra again, but now not in log scale. Make the same selection
on exposure time.

>>> mydata = search('HD170580',data_type='cosmicsremoved_wavelength')
>>> myselection = mydata[mydata['exptime']>500]

Then read them all in and shift according to the barycentric velocity. Also,
convert the spectra to velocity space afterwards for plotting purposes.

>>> for rv,fname in zip(myselection['bvcor'],myselection['filename']):
...     wave,flux = fits.read_spectrum(fname)
...     wave_shifted = model.doppler_shift(wave,rv)
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
>>> print [(i,list(data['observer']).count(i)) for i in set(data['observer'])]
[('Robin Lombaert', 4), ('Steven Bloemen', 6), ('Pieter Degroote', 25), ('Michel Hillen', 3)]



"""
import re
import sys
import glob
import os
import logging
import datetime
import subprocess
import numpy as np
import pyfits

from ivs.catalogs import sesame
from ivs.io import ascii
from ivs.aux import loggers
from ivs.observations import airmass
from ivs.observations.barycentric_correction import helcorr
from ivs.units import conversions
from ivs import config

logger = logging.getLogger("CAT.HERMES")
logger.addHandler(loggers.NullHandler)

#{ User functions

def search(ID=None,time_range=None,data_type='cosmicsremoved_log',radius=1.,filename=None):
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
    fits header, small letters are calculated values):
        
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
    
    >>> telalts = [pyfits.getheader(fname)['telalt'] for fname in myselection['filename']]
    
    Search for all data of HD50230 taken in the night of 22 September 2009:
    
    >>> data = hermes.search('HD50230',time_range='2009-9-22')
    
    Or within an interval of a few days:
    
    >>> data = hermes.search('HD50230',time_range=('2009-9-23','2009-9-30'))
    
    Search for all data observed in a given night:
    
    >>> data = hermes.search(time_range='2009-9-22')
    
    B{Warning:} the heliocentric correction is not calculated when no ID is given,
    so make sure it is present in the header if you need it, or calculate it yourself.
    
    @param ID: ID of the star, understandable by SIMBAD
    @type ID: str
    @param time_range: range of dates to confine the search to
    @type time_range: tuple strings of type '2009-09-23T04:24:35.712556' or '2009-09-23'
    @type time_range: tuple strings
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
    data = ascii.read2recarray(config.get_datafile(os.path.join('catalogs','hermes'),'HermesFullDataOverview.tsv'),splitchar='\t')
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
        
    #-- if some data is found, we check if the C{data_type} string is contained
    #   with the file's name. If not, we remove it.
    if np.any(keep):
        data = data[keep]
    
        if data_type is not None:
            data_type == data_type.lower()
            keep = np.array([(data_type in ff.lower() and True or False) for ff in data['filename']])
            data = data[keep]
            seqs = sorted(set(data['unseq']))
        logger.info('%s: Found %d spectra (data type=%s with unique unseqs)'%(ID,len(seqs),data_type))
    else:
        data = data[:0]
        logger.info('%s: Found no spectra'%(ID))
    
    #-- we now check if the barycentric correction was calculated properly.
    #   If not, we calculate it here, but only if the object was found in
    #   SIMBAD. Else, we have no information on the ra and dec (if bvcorr was
    #   not calculated, ra and dec are not in the header).
    for obs in data:
        if info:
            jd  = _timestamp2jd(obs['date-avg'])
            bvcorr, hjd = helcorr(ra, dec, jd)
        else:
            break
        if np.isnan(obs['bvcor']): obs['bvcor'] = bvcorr
        if np.isnan(obs['bjd']):   obs['bjd'] = hjd
    
    #-- do we need the information as a file, or as a numpy array?
    if filename is not None:
        ascii.write_array(data,filename,auto_width=True,header=True)
    else:
        return data
    
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
    @param strength: strenght of the line (1-normalised flux minimum)
    @type strength: 1D numpy array
    """
    data = np.array([wavelength,depth]).T
    hdulist = pyfits.PrimaryHDU(data1)
    hdulist.writeto(filename)
    hdulist.close()


#}

#{ Hermes DRS wrappers

def CCFList(ID,config_dir=None,out_dir=None,mask_file=None,cosmic_clipping=True,flatfield=False):
    """
    Calculate radial velocities for Hermes stars.
    
    WARNING: you have to have the DRS installed locally!
    WARNING: this function changes your hermesConfig.xml file, you'd better
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
    except OSError, e:
        logger.error("Execution (%s) failed: %s"%(cmd,e))
    
    

    
    


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
    
    This file can most easily be read with the L{ivs.io.ascii} module and the
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
        obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*wavelength_merged.fits')))
        obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*wavelength_merged_c.fits')))
        obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*log_merged.fits')))
        obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*log_merged_c.fits')))
    
    #-- keep track of what is already in the file, if it exists:
    try:
        overview_file = os.path.join('/STER/100/pieterd/IVSDATA/catalogs/hermes','HermesFullDataOverview.tsv')
        #overview_file = config.get_datafile(os.path.join('catalogs','hermes'),'HermesFullDataOverview.tsv')
        overview_data = ascii.read2recarray(overview_file,splitchar='\t')
        outfile = open(overview_file,'a')
        logger.info('Found %d FITS files: appending to overview file %s'%(len(obj_files),overview_file))
    #   if not, begin a new file
    except IOError:
        overview_file = os.path.join('/STER/100/pieterd/IVSDATA/catalogs/hermes','HermesFullDataOverview.tsv')
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
        header = pyfits.getheader(obj_file)
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
    timestamp = [float(i) for i in ' '.join(' '.join(' '.join(timestamp.split('-')).split('T')).split(':')).split()]
    #-- only the day is given, make sure to switch it to mid-day
    if len(timestamp)==3:
        timestamp += [12,0,0]
    return datetime.datetime(*timestamp)

#}

if __name__=="__main__":
    import time
    import sys
    import doctest
    import shutil
    import pylab as pl
    
    if len(sys.argv[1:])==0:
        doctest.testmod()
        pl.show()
    
    elif sys.argv[1].lower()=='update':
        logger = loggers.get_basic_logger()
        
        while 1:
            make_data_overview()
            
            source = '/STER/100/pieterd/IVSDATA/catalogs/hermes/HermesFullDataOverview.tsv'
            destination = '/STER/mercator/hermes/HermesFullDataOverview.tsv'
            if os.path.isfile(destination):
                original_size = os.path.getsize(destination)
                logger.info("Original file size: %.6f MB"%(original_size/1.0e6))
            else:
                logger.info('New file will be created')
            new_size = os.path.getsize(source)
            logger.info("New file size: %.6f MB"%(new_size/1.0e6))
            os.system('cp %s %s'%(source,destination))
            logger.info('Copied %s to %s'%(source,destination))
            
            logger.info('Going to bed know... see you tomorrow!')
            time.sleep(24*3600)
            logger.info('Rise and shine!')
            
            
    elif sys.argv[1].lower()=='copy':
        while 1:
            source = '/STER/100/pieterd/IVSDATA/catalogs/hermes/HermesFullDataOverview.tsv'
            destination = '/STER/mercator/hermes/HermesFullDataOverview.tsv'
            if os.path.isfile(destination):
                original_size = os.path.getsize(destination)
                logger.info("Original file size: %.5f kB"%(original_size/1000.))
            else:
                logger.info('New file will be created')
            new_size = os.path.getsize(source)
            logger.info("New file size: %.5f kB"%(new_size/1000.))
            shutil.copy(source,destination)
            logger.info('Copied %s to %s'%(source,destination))
            time.sleep(24*3600)
            
    else:
        logger = loggers.get_basic_logger()
        for target in sys.argv[1:]:
            make_list_star(target)
    
