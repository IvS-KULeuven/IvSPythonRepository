"""
Retrieve light curves from the Kepler satellite mission.
"""
import urllib
import os
import astropy.io.fits as pf
import logging
import numpy as np
from ivs.catalogs import mast
from ivs.units import conversions
from ivs.inout import ascii
from ivs.inout import fits
from ivs import config

logger = logging.getLogger("CAT.KEPLER")

def download_light_curve(KIC,directory=''):
    """
    Download a light curve file from the data archive.
    
    @param KIC: kic number
    @type KIC: integer
    @param directory: directory to save data to, defaults to cwd
    @type directory: string
    @return: list of filenames
    @rtype: list of str
    """
    #-- get info on what is available
    results,units,comms = mast.search('kepler/data_search',ktc_kepler_id=KIC)
    #-- reconstruct the links
    links = []
    for name in results['Dataset Name']:
        num = '%09d'%(int(KIC))
        folder = num[:4]
        links.append("http://archive.stsci.edu/pub/kepler/lightcurves/%s/%s/%s_llc.fits"%(folder,num,name.lower()))
    filenames = []
    #-- and download the files
    logger.info("Found %d public/proprietary light curves for KIC%s"%(len(links),KIC))
    for base_url in links:
        url = urllib.URLopener()
        filename = os.path.basename(base_url)
        if directory:
            filename = os.path.join(directory,filename)
        try:
            filen,msg = url.retrieve(base_url,filename=filename)
        except IOError:
            logger.info('... %s is proprietary'%(base_url))
        else:
            logger.info('... downloaded %s'%(filename))
            filenames.append(filen)
        url.close()
        
    return filenames

def get_data(KIC):
    """
    Retrieve Kepler timeseries from a remote data repository.
    
    Fields are 'HJD','flux','e_flux','bkg','quarter'.
    
    @param KIC: kic number or list of filenames
    @type KIC: integer or list
    @return: data, header
    @rtype: recarray, dict
    """
    times = []
    flux = []
    e_flux = []
    background = []
    quarter = []
    if isinstance(KIC,str) or isinstance(KIC,int):
        filenames = download_light_curve(KIC)
    else:
        filenames = KIC
    for filename in filenames:
        header = pf.getheader(filename)
        data = fits.read2recarray(filename,ext='LIGHTCURVE')
        times.append(data['TIME']+2454833.)
        flux.append(data['SAP_FLUX'])
        e_flux.append(data['SAP_FLUX_ERR'])
        background.append(data['SAP_BKG'])
        quarter.append(np.ones(len(data))*header['quarter'])
    data = np.rec.fromarrays([np.hstack(times),np.hstack(flux),np.hstack(e_flux),
                              np.hstack(background),np.hstack(quarter)],
                              names=['HJD','flux','e_flux','bkg','quarter'])
    return data,header
   

def systematics(units='muHz'):
    """
    Return a list of known systematic effects from the Kepler satellite.
    """
    ffile = config.get_datafile('catalogs/kepler','systematics.dat')
    systems = ascii.read2recarray(ffile)
    systems['frequency'] = conversions.nconvert(systems['unit'],units,systems['frequency'])
    systems['e_frequency'] = conversions.nconvert(systems['unit'],units,systems['e_frequency'])
    systems['w_frequency'] = conversions.nconvert(systems['unit'],units,systems['w_frequency'])
    return systems
    
