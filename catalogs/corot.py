"""
Retrieve CoRoT data from a local data repository.

check out http://nsted.ipac.caltech.edu/
"""

import logging
import os
import urllib
import astropy.io.fits as pf
import numpy as np
from ivs.aux import loggers
from ivs.catalogs import sesame
from ivs.catalogs import vizier
from ivs.inout import fits
from ivs import config

logger = logging.getLogger("CAT.COROT")
logger.addHandler(loggers.NullHandler())

def get_sismo_data(ID):
    """
    Retrieve CoRoT timeseries from a local data repository.
    
    The output record array has fields 'HJD', 'flux', 'e_flux', 'flag'.
    
    @param ID: ID of the target: either an integer (CoRoT ID), an SIMBAD-recognised
    target name, or a valid CoRoT FITS file
    @type ID: int or str
    @return: data, header
    @rtype: numpy recarray, dict
    """
    #-- data on one target can be spread over multiple files: collect the
    #   data
    data = []
        
    if isinstance(ID,str) and os.path.isfile(ID):
        header = pf.getheader(ID)
        times,flux,error,flags = fits.read_corot(ID)
        data.append([times,flux,error,flags])
    else:
        #-- resolve the target's name: it's either a target name or CoRoT ID.
        try:
            ID = int(ID)
        except ValueError:
            info = sesame.search(ID,db='S')
            IDs = [alias for alias in info['alias'] if 'HD' in alias]
            if len(IDs)!=1:
                logger.error("Data retrieval for %s not possible. Reason: no HD number resolved" % (ID))
                return
            ID = IDs[0]
        #-- collect the files containing data on the target
        catfiles = config.glob((os.sep).join(['catalogs','corot','sismo']),'*.fits')
        for catfile in catfiles:
            try:
                header = pf.getheader(catfile)
            except IOError:
                continue
            if header['starname']==ID or header['corotid'].replace(' ','')=='%s'%(ID):
                times,flux,error,flags = fits.read_corot(catfile)
                data.append([times,flux,error,flags])
    #-- now make a record array and sort according to times
    if not data:
        raise ValueError('target {0} not in offline CoRoT data repository'.format(ID))
    data = np.hstack(data)        
    data = np.rec.fromarrays(data,dtype=[('HJD','>f8'),('flux','>f8'),('e_flux','>f8'),('flag','i')])
    sa = np.argsort(data['HJD'])
    return data[sa],header

def get_exo_data(ID,type_data='white'):
    """
    Retrieve CoRoT timeseries from a remote data repository.
    
    The output record array has fields 'HJD', 'flux', 'e_flux', 'flag'.
    
    @param ID: ID of the target: either an integer (CoRoT ID), an SIMBAD-recognised
    target name, or a valid CoRoT FITS file
    @type ID: int or str
    @param type_data: 'white' or 'colors' (if available!)
    @type type_data: str
    @return: data, header
    @rtype: numpy recarray, dict
    """
    cat,units,comms = get_exo_catalog()
    header = None
    data = []
    if isinstance(ID,str) and os.path.isfile(ID):
        header = pf.getheader(ID)
        times,flux,error,flags = fits.read_corot(ID)
        data.append([times,flux,error,flags])
    else:
        #-- resolve the target's name: it's either a target name or CoRoT ID.
        try:
            ID = int(ID)
        except ValueError:
            logger.error("Data retrieval for %s not possible. Reason: ID not resolved" % (ID))
            return None
        #-- collect the files containing data on the target
        for filename in cat['FileName'][cat['CoRoT']==ID]:
            url = urllib.URLopener()
            filen,msg = url.retrieve(filename)
            try:
                header = pf.getheader(filen)
            except IOError:
                continue
            times,flux,error,flags = fits.read_corot(filen,type_data=type_data)
            url.close()
            data.append([times,flux,error,flags])
        
    #-- now make a record array and sort according to times
    if not data:
        raise ValueError('target {0} not in online CoRoT data repository'.format(ID))
    data = np.hstack(data)        
    data = np.rec.fromarrays(data,dtype=[('HJD','>f8'),('flux','>f8'),('e_flux','>f8'),('flag','i')])
    sa = np.argsort(data['HJD'])
    return data[sa],header


def get_exo_catalog():
    """
    Retrieve the locally saved CoRoT exoplanet database.
    """
    exofile = config.get_datafile('catalogs/corot/exo','exo.tsv')
    data,units,comms = vizier.tsv2recarray(exofile)
    return data,units,comms
    
def resolve(corot_id):
    """
    Convert a CoRoT ID to ra,dec.
    
    @param corot_id: CoRoT exoplanet identification number
    @type corot_id: int
    @return: RA, DEC (degrees)
    @rtype: float,float
    """
    exofile = config.get_datafile('catalogs/corot/exo','exo.tsv')
    data,units,comms = vizier.tsv2recarray(exofile)
    index = np.argmin(np.abs(data['CoRoT']-corot_id))
    if data['CoRoT'][index]-corot_id==0:
        logger.info('CoRoT %s resolved to RA=%s, DEC=%s'%(corot_id,data[index]['_RAJ2000'],data[index]['_DEJ2000']))
        return data[index]['_RAJ2000'],data[index]['_DEJ2000']
    else:
        logger.info('CoRoT %s not resolved by CoRoT catalog')
