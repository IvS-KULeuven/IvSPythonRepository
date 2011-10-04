"""
Retrieve CoRoT data from a local data repository.

check out http://nsted.ipac.caltech.edu/
"""

import logging
import os
import pyfits
import numpy as np
from ivs.aux import loggers
from ivs.catalogs import sesame
from ivs.io import fits
from ivs import config

logger = logging.getLogger("CAT.COROT")
logger.addHandler(loggers.NullHandler)

def get_sismo_data(ID):
    """
    Retrieve CoRoT timeseries from a local data repository.
    
    The output record array has fields 'HJD', 'flux', 'e_flux', 'flag'.
    
    @param ID: ID of the target: either an integer (CoRoT ID) or an official
    target name
    @type ID: int or str
    @param channel: type of CoRoT data (sismo or exo)
    @type channel: str ('sismo' or 'exo')
    @return: data, header
    @rtype: numpy recarray, dict
    """
    
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
    
    #-- data on one target can be spread over multiple files: collect the
    #   data
    data = []
    #-- collect the files containing data on the target
    catfiles = config.glob((os.sep).join(['catalogs','corot','sismo']),'*.fits')
    for catfile in catfiles:
        try:
            header = pyfits.getheader(catfile)
        except IOError:
            continue
        if header['starname']==ID or header['corotid'].replace(' ','')=='%s'%(ID):
            times,flux,error,flags = fits.read_corot(catfile)
            data.append([times,flux,error,flags])
    #-- now make a record array and sort according to times
    data = np.hstack(data)        
    data = np.rec.fromarrays(data,dtype=[('HJD','>f8'),('flux','>f8'),('e_flux','>f8'),('flag','i')])
    sa = np.argsort(data['HJD'])
    return data[sa],header
        
        