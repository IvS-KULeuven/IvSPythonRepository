"""
Retrieve data from the P7 photometer.
"""
import logging
import os
import astropy.io.fits as pf
import numpy as np
from ivs.aux import loggers
from ivs.catalogs import sesame
from ivs.inout import fits
from ivs import config
        
logger = logging.getLogger("")
logger.addHandler(loggers.NullHandler)



def getP7Data(ID=None,code=None,include_nans=True):
    """
    Extract P7 timeseries from the catalog.
    
    WARNING: only B{very} few target ID's can be resolved (HD,HIP,SAO and that's
    about it)
    
    WARNING: there could be nan's in the data somewhere. If you don't want
    nan's anywhere, set 'include_nans' to False.
    
    @param ID: target ID (limited!)
    @type ID: str
    @param code: target's GENEVA code (e.g. 100180642 for HD180642)
    @type code: int
    @return: record array containing times (HJD) and corresponding GENEVA mags,
    and a dictionary with header information (only source=P7)
    @rtype: np record array,dict
    """
    if ID is not None:
        if not 'HD' in ID or not 'SAO' in ID or not 'HIC' in ID:
            info = sesame.search(ID)
            print info
            if 'alias' in info:
                for alias in info['alias']:
                    if 'HD' in alias:
                        ID = alias
                        break
                    if 'SAO' in alias:
                        ID = alias
                        break
                    if 'HIC' in alias:
                        ID = alias
                        break
        # this should resolve the GENEVA name
        code = _geneva_name_resolver(ID=ID)
        
    catfile = config.get_datafile('catalogs/p7','p7photometry.fits')
    ff = pf.open(catfile)

    valid = ff[1].data.field('CODE')==code
    hjd = ff[1].data.field('HJD')[valid]
    U = ff[1].data.field('U')[valid]
    B = ff[1].data.field('B')[valid]
    B1 = ff[1].data.field('B1')[valid]
    B2 = ff[1].data.field('B2')[valid]
    V = ff[1].data.field('V')[valid]
    V1 = ff[1].data.field('V1')[valid]
    G = ff[1].data.field('G')[valid]
    ff.close()
    
    data = np.rec.fromarrays([hjd,U,B,B1,B2,V,V1,G],names='HJD,U,B,B1,B2,V,V1,G')
    
    logger.info('Retrieved %d photometric points from P7'%(len(data)))
    
    if not include_nans:
        nans = np.isnan(data['HJD'])
        for name in data.dtype.names:
            nans = nans | np.isnan(data[name])
        data = data[-nans]
        logger.info('Keeping %d photometric points without "NaN" from P7'%(len(data)))
    
    return data,{'source':'P7'}


def _geneva_name_resolver(code=None,ID=None):
    """
    Resolve the GENEVA object codes.
    
    @param ID: target ID (not implemented yet)
    @type ID: str
    @param code: target's GENEVA code (e.g. 100180642 for HD180642)
    @type code: int
    @return: Geneva code or target ID (dependent on input)
    @rtype: str
    """
    if code is not None:
        if not isinstance(code,str):
            code = '%+09d'%(code)
        if code[1:3]=='10':   ID = 'HD'+code[3:]
        elif code[1:3]=='15': ID = 'SAO'+code[3:]
        elif code[1:3]=='16': ID = 'HIC'+code[3:]
        elif code[1:3]=='17': ID = 'PPM'+code[3:]
        elif code[:3]=='-00': ID = 'BD%s%s %s'%(code[0],code[3:5],code[5:])
        elif code[:3]=='+00': ID = 'BD%s%s %s'%(code[0],code[3:5],code[5:])
        return ID
    
    elif ID is not None:
        if    'HD' in ID: code = '10%06d'%(int(ID.split('HD')[1]))
        elif 'SAO' in ID: code = '15%06d'%(int(ID.split('SAO')[1]))
        elif 'HIC' in ID: code = '16%06d'%(int(ID.split('HIC')[1]))
        elif 'PPM' in ID: code = '17%06d'%(int(ID.split('PPM')[1]))
        else: code = ID
        return  int(code)
        




