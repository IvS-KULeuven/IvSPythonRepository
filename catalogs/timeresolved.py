# -*- coding: utf-8 -*-
"""
Retrieve epoch photometry from the internet.

Author: Joris De Ridder & Pieter Degroote

Error messages are written to the logger "timeresolved".
"""

from __future__ import with_statement
import httplib
import logging
import os
import pyfits
import numpy as np
from ivs.misc import loggers
from ivs.catalogs import sesame
from ivs.io import fits
from ivs import config
        
logger = logging.getLogger("timeresolved")
logger.addHandler(loggers.NullHandler)

def getHipData(ID,outputFileName=None):

    """
    Retrieve Hipparcos epoch photometry from the ESA website.
    The time series together with some header information is stored in a record
    array and dictionary, and optionally written in a specified file.
    
    The time points are given in barycentric Julian Date and are corrected
    for the offset of 2440000.0 in the original files, but B{only in the output
    record array}. The output files display the B{original contents}.
       
    For more information:
    C{http://www.rssd.esa.int/SA-general/Projects/Hipparcos/CATALOGUE_VOL1/sect2_05.ps.gz}
    
    Example:
    
    >>> data,header = getHipData(1234)
    >>> data = data[data['q_mag'] <= 2]         # keep only the good points
    
    To write the retrieved data to a file:
    
    >>> data, header = getHipData(1234 , "myfile.txt")
    
    To store the different columns in separate arrays:
    
    >>> data, header = getHipData(1234)
    >>> time = data['time']
    >>> magnitude = data['mag']
    >>> errorbar = data['e_mag']
    >>> qualityflag = data['q_mag']
    
    
    @param ID: identification of the star: if you give an integer or string that
    can be converted to an integer, it is assumed to be the hipparcos number of
    the star.  E.g. 1234 or "1234". If it is not an integer, the star will
    be resolved via sesame to get the HIP number if possible
    @type ID: integer or string
    @param outputFileName: the name of the file that will be created
                           to save the Hipparcos time series
    @type outputFileName: string
    @return: record array with fields time, mag, e_mag (errorbar), 
             q_mag (quality flag), and a dictionary containing the 
             header information. The header dictionary is of style
             {'HH14': ('A', 'Annex flag (light curves)'), ...}
    @rtype: rec array, dict
    """

    server = "www.rssd.esa.int"
    webpage = "/hipparcos_scripts/HIPcatalogueSearch.pl?hipepId="
    
    # Resolve the name if necessary (i.e., if it's not a HIP number). If the 
    # star has no HIP number, log an error and return None
    try:
        ID = int(ID)
    except ValueError:
        info = sesame.search(ID,db='S')
        IDs = [alias for alias in info['alias'] if 'HIP' in alias]
        if len(IDs)!=1:
            logger.error("Data retrieval for %s not possible. Reason: no HIP number resolved" % (ID))
            return
        hipnr = IDs[0].split(' ')[1]
    
    # Connect to the website, en retrieve the wanted webpage
    
    conn = httplib.HTTPConnection(server)
    conn.request("GET", webpage + str(hipnr))
    response = conn.getresponse()
    if response.reason != "OK":
        logger.error("Data retrieval for HIP%s not possible. Reason: %s" % (str(hipnr), response.reason))
        return
    else:
        logger.info("Data retrieval for HIP%s: OK" % str(hipnr))
        
    contents = response.read()
    conn.close()
    
    # Parse the webpage, to remove the html codes (line starts with <").
    # Put a "#" in front of the header information, and format nicely.
    # Write to the output file if asked for.
    
    data = []
    header = {}
    
    if outputFileName:
        outputFile = open(outputFileName,'w')
    
    for line in contents.split('\n'):
        if line == "": continue
        if not line.startswith("<"):
            line = line.replace("|", " ").replace("\r", "")
            
            # This is the header
            
            if not line[0].isdigit():
                sline = line.split(':')
                
                # Only keep header entries of the style 
                # "key: value information" in the dictionary
                
                if len(sline)==2:
                    key,info = sline
                    info = info.split()
                    header[key] = (info[0]," ".join(info[1:]))
                if outputFileName:
                    line = "# " + line
                    
            # This is the real contents
            
            else:
                data.append(tuple(line.split()))
            if outputFileName:
                outputFile.write(line + "\n")
    if outputFileName:
        outputFile.close()
    
    # Make a record array.
    # Choose the header names to be in the VizieR style.
    
    dtypes = [('time','>f8'),('mag','>f8'),('e_mag','>f8'),('q_mag','i')]
    data = np.rec.array(data,dtype=dtypes)
    
    # Fix the time offset
    
    data['time'] += 2440000.0
    
    return data,header





def getCoRoTData(ID,channel='sismo'):
    """
    Retrieve CoRoT timeseries from a local data repository.
    
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
    
    if channel=='sismo':
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
            if header['starname']==ID or header['corotid']=='%s'%(ID):
                times,flux,error,flags = fits.read_corot(catfile)
                data.append([times,flux,error,flags])
        #-- now make a record array and sort according to times
        data = np.hstack(data)        
        data = np.rec.fromarrays(data,dtype=[('HJD','>f8'),('flux','>f8'),('e_flux','>f8'),('flag','i')])
        sa = np.argsort(data['HJD'])
        return data[sa],header
        
        













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
        # this should resolve the GENEVA name
        code = _geneva_name_resolver(ID=ID)
        
    catfile = config.get_datafile('catalogs/p7','p7photometry.fits')
    ff = pyfits.open(catfile)
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
        











if __name__=="__main__":
    import doctest
    doctest.testmod()
    