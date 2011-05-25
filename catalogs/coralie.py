# -*- coding: utf-8 -*-
"""
Interface to s1dA spectra from the Coralie spectrograph.
"""
import re
import sys
import glob
import os
import logging
import numpy as np
import pyfits
from ivs.catalogs import sesame
from ivs.io import ascii
from ivs.misc import loggers
from ivs import config

logger = logging.getLogger("CAT.CORALIE")
logger.addHandler(loggers.NullHandler)

def search(ID,radius=1.,filename=None):
    """
    Retrieve datafiles from the Coralie catalogue.
    
    We search on coordinates, pulled from SIMBAD. If the star ID is not
    recognised, a string search is performed to match the 'targ name' field in the
    FITS headers.
    
    Only the s1d_A data are searched.
    
    @param ID: ID of the star, understandable by SIMBAD
    @type ID: str
    @param radius: search radius around the coordinates
    @type radius: 1
    @param filename: write summary to outputfile if not None
    @type filename: str
    @return: record array with summary information on the observations, as well
    as their location (column 'filename')
    @rtype: numpy rec array
    """
    data = ascii.read2recarray(config.get_datafile(os.path.join('catalogs','coralie'),'CoralieFullDataOverview.tsv'),splitchar='\t')
    info = sesame.search(ID)
    if info:
        ra,dec = info['jradeg'],info['jdedeg']
        keep = np.sqrt((data['ra']-ra)**2 + (data['dec']-dec)**2) < radius/60.
    else:
        keep = [((re.compile(ID).search(objectn) is not None) and True or False) for objectn in data['object']]
        keep = np.array(keep)
    
    data = data[keep]
    
    logger.info('Found %d spectra'%(len(data)))
    
    if filename is not None:
        ascii.write_array(data,filename,auto_width=True,header=True)
    else:
        return data
    
    


def make_data_overview():
    """
    Summarize all COralie data in a file for easy data retrieval.
    """
    #-- all Coralie data directories
    obj_files = []
    for root,dirs,files in os.walk(config.ivs_dirs['coralie']):
        for name in files:
            if 's1d_A' in name:
                obj_files.append(os.path.join(root,name))
    
    #-- and summarize the contents in a tab separated file (some columns contain spaces)
    outfile = open('CoralieFullDataOverview.tsv','w')
    outfile.write('#unseq prog_id obsmode bvcor observer object ra dec bjd exptime date-avg filename\n')
    outfile.write('#i i a20 >f8 a50 a50 >f8 >f8 >f8 >f8 a30 a200\n')
    for i,obj_file in enumerate(obj_files):
        print i,len(obj_files)
        #-- keep track of: UNSEQ, BJD, BVCOR, OBSERVER, RA, DEC , PROG_ID, OBSMODE, EXPTIME, DATE-AVG, OBJECT and filename
        contents = dict(unseq=-1,prog_id=-1,obsmode='CORALIE',bvcor=0,observer='nan',
                        object='nan',ra=np.nan,dec=np.nan,
                        bjd=np.nan,exptime=np.nan,
                        filename=os.path.realpath(obj_file))
        contents['date-avg'] = 'nan'
        try:
            header = pyfits.getheader(obj_file)
        except:
            continue
        if 'ESO DRS TEXP' in header:     contents['exptime'] = float(header['ESO DRS TEXP'])
        if 'ESO DRS OBSERVER' in header: contents['observer'] = header['ESO DRS OBSERVER']
        if 'ESO DRS BJD' in header:      contents['bjd'] = float(header['ESO DRS BJD'])
        if 'ESO DRS FDATE' in header:    contents['date-avg'] = header['ESO DRS FDATE']
        if 'ESO TEL TARG ALPHA' in header:    contents['ra']  = float(header['ESO TEL TARG ALPHA'])
        if 'ESO TEL TARG DELTA' in header:    contents['dec'] = float(header['ESO TEL TARG DELTA'])
        if 'ESO OBS TARG NAME' in header:     contents['object'] = header['ESO OBS TARG NAME']
            
        outfile.write('%(unseq)d\t%(prog_id)d\t%(obsmode)s\t%(bvcor)f\t%(observer)s\t%(object)s\t%(ra)f\t%(dec)f\t%(bjd)f\t%(exptime)f\t%(date-avg)s\t%(filename)s\n'%contents)
        outfile.flush()
    outfile.close()
