# -*- coding: utf-8 -*-
"""
Interface the spectra from the Hermes spectrograph.
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

logger = logging.getLogger("CAT.HERMES")
logger.addHandler(loggers.NullHandler)

def search(ID,data_type='cosmicsremoved',radius=1.,filename=None):
    """
    Retrieve datafiles from the Hermes catalogue.
    
    We search on coordinates, pulled from SIMBAD. If the star ID is not
    recognised, a string search is performed to match the 'object' field in the
    FITS headers.
    
    @param ID: ID of the star, understandable by SIMBAD
    @type ID: str
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
    data = ascii.read2recarray(config.get_datafile(os.path.join('catalogs','hermes'),'HermesFullDataOverview.tsv'),splitchar='\t')
    info = sesame.search(ID)
    
    if info:
        ra,dec = info['jradeg'],info['jdedeg']
        keep = np.sqrt((data['ra']-ra)**2 + (data['dec']-dec)**2) < radius/60.
    else:
        keep = [((re.compile(ID).search(objectn) is not None) and True or False) for objectn in data['object']]
        keep = np.array(keep)
    
    if np.any(keep):
        data = data[keep]
    
        if data_type is not None:
            data_type == data_type.lower()
            keep = np.array([(data_type in ff.lower() and True or False) for ff in data['filename']])
            data = data[keep]
            seqs = sorted(set(data['unseq']))
        logger.info('Found %d spectra (data type=%s with unique unseqs)'%(len(seqs),data_type))
    else:
        data = data[:0]
        logger.info('Found no spectra')
    
    
    
    if filename is not None:
        ascii.write_array(data,filename,auto_width=True,header=True)
    else:
        return data
    
def make_list_star(ID):
    """
    Mimics HermesTool MakeListStar without airmass and pm column.
    
    This works as input for HermesTool CCFList.py
    """
    data = search(ID)
    ascii.write_array([data['unseq'],data['date-avg'],[ID for i in data['unseq']],data['bjd'],data['bvcor'],data['prog_id'],data['exptime']],'%s.list'%(ID),axis0='cols')
    


def make_data_overview():
    """
    Summarize all Hermes data in a file for easy data retrieval.
    """
    #-- all hermes data directories
    dirs = sorted(glob.glob(os.path.join(config.ivs_dirs['hermes'],'20??????')))
    dirs = [idir for idir in dirs if os.path.isdir(idir)]
    obj_files = []
    #-- collect in those directories the raw and reduced files
    for idir in dirs:
        obj_files += sorted(glob.glob(os.path.join(idir,'raw','*OBJ*.fits')))
        obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*wavelength_merged.fits')))
        obj_files += sorted(glob.glob(os.path.join(idir,'reduced','*OBJ*wavelength_merged_c.fits')))
    
    #-- and summarize the contents in a tab separated file (some columns contain spaces)
    outfile = open('HermesFullDataOverview.tsv','w')
    outfile.write('#unseq prog_id obsmode bvcor observer object ra dec bjd exptime date-avg filename\n')
    outfile.write('#i i a20 >f8 a50 a50 >f8 >f8 >f8 >f8 a30 a200\n')
    for i,obj_file in enumerate(obj_files):
        sys.stdout.write(chr(27)+'[s') # save cursor
        sys.stdout.write(chr(27)+'[2K') # remove line
        sys.stdout.write('Scanning %5d / %5d FITS files'%(i+1,len(obj_files)))
        sys.stdout.flush() # flush to screen
        #-- keep track of: UNSEQ, BJD, BVCOR, OBSERVER, RA, DEC , PROG_ID, OBSMODE, EXPTIME, DATE-AVG, OBJECT and filename
        contents = dict(unseq=-1,prog_id=-1,obsmode='nan',bvcor=np.nan,observer='nan',
                        object='nan',ra=np.nan,dec=np.nan,
                        bjd=np.nan,exptime=np.nan,
                        filename=os.path.realpath(obj_file))
        contents['date-avg'] = 'nan'
        header = pyfits.getheader(obj_file)
        for key in contents:
            if key in header and key in ['unseq','prog_id']:
                try: contents[key] = int(header[key])
                except: pass
            elif key in header and key in ['obsmode','observer','object','date-avg']:
                contents[key] = str(header[key])
            elif key in header and key in ['bvcor','ra','dec','bjd','exptime']:
                contents[key] = float(header[key])
        outfile.write('%(unseq)d\t%(prog_id)d\t%(obsmode)s\t%(bvcor)f\t%(observer)s\t%(object)s\t%(ra)f\t%(dec)f\t%(bjd)f\t%(exptime)f\t%(date-avg)s\t%(filename)s\n'%contents)
        outfile.flush()
        sys.stdout.write(chr(27)+'[u') # reset cursor
    outfile.close()

if __name__=="__main__":
    import shutil
    make_data_overview()
    shutil.move('HermesFullDataOverview.tsv','/STER/pieterd/IVSDATA/catalogs/hermes/HermesFullDataOverview.tsv')
    