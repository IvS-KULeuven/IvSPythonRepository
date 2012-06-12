# -*- coding: utf-8 -*-
"""
Crossmatch catalogs for data on one target, or crossmatch entire catalogs.

This module is callable from the command line. E.g., simply type::
    
    $:> python crossmatch.py get_photometry ID=HD180642

or ask for help on any function::
    
    $:> python crossmatch.py get_photometry --help
"""
import logging
import itertools
import pylab as pl
import numpy as np

from ivs.catalogs import vizier
from ivs.catalogs import gator
from ivs.catalogs import gcpd
from ivs.catalogs import mast
from ivs.catalogs import sesame
from ivs.units import conversions
from ivs.aux import numpy_ext
from ivs.aux import progressMeter
from ivs.aux import loggers
from ivs.aux import argkwargparser
from ivs.io import ascii

from scipy.spatial import KDTree

logger = logging.getLogger("CAT.XMATCH")

def get_photometry(ID=None,to_units='erg/s/cm2/A',extra_fields=[],include=None,
         exclude=None,**kwargs):
    """
    Collect photometry from different sources.
    
    The output consists of a record array containing the following keys:
    
    'meas': the measurement's value directly from the catalog in original units
    'e_meas': the error on the measurements
    'flag': any flags that are associated with the measurement in a catalog
    'unit': the unit of the original measurement
    'source' the source catalog's name
    'photband': photometric pass bands' name
    'cwave': the effective wavelength of the passband
    'cmeas': converted measurement (to C{to_units})
    'e_cmeas': error on converted measurment
    'cunit': converted unit
    
    Be aware that some of the values can be 'nan': e.g. sometimes no error is
    listed in the catalog, or no flag. Also the 'cwave' column will be nan for
    all photometric colours (e.g. B-V)
    
    If you define C{extra_fields}, make sure all the {get_photometry} know how
    to handle it: probably some default values need to be inserted if these
    extra columns are not available in some catalog. It is safest just to leave
    it blank.
    
    You can include or exclude search sources via C{include} and C{exclude}.
    When given, these should be a list containing strings. The default is to
    include C{gator}, C{vizier} and C{gcpd}.
    
    Extra keyword arguments are passed to each C{get_photometry} functions in
    this package's modules.
    
    Example usage:
    
        1. You want to download all available photometry and write the results to
        an ASCII file for later reference.
    
        >>> master = get_photometry('vega')
        >>> ascii.write_array(master,header=True,auto_width=True)
    
        2. You want to plot the raw, unmodelled SED of an object:
    
        >>> master = get_photometry('vega')
        >>> pl.errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ko')
        >>> pl.gca().set_xscale('log',nonposx='clip')
        >>> pl.gca().set_yscale('log',nonposy='clip')
        
    We made no difference between colors (B-V) and magnitude (V), because the
    'cwave' for colors is 'nan', so they will not be plotted anyway.The final
    two lines are just to correct errorbars that go below zero in a logarithmic
    plot.
    
    @param ID: the target's name, understandable by SIMBAD
    @type ID: str
    @param to_units: units to convert everything to.
    @type to_units:
    @param include: sources to include
    @type include: list of strings (from C{gator}, C{vizier} or C{gcpd})
    @param exclude: sources to include
    @type exclude: list of strings (from C{gator}, C{vizier} or C{gcpd})
    @return: record array where eacht entry is a photometric measurement
    @rtype: record array
    """
    #-- make sure all catalog names are lower case
    if include is not None: include = [i.lower() for i in include]
    if exclude is not None: exclude = [i.lower() for i in exclude]
    
    #-- check which sources to include/exclude
    searchables = ['gator','vizier','gcpd','mast']
    if include is not None:
        searchables = include
    if exclude is not None:
        searchables = list( set(searchables)- set(exclude))
    
    #-- and search photometry
    if 'mast' in searchables:
        kwargs['master'] = mast.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    if 'gator' in searchables:
        kwargs['master'] = gator.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    if 'vizier' in searchables:
        #-- first query catalogs that can only be queried via HD number
        info = sesame.search(ID=ID,fix=True)
        if 'alias' in info:
            HDnumber = [name for name in info['alias'] if name[:2]=='HD']
            if HDnumber:
                kwargs['master'] = vizier.get_photometry(extra_fields=extra_fields,constraints=['HD=%s'%(HDnumber[0][3:])],sources=['II/83/catalog','V/33/phot'],sort=None,**kwargs)
        #-- then query catalogs that can only be queried via another catalog
        results,units,comms = vizier.search('J/A+A/380/609/table1',ID=ID)
        if results is not None:
            catname = results[0]['Name'].strip()
            kwargs['master'] = vizier.get_photometry(extra_fields=extra_fields,constraints=['Name={0}'.format(catname)],sources=['J/A+A/380/609/table{0}'.format(tnr) for tnr in range(2,5)],sort=None,**kwargs)
        kwargs['master'] = vizier.get_photometry(extra_fields=extra_fields,constraints=['HD=%s'%(HDnumber[0][3:])],sources=['II/83/catalog','V/33/phot'],sort=None,**kwargs)
        #-- then query normal catalogs
        kwargs['master'] = vizier.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
        
        
    if 'gcpd' in searchables:
        kwargs['master'] = gcpd.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    master = kwargs['master']
    
    #-- now make a summary of the contents:
    photbands = [phot.split('.')[0]  for phot in master['photband']]
    contents = [(i,photbands.count(i)) for i in sorted(list(set(photbands)))]
    for phot in contents:
        logger.info('%10s: found %d measurements'%phot)
    return master

def add_bibcodes(master):
    """
    Add bibcodes to a master record.
    
    @param master: master record of photometry
    @type master: numpy record array
    @return: master record of photometry with additional column
    @rtype: numpy record array
    """
    bibcodes = []
    for source in master['source']:
        bibcodes.append('None')
        #-- check if it's a ViZieR catalog
        if source.count('/')>=2:
            #-- some catalogs have no ADS links but do have a bibcode, or the
            #   bibcode is manually put in the *cats_phot.cfg file
            if 'bibcode' in vizier.cat_info.options(source):
                bibcodes[-1] = vizier.cat_info.get(source,'bibcode')
            #-- others may be retrieved immediately from the ViZieR site
            else:
                bibcodes[-1] = vizier.catalog2bibcode(source)
        elif source=='GCPD':
            bibcodes[-1] = '1997A&AS..124..349M'
        elif source in gator.cat_info.sections():
            if 'bibcode' in gator.cat_info.options(source):
                bibcodes[-1] = gator.cat_info.get(source,'bibcode')
        elif source in mast.cat_info.sections():
            if 'bibcode' in mast.cat_info.options(source):
                bibcodes[-1] = mast.cat_info.get(source,'bibcode')
        if bibcodes[-1]=='None' or bibcodes[-1] is None:
            logger.info('Bibcode of source {0:25s} not found'.format(source))
    bibcodes = np.array(bibcodes,str)
    bibcodes = np.rec.fromarrays([bibcodes],names=['bibcode'])
    return numpy_ext.recarr_join(master,bibcodes)


def make_bibtex(master,ID):
    """
    Make a bib file from a master record
    """
    bibcodes = sorted(set(list(master['bibcode'])))
    with open('{0}.bib'.format(ID),'w') as ff:
        ff.write('% \citet{{{0}}}\n\n\n'.format(','.join(bibcodes)))
        for bibcode in bibcodes:
            try:
                ff.write(vizier.bibcode2bibtex(bibcode)+'\n')
            except IOError:
                logger.info('Error retrieving bibtex for {0}'.format(bibcode))
        
    


def xmatch(coords1,coords2):
    """
    Crossmatch two sets of coordinates.
    """
    tree = KDTree(coords1)
    distance,order = tree.query(coords2)
    
    match_success = distance<(tol/(60.))
    
    
    #-- first tuple contains matches (cat1,cat2)
    #   second tuple contains mismatches (cat1,cat2)
    #   so you if you call the function via 
    #   >>> match,mismatch = xmatch(coord1,coord2)
    #   you probably want to do something like
    #   >>> matched_cat1,matched_cat2 = cat1[match[0]],cat2[match[1]]
    #   >>> mismatched_cat1,mismatched_cat2 = cat1[mismatch[0]],cat2[mismatch[1]]
    #   Then, you can just glue together these four parts:
    #   [    matched_cat1        matched_cat2    ]
    #   [ mismatched_cat1          nan           ]
    #   [       nan              mismatched_cat2 ]
    return (order[match_success],match_success),(order[-match_success],-match_success)

def photometry2str(master):
    """
    String representation of master record array
    
    @param master: master record array containing photometry
    @type master: numpy record array
    """
    master = master[np.argsort(master['photband'])]
    txt = '#%19s %12s %12s %12s %10s %12s %12s %11s %s\n'%('PHOTBAND','MEAS','E_MEAS','UNIT','CWAVE','CMEAS','E_CMEAS','UNIT','SOURCE')
    txt+= '#=========================================================================================================================\n'
    for i,j,k,l,m,n,o,p in zip(master['photband'],master['meas'],master['e_meas'],master['unit'],master['cwave'],master['cmeas'],master['e_cmeas'],master['source']):
        txt+='%20s %12g %12g %12s %10.0f %12g %12g erg/s/cm2/A %s\n'%(i,j,k,l,m,n,o,p)
    return txt    

if __name__=="__main__":
    logger = loggers.get_basic_logger()
    
    method,args,kwargs = argkwargparser.parse()
    print_help = '--help' in args or '-h' in args
    if print_help:
        help(globals()[method])
    else:
        master = globals()[method](*args,**kwargs)
        print photometry2str(master)
        
    

    #xmatch_vizier('II/169/main','B/pastel/pastel')