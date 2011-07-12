# -*- coding: utf-8 -*-
"""
Crossmatch catalogs for data on one target, or crossmatch entire catalogs.
"""
import logging
import itertools
import pylab as pl
import numpy as np

from ivs.catalogs import vizier
from ivs.catalogs import gator
from ivs.catalogs import gcpd
from ivs.catalogs import mast
from ivs.units import conversions
from ivs.aux import numpy_ext
from ivs.aux import progressMeter
from ivs.aux import loggers
from ivs.io import ascii

from scipy.spatial import KDTree

logger = logging.getLogger("CAT.XMATCH")

def get_photometry(ID,to_units='erg/s/cm2/A',extra_fields=[],include=None,
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
    searchables = ['gator','vizier','gcpd']
    if include is not None:
        searchables = include
    if exclude is not None:
        for exclude_ in exclude:
            searchables.remove(exclude_)
    
    #-- and search photometry
    if 'gator' in searchables:
        kwargs['master'] = gator.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    if 'vizier' in searchables:
        kwargs['master'] = vizier.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    if 'gcpd' in searchables:
        kwargs['master'] = gcpd.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    #master = mast.get_photometry(ID=ID,to_units=to_units,extra_fields=extra_fields,**kwargs)
    master = kwargs['master']
    
    #-- now make a summary of the contents:
    photbands = [phot.split('.')[0]  for phot in master['photband']]
    contents = [(i,photbands.count(i)) for i in sorted(list(set(photbands)))]
    for phot in contents:
        logger.info('%10s: found %d measurements'%phot)
    return master

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
    

if __name__=="__main__":
    logger = loggers.get_basic_logger()
    #xmatch_vizier('II/169/main','B/pastel/pastel')