# -*- coding: utf-8 -*-
"""
Interface to Sesame for general information on a star (SIMBAD)
"""
import urllib
import logging

import numpy as np
from ivs.units import conversions
from ivs.aux import xmlparser
from ivs.catalogs import vizier

logger = logging.getLogger("CAT.SESAME")

def get_URI(ID,db='S'):
    """
    Build Sesame URI from available options.
    
    @param ID: name of the star
    @type ID: str
    @keyword db: database (one of 'S','N', or 'A')
    @type db: str
    @return: uri name
    @rtype: str
    """
    #mirrors:
    # http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oxpsIF/~%s?%s'
    ID = urllib.quote(ID)
    return 'http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxpsIF/%s?%s'%(db,ID)






def search(ID,db='S',fix=False):
    """
    Query Simbad, NED and/or Vizier for information on an identifier.
    
    This retrieves basic information on a star, e.g. as shown in a typical
    Simbad page: coordinates, spectral type, fluxes, aliases, references...
    
    Database C{db} is one of 'S' (SIMBAD), 'N' (NED), 'V' (Vizier)
    or 'A' (all).
    
    This function returns a (sometimes nested) dictionary. Example output is
    given below, where nested dictionaries are shown with the separator '.'
    between the keys.
    
    If you set C{fix} to C{False}, following values will be updated:
    
        1. the spectral type will be replaced by the one from the Skiff (2010)
        catalog if possible.
        2. The parallax will be replaced with the value from the new Van Leeuwen
        reduction.
        3. The galactic coordinates will be added (converted from RA and DEC)
        4. The proper motions will be taken from the PPMXL catalog from Roeser
        2010.
    
    Example usage:
    
    >>> info = search('vega',db='S')
    >>> print info['jpos']
    18:36:56.33 +38:47:01.2
    >>> print info['jdedeg']
    38.78369194
    >>> print info['alias'][1]
    * alf Lyr
    >>> print info['plx']['v']
    128.93
    >>> print info['mag']['B']['v']
    0.03
    
    This is an exhaustive list of example contents::
        Vel.e    = 0.9
        Vel.q    = A
        Vel.r    = 1979IAUS...30...57E
        Vel.v    = -13.9
        alias    = [u'V* alf Lyr', u'* alf Lyr', u'* 3 Lyr', u'ADS 11510 A', u'AG+38 1711', u'ASCC 507896', u'BD+38 3238', u'CCDM J18369+3847A', u'CEL 4636', u'CSI+38 3238 1', u'CSV 101745', u'1E 183515+3844.3', u'EUVE J1836+38.7', u'FK5 699', u'GC 25466', u'GCRV 11085', u'GEN# +1.00172167', u'GJ 721', u'HD 172167', u'HGAM 706', u'HIC 91262', u'HIP 91262', u'HR 7001', u'IDS 18336+3841 A', u'IRAS 18352+3844', u'IRC +40322', u'JP11 2999', u'LSPM J1836+3847', u'LTT 15486', u'2MASS J18365633+3847012', u'N30 4138', u'NAME VEGA', u'NLTT 46746', u'NSV 11128', u'8pc 128.93', u'PLX 4293.00', u'PLX 4293', u'PMC 90-93 496', u'PPM 81558', u'RAFGL 2208', u'ROT 2633', u'SAO 67174', u'SKY# 34103', u'TD1 22883', u'TYC 3105-2070-1', u'UBV 15842', u'UBV M 23118', u'USNO-B1.0 1287-00305764', u'USNO 882', u'uvby98 100172167 V', u'WDS J18369+3846A', u'Zkh 277', u'[HFE83] 1223']
        errDEmas = 5.4
        errRAmas = 5.16
        jdedeg   = 38.78369194
        jpos     = 18:36:56.33 +38:47:01.2
        jradeg   = 279.234735
        mag.B.q  = C
        mag.B.v  = 0.03
        mag.H.q  = C
        mag.H.r  = 2003yCat.2246....0C
        mag.H.v  = -0.03
        mag.I.q  = E
        mag.I.r  = 2003AJ....125..984M
        mag.I.v  = 0.2
        mag.J.q  = C
        mag.J.r  = 2003yCat.2246....0C
        mag.J.v  = -0.18
        mag.K.q  = C
        mag.K.r  = 2003yCat.2246....0C
        mag.K.v  = 0.13
        mag.R.q  = E
        mag.R.r  = 2003AJ....125..984M
        mag.R.v  = 0.1
        mag.V.q  = C
        mag.V.v  = 0.03
        nrefs    = 1860.0
        oid      = @2900336
        oname    = NAME VEGA
        otype    = V*
        plx.e    = 0.55
        plx.q    = A
        plx.r    = [u'1997A', u'&', u'A...323L..49P']
        plx.v    = 128.93
        pm.e     = 0.83
        pm.epmDE = 0.6
        pm.epmRA = 0.57
        pm.pa    = 35.0
        pm.pmDE  = 287.47
        pm.pmRA  = 201.03
        pm.q     = A
        pm.r     = [u'1997A', u'&', u'A...323L..49P']
        pm.v     = 350.79
        refPos   = [u'1997A', u'&', u'A...323L..49P']
        spNum    = 0.0000C800.0030.0000000000000000
        spType   = A0V
            
            
    >>> info = search('vega',db='N')
    >>> for key1 in sorted(info.keys()):
    ...    print '%s = %s'%(key1.ljust(8),info[key1])
    INFO     = from cache
    alias    = [u'alpha Lyr', u'HR 7001', u'HD 172167', u'IRAS  18352+3844', u'IRAS F18352+3844']
    errDEmas = 4824.0
    errRAmas = 19570.0
    jdedeg   = 38.782316
    jpos     = 18:36:55.70 +38:46:56.3
    jradeg   = 279.2321017
    oname    = VEGA 
    otype    = !*
    refPos   = 1990IRASF.C...0000M
    
    @param ID: name of the source
    @type ID: str
    @param db: database to use
    @type db: str ('N','S','V','A')
    @return: (nested) dictionary containing information on star
    @rtype: dictionary
    """
    base_url = get_URI(ID,db=db)
    ff = urllib.urlopen(base_url)
    xmlpage = ""
    for line in ff.readlines():
        line_ = line[::-1].strip(' ')[::-1]
        if line_[0]=='<':
            line = line_
        xmlpage+=line.strip('\n')
    database = xmlparser.XMLParser(xmlpage).content
    try:
        database = database['Sesame']['Target']['%s'%(db)]['Resolver']
        database = database[database.keys()[0]]
    except KeyError,IndexError:
        #-- we found nothing!
        database = {}
    ff.close()
    
    if fix:
        #-- fix the parallax: make sure we have the Van Leeuwen 2007 value.
        #   simbad seems to have changed to old values to the new ones somewhere
        #   in 2011. We check if this is the case for all stars:
        if 'plx' in database and not ('2007' in database['plx']['r']):
            data,units,comms = vizier.search('I/311/hip2',ID=ID)
            if data is not None and len(data):
                if not 'plx' in database:
                    database['plx'] = {}
                database['plx']['v'] = data['Plx'][0]
                database['plx']['e'] = data['e_Plx'][0]
                database['plx']['r'] = 'I/311/hip2'
        #-- fix the spectral type
        data,units,comms = vizier.search('B/mk/mktypes',ID=ID)
        if data is not None and len(data):
            database['spType'] = data['SpType'][0]
        #-- add galactic coordinates (in degrees)
        ra,dec = database['jpos'].split()
        gal = conversions.convert('equ','gal',(str(ra),str(dec)),epoch='2000')
        gal = float(gal[0])/np.pi*180,float(gal[1])/np.pi*180
        database['galpos'] = gal
        #-- fix the proper motions
        data,units,comms = vizier.search('I/317/sample',ID=ID)
        if data is not None and len(data):
            if not 'pm' in database:
                database['pm'] = {}
            database['pm']['pmRA'] = data['pmRA'][0]
            database['pm']['pmDE'] = data['pmDE'][0]
            database['pm']['epmRA'] = data['e_pmRA'][0]
            database['pm']['epmDE'] = data['e_pmDE'][0]
            database['pm']['r'] = 'I/317/sample'
    return database
    
if __name__=="__main__":
    import doctest
    doctest.testmod()