# -*- coding: utf-8 -*-
"""
Interface to Geneva's Genaral Catalogue of Photometric Data
"""
import urllib
import logging
import os
import ConfigParser

import numpy as np

from ivs.aux import loggers
from ivs.aux import numpy_ext
from ivs.catalogs import sesame
from ivs.catalogs import vizier
from ivs.sed import filters
from ivs.units import conversions

#-- the photometric systems in GCPD are represented by a number
systems = {'JOHNSON':1,
           'STROMGREN':4,
           'ARGUE':10,
           'GENEVA':13,
           'KRON':19,
           'VILNIUS':21,
           'WOOD':22,
           'WBVR':78,
           'STRAIZYS':80,
           'DDO':12}

logger = logging.getLogger("CAT.VIZIER")
logger.addHandler(loggers.NullHandler())

basedir = os.path.dirname(os.path.abspath(__file__))

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open(os.path.join(basedir,'gcpd_cats_phot.cfg')))


#{ Basic interfaces

def search(name,**kwargs):
    """
    Search and retrieve information from the GCPD catalog.

    @param name: name of photometric system
    @type name: string
    """
    base_url = _get_URI(name,**kwargs)

    #-- the data is listed in two lines: one with the header, one with
    #   the values
    webpage = urllib.urlopen(base_url)
    entries,values = None,None
    log_message = '0'
    start = -1
    for line in webpage:
        line = line.replace('\n','')
        #-- change this marker to be much more general: now it only works for
        #   geneva photometry
        if r'<HR>' in line: start = 0
        if start==1 and r'<B>' in line and not 'No values' in line:
            entries = line.split('<B>')[1].split('</B>')[0].split('\t')
            #-- sometimes there's two columns with the name N: make them unique
            entries = [(entry=='N' and 'N%d'%(i) or entry) for i,entry in enumerate(entries)]
            log_message = '1'
            continue
        if entries is not None:
            values = line.split('\t')
            break
        start += 1
    logger.info('Querying GCPD for %s (%s)'%(name,log_message))
    webpage.close()

    #-- assume that all entries are floats, and are values are magnitudes
    if entries is not None:
        if len(values) < len(entries):
            values += ['nan']*(len(entries)-len(values))
        #-- in some columns, there are '*' or 'STD' signs
        values = [value.replace('*','').replace('STD','').replace('/','') for value in values]
        #-- some columns have no values
        values = [(value and value or 'nan') for value in values]
        #-- some columns are upper limits (or whatever ':' means)
        values = [(':' in value and 'nan' or value) for value in values]
        dtypes = np.dtype([(i,'f8') for i in entries])
        units = {}
        for entry in entries:
            units[entry] = 'mag'
        cols = [np.cast[dtypes[i]](np.array([values[i]])) for i in range(len(values))]
        results = np.rec.array(cols, dtype=dtypes)
    else:
        results,units,comms = None,None,None

    return results, units, None


def gcpd2phot(source,results,units,master=None,e_flag='e_',q_flag='q_',extra_fields=None):
    """
    Convert/combine GCPD record arrays to measurement record arrays.

    Every line in the combined array represents a measurement in a certain band.

    The standard columns are:

        1. C{meas}: containing the photometric measurement
        2. C{e_meas}: the error on the photometric measurement
        3. C{flag}: an optional quality flag
        4. C{unit}: the unit of the measurement
        5. C{photband}: the photometric passband (FILTER.BAND)
        6. C{source}: name of the source catalog

    If you give a C{master}, the information will be added to a previous
    record array. If not, a new master will be created.

    Colors will be expanded, derived from the other columns and added to the
    master.

    The result is a record array with each row a measurement.

    Extra fields are not available for the GCPD, they will be filled in with
    nans.

    @param source: name of the VizieR source
    @type source: str
    @param results: results from VizieR C{search}
    @type results: record array
    @param units: header of Vizier catalog with key name the column name and
    key value the units of that column
    @type units: dict
    @param master: master record array to add information to
    @type master: record array
    @param e_flag: flag denoting the error on a column
    @type e_flag: str
    @param q_flag: flag denoting the quality of a measurement
    @type q_flag: str
    @param extra_fields: any extra columns you want to add information from
    @type extra_fields: list of str
    @return: array with each row a measurement
    @rtype: record array
    """
    if cat_info.has_option(source,'e_flag'):
        e_flag = cat_info.get(source,'e_flag')

    #-- basic dtypes
    dtypes = [('meas','f8'),('e_meas','f8'),('flag','a20'),
                  ('unit','a30'),('photband','a30'),('source','a50')]

    #-- extra can be added, but only if a master is already given!! The reason
    #   is that thre GCPD actually does not contain any extra information, so
    #   we will never be able to add it and will not know what dtype the extra
    #   columns should be
    #-- extra can be added:
    names = list(results.dtype.names)
    if extra_fields is not None:
        for e_dtype in extra_fields:
            dtypes.append((e_dtype,'f8'))

    #-- create empty master if not given
    newmaster = False
    if master is None or len(master)==0:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
        newmaster = True

    #-- add fluxes and magnitudes to the record array
    cols_added = 0
    for key in cat_info.options(source):
        if key[:2] =='e_' or key=='bibcode':
            continue
        photband = cat_info.get(source,key)
        #-- contains measurement, error, quality, units, photometric band and
        #   source
        cols = [results[key][:1],
                (e_flag+key in cat_info.options(source) and results[cat_info.get(source,e_flag+key)][:1] or np.ones(len(results[:1]))*np.nan),
                (q_flag+key in results.dtype.names and results[q_flag+key][:1] or np.ones(len(results[:1]))*np.nan),
                np.array(len(results[:1])*[units[key]],str),
                np.array(len(results[:1])*[photband],str),
                np.array(len(results[:1])*['GCPD'],str)]
        #-- add any extra fields if desired.
        if extra_fields is not None:
            for e_dtype in extra_fields:
                cols += [e_dtype in results.dtype.names and results[e_dtype][:1] or np.ones(len(results[:1]))*np.nan]
        #-- add to the master
        rows = []
        for i in range(len(cols[0])):
            rows.append(tuple([col[i] for col in cols]))
        master = np.core.records.fromrecords(master.tolist()+rows,dtype=dtypes)
        cols_added += 1

    #-- fix colours: we have to run through it two times to be sure to have
    #   all the colours
    N = len(master)-cols_added
    master_ = vizier._breakup_colours(master[N:])
    for i in range(5):
        master_ = vizier._breakup_colours(master_)
    #-- combine and return
    master = np.core.records.fromrecords(master.tolist()[:N]+master_.tolist(),dtype=dtypes)

    #-- skip first line from building
    if newmaster: master = master[1:]
    return master

def get_photometry(ID=None,extra_fields=[],**kwargs):
    """
    Download all available photometry from a star to a record array.

    Extra fields will not be useful probably.

    For extra kwargs, see L{_get_URI} and L{gcpd2phot}
    """
    to_units = kwargs.pop('to_units','erg/s/cm2/AA')
    master_ = kwargs.get('master',None)
    master = None
    #-- retrieve all measurements
    for source in cat_info.sections():
        results,units,comms = search(source,ID=ID,**kwargs)
        if results is not None:
            master = gcpd2phot(source,results,units,master,extra_fields=extra_fields)

    #-- convert the measurement to a common unit.
    if to_units and master is not None:
        #-- prepare columns to extend to basic master
        dtypes = [('cwave','f8'),('cmeas','f8'),('e_cmeas','f8'),('cunit','a50')]
        cols = [[],[],[],[]]
        #-- forget about 'nan' errors for the moment
        no_errors = np.isnan(master['e_meas'])
        master['e_meas'][no_errors] = 0.
        #-- extend basic master
        try:
            zp = filters.get_info(master['photband'])
        except:
            print master['photband']
            raise
        for i in range(len(master)):
            to_units_ = to_units+''
            try:
                value,e_value = conversions.convert(master['unit'][i],to_units,master['meas'][i],master['e_meas'][i],photband=master['photband'][i])
            except ValueError: # calibrations not available
                # if it is a magnitude color, try converting it to a flux ratio
                if 'mag' in master['unit'][i]:
                    try:
                        value,e_value = conversions.convert('mag_color','flux_ratio',master['meas'][i],master['e_meas'][i],photband=master['photband'][i])
                        to_units_ = 'flux_ratio'
                    except ValueError:
                        value,e_value = np.nan,np.nan
                # else, we are powerless...
                else:
                    value,e_value = np.nan,np.nan
            try:
                eff_wave = filters.eff_wave(master['photband'][i])
            except IOError:
                eff_wave = np.nan
            cols[0].append(eff_wave)
            cols[1].append(value)
            cols[2].append(e_value)
            cols[3].append(to_units_)
        master = numpy_ext.recarr_addcols(master,cols,dtypes)
        #-- reset errors
        master['e_meas'][no_errors] = np.nan
        master['e_cmeas'][no_errors] = np.nan

    #-- if a master is given as a keyword, and data is found in this module,
    #   append the two
    if master_ is not None and master is not None:
        master = numpy_ext.recarr_addrows(master_,master.tolist())
    elif master is None:
        master = master_
    #-- and return the results
    return master

#}
#{ Internal helper functions

def _get_URI(name='GENEVA',ID=None,**kwargs):
    """
    Build GCPD URI from available options.

    kwargs are to catch unused arguments.

    @param name: photometric system name (E.g. JOHNSON, STROMGREN, GENEVA...)
    @type name: string (automatically uppercased)
    @keyword ID: star name (should be resolvable by SIMBAD)
    @type ID: string
    @return:
    """
    #-- GCPD is poor at recognising aliases: therefore we try different
    #   identifiers retrieved from Sesame that GCPD understands
    recognized_alias = ['HD','BD',"CD"]

    try:
        aliases = sesame.search(ID)['alias']
        for alias in aliases:
            if alias[:2] in recognized_alias:
                ID = alias[:2]+' '+alias[2:]
                break
        else:
            logger.error('Star %s has no aliases recognised by GCPD: query will not return results'%(ID))
    except KeyError:
        logger.error('Unknown star %s: GCPD query will not return results'%(ID))


    base_url = 'http://obswww.unige.ch/gcpd/cgi-bin/photoSys.cgi?phot=%02d&type=original&refer=with&mode=starno&ident=%s'%(systems[name],urllib.quote(ID))
    logger.debug(base_url)
    return base_url

#}

if __name__=="__main__":
    results,units,comms = search('GENEVA',ID='HD180642')
    master = gcpd2phot('GENEVA',results,units)
    print master
    print ""
    master = get_photometry(ID='vega')
    print master
