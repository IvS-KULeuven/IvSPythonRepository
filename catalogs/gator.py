# -*- coding: utf-8 -*-
"""
Interface to the GATOR search engine
"""
import os
import urllib
import logging
import ConfigParser
import numpy as np

from ivs.misc import loggers
from ivs.misc import numpy_ext
from ivs.io import ascii
from ivs.sed import filters
from ivs.units import conversions


logger = logging.getLogger("CAT.GATOR")
logger.addHandler(loggers.NullHandler)


basedir = os.path.dirname(os.path.abspath(__file__))

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open(os.path.join(basedir,'gator_cats_phot.cfg')))

def _get_URI(catalog,**kwargs):
    """
    Build GATOR URI from available options.
    
    @param catalog: name of a GATOR catalog (e.g. 'wise_prelim_p3as_psd')
    @type catalog: str
    @keyword objstr: target name
    @type objstr: str
    @keyword outfmt: type of output format
    @type outfmt: str ('1','2',...)
    @keyword radius: search radius (arcseconds)
    @type radius: float
    @return: url
    @rtype: str
    """
    base_url = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'
    base_url += 'catalog=%s'%(catalog)
    
    if 'objstr' in kwargs: base_url += '&objstr=%s'%(urllib.quote(kwargs['objstr']))
    if 'outfmt' in kwargs: base_url += '&outfmt=%s'%(kwargs['outfmt'])
    if 'radius' in kwargs: base_url += '&radius=%s'%(kwargs['radius'])
    return base_url


def search(catalog,**kwargs):
    """
    Search and retrieve information from a Gator catalog.
    
    Two ways to search for data within a catalog C{name}:
        
        1. You're looking for info on B{one target}, then give the target's
        C{ID} or coordinates (C{ra} and C{dec}), and a search C{radius}.
    
        2. You're looking for information of B{a whole field}, then give the
        field's coordinates (C{ra} and C{dec}), and C{radius}.
    
    If you have a list of targets, you need to loop this function.
    
    If you supply a filename, the results will be saved to that path, and you
    will get the filename back as received from urllib.URLopener (should be the
    same as the input name, unless something went wrong).
    
    If you don't supply a filename, you should leave C{filetype} to the default
    C{tsv}, and the results will be saved to a temporary
    file and deleted after the function is finished. The content of the file
    will be read into a dictionary, as well as the units (two separate
    dictionaries with the same keys, depending on the colum names in the
    catalog). The entries in the dictionary are of type C{ndarray}, and will
    be converted to a float-array if possible. If not, the array will consist
    of strings. The comments are also returned as a list of strings.
        
    
    @param catalog: name of a GATOR catalog (e.g. 'II/246/out')
    @type catalog: str
    @keyword filename: name of the file to write the results to (no extension)
    @type filename: str
    @return: filename / catalog data columns, units, comments
    @rtype: str/ record array, dict, list of str
    """
    filename = kwargs.pop('filename',None) # remove filename from kwargs
    filetype = kwargs.get('filetype','1')
    outfmt = kwargs.setdefault('outfmt','1')
    if filename is not None and '.' in os.path.basename(filename):
        filetype = os.path.splitext(filename)[1][1:]
    elif filename is not None:
        filename = '%s.%s'%(filename,filetype)
    
    #-- gradually build URI
    base_url = _get_URI(catalog,**kwargs)
    print base_url
    #-- prepare to open URI
    logger.info('Querying GATOR source %s'%(catalog))
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename=filename)
    #   maybe we are just interest in the file, not immediately in the content
    if filename is not None:
        logger.info('Downloaded results to file %s'%(filen))
        url.close()
        return filen
    
    #   otherwise, we read everything into a dictionary
    if filetype=='1':
        try:
            results,units,comms = txt2recarray(filen)
        #-- raise an exception when multiple catalogs were specified
        except ValueError:
            raise ValueError, "failed to read %s, perhaps multiple catalogs specified (e.g. III/168 instead of III/168/catalog)"%(catalog)
        url.close()
        logger.debug('Results converted to record array (found %d targets)'%(results is not None and len(results) or -1))
        return results,units,comms


def txt2recarray(filename):
    """
    Read a GATOR outfmt=1__ (ASCII) file into a record array.
    
    @param filename: name of the TSV file
    @type filename: str
    @return: catalog data columns, units, comments
    @rtype: record array, dict, list of str
    """
    ff = open(filename,'r')
    data = []
    comms = []
    while 1:  # might call read several times for a file
        line = ff.readline()
        if not line: break  # end of file
        
        #-- strip return character from line
        if line.isspace():
            continue # empty line
        
        #-- remove return characters
        line = line.replace('\n','')
        #-- when reading a comment line
        if line[0] =='\\':
            continue # treat next line
        if line[0] == '|':
            comms.append(line)
            indices = []
            for i in range(len(line)):
                if line[i]=='|': indices.append(i)
            continue
        data.append([])
        for i in range(len(indices)-1):
            data[-1].append(line[indices[i]:indices[i+1]].strip())
    
    results = None
    units = {}
    #-- retrieve the data and put it into a record array
    if len(data)>0:
        #-- now convert this thing into a nice dictionary
        data = np.array(data)
        #-- retrieve the format of the columns. They are given in the
        #   Fortran format. In rare cases, columns contain multiple values
        #   themselves (so called vectors). In those cases, we interpret
        #   the contents as a long string
        names = [head.strip() for head in comms[0].split('|')[1:-1]]
        formats = comms[1]
        formats = formats.replace('double','>f8')
        formats = formats.replace('char','a100')
        formats = formats.replace('int','>f8')
        formats = formats.replace('long','>f8')
        formats = [head.strip() for head in formats.split('|')[1:-1]]
        units_ = [head.strip() for head in comms[2].split('|')[1:-1]]
        #-- define dtypes for record array
        dtypes = np.dtype([(i,j) for i,j in zip(names,formats)])
        #-- remove spaces or empty values
        cols = []
        for i,key in enumerate(names):
             col = data[:,i]
             #-- fill empty values with nan
             cols.append([(row.isspace() or not row) and np.nan or row for row in col])
             ##-- fix unit name
             units[key] = units_[i]
        #-- define columns for record array and construct record array
        cols = data.T
        cols = [np.cast[dtypes[i]](cols[i]) for i in range(len(cols))]
        results = np.rec.array(cols, dtype=dtypes)
    return results,units,comms


def list_catalogs():
    """
    Return a list of all availabe GATOR catalogues.
    
    @return: list of gator catalogs and discriptions
    @rtype: list of string tuples
    """
    url = urllib.URLopener()
    filen,msg = url.retrieve('http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mode=ascii')
    results,units,comms = txt2recarray(filen)
    cats = []
    for i in range(len(results)):
        cats += [(results['catname'][i],results['description'][i])]
        logger.info('%-20s: %s'%cats[-1])
    url.close()
    return cats


def gator2phot(source,results,units,master=None,extra_fields=None):
    """
    Convert/combine Gator record arrays to measurement record arrays.
    
    Every line in the combined array represents a measurement in a certain band.
    
    This is probably only useful if C{results} contains only information on
    one target (or you have to give 'ID' as an extra field, maybe).
    
    The standard columns are:
    
        1. C{meas}: containing the photometric measurement
        2. C{e_meas}: the error on the photometric measurement
        3. C{flag}: an optional quality flag
        4. C{unit}: the unit of the measurement
        5. C{photband}: the photometric passband (FILTER.BAND)
        6. C{source}: name of the source catalog
    
    You can add extra information from the Gator catalog via the list of keys
    C{extra_fields}.
    
    If you give a C{master}, the information will be added to a previous
    record array. If not, a new master will be created.
    
    The result is a record array with each row a measurement.
    
    @param source: name of the Gator source
    @type source: str
    @param results: results from Gator C{search}
    @type results: record array
    @param units: header of Gator catalog with key name the column name and
    key value the units of that column
    @type units: dict
    @param master: master record array to add information to
    @type master: record array
    @param extra_fields: any extra columns you want to add information from
    @type extra_fields: list of str
    @return: array with each row a measurement
    @rtype: record array
    """
    e_flag = 'e_'
    q_flag = 'q_'
    #-- basic dtypes
    dtypes = [('meas','>f4'),('e_meas','>f4'),('flag','a20'),
                  ('unit','a30'),('photband','a30'),('source','a50')]
    
    #-- extra can be added:
    names = list(results.dtype.names)
    if extra_fields is not None:
        translation = dict(dist='_r',ra='_RAJ2000',dec='_DEJ2000')
        for e_dtype in extra_fields:
            dtypes.append((translation[e_dtype],results.dtype[names.index(e_dtype)].str))
    
    #-- create empty master if not given
    newmaster = False
    if master is None or len(master)==0:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
        newmaster = True
    
    #-- add fluxes and magnitudes to the record array
    cols_added = 0
    for key in cat_info.options(source):
        if key[:2] in [e_flag,q_flag]:
            continue
        photband = cat_info.get(source,key)
        #-- contains measurement, error, quality, units, photometric band and
        #   source
        cols = [results[key][:1],
                (e_flag+key in cat_info.options(source) and results[cat_info.get(source,e_flag+key)][:1] or np.ones(len(results[:1]))*np.nan),
                (q_flag+key in cat_info.options(source) and results[cat_info.get(source,q_flag+key)][:1] or np.ones(len(results[:1]))*np.nan),
                np.array(len(results[:1])*[units[key]],str),
                np.array(len(results[:1])*[photband],str),
                np.array(len(results[:1])*[source],str)]
        #-- add any extra fields if desired.
        if extra_fields is not None:
            for e_dtype in extra_fields:
                cols.append(results[:1][e_dtype])
        #-- add to the master
        rows = []
        for i in range(len(cols[0])):
            rows.append(tuple([col[i] for col in cols]))
        master = np.core.records.fromrecords(master.tolist()+rows,dtype=dtypes)
        cols_added += 1
    
    #-- skip first line from building 
    if newmaster: master = master[1:]
    return master


def get_photometry(ID,extra_fields=['dist','ra','dec'],**kwargs):
    """
    Download all available photometry from a star to a record array.
    
    For extra kwargs, see L{_get_URI} and L{gator2phot}
    
    Example usage:
    
    >>> import pylab
    >>> import vizier
    >>> name = 'kr cam'
    >>> master = vizier.get_photometry(name,to_units='erg/s/cm2/A',extra_fields=[])
    >>> master = get_photometry(name,to_units='erg/s/cm2/A',extra_fields=[],master=master)
    >>> p = pylab.figure()
    >>> wise = np.array(['WISE' in photband and True or False for photband in master['photband']])
    >>> p = pylab.errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ko')
    >>> p = pylab.errorbar(master['cwave'][wise],master['cmeas'][wise],yerr=master['e_cmeas'][wise],fmt='ro',ms=8)
    >>> p = pylab.gca().set_xscale('log')
    >>> p = pylab.gca().set_yscale('log')
    >>> p = pylab.show()
    """
    to_units = kwargs.pop('to_units','erg/s/cm2/A')
    master_ = kwargs.get('master',None)
    master = None
    #-- retrieve all measurements
    for source in cat_info.sections():
        results,units,comms = search(source,objstr=ID,**kwargs)
        if results is not None:
            master = gator2phot(source,results,units,master,extra_fields=extra_fields)
    
    #-- convert the measurement to a common unit.
    if to_units and master is not None:
        #-- prepare columns to extend to basic master
        dtypes = [('cwave','>f4'),('cmeas','>f4'),('e_cmeas','>f4'),('cunit','a50')]
        cols = [[],[],[],[]]
        #-- forget about 'nan' errors for the moment
        no_errors = np.isnan(master['e_meas'])
        master['e_meas'][no_errors] = 0.
        #-- extend basic master
        zp = filters.get_info(master['photband'])
        for i in range(len(master)):
            try:
                value,e_value = conversions.convert(master['unit'][i],to_units,master['meas'][i],master['e_meas'][i],photband=master['photband'][i])
            except ValueError: # calibrations not available
                value,e_value = np.nan,np.nan
            try:
                eff_wave = filters.eff_wave(master['photband'][i])
            except IOError:
                eff_wave = np.nan
            cols[0].append(eff_wave)
            cols[1].append(value)
            cols[2].append(e_value)
            cols[3].append(to_units)
        master = numpy_ext.recarr_addcols(master,cols,dtypes)
        #-- reset errors
        master['e_meas'][no_errors] = np.nan
        master['e_cmeas'][no_errors] = np.nan
    
    if master_ is not None and master is not None:
        master = numpy_ext.recarr_addrows(master_,master.tolist())
    elif master_ is not None:
        master = master_
    
    #-- and return the results
    return master



if __name__=="__main__":
    import doctest
    doctest.testmod()
    #logger = loggers.get_basic_logger()
    #catalogs = ['wise_prelim_p3as_psd']
    #import vizier
    #import pylab
    ##list_gator_catalogs()
    ##sys.exit()
    #name = 'kr cam'
    #master = vizier.get_photometry(name,to_units='erg/s/cm2/A',extra_fields=[])
    #master = get_photometry(name,to_units='erg/s/cm2/A',extra_fields=[],master=master)
    
    #pylab.figure()
    #wise = np.array(['WISE' in photband and True or False for photband in master['photband']])
    #pylab.errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ko')
    #pylab.errorbar(master['cwave'][wise],master['cmeas'][wise],yerr=master['e_cmeas'][wise],fmt='ro',ms=8)
    #pylab.gca().set_xscale('log')
    #pylab.gca().set_yscale('log')
    #pylab.show()