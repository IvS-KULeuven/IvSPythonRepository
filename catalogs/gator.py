# -*- coding: utf-8 -*-
"""
Interface to the GATOR search engine
"""
import os
import urllib
import logging
import ConfigParser
import numpy as np

from ivs.aux import loggers
from ivs.aux import numpy_ext
from ivs.io import ascii
from ivs.sed import filters
from ivs.units import conversions


logger = logging.getLogger("CAT.GATOR")
logger.addHandler(loggers.NullHandler())


basedir = os.path.dirname(os.path.abspath(__file__))

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open(os.path.join(basedir,'gator_cats_phot.cfg')))


#{ Basic interfaces

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
    filetype = kwargs.setdefault('filetype','1')
    if filename is not None and '.' in os.path.basename(filename):
        filetype = os.path.splitext(filename)[1][1:]
    elif filename is not None:
        filename = '%s.%s'%(filename,filetype)
    
    #-- gradually build URI
    base_url = _get_URI(catalog,**kwargs)
    #-- prepare to open URI
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename=filename)
    #   maybe we are just interest in the file, not immediately in the content
    if filename is not None:
        logger.info('Querying GATOR source %s and downloading to %s'%(catalog,filen))
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
        logger.info('Querying GATOR source %s (%d)'%(catalog,(results is not None and len(results) or 0)))
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





def get_photometry(ID=None,extra_fields=['dist','ra','dec'],**kwargs):
    """
    Download all available photometry from a star to a record array.
    
    For extra kwargs, see L{_get_URI} and L{gator2phot}
    
    Example usage:
    
    >>> import pylab
    >>> import vizier
    >>> name = 'kr cam'
    >>> master = vizier.get_photometry(name,to_units='erg/s/cm2/AA',extra_fields=[])
    >>> master = get_photometry(name,to_units='erg/s/cm2/AA',extra_fields=[],master=master)
    >>> p = pylab.figure()
    >>> wise = np.array(['WISE' in photband and True or False for photband in master['photband']])
    >>> p = pylab.errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ko')
    >>> p = pylab.errorbar(master['cwave'][wise],master['cmeas'][wise],yerr=master['e_cmeas'][wise],fmt='ro',ms=8)
    >>> p = pylab.gca().set_xscale('log')
    >>> p = pylab.gca().set_yscale('log')
    >>> p = pylab.show()
    
    Other examples:
    >>> master = get_photometry(ra=71.239527,dec=-70.589427,to_units='erg/s/cm2/AA',extra_fields=[],radius=1.)
    >>> master = get_photometry(ID='J044458.39-703522.6',to_units='W/m2',extra_fields=[],radius=1.)
    """
    kwargs['ID'] = ID
    to_units = kwargs.pop('to_units','erg/s/cm2/AA')
    master_ = kwargs.get('master',None)
    master = None
    #-- retrieve all measurements
    for source in cat_info.sections():
        results,units,comms = search(source,**kwargs)
        if results is not None:
            master = gator2phot(source,results,units,master,extra_fields=extra_fields)
    
    #-- convert the measurement to a common unit.
    if to_units and master is not None:
        #-- prepare columns to extend to basic master
        dtypes = [('cwave','f8'),('cmeas','f8'),('e_cmeas','f8'),('cunit','a50')]
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
            except AssertionError: # the error or flux must be positive number
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


#}

#{ Convenience functions


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
    indices = None
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
        if indices is None: break
        data.append([])
        for i in range(len(indices)-1):
            data[-1].append(line[indices[i]:indices[i+1]].strip())
    results = None
    units = {}
    #-- retrieve the data and put it into a record array
    if comms:
        #-- now convert this thing into a nice dictionary
        data = np.array(data)
        #-- retrieve the format of the columns. They are given in the
        #   Fortran format. In rare cases, columns contain multiple values
        #   themselves (so called vectors). In those cases, we interpret
        #   the contents as a long string
        names = [head.strip() for head in comms[0].split('|')[1:-1]]
        formats = comms[1]
        formats = formats.replace('double','f8')
        formats = formats.replace('char','a100')
        formats = formats.replace('int','f8')
        formats = formats.replace('long','f8')
        formats = [head.strip() for head in formats.split('|')[1:-1]]
        units_ = [head.strip() for head in comms[2].split('|')[1:-1]]
        #-- define dtypes for record array
        dtypes = np.dtype([(i,j) for i,j in zip(names,formats)])
        #-- remove spaces or empty values
        cols = []
        for i,key in enumerate(names):
            units[key] = units_[i]
            if len(data)==0: break
            col = data[:,i]
            #-- fill empty or null values with nan
            col = [(row.isspace() or not row or row=='null') and 'nan' or row for row in col]
            cols.append(col)
            
        #-- define columns for record array and construct record array
        if len(data)==0:
            results = None
        else:
            cols = np.array(cols)
            cols = [np.cast[dtypes[i]](cols[i]) for i in range(len(cols))]
            results = np.rec.array(cols, dtype=dtypes)
    return results,units,comms



def gator2phot(source,results,units,master=None,extra_fields=['_r','_RAJ2000','_DEJ2000']):
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
    dtypes = [('meas','f8'),('e_meas','f8'),('flag','a20'),
                  ('unit','a30'),('photband','a30'),('source','a50')]
    
    #-- extra can be added:
    names = list(results.dtype.names)
    translation = {'_r':'dist','_RAJ2000':'ra','_DEJ2000':'dec'}
    if extra_fields is not None:
        for e_dtype in extra_fields:
            dtypes.append((e_dtype,results.dtype[names.index(translation[e_dtype])].str))
    
    #-- create empty master if not given
    newmaster = False
    if master is None or len(master)==0:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
        newmaster = True
    
    #-- add fluxes and magnitudes to the record array
    cols_added = 0
    for key in cat_info.options(source):
        if key[:2] in [e_flag,q_flag] or key =='bibcode':
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
                cols += [translation[e_dtype] in results.dtype.names and results[translation[e_dtype]][:1] or np.ones(len(results[:1]))*np.nan]
        #-- add to the master
        rows = []
        for i in range(len(cols[0])):
            rows.append(tuple([col[i] for col in cols]))
        master = np.core.records.fromrecords(master.tolist()+rows,dtype=dtypes)
        cols_added += 1
    
    #-- skip first line from building 
    if newmaster: master = master[1:]
    return master

#}

#{ Internal helper functions

def _get_URI(name,ID=None,ra=None,dec=None,radius=1.,filetype='1',spatial='cone',**kwargs):
    """
    Build GATOR URI from available options.
    
    Filetype should be one of:
    
    6___ returns a program interface in XML 
    3___ returns a VO Table (XML) 
    2___ returns SVC (Software handshaking structure) message 
    1___ returns an ASCII table 
    0___ returns Gator Status Page in HTML (default)
    
    kwargs are to catch unused arguments.
    
    @param name: name of a GATOR catalog (e.g. 'wise_prelim_p3as_psd')
    @type name: str
    @keyword ID: target name
    @type ID: str
    @param filetype: type of the retrieved file 
    @type filetype: str (one of '0','1','2','3','6'... see GATOR site)
    @keyword radius: search radius (arcseconds)
    @type radius: float
    @param ra: target's right ascension
    @type ra: float
    @param dec: target's declination
    @type dec: float
    @param spatial: type of spatial query
    @type spatial: str (one of Cone, Box (NIY), Polygon (NIY), NONE (NIY))
    @return: url
    @rtype: str
    """
    base_url = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'
    base_url += 'catalog=%s'%(name)
    #base_url += '&spatial=cone'
    
    #-- in GATOR, many things should be given via the keyword 'objstr': right
    #   ascension, declination, target name...
    objstr = None
    
    #-- transform input to the 'objstr' paradigm
    if ID is not None:
        #-- if the ID is given in the form 'J??????+??????', derive the
        #   coordinates of the target from the name.
        if ID[0]=='J':
            ra = int(ID[1:3]),int(ID[3:5]),float(ID[5:10])
            dec = int(ID[10:13]),int(ID[13:15]),float(ID[15:])
            ra = '%02d+%02d+%.2f'%ra
            dec = '+%+02d+%02d+%.2f'%dec
            objstr = ra+dec
        else:
            objstr = ID
            
    #-- this is when ra and dec are given
    if ra is not None and dec is not None:
        ra = str(ra)
        dec = str(dec)
        ra = ra.replace(' ','+').replace(':','+')
        dec = dec.replace(' ','+').replace(':','+')
        objstr = '+'.join([ra,dec])
    
    #-- build up the URI
    if 'objstr' is not None:   base_url += '&objstr=%s'%(objstr)
    if 'filetype' is not None: base_url += '&outfmt=%s'%(filetype)
    if 'radius' is not None:   base_url += '&radius=%s'%(radius)
    
    logger.debug(base_url)
    
    return base_url


#}

if __name__=="__main__":
    import vizier
    from ivs.misc import loggers
    logger = loggers.get_basic_logger("")
    #-- example 1
    #master = get_photometry(ra=71.239527,dec=-70.589427,to_units='erg/s/cm2/AA',extra_fields=[],radius=1.)
    #master = vizier.get_photometry(ra=71.239527,dec=-70.589427,to_units='erg/s/cm2/AA',extra_fields=[],radius=5.,master=master)
    
    #-- example 2
    #master = get_photometry(ID='J044458.39-703522.6',to_units='W/m2',extra_fields=[],radius=1.)
    #master = vizier.get_photometry(ID='J044458.39-703522.6',to_units='W/m2',extra_fields=[],radius=5.,master=master)
    
    
    #-- example 3
    #master = get_photometry(ID='HD43317',to_units='W/m2',extra_fields=[],radius=1.)
    #master = vizier.get_photometry(ID='HD43317',to_units='W/m2',extra_fields=[],radius=5.,master=master)
    
    #-- example 4
    #master = get_photometry(ID='HD143454',to_units='erg/s/cm2/AA',extra_fields=[],radius=1.)
    #master = vizier.get_photometry(ID='HD143454',to_units='erg/s/cm2/AA',extra_fields=[],radius=30.,master=master)
    #-- example 5
    master = get_photometry(ID='RR Aql',to_units='erg/s/cm2/AA',extra_fields=[],radius=1.)
    master = vizier.get_photometry(ID='RR Aql',to_units='erg/s/cm2/AA',extra_fields=[],radius=30.,master=master)
    
    
    print master
    
    from pylab import *
    figure()
    #loglog(master['cwave'],master['cmeas'],'ko')
    #errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ro')
    plot(np.log10(master['cwave']),np.log10(master['cmeas']),'ko')
    #errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ro')
    #plot(np.log10(1538.62),np.log10(conversions.convert('muJy','W/m2',1202.758,photband='GALEX.FUV'))+13,'bo')
    #plot(np.log10(2315.66),np.log10(conversions.convert('muJy','W/m2',2896.989,photband='GALEX.NUV'))+13,'bo')
    #plot(np.log10(1538.62),np.log10(conversions.convert('muJy','W/m2',5025.582,photband='GALEX.FUV'))+13,'go')
    #plot(np.log10(2315.66),np.log10(conversions.convert('muJy','W/m2',7106.731,photband='GALEX.NUV'))+13,'go')
    #ylim(-1.5,2.5)
    #xlim(3.1,4.7)
    
    
    show()
    
    sys.exit()
    import doctest
    doctest.testmod()
    #logger = loggers.get_basic_logger()
    #catalogs = ['wise_prelim_p3as_psd']
    #import vizier
    #import pylab
    ##list_gator_catalogs()
    ##sys.exit()
    #name = 'kr cam'
    #master = vizier.get_photometry(name,to_units='erg/s/cm2/AA',extra_fields=[])
    #master = get_photometry(name,to_units='erg/s/cm2/AA',extra_fields=[],master=master)
    
    #pylab.figure()
    #wise = np.array(['WISE' in photband and True or False for photband in master['photband']])
    #pylab.errorbar(master['cwave'],master['cmeas'],yerr=master['e_cmeas'],fmt='ko')
    #pylab.errorbar(master['cwave'][wise],master['cmeas'][wise],yerr=master['e_cmeas'][wise],fmt='ro',ms=8)
    #pylab.gca().set_xscale('log')
    #pylab.gca().set_yscale('log')
    #pylab.show()