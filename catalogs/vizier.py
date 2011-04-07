# -*- coding: utf-8 -*-
"""
Interface to the ViZieR website.

Download or retrieve ViZieR catalogs.

You can define 'standard' photometric catalogs in the C{vizier_cats.cfg} file.
This files is basically a translator for VizieR column headers to photometric
passbands (and colors).
"""
#-- standard libraries
import numpy as np
import urllib
import logging
import os
import ConfigParser

#-- IvS repository
from ivs.io import ascii
from ivs.units import conversions

logger = logging.getLogger("CAT.VIZIER")

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'vizier_cats.cfg')))


def search(source,**kwargs):
    """
    Search and retrieve information from a VizieR catalog.
    
    Two ways to search for data within a catalog C{source}:
        
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
    
    WARNING: when retrieving a FITS file, ViZieR sometimes puts weird formats
    into the header ('3F10.6E' in the 2MASS catalog), which cannot be read by
    the C{pyfits} module. One option is to download to another format, or to
    restrict the columns with C{out_all=None}.
    
    Example usage:
    
        1. Look for the Geneva V magnitude and coordinates of Vega in the GENEVA
        catalog of Rufener.
    
        >>> results,units,comms = search('II/169/main',ID='vega',radius=60.)
        >>> print "Vega: Vmag = %.3f %s, RA = %.2f %s, DEC = %.2f %s"%(results['Vmag'],units['Vmag'],results['_RAJ2000'],units['_RAJ2000'],results['_DEJ2000'],units['_DEJ2000'])
        Vega: Vmag = 0.061 mag, RA = 279.24 deg, DEC = 38.77 deg
    
        2. Search for all targets in the 2MASS catalog in a particular field.
        Download the results to a FITS file, read the file, and plot the results
        to the screen.
    
        >>> filename = search('II/246/out',ra=100.79,dec=0.70,radius=1000.,filetype='fits',filename='2mass_test',out_all=None)
    
        Now read in the FITS-file and plot the contents
        
        >>> import pyfits,pylab,os
        >>> ff = pyfits.open('2mass_test.fits')
        >>> p = pylab.gcf().canvas.set_window_title('test of <search>')
        >>> p = pylab.scatter(ff[1].data.field('_RAJ2000'),ff[1].data.field('_DEJ2000'),c=ff[1].data.field('Jmag'),s=(20-ff[1].data.field('Jmag'))**2,cmap=pylab.cm.hot_r,edgecolors='none')
        >>> p = pylab.colorbar()
        >>> p = p.set_label('Jmag')
        >>> p,q = pylab.xlabel('RA [deg]'),pylab.ylabel('DEC [deg]')
        >>> ff.close()
        >>> os.remove('2mass_test.fits')
        
    
    @param source: name of a ViZieR catalog (e.g. 'II/246/out')
    @type source: str
    @keyword filetype: type of the retrieved file 
    @type filetype: str (one of 'tsv','csv','ascii'... see ViZieR site)
    @keyword oc: coordinates
    @type oc: str (one of 'deg'...)
    @keyword out_all: retrieve all or basic information
    @type out_all: boolean (True = all, None = basic)
    @keyword filename: name of the file to write the results to (no extension)
    @type filename: str
    @keyword ID: target name
    @type ID: str
    @keyword ra: target's right ascension
    @type ra: float
    @keyword dec: target's declination
    @type dec: float
    @keyword radius: search radius (arcseconds)
    @type radius: float
    @return: filename / catalog data columns, units, comments
    @rtype: str/ dict, dict, list of str
    """
    #-- options for output format
    filetype = kwargs.get('filetype','tsv')
    oc = kwargs.get('oc','deg')
    oc_eq = kwargs.get('oc_eq','J2000')
    out_all = kwargs.get('out_all',True) # if None, basic
    out_max = kwargs.get('out_max',1000000) # not to unlimited
    filename = kwargs.get('filename',None)
    if filename is not None:
        filename = '%s.%s'%(filename,filetype)
    
    #-- options for quering
    ID = kwargs.get('ID',None)
    ra = kwargs.get('ra',None)
    dec = kwargs.get('dec',None)
    radius = kwargs.get('radius',5.)
    
    #-- gradually build URI
    source_ = urllib.quote(source)
    source_ = source_.replace('/','2%F')
    base_url = 'http://vizier.u-strasbg.fr/viz-bin/asu-%s/VizieR?&-source=%s&-oc=%s,eq=%s'%(filetype,source,oc,oc_eq)
    if out_all: base_url += '&-out.all'
    if out_max: base_url += '&-out.max=%s'%(out_max)
    if radius:  base_url += '&-c.rs=%s'%(radius)
    if ID is not None: base_url += '&-c=%s'%(ID)
    if ra is not None: base_url += '&-c.ra=%s&-c.dec=%s'%(ra,dec)
    
    #-- prepare to open URI
    logger.info('Querying ViZieR source %s'%(source))
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename=filename)
    #   maybe we are just interest in the file, not immediately in the content
    if filename is not None:
        url.close()
        logger.info('Downloaded results to file %s'%(filen))
        return filen
    
    #   otherwise, we read everything into a dictionary
    results = None
    units = {}
    if filetype=='tsv':
        import os
        os.system('cp %s test.csv'%(filen))
        data,comms = ascii.read2array(filen,dtype=np.str,splitchar='\t',return_comments=True)
        
        #-- retrieve the data and put it into a record array
        if len(data)>0:
            #-- now convert this thing into a nice dictionary
            data = np.array(data)
            #-- retrieve the format of the columns. They are given in the
            #   Fortran format. In rare cases, columns contain multiple values
            #   themselves (so called vectors). In those cases, we interpret
            #   the contents as a long string
            formats = np.zeros_like(data[0])
            for line in comms:                  
                line = line.split('\t')
                if len(line)<3: continue
                for i,key in enumerate(data[0]):
                    if key == line[1] and line[0]=='Column': # this is the line with information
                        formats[i] = line[2].replace('(','').replace(')','').lower()
                        if formats[i][0].isdigit(): formats[i] = 'a100'
                        elif 'f' in formats[i]: formats[i] = '>f4' # floating point
                        elif 'i' in formats[i]: formats[i] = '>f4' # integer, but make it float to contain nans
                        elif 'e' in formats[i]: formats[i] = '>f4' # exponential
            #-- define dtypes for record array
            dtypes = np.dtype([(i,j) for i,j in zip(data[0],formats)])
            #-- remove spaces or empty values
            cols = []
            for i,key in enumerate(data[0]):
                 col = data[3:,i]
                 #-- fill empty values with nan
                 cols.append([(row.isspace() or not row) and np.nan or row for row in col])
                 units[key] = data[1,i]
            #-- define columns for record array and construct record array
            cols = [np.cast[dtypes[i]](cols[i]) for i in range(len(cols))]
            results = np.rec.array(cols, dtype=dtypes)
    url.close()
    logger.debug('Results converted to record array')
    return results,units,comms

def vizier2phot(source,results,units,master=None,e_flag='e_',q_flag='q_',extra_fields=None):
    """
    Convert/combine ViZieR record arrays to measurement record arrays.
    
    Every line in the combined array represent a measurement in a certain band.
    
    This is probably only useful if C{results} contains only information on
    one target (or you have the give 'ID' as an extra column).
    
    The standard columns are:
    
        1. C{meas}: containing the photometric measurement
        2. C{e_meas}: the error on the photometric measurement
        3. C{flag}: an optional quality flag
        4. C{unit}: the unit of the measurement
        5. C{photband}: the photometric passband (FILTER.BAND)
        6. C{source}: name of the source catalog
    
    You can add extra information from the VizieR catalog via the list of keys
    C{extra_fields}.
    
    If you give a C{master}, the information will be added to a previous
    record array. If not, a new master will be created.
    
    Colors will be expanded, derived from the other columns and added to the
    master.
    
    The result is a record array with each row a measurement.
    
    Example usage:
    
    First look for all photometry of Vega in all VizieR catalogs:
    
    >>> from ivs.reduction.photometry import calibration
    >>> master = None
    >>> for source in cat_info.sections():
    ...     results,units,comms = search(source,ID='vega',radius=60.)
    ...     if results is not None:
    ...         master = vizier2phot(source,results,units,master,extra_fields=['_r','_RAJ2000','_DEJ2000'])
    
    Keep only observations we have an error of, and convert every observation to 'Jy'
    
    >>> master = master[-np.isnan(master['e_meas'])]
    >>> for i in range(len(master)):
    ...     try:
    ...         myvalue = conversions.convert(master[i]['unit'],'Jy',master[i]['meas'],photband=master[i]['photband'])
    ...     except ValueError:
    ...         myvalue = np.nan
    ...     print '%15s %10.3e+/-%10.3e %15s %10.3e %15s %6.2f %6.2f %6.3f %25s'%(master[i]['photband'],master[i]['meas'],master[i]['e_meas'],master[i]['unit'],myvalue,'Jy',master[i]['_RAJ2000'],master[i]['_DEJ2000'],master[i]['_r'],master[i]['source'])
          JOHNSON.V  3.300e-02+/- 1.200e-02             mag  3.598e+03              Jy 279.23  38.78  0.000           II/168/ubvmeans
        JOHNSON.B-V -1.000e-03+/- 5.000e-03             mag        nan              Jy 279.23  38.78  0.000           II/168/ubvmeans
        JOHNSON.U-B -6.000e-03+/- 6.000e-03             mag        nan              Jy 279.23  38.78  0.000           II/168/ubvmeans
          JOHNSON.B  3.200e-02+/- 1.300e-02             mag  4.054e+03              Jy 279.23  38.78  0.000           II/168/ubvmeans
          JOHNSON.U  2.600e-02+/- 1.432e-02             mag  1.807e+03              Jy 279.23  38.78  0.000           II/168/ubvmeans
           TD1.1965  4.928e-09+/- 1.300e-11      10mW/m2/nm  6.347e+02              Jy 279.23  38.78 18.520            II/59B/catalog
           TD1.1565  5.689e-09+/- 1.700e-11      10mW/m2/nm  4.648e+02              Jy 279.23  38.78 18.520            II/59B/catalog
           TD1.2365  3.700e-09+/- 1.000e-11      10mW/m2/nm  6.903e+02              Jy 279.23  38.78 18.520            II/59B/catalog
           TD1.2740  3.123e-09+/- 9.000e-12      10mW/m2/nm  7.821e+02              Jy 279.23  38.78 18.520            II/59B/catalog
        AKARI.WIDEL  4.047e+00+/- 3.500e-01              Jy  4.047e+00              Jy 279.23  38.78  3.400                II/298/fis
        AKARI.WIDES  6.201e+00+/- 1.650e-01              Jy  6.201e+00              Jy 279.23  38.78  3.400                II/298/fis
         AKARI.N160  3.221e+00+/- 2.550e-01              Jy  3.221e+00              Jy 279.23  38.78  3.400                II/298/fis
          AKARI.N60  6.582e+00+/- 2.090e-01              Jy  6.582e+00              Jy 279.23  38.78  3.400                II/298/fis
         DIRBE.F140  2.557e+02+/- 5.223e+03              Jy  2.557e+02              Jy 279.23  38.78  0.120            J/ApJS/154/673
         DIRBE.F240  8.290e+01+/- 2.881e+03              Jy  8.290e+01              Jy 279.23  38.78  0.120            J/ApJS/154/673
         DIRBE.F4_9  1.504e+02+/- 6.200e+00              Jy  1.504e+02              Jy 279.23  38.78  0.120            J/ApJS/154/673
         DIRBE.F2_2  6.217e+02+/- 9.500e+00              Jy  6.217e+02              Jy 279.23  38.78  0.120            J/ApJS/154/673
         DIRBE.F3_5  2.704e+02+/- 1.400e+01              Jy  2.704e+02              Jy 279.23  38.78  0.120            J/ApJS/154/673
          DIRBE.F12  2.910e+01+/- 1.570e+01              Jy  2.910e+01              Jy 279.23  38.78  0.120            J/ApJS/154/673
          DIRBE.F25  1.630e+01+/- 3.120e+01              Jy  1.630e+01              Jy 279.23  38.78  0.120            J/ApJS/154/673
          DIRBE.F60  1.220e+01+/- 5.610e+01              Jy  1.220e+01              Jy 279.23  38.78  0.120            J/ApJS/154/673
         DIRBE.F100 -7.000e-01+/- 7.790e+01              Jy -7.000e-01              Jy 279.23  38.78  0.120            J/ApJS/154/673
           GENEVA.V  6.100e-02+/- 2.500e-02             mag  3.629e+03              Jy 279.24  38.77 56.000               II/169/main
           GENEVA.B -8.980e-01+/- 2.500e-02             mag  3.837e+03              Jy 279.24  38.77 56.000               II/169/main
           GENEVA.G  1.270e+00+/- 2.500e-02             mag  3.450e+03              Jy 279.24  38.77 56.000               II/169/main
          GENEVA.B2  6.120e-01+/- 2.500e-02             mag  4.152e+03              Jy 279.24  38.77 56.000               II/169/main
           GENEVA.U  6.070e-01+/- 2.500e-02             mag  1.259e+03              Jy 279.24  38.77 56.000               II/169/main
          GENEVA.V1  7.640e-01+/- 2.500e-02             mag  3.673e+03              Jy 279.24  38.77 56.000               II/169/main
          GENEVA.B1  2.000e-03+/- 2.500e-02             mag  3.521e+03              Jy 279.24  38.77 56.000               II/169/main
           IRAS.F60  9.510e+00+/- 8.000e+00              Jy  9.510e+00              Jy 279.23  38.78  9.400               II/125/main
          IRAS.F100  7.760e+00+/- 9.000e+00              Jy  7.760e+00              Jy 279.23  38.78  9.400               II/125/main
           IRAS.F12  4.160e+01+/- 4.000e+00              Jy  4.160e+01              Jy 279.23  38.78  9.400               II/125/main
           IRAS.F25  1.100e+01+/- 5.000e+00              Jy  1.100e+01              Jy 279.23  38.78  9.400               II/125/main
         AKARI.L18W  1.254e+01+/- 7.770e-02              Jy  1.254e+01              Jy 279.24  38.78  2.540                II/297/irc
          AKARI.S9W  5.670e+01+/- 4.010e-01              Jy  5.670e+01              Jy 279.24  38.78  2.540                II/297/irc
          JOHNSON.J        nan+/- 2.100e-01             mag       -nan              Jy 279.23  38.78  0.000   J/PASP/120/1128/catalog
          JOHNSON.K  1.300e-01+/- 1.900e-01             mag  6.051e+02              Jy 279.23  38.78  0.000   J/PASP/120/1128/catalog
          JOHNSON.H        nan+/- 1.500e-01             mag       -nan              Jy 279.23  38.78  0.000   J/PASP/120/1128/catalog
             ANS.25  4.600e-02+/- 4.000e-03             mag        nan              Jy 279.23  38.78 14.000                 II/97/ans
            ANS.15W -4.410e-01+/- 1.200e-02             mag        nan              Jy 279.23  38.78 14.000                 II/97/ans
             ANS.33  1.910e-01+/- 3.000e-03             mag        nan              Jy 279.23  38.78 14.000                 II/97/ans
             ANS.18 -4.620e-01+/- 3.000e-03             mag        nan              Jy 279.23  38.78 14.000                 II/97/ans
            ANS.15N -4.910e-01+/- 1.000e-03             mag        nan              Jy 279.23  38.78 14.000                 II/97/ans
          JOHNSON.B  1.900e-02+/- 1.000e-02             mag  4.102e+03              Jy 279.23  38.78  3.070               I/280B/ascc
          JOHNSON.V  7.400e-02+/- 2.000e-03             mag  3.465e+03              Jy 279.23  38.78  3.070               I/280B/ascc
            2MASS.J -1.770e-01+/- 2.060e-01             mag  1.845e+03              Jy 279.23  38.78  0.034                II/246/out
           2MASS.KS  1.290e-01+/- 1.860e-01             mag  6.083e+02              Jy 279.23  38.78  0.034                II/246/out
            2MASS.H -2.900e-02+/- 1.460e-01             mag  1.067e+03              Jy 279.23  38.78  0.034                II/246/out
                
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
    if cat_info.has_option(source,'q_flag'):
        q_flag = cat_info.get(source,'q_flag')
    if cat_info.has_option(source,'e_flag'):
        e_flag = cat_info.get(source,'e_flag')
    
    #-- basic dtypes
    dtypes = [('meas','>f4'),('e_meas','>f4'),('flag','a20'),
                  ('unit','a30'),('photband','a30'),('source','a50')]
    
    #-- extra can be added:
    names = list(results.dtype.names)
    if extra_fields is not None:
        for e_dtype in extra_fields:
            dtypes.append((e_dtype,results.dtype[names.index(e_dtype)].str))
    
    #-- create empty master if not given
    if master is None:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
    
    #-- add fluxes and magnitudes to the record array
    cols_added = 0
    for key in cat_info.options(source):
        if key in ['e_flag','q_flag']:
            continue
        photband = cat_info.get(source,key)
        #-- contains measurement, error, quality, units, photometric band and
        #   source
        cols = [results[key],
                (e_flag+key in results.dtype.names and results[e_flag+key] or np.ones(len(results))*np.nan),
                (q_flag+key in results.dtype.names and results[q_flag+key] or np.ones(len(results))*np.nan),
                np.array(len(results)*[units[key]],str),
                np.array(len(results)*[photband],str),
                np.array(len(results)*[source],str)]
        #-- add any extra fields if desired.
        if extra_fields is not None:
            for e_dtype in extra_fields:
                cols.append(results[e_dtype])
        #-- add to the master
        rows = []
        for i in range(len(cols[0])):
            rows.append(tuple([col[i] for col in cols]))
        master = np.core.records.fromrecords(master.tolist()+rows,dtype=dtypes)
        cols_added += 1
    
    #-- fix colours: we have to run through it two times to be sure to have
    #   all the colours
    N = len(master)-cols_added
    master_ = breakup_colours(master[N:])
    master_ = breakup_colours(master_)
    #-- combine and return
    master = np.core.records.fromrecords(master.tolist()[:N]+master_.tolist(),dtype=dtypes)
    return master
    
def breakup_colours(master):
    """
    From colors and one magnitude measurement, derive the other magnitudes.
    
    @param master: master record array from vizier2phot.
    @type master: record array
    @return: master with added magnitudes
    @rtype: record array
    """
    names = list(master.dtype.names)
    photbands = list(master['photband'])
    for i,photband in enumerate(photbands):
        system,color = photband.split('.')
        if '-' in color: # we have a colour            
            #-- in which column are the measurements (and error) located?
            index_meas, index_emeas = names.index('meas'),names.index('e_meas')
            index_band = names.index('photband')
            row = list(master[i])
            meas,e_meas = row[index_meas],row[index_emeas]
            
            band1,band2 = ['%s.%s'%(system,band) for band in color.split('-')]
            band1_present = band1 in photbands
            band2_present = band2 in photbands
            
            if band1_present and not band2_present:
                #-- which row do we need to compute the component of the colour?
                index1 = photbands.index(band1)    
                row1 = list(master[index1])
                row1[index_meas] = row1[index_meas] - row[index_meas]
                errs = np.array([row[index_emeas],row1[index_emeas]],float)
                #-- allow for nan errors if all errors on the photbands are nan
                if sum(np.isnan(errs))<len(errs):
                    errs[np.isnan(errs)] = 0.
                row1[index_emeas]= np.sqrt(np.sum(errs**2))
                row1[index_band] = band2
                logger.debug("Added band %s = %s - %s (a)"%(band2,band1,photband))
                master = np.core.records.fromrecords(master.tolist()+[tuple(row1)],dtype=master.dtype)
            elif band2_present and not band1_present:
                #-- which row do we need to compute the component of the colour?
                index1 = photbands.index(band2)    
                row1 = list(master[index1])
                row1[index_meas] = row[index_meas] + row1[index_meas]
                errs = np.array([row[index_emeas],row1[index_emeas]],float)
                #-- allow for nan errors if all errors on the photbands are nan
                if sum(np.isnan(errs))<len(errs):
                    errs[np.isnan(errs)] = 0.
                row1[index_emeas]= np.sqrt(np.sum(errs**2))
                row1[index_band] = band1
                master = np.core.records.fromrecords(master.tolist()+[tuple(row1)],dtype=master.dtype)
                logger.debug("Added band %s = %s + %s (b)"%(band1,band2,photband))
                
    return master
    

def test():
    """
    Execute all docstrings.
    
    >>> import pylab
    >>> p = pylab.show()
    """
    import doctest
    doctest.testmod()

if __name__=="__main__":
    test()