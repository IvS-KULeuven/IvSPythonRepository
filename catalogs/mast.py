# -*- coding: utf-8 -*-
"""
Interface to the MAST archive.

Because the MAST archive is very inhomegeneous, this module is very limited in
use, and sometimes confusing. It is probably best to check the data or retrieve
the data manually from the archive.
"""
import urllib
import socket
import logging
import os
import ConfigParser

import numpy as np
import astropy.io.fits as pf

from ivs.sed import filters
from ivs.aux import loggers
from ivs.aux import numpy_ext
from ivs.inout import ascii
from ivs.inout import http
from ivs.units import conversions
from ivs.catalogs import vizier
from ivs.catalogs import sesame


logger = logging.getLogger("CAT.MAST")
logger.addHandler(loggers.NullHandler())

basedir = os.path.dirname(os.path.abspath(__file__))

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open(os.path.join(basedir,'mast_cats_phot.cfg')))



def _get_URI(name,ID=None,ra=None,dec=None,radius=5.,filetype='CSV',
             out_max=100000,resolver='SIMBAD',coord='dec',**kwargs):
    """
    Build MAST URI from available options.

    Filetype should be one of:

    CSV,SSV,PSV...

    kwargs are to catch unused arguments.

    @param name: name of a MAST mission catalog (e.g. 'hst') or search type
    ('images','spectra')
    @type name: str
    @keyword ID: target name
    @type ID: str
    @param filetype: type of the retrieved file
    @type filetype: str (one of 'CSV','PSV'... see MAST site)
    @keyword radius: search radius (arcseconds)
    @type radius: float
    @param ra: target's right ascension
    @type ra: float
    @param dec: target's declination
    @type dec: float
    @param coord: coordinate output format
    @type coord: str: one of 'sex',dec','dechr'
    @param out_max: maximum number of rows
    @type out_max: integer
    @return: url
    @rtype: str
    """
    base_url = 'http://archive.stsci.edu/'

    if name == 'images':
        base_url += 'siap/search2.php?'
    elif name == 'spectra':
        base_url += 'ssap/search2.php?'
    #elif name == 'galex':
    #    base_url = 'http://galex.stsci.edu/search.php?action=Search'
    else:
        base_url += '%s/search.php?action=Search'%(name)

    #-- only use the ID if ra and dec are not given
    if ID is not None and (ra is None and dec is None):
        base_url += '&target=%s'%(ID)

    if ra is not None and dec is not None:
        base_url += '&radius=%s'%(radius/60.)
        if ra is not None: base_url += '&RA=%s'%(ra)
        if dec is not None: base_url += '&DEC=%s'%(dec)
    elif radius is not None:
        base_url += '&SIZE=%s'%(radius/60.)

    for kwarg in kwargs:
        base_url += '&%s=%s'%(kwarg,kwargs[kwarg])

    if name not in ['spectra','images']:
        base_url += '&outputformat=%s'%(filetype)
        base_url += '&resolver=%s'%(resolver)
        #base_url += '&max_records=%d'%(out_max)
        base_url += '&coordformat=%s'%(coord)
    logger.debug('url: %s'%(base_url))
    return base_url



def galex(**kwargs):
    """
    Cone search for Galex targets.
    """
    ID = kwargs.pop('ID',None)
    if ID is None:
        ra,dec = kwargs.pop('ra'),kwargs.pop('dec')
    else:
        info = sesame.search(ID,db='A')
        if not 'jradeg' in info: return None,None,None
        ra,dec = info['jradeg'],info['jdedeg']
    radius = 0.1
    #radius = radius/60.
    base_url = 'http://galex.stsci.edu/gxws/conesearch/conesearch.asmx/ConeSearchToXml?ra={0:f}&dec={1:f}&sr={2:f}&verb=1'.format(ra,dec,radius)
    #base_url = 'http://galex.stsci.edu/GR4/?page=searchresults&RA={ra:f}&DEC={dec:f}&query=no'.format(ra=ra,dec=dec)
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename=None)
    fuv_flux,e_fuv_flux = None,None
    columns = ['_r','ra','dec','fuv_flux','fuv_fluxerr','nuv_flux','nuv_fluxerr']
    values = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

    #flux in microJy
    units = dict(fuv_flux='muJy',nuv_flux='muJy')

    got_target = False
    with open(filen,'r') as ff:
        for line in ff.readlines():
            for i,col in enumerate(columns):
                if col+'>' in line:
                    values[i] = np.float(line.split('>')[1].split('<')[0])
                    got_target = (col=='fuv_fluxerr')
            if got_target:
                break

    values[0] = np.sqrt( (values[1]-ra)**2 + (values[2]-dec)**2)*3600
    columns[1] = '_RAJ2000'
    columns[2] = '_DEJ2000'

    results = np.rec.fromarrays(np.array([values]).T,names=columns)
    if np.all(np.isnan(np.array(values))):
        results = None

    return results,units,None




def search(catalog,**kwargs):
    """
    Search and retrieve information from a MAST catalog.

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


    @param catalog: name of a MAST mission catalog
    @type catalog: str
    @keyword filename: name of the file to write the results to (no extension)
    @type filename: str
    @return: filename / catalog data columns, units, comments
    @rtype: str/ record array, dict, list of str
    """
    filename = kwargs.pop('filename',None) # remove filename from kwargs
    filetype = kwargs.setdefault('filetype','CSV')
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
        logger.info('Querying MAST source %s and downloading to %s'%(catalog,filename))
        url.close()
        return filen

    #   otherwise, we read everything into a dictionary

    if filetype=='CSV' and not filename:
        try:
            results,units,comms = csv2recarray(filen)
        #-- raise an exception when multiple catalogs were specified
        except ValueError:
            url.close()
            #raise ValueError, "failed to read %s, perhaps multiple catalogs specified (e.g. III/168 instead of III/168/catalog)"%(catalog)
            results,units,comms = None,None,None
        url.close()
        logger.info('Querying MAST source %s (%d)'%(catalog,(results is not None and len(results) or 0)))
        return results,units,comms
    else:
        return filename


def mast2phot(source,results,units,master=None,extra_fields=None):
    """
    Convert/combine MAST record arrays to measurement record arrays.

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

    You can add extra information from the Mast catalog via the list of keys
    C{extra_fields}.

    If you give a C{master}, the information will be added to a previous
    record array. If not, a new master will be created.

    The result is a record array with each row a measurement.

    @param source: name of the Mast source
    @type source: str
    @param results: results from Mast C{search}
    @type results: record array
    @param units: header of Mast catalog with key name the column name and
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

    #-- MAST has no unified terminology:
    translations = {'kepler/kgmatch':{'_r':"Ang Sep (')",
                                      '_RAJ2000':'KIC RA (J2000)',
                                      '_DEJ2000':'KIC Dec (J2000)'}}

    #-- extra can be added:
    names = list(results.dtype.names)
    if extra_fields is not None:
        if source in translations:
            translation = translations[source]
        else:
            translation = {'_r':'_r','_RAJ2000':'_RAJ2000','_DEJ2000':'_DEJ2000'}
        for e_dtype in extra_fields:
            try:
                dtypes.append((e_dtype,results.dtype[names.index(translation[e_dtype])].str))
            except ValueError:
                if e_dtype != '_r': raise
                dtypes.append((e_dtype,results.dtype[names.index('AngSep')].str))


    #-- create empty master if not given
    newmaster = False
    if master is None or len(master)==0:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
        newmaster = True

    #-- add fluxes and magnitudes to the record array
    cols_added = 0
    for key in cat_info.options(source):
        if key[:2] in [e_flag,q_flag] or '_unit' in key or key=='bibcode':
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
    master_ = vizier._breakup_colours(master_)
    #-- combine and return
    master = np.core.records.fromrecords(master.tolist()[:N]+master_.tolist(),dtype=dtypes)

    #-- skip first line from building
    if newmaster: master = master[1:]
    return master



def csv2recarray(filename):
    """
    Read a MAST csv (comma-sep) file into a record array.

    @param filename: name of the TCSV file
    @type filename: str
    @return: catalog data columns, units, comments
    @rtype: record array, dict, list of str
    """
    data,comms = ascii.read2array(filename,dtype=np.str,splitchar=',',return_comments=True)
    results = None
    units = {}
    #-- retrieve the data and put it into a record array
    if len(data)>1:
        #-- now convert this thing into a nice dictionary
        data = np.array(data)
        #-- retrieve the format of the columns. They are given in the
        #   Fortran format. In rare cases, columns contain multiple values
        #   themselves (so called vectors). In those cases, we interpret
        #   the contents as a long string
        formats = np.zeros_like(data[0])
        for i,fmt in enumerate(data[1]):
            if 'string' in fmt or fmt=='datetime': formats[i] = 'a100'
            if fmt=='integer': formats[i] = 'f8'
            if fmt=='ra': formats[i] = 'f8'
            if fmt=='dec': formats[i] = 'f8'
            if fmt=='float': formats[i] = 'f8'
        #-- define dtypes for record array
        dtypes = np.dtype([(i,j) for i,j in zip(data[0],formats)])
        #-- remove spaces or empty values
        cols = []
        for i,key in enumerate(data[0]):
             col = data[2:,i]
             #-- fill empty values with nan
             cols.append([(row.isspace() or not row) and np.nan or row for row in col])
             #-- fix unit name
             for source in cat_info.sections():
                if cat_info.has_option(source,data[0,i]+'_unit'):
                    units[key] = cat_info.get(source,data[0,i]+'_unit')
                    break
             else:
                units[key] = 'nan'
        #-- define columns for record array and construct record array
        cols = [np.cast[dtypes[i]](cols[i]) for i in range(len(cols))]
        results = np.rec.array(cols, dtype=dtypes)
    else:
        results = None
        units = {}
    return results,units,comms



def get_photometry(ID=None,extra_fields=['_r','_RAJ2000','_DEJ2000'],**kwargs):
    """
    Download all available photometry from a star to a record array.

    For extra kwargs, see L{_get_URI} and L{mast2phot}
    """
    to_units = kwargs.pop('to_units','erg/s/cm2/AA')
    master_ = kwargs.get('master',None)
    master = None
    #-- retrieve all measurements
    for source in cat_info.sections():
        if source=='galex':
            results,units,comms = galex(ID=ID,**kwargs)
        else:
            results,units,comms = search(source,ID=ID,**kwargs)

        if results is not None:
            master = mast2phot(source,results,units,master,extra_fields=extra_fields)

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
            except AssertionError: # postive flux and errors!
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


#@retry_http(3)
def get_dss_image(ID,ra=None,dec=None,width=5,height=5):
    """
    Retrieve an image from DSS

    plot with

    >>> data,coords,size = mast.get_image('HD21389')
    >>> pl.imshow(data[::-1],extent=[coords[0]-size[0]/2,coords[0]+size[0]/2,
                                    coords[1]-size[1]/2,coords[1]+size[1]/2])
    """
    #-- set a reasonable timeout
    timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30.)
    if ra is None or dec is None:
        info = sesame.search(ID)
        ra,dec = info['jradeg'],info['jdedeg']
    url  = urllib.URLopener()
    myurl = "http://archive.stsci.edu/cgi-bin/dss_search?ra=%s&dec=%s&equinox=J2000&height=%s&generation=%s&width=%s&format=FITS"%(ra,dec,height,'2i',width)
    out = url.retrieve(myurl)
    data1 = pf.getdata(out[0])
    #-- reset timeout to original value
    socket.setdefaulttimeout(timeout)
    return data1,(ra,dec),(width/60.,height/60.)



def get_FUSE_spectra(ID=None,directory=None,cat_info=False,select=['ano']):
    """
    Get Fuse spectrum.

    select can be ['ano','all']

    # what about MDRS aperture?
    """
    direc = (directory is None) and os.getcwd() or directory
    if not os.path.isdir(direc):
        os.mkdir(direc)
    data = search('fuse',ID=ID)
    download_link = 'https://archive.stsci.edu/pub/fuse/data/{folder}/{filename}'
    #-- maybe no information was found
    if data[0] is None:
        return None
    elif cat_info:
        return data
    output = []
    for spectrum in data[0]:
        folder = spectrum['Data ID'][:8]
        aperture = spectrum['Aperture']
        noexp = int(spectrum['No. Exps'])
        an = dict(MDWRS=2,LWRS=4)
        if not aperture in an:
            continue
        an = an[aperture]
        for exp in range(noexp):
            data_id = folder+'{0:03d}00'.format(exp)
            for type in select:
                filename = '{data_id}{type}{an}ttagfcal.fit.gz'.format(data_id=data_id,type=type,an=an)
                link = download_link.format(folder=folder,filename=filename)
                link = link.lower()
                filename = os.path.join(direc,filename)
                myfile = http.download(link,filename=filename)
                #-- if the file is too small to be a science file, it is
                #   actually an HTML file stating that the URL does not exist.
                if os.path.getsize(myfile)<1000:
                    os.unlink(myfile)
                    continue
                logger.info('FUSE spectrum %s to %s'%(link,direc))
                output.append(myfile)
    return output
    #https://archive.stsci.edu/pub/fuse/data/a0010101/a001010100000all4ttagfcal.fit.gz

if __name__=="__main__":
    #get_FUSE_spectrum(ID='hd163296')
    print get_FUSE_spectrum(ID='hd163296').dtype.names

    raise SystemExit
    mission = 'fuse'
    base_url = _get_URI(mission,ID='hd163296')
    print base_url
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename='%s.test'%(mission))
    url.close()
    raise SystemExit

    missions = ['hst',# - Hubble Space Telescope
                'kepler',# - Kepler Data search form
                'kepler/kepler_fov',# - Kepler Target search form
                'kepler/kic10',# - Kepler Input Catalog search form
                'kepler/kgmatch',# - Kepler/Galex Cross match
                'iue',# - International Ultraviolet Explorer
                'hut',# - Hopkins Ultraviolet Telescope
                'euve',# - Extreme Ultraviolet Explorer
                'fuse',# - Far-UV Spectroscopic Explorer
                'uit',# - Ultraviolet Imageing Telescope
                'wuppe',# - Wisconsin UV Photo-Polarimeter Explorer
                'befs',# - Berkeley Extreme and Far-UV Spectrometer
                'tues',# - Tübingen Echelle Spectrograph
                'imaps',# - Interstellar Medium Absorption Profile Spectrograph
                'hlsp',# - High Level Science Products
                'pointings',# - HST Image Data grouped by position
                'copernicus',# - Copernicus Satellite
                'hpol',# - ground based spetropolarimeter
                'vlafirst',# - VLA Faint Images of the Radio Sky (21-cm)
                'xmm-om']# - X-ray Multi-Mirror Telescope Optical Monitor

    #search('euve',ID='hd46149')
    #filename = search('kepler/kgmatch',ID='KIC3955868',filename='kepler.test')

    #results,units,comms = search('kepler/kgmatch',ID='T CrB')
    #sys.exit()
    out = search('kepler/kepler_fov',ID='TYC3134-165-1',filename='kepler_fov.test',radius=1.)
    out = search('kepler/kgmatch',ID='3749404',filename='kgmatch.test',radius=1.)
    #results,units,comms = search('kepler/kgmatch',ID='3749404')
    master = mast2phot('kepler/kgmatch',results,units,master=None,extra_fields=None)
    print master
    #data,units,comms = search('spectra',ID='hd46149',filename='ssap.test')
    sys.exit()

    for mission in missions:
        base_url = _get_URI(mission,ID='hd46149')
        #ff = urllib.urlopen(base_url)
        try:
            url = urllib.URLopener()
            filen,msg = url.retrieve(base_url,filename='%s.test'%(mission))
            url.close()
        except IOError:
            print 'failed'
            continue


