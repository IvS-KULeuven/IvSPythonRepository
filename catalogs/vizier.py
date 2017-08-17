# -*- coding: utf-8 -*-
"""
Interface to the VizieR website.

Download or retrieve VizieR catalogs.

The basic interface C{search} lets you download B{entire catalogs or parts} of
them. The structure array containts then one row per target, the columns
denoting the different columns of the catalogs. You can also specify two
catalogs in C{xmatch}; the second will then be cross-matched against the first.

The convenience functions (C{get_photometry},...) aim to collect information
from different catalogs on B{one target} in one array. Each row represents
one measurement from one catalog (typically you'll have many rows both from one
catalog but also from different catalogs). The columns then denote the contents
of each row (e.g. the magnitude, photometric passband etc).

Section 1. Download catalogs
============================

Section 1.1. To a file
----------------------

Download the entire Van Leeuwen Hipparcos New Reduction catalog to a file. The
filename is returned as a check for success.

>>> filename = search('I/311/hip2',filename='vanleeuwen.tsv')

Download a 60 arcsec circular area from the catalog around the coordinates
ra=237.1, dec=-10.10

>>> filename = search('I/311/hip2',ra=237.1,dec=-10.10,radius=60.,filename='vanleeuwen.tsv')

Search for the presence of a target in the catalog. The downloaded file will 
contain no rows if the target is not in the catalog. If more than one target are
in the search radius around the target, there will be more than one row. They
are ordered via distance to the target, so it's probably the first one you need.

>>> filename = search('I/311/hip2',ID='vega',filename='vanleeuwen.tsv')
>>> filename = search('I/311/hip2',ID='vega',filename='vanleeuwen.tsv',radius=60.)

Section 1.2 To a RecordArray
----------------------------

Instead of downloading to a file and then reading in the file for further
analysis, you can download the contents of the file directly to a record array,
retrieving the units and comments from the catalog in one go. The equivalent of
the third example above becomes

>>> rec_arr,unit_dict,comment_str = search('I/311/hip2',ID='vega')

With these record arrays, its very easy to select targets based on criteria.
For example, if you want to extract 2MASS targets in a certain area with a
negative H-K index, you can do

>>> data,units,comms = search('II/246/out',ra=237.1,dec=-10.10,radius=600.)
>>> selection = (data['Hmag'] - data['Kmag']) < 0
>>> data = data[selection] 

You can also read in a data file you've previously downloaded via

>>> data,units,comms = tsv2recarray('vanleeuwen.tsv')

Section 1.3 List relevant catalogs
----------------------------------

To know in which catalogs your target is present, list them all via

>>> my_cats = list_catalogs('vega')

Now you could iterate over them and download them to a file or whatever.

Section 2. Convenience functions
================================

You can define 'standard' photometric catalogs in the C{vizier_cats.cfg} file.
This file is basically a translator for VizieR column headers to photometric
passbands (and colors). For examples, see the file itself.

You can add catalogs on the fly via

>>> cat_info.add_section('my_new_catalog')
>>> cat_info.set('my_new_catalog','Bmag','JOHNSON.B')
"""
#-- standard libraries
import numpy as np
import urllib
import logging
import os
import itertools
import astropy.io.fits as pf
import tarfile
import tempfile
import shutil
import ConfigParser
from scipy.spatial import KDTree

#-- IvS repository
from ivs.inout import ascii
from ivs.inout import fits
from ivs.inout import http
from ivs.units import conversions
from ivs.aux import loggers
from ivs.aux import numpy_ext
from ivs.sed import filters
from ivs import config

logger = logging.getLogger("CAT.VIZIER")
logger.addHandler(loggers.NullHandler())

basedir = os.path.dirname(os.path.abspath(__file__))

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open(os.path.join(basedir,'vizier_cats_phot.cfg')))

cat_info_fund = ConfigParser.ConfigParser()
cat_info_fund.optionxform = str # make sure the options are case sensitive
cat_info_fund.readfp(open(os.path.join(basedir,'vizier_cats_fund.cfg')))

mirrors = {'cycle': itertools.cycle(['vizier.u-strasbg.fr',      # France
                                     'vizier.nao.ac.jp',         # Japan
                                     'vizier.hia.nrc.ca',        # Canada
                                     'vizier.ast.cam.ac.uk',     # UK
                                     'vizier.cfa.harvard.edu',   # USA CFA/harvard
                                     'www.ukirt.jach.hawaii.edu',# USA Ukirt
                                     'vizier.inasan.ru',         # Russia
                                     'vizier.iucaa.ernet.in',    # India
                                     'data.bao.ac.cn'])}        # China
mirrors['current'] = mirrors['cycle'].next()

#{ Basic interfaces

def change_mirror():
    """
    Cycle through the mirrors of ViZieR.
    """
    mirrors['current'] = mirrors['cycle'].next()
    logger.info("Changed cycle to {}".format(mirrors['current']))

def search(name,filetype='tsv',filename=None,**kwargs):
    """
    Search and retrieve information from a VizieR catalog.
    
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
    dictionaries with the same keys, depending on the column names in the
    catalog). The entries in the dictionary are of type C{ndarray}, and will
    be converted to a float-array (no integers, we need to support nans) if
    possible. If not, the array will consist of strings. The comments are also
    returned as a list of strings.
    
    WARNING: when retrieving a FITS file, ViZieR sometimes puts weird formats
    into the header ('3F10.6E' in the 2MASS catalog), which cannot be read by
    the C{astropy.io.fits} module. These columns are actually multi-dimensional vectors.
    One option is to download to another format, or to restrict the columns with
    C{out_all=None}.
    
    Example usage:
    
        1. Look for the Geneva V magnitude and coordinates of Vega in the GENEVA
        catalog of Rufener.
    
        >>> results,units,comms = search('II/169/main',ID='vega',radius=60.)
        >>> print "Vega: Vmag = %.3f %s, RA = %.2f %s, DEC = %.2f %s"%(results['Vmag'],units['Vmag'],results['_RAJ2000'],units['_RAJ2000'],results['_DEJ2000'],units['_DEJ2000'])
        Vega: Vmag = 0.061 mag, RA = 279.24 deg, DEC = 38.77 deg
    
        2. Search for all targets in the 2MASS catalog in a particular field.
        Download the results to a FITS file, read the file, and plot the results
        to the screen.
    
        >>> #filename = search('II/246/out',ra=100.79,dec=0.70,radius=1000.,filetype='fits',filename='2mass_test',out_all=None)
    
        Now read in the FITS-file and plot the contents
        
        >>> #import pylab
        >>> #import astropy.io.fits as pf
        >>> #ff = pf.open('2mass_test.fits')
        >>> #p = pylab.gcf().canvas.set_window_title('test of <search>')
        >>> #p = pylab.scatter(ff[1].data.field('_RAJ2000'),ff[1].data.field('_DEJ2000'),c=ff[1].data.field('Jmag'),s=(20-ff[1].data.field('Jmag'))**2,cmap=pylab.cm.hot_r,edgecolors='none')
        >>> #p = pylab.colorbar()
        >>> #p = p.set_label('Jmag')
        >>> #p,q = pylab.xlabel('RA [deg]'),pylab.ylabel('DEC [deg]')
        >>> #ff.close()
        >>> #os.remove('2mass_test.fits')
        
    
    @param name: name of a ViZieR catalog (e.g. 'II/246/out')
    @type name: str
    @param filetype: type of the file to write the results to ('tsv' if no file desired)
    @type filetype: string (one of 'tsv','fits','ascii','csv'...
    @param filename: name of the file to write the results to (no extension)
    @type filename: str
    @return: filename / catalog data columns, units, comments
    @rtype: str/ record array, dict, list of str
    """
    #-- two ways of giving filename: either complete filename with extension,
    #   or filename without extension, but C{filetype} speficied.
    if filename is not None and '.' in os.path.basename(filename):
        filetype = os.path.splitext(filename)[1][1:]
    elif filename is not None:
        filename = '%s.%s'%(filename,filetype)
    
    #-- gradually build URI
    base_url = _get_URI(name=name,**kwargs)
    
    #-- prepare to open URI
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename=filename)
    #   maybe we are just interest in the file, not immediately in the content
    if filename is not None:
        logger.info('Querying ViZieR source %s and downloading to %s'%(name,filen))
        url.close()
        return filen
    
    #   otherwise, we read everything into a dictionary
    if filetype=='tsv':
        try:
            results,units,comms = tsv2recarray(filen)
        #-- raise an exception when multiple catalogs were specified
        except ValueError:
            raise ValueError, "failed to read %s, perhaps multiple catalogs specified (e.g. III/168 instead of III/168/catalog)"%(name)
        url.close()
        logger.info('Querying ViZieR source %s (%d)'%(name,(results is not None and len(results) or 0)))
        return results,units,comms
    

def list_catalogs(ID,filename=None,filetype='tsv',**kwargs):
    """
    Print and return all catalogs containing information on the star.
    
    If you give C{filetype} and C{filename}, all information will be downloaded
    to that file.
    
    Extra kwargs: see L{_get_URI}.
    
    @param ID: identification of the star
    @type ID: str
    @keyword filetype: type of the output file ('fits','tsv','csv'...)
    @type filetype: str
    @keyword filename: name of the output file
    @type filename: str
    @return: dictionary with keys the catalog ID from VizieR and entries the
    names of the catalogs
    @rtype: dictionary
    """
    base_url = _get_URI(ID=ID,filetype='fits',**kwargs)
    
    #-- download the file
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url,filename=filename)
    
    #-- if it is a FITS file, we extract all catalog IDs. We download the
    #   individual catalogs to retrieve their title.
    if filetype=='fits':
        mycats = {}
        ff = pf.open(filen)
        for ext in range(1,len(ff)):            
            name = ff[ext].header['CDS-name']
            results,units,comms = search(name,ID=ID,**kwargs)
            for line in comms:
                if 'Title:' in line:
                    title = line.strip().split(':')[1]
                    break
            mycats[name] = title
            logger.info('%25s %s'%(name,title))
            
            photometry = [col for col in units.keys() if 'mag' in units[col]]
            rv = [col for col in units.keys() if 'rv' in col.lower()]
            vsini = [col for col in units.keys() if 'sin' in col.lower()]
            sptype = [col for col in units.keys() if col.lower()=='sp' or col.lower()=='sptype']
            fund = [col for col in units.keys() if 'logg' in col.lower() or 'teff' in col.lower()]
            
        ff.close()
        url.close()
        return mycats
    else:
        url.close()
        return filen        


def xmatch(source1,source2,output_file=None,tol=1.,**kwargs):
    """
    Crossmatch two vizier catalogs via a fast KDTree.
    
    The limit for these catalogs is probably somewhere between ~100000 entries,
    so make sure your catalogs do not contain to many targets. You can always
    do a subselection via the keyword arguments (e.g. give ra, dec and radius).
    
    An output tsv file will be written (by default named 'source1__source2',
    which can be read in via C{tsv2recarray} for further analysis.
    
    tolerance is in arcseconds.
    
    Extra keyword arguments are passed to C{search}. Column names of second
    source will be appended with postfix '_2', to avoid clashes of double-defined
    column headers.
    """
    #-- construct default filename.
    if output_file is None:
        output_file = "__".join([source1,source2]).replace('/','_').replace('+','')+'.tsv'
    
    #-- download the two catalogs
    cat1,units1,comms1 = search(source1,**kwargs)
    cat2,units2,comms2 = search(source2,**kwargs)
    
    logger.info('Start Vizier Xmatch')
    coords1 = np.array([cat1['_RAJ2000'],cat1['_DEJ2000']]).T
    coords2 = np.array([cat2['_RAJ2000'],cat2['_DEJ2000']]).T
    
    logger.info('Building KDTree of shape %d,%d'%coords1.shape)
    tree = KDTree(coords1)
    
    logger.info('Querying KDTree with %d entries'%(len(coords2)))
    distance,order = tree.query(coords2)
    
    keep = distance<(tol/(60.))
    
    logger.info('Matched %d points (tol<%.3g arcsec)'%(sum(keep),tol))

    #-- this the subset matching both catalogues
    cat1 = cat1[order[keep]]
    cat2 = cat2[keep]
    
    #-- now write it to a vizier-like file...
    #-- first append '2' to each column name of the second source in the
    #   comments, to make sure there are no doubles.
    for i,line in enumerate(comms2):
        line = line.split('\t')
        if line[0]=='Column':
            line[1] = line[1]+'_2'
            comms2[i] = '\t'.join(line)
    #-- change unit dictionaries and dtypes
    units2_ = {}
    for key in units2: units2_[key+'_2'] = units2[key]
    cat2.dtype.names = [name+'_2' for name in cat2.dtype.names]
    
    
    ff = open(output_file,'w')
    ff.write('\n#'.join(comms1))
    ff.write('\n#'.join(comms2))
    
    names1 = list(cat1.dtype.names)
    names2 = list(cat2.dtype.names)
    dtypes = [(name,cat1.dtype[names1.index(name)].str) for name in names1]
    dtypes += [(name,cat2.dtype[names2.index(name)].str) for name in names2]
    
    ff.write('\n')
    for nr,i in enumerate(dtypes):
        ff.write(str(i[0]))
        if nr<(len(dtypes)-1): ff.write('\t')
    ff.write('\n')
    for nr,i in enumerate(dtypes):
        if i[0]   in units1:   ff.write(units1[i[0]])
        elif i[0] in units2_:  ff.write(units2_[i[0]])
        else:
            raise ValueError,'this cannot be'
        if nr<(len(dtypes)-1): ff.write('\t')
        
    ff.write('\n')
    ff.write('\t'.join(['---']*len(dtypes)))
    ff.write('\n')
    
    for row1,row2 in itertools.izip(cat1,cat2):
        ff.write('\t'.join([str(x) for x in row1])+'\t')
        ff.write('\t'.join([str(x) for x in row2])+'\n')
    
    ff.close()    




#}

#{ Interface to specific catalogs

def get_IUE_spectra(ID=None,directory=None,unzip=True,cat_info=False,select='low',**kwargs):
    """
    Download IUE spectra.
    
    If you want to download all the spectral files, set C{directory='/home/user/'}
    or whatever. All the tarfiles will be downloaded to this directory, they
    will be untarred, science data extracted and all unnecessary files and
    directories will be deleted. If you don't set a directory, it will default
    to the CWD.
    
    If you don't wish to unzip them, set unzip=False
    
    DEPRECATED: If you don't give a directory, the function will return a list
    of all extracted spectra (no data files are kept).
    
    You can retrieve the contents of the vizier catalog via {cat_info=True}. The
    files will not be downloaded in this case.
    """
    if directory is None:
        direc = os.getcwd()
        directory = os.getcwd()
        filename = None
    else:
        direc = directory
        if not os.path.isdir(direc):
            os.mkdir(direc)
        
    output = []
    #-- construct the download link form the camera and image data
    data,units,comments = search('VI/110/inescat/',ID=ID,**kwargs)
    
    if cat_info:
        return data,units,comments
    
    if data is None:
        return output
    
    for spectrum in data:
        download_link = "http://archive.stsci.edu/cgi-bin/iue_retrieve?iue_mark=%s%05d&mission=iue&action=Download_MX"%(spectrum['Camera'].strip(),int(spectrum['Image']))
        logger.info('IUE spectrum %s/%s: %s'%(spectrum['Camera'],spectrum['Image'],download_link))
    
        #-- prepare to download the spectra to a temparorary file
        if directory is not None:
            filename = download_link.split('iue_mark=')[1].split('&')[0]
            filename = os.path.join(direc,filename)
        #-- download file and retrieve the path to the downloaded file
        mytarfile = http.download(download_link,filename=filename)
        if filename is None:
            mytarfile,url = mytarfile
        #-- perhaps unzipping is not necessary
        if directory is not None and not unzip:
            output.append(mytarfile)
            url.close()
            continue
        #-- else unpack tar files and copy spectra to the directory,
        #   remember the location of the files:
        if not tarfile.is_tarfile(mytarfile):
            logger.info("Not a tarfile: %s"%(mytarfile))
            continue
        tarf = tarfile.open(mytarfile)
        files = tarf.getmembers(),tarf.getnames()
        deldirs = []
        outfile = None
        for mem,name in zip(*files):
            #-- first check if it is a spectrum, but skip low or high res if needed.
            myname = name.lower()
            skip = False
            if 'lo' in select and not 'lo' in os.path.basename(myname):
                skip = True
            if 'hi' in select and not 'hi' in os.path.basename(myname):
                skip = true
            if ('lwp' in name or 'swp' in name or 'lwr' in name) and not skip:
                #-- first unpack this file
                tarf.extract(mem,path=direc)
                #-- move it to the direc directory
                outfile = os.path.join(direc,os.path.basename(name))
                dirname = os.path.dirname(name)
                shutil.move(os.path.join(direc,name),outfile)
                logger.info("Extracted %s to %s"%(name,outfile))
                #-- remove any left over empty directory
                deldirs.append(os.path.join(direc,dirname))
            else:
                logger.debug("Did not extract %s (probably not a spectrum)"%(name))
        #-- remove the tar file:
        tarf.close()
        os.unlink(mytarfile)
        for dirname in deldirs:
            if dirname and os.path.isdir(dirname) and not os.listdir(dirname):
                os.rmdir(dirname)
                logger.debug("Deleted left over (empty) directory %s"%(dirname))
        if filename is None: url.close()
        
        #-- only read in the data if they need to be extracted
        if directory is not None and outfile:
            output.append(outfile)
            continue
        if outfile and os.path.isfile(outfile):
            wavelength,flux,error,header = fits.read_iue(outfile,return_header=True)
            os.unlink(outfile)
            output.append([wavelength,flux,error,header])
        else:
            logger.info('Unsuccesfull extraction of %s'%(outfile))
        
    return output
    
def get_UVSST_spectrum(units='erg/s/cm2/AA',**kwargs):
    """
    Get a spectrum from the UVSST spectrograph onboard TD1.
    
    From vizier catalog III/39A.
    
    Also have a look at II/86/suppl.
    """
    kwargs.setdefault('out_max',10)
    kwargs.setdefault('radius',60.)
    data,units_o,comments = search('III/39A/catalog',**kwargs)
    if data is None: return None,None
    fields = sorted([field for field in data.dtype.names if (field[0]=='F' and len(field)==5)])
    spectrum = np.zeros((len(fields),3))
    units_ = []
    for i,field in enumerate(fields):
        wave,flux = float(field[1:]),data[field][0]
        if flux==0: continue
        #-- there is an error in the units ***in vizier***!
        flux = conversions.convert(units_o[field],units,100*flux,wave=(wave,'AA'),unpack=True)
        e_flux = data['e_'+field][0]/100.*flux
        spectrum[i][0] = wave
        spectrum[i][1] = flux
        spectrum[i][2] = e_flux
        units_.append(units_o[field])
    return spectrum[spectrum[:,0]>0]
#}


#{ Convenience functions
def get_photometry(ID=None,extra_fields=['_r','_RAJ2000','_DEJ2000'],take_mean=False,**kwargs):
    """
    Download all available photometry from a star to a record array.
    
    For extra kwargs, see L{_get_URI} and L{vizier2phot}
    
    """
    kwargs['ID'] = ID
    to_units = kwargs.pop('to_units','erg/s/cm2/AA')
    sources = kwargs.get('sources',cat_info.sections())
    master_ = kwargs.get('master',None)
    master = None
    #-- retrieve all measurements
    for source in sources:
        results,units,comms = search(source,**kwargs)
        if results is None: continue
        master = vizier2phot(source,results,units,master,extra_fields=extra_fields,take_mean=take_mean)
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
            to_units_ = to_units+''
            try:
                value,e_value = conversions.convert(master['unit'][i],to_units,master['meas'][i],master['e_meas'][i],photband=master['photband'][i])
            except ValueError: # calibrations not available, or its a color
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
    
    if master_ is not None and master is not None:
        master = numpy_ext.recarr_addrows(master_,master.tolist())
    elif master is None:
        master = master_
    
    #-- and return the results
    return master


def quality_check(master,ID=None,return_master=True,**kwargs):
    """
    Perform quality checks on downloaded data.
    
    This function translates flags in to words, and looks up additional
    information in selected catalogs.
    """
    messages = ['' for i in range(len(master))]
    #-- for some, we can do a quality check given the information that is
    #   already available in the master record
    #-- IRAS
    logger.info('Checking flags')
    for i,entry in enumerate(master):
        if 'USNO' in str(entry['photband']):
            messages[i] = '; '.join([messages[i],'unreliable zeropoint and transmission profile'])
            continue
        if 'COUSINS' in str(entry['photband']):
            messages[i] = '; '.join([messages[i],'unreliable zeropoint'])
            continue
        flag = entry['flag']
        try:
            flag = float(flag)
        except:
            continue
        if np.isnan(flag): continue
        if entry['source'].strip() in ['II/125/main','II/275/fsr']:
            flag = float(int(flag))
            if flag==3: messages[i] = '; '.join([messages[i],'high quality'])
            if flag==2: messages[i] = '; '.join([messages[i],'moderate quality'])
            if flag==1: messages[i] = '; '.join([messages[i],'upper limit'])
        if entry['source'].strip() in ['II/298/fis','II/297/irc']:
            flag = float(int(flag))
            if flag==3: messages[i] = '; '.join([messages[i],'high quality'])
            if flag==2: messages[i] = '; '.join([messages[i],'source is confirmed but the flux is not reliable'])
            if flag==1: messages[i] = '; '.join([messages[i],'the source is not confirmed'])
            if flag==0: messages[i] = '; '.join([messages[i],'not observed'])
        if 'DIRBE' in entry['source'].strip():
            flag = '{0:03d}'.format(int(float(flag)))
            if flag[0]=='1': messages[i] = '; '.join([messages[i],'IRAS/2MASS companion greater than DIRBE noise level'])
            if flag[1]=='1': messages[i] = '; '.join([messages[i],'possibly extended emission or highly variable source']) # discrepancy between DIRBE and IRAS/2MASS flux density
            if flag[2]=='1': messages[i] = '; '.join([messages[i],'possibly affected by nearby companion'])
        if entry['source'].strip() in ['V/114/msx6_main']:
            if flag==4: messages[i] = '; '.join([messages[i],'excellent'])
            if flag==3: messages[i] = '; '.join([messages[i],'good'])
            if flag==2: messages[i] = '; '.join([messages[i],'fair'])
            if flag==1: messages[i] = '; '.join([messages[i],'on limit'])
            if flag==0: messages[i] = '; '.join([messages[i],'not detected'])
    
    #-- for other targets, we need to query additional information
    sources_with_quality_check = ['B/denis/denis','II/311/wise','II/246/out']
    sources = set(list(master['source'])) & set(sources_with_quality_check)
    
    denis_image_flags = {'01':'clouds during observation',
                         '02':'electronic Read-out problem',
                         '04':'internal temperature problem',
                         '08':'very bright star',
                         '10':'bright star',
                         '20':'stray light',
                         '40':'unknown problem'}
    denis_source_flags = {'01':'source might be a dust on mirror',
                          '02':'source is a ghost detection of a bright star',
                          '04':'source is saturated',
                          '08':'source is multiple detect',
                          '10':'reserved'}
    wise_conf_flag = {'d':'diffraction spike',
                      'p':'contaminated by latent image left by bright star',
                      'h':'halo of nearby bright source',
                      'o':'optical ghost'}
    wise_var_flag = {'n':'too few measurements to decide if variable',
                     '0':'most likely not variable (0, 0=certainly not-5=probably not)',
                     '1':'most likely not variable (1, 0=certainly not-5=probably not)',
                     '2':'most likely not variable (2, 0=certainly not-5=probably not)',
                     '3':'most likely not variable (3, 0=certainly not-5=probably not)',
                     '4':'most likely not variable (4, 0=certainly not-5=probably not)',
                     '5':'most likely not variable (5, 0=certainly not-5=probably not)',
                     '6':'likely variable (6, 6=likely-7=more likely)',
                     '7':'likely variable (7, 6=likely-7=more likely)',
                     '8':'most likely variable (8, 8=most likely-9=almost certain)',
                     '9':'most likely variable (9, 8=most likely-9=almost certain)'}
    twomass_qual_flag = {'X':'detection, but no valid brightness estimate (X)',
                         'U':'upper limit (U)',
                         'F':'error estimate not reliable (F)',
                         'E':'poor PSF fit (E)',
                         'A':'high quality (A)',
                         'B':'high quality (B)',
                         'C':'high quality (C)',
                         'D':'high quality (D)'}
    
    
    indices = np.arange(len(master))
    logger.info('Checking source catalogs for additional information')
    for source in sorted(sources):
        results,units,comms = search(source,ID=ID,**kwargs)
        if results is None: continue
        #-- the DENIS database
        if source=='B/denis/denis':
            for iflag,photband in zip(['Iflg','Jflg','Kflg'],['DENIS.I','DENIS.J','DENIS.KS']):
                index = indices[(master['source']==source) & (master['photband']==photband)]
                if not len(index):
                    continue
                flag = float(results[0][iflag])
                if np.isnan(flag):
                    messages[index] = '; '.join([messages[index],'high quality'])
                    continue
                flag = '{0:04d}'.format(int(flag))
                image_flag = flag[:2]
                source_flag = flag[2:]
                #-- keep track for output
                if image_flag in denis_image_flags:
                    messages[index] = '; '.join([messages[index],denis_image_flags[image_flag]])
                if source_flag in denis_source_flags:
                    messages[index] = '; '.join([messages[index],denis_source_flags[source_flag]]) 
        if source=='II/311/wise':
            conf = results[0]['ccf']
            var = results[0]['var']
            ex = int(results[0]['ex'])
            for i,photband in enumerate(['WISE.W1','WISE.W2','WISE.W3','WISE.W4']):
                index = indices[(master['source']==source) & (master['photband']==photband)]
                if len(index)!=1:
                    logger.warning("Skipping WISE flags, don't know what to do with {}".format(index))
                    continue
                if conf[i]!=' ' and conf[i] in wise_conf_flag:
                    messages[index] = '; '.join([messages[index],wise_conf_flag[conf[i].lower()]])
                if var[i]!=' ':
                    messages[index] = '; '.join([messages[index],wise_var_flag[var[i].lower()]])
                messages[index] = '; '.join([messages[index],(ex==0 and 'point source' or 'extended source')])
        if source=='II/246/out':
            flag = results[0]['Qflg'].strip()
            for i,photband in enumerate(['2MASS.J','2MASS.H','2MASS.KS']):
                if flag[i] in twomass_qual_flag:
                    index = indices[(master['source']==source) & (master['photband']==photband)]
                    if not len(index):
                        continue
                    messages[index] = '; '.join([messages[index],twomass_qual_flag[flag[i]]])
            
    
    #-- strip first '; ':
    for i in range(len(messages)):
        if messages[i]:
            messages[i] = messages[i][2:]
        else:
            messages[i] = '-'
        messages[i] = messages[i].replace(' ','_')
    
    #-- perhaps we want to return the master record array with an extra column
    messages = np.rec.fromarrays([messages],names=['comments'])
    if return_master and not 'comments' in master.dtype.names:
        return numpy_ext.recarr_join(master,messages)
    elif return_master:
        master['comments'] = messages
        return master
    else:
        return messages
                            
            
    




def tsv2recarray(filename):
    """
    Read a Vizier tsv (tab-sep) file into a record array.
    
    @param filename: name of the TSV file
    @type filename: str
    @return: catalog data columns, units, comments
    @rtype: record array, dict, list of str
    """
    data,comms = ascii.read2array(filename,dtype=np.str,splitchar='\t',return_comments=True)
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
        formats = np.zeros_like(data[0])
        for line in comms:                  
            line = line.split('\t')
            if len(line)<3: continue
            for i,key in enumerate(data[0]):
                if key == line[1] and line[0]=='Column': # this is the line with information
                    formats[i] = line[2].replace('(','').replace(')','').lower()
                    if formats[i][0].isdigit(): formats[i] = 'a100'
                    elif 'f' in formats[i]: formats[i] = 'f8' # floating point
                    elif 'i' in formats[i]: formats[i] = 'f8' # integer, but make it float to contain nans
                    elif 'e' in formats[i]: formats[i] = 'f8' # exponential
                    #-- see remark about the nans a few lines down
                    if formats[i][0]=='a':
                        formats[i] = 'a'+str(int(formats[i][1:])+3)
        #-- define dtypes for record array
        dtypes = np.dtype([(i,j) for i,j in zip(data[0],formats)])
        #-- remove spaces or empty values
        cols = []
        for i,key in enumerate(data[0]):
             col = data[3:,i]
             #-- fill empty values with nan, and make sure each string in the
             #   array is at least three spaces long (otherwise we cannot fit
             #   'nan' in the row and the casting will not work properly in a
             #   few lines
             cols.append([(row.isspace() or not row) and np.nan or 3*' '+row for row in col])
             #-- fix unit name
             #if source in cat_info.sections() and cat_info.has_option(source,data[1,i]):
             #   units[key] = cat_info.get(source,data[1,i])
             #else:  
             units[key] = data[1,i]
        #-- define columns for record array and construct record array
        cols = [np.cast[dtypes[i]](cols[i]) for i in range(len(cols))]
        results = np.rec.array(cols, dtype=dtypes)
    return results,units,comms

def vizier2phot(source,results,units,master=None,e_flag='e_',q_flag='q_',extra_fields=None,take_mean=False):
    """
    Convert/combine VizieR record arrays to measurement record arrays.
    
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
    
    You can add extra information from the VizieR catalog via the list of keys
    C{extra_fields}.
    
    If you give a C{master}, the information will be added to a previous
    record array. If not, a new master will be created.
    
    Colors will be expanded, derived from the other columns and added to the
    master.
    
    The result is a record array with each row a measurement.
    
    Example usage:
    
    First look for all photometry of Vega in all VizieR catalogs:
    
    >>> from ivs.sed import filters
    >>> import pylab
    >>> master = None
    >>> for source in cat_info.sections():
    ...     results,units,comms = search(source,ID='vega',radius=60.)
    ...     if results is not None:
    ...         master = vizier2phot(source,results,units,master,extra_fields=['_r','_RAJ2000','_DEJ2000'])
    
    Keep only observations we have an measurement and error of, convert every
    observation to 'Jy' and keep track of the results to plot.
    
    >>> master = master[(-np.isnan(master['e_meas'])) & (-np.isnan(master['meas']))]
    >>> eff_waves = filters.eff_wave(master['photband'])
    >>> myvalue,e_myvalue = conversions.nconvert(master['unit'],'erg/s/cm2/AA',master['meas'],master['e_meas'],photband=master['photband'])
    >>> for i in range(len(master)):
    ...    print '%15s %10.3e+/-%10.3e %11s %10.3e %3s %6.2f %6.2f %6.3f %23s'%(master[i]['photband'],master[i]['meas'],master[i]['e_meas'],master[i]['unit'],myvalue[i],'Jy',master[i]['_RAJ2000'],master[i]['_DEJ2000'],master[i]['_r'],master[i]['source'])
           GENEVA.V  6.100e-02+/- 2.500e-02         mag  3.620e-09  Jy 279.24  38.77 56.000             II/169/main
           GENEVA.B -8.980e-01+/- 2.500e-02         mag  6.518e-09  Jy 279.24  38.77 56.000             II/169/main
           GENEVA.U  6.070e-01+/- 2.500e-02         mag  3.223e-09  Jy 279.24  38.77 56.000             II/169/main
          GENEVA.V1  7.640e-01+/- 2.500e-02         mag  3.782e-09  Jy 279.24  38.77 56.000             II/169/main
          GENEVA.B1  2.000e-03+/- 2.500e-02         mag  6.584e-09  Jy 279.24  38.77 56.000             II/169/main
          GENEVA.B2  6.120e-01+/- 2.500e-02         mag  6.208e-09  Jy 279.24  38.77 56.000             II/169/main
           GENEVA.G  1.270e+00+/- 2.500e-02         mag  3.111e-09  Jy 279.24  38.77 56.000             II/169/main
            2MASS.J -1.770e-01+/- 2.060e-01         mag  3.591e-10  Jy 279.23  38.78  0.034              II/246/out
            2MASS.H -2.900e-02+/- 1.460e-01         mag  1.176e-10  Jy 279.23  38.78  0.034              II/246/out
           2MASS.KS  1.290e-01+/- 1.860e-01         mag  3.799e-11  Jy 279.23  38.78  0.034              II/246/out
           IRAS.F12  4.160e+01+/- 1.664e+00          Jy  1.024e-13  Jy 279.23  38.78  9.400             II/125/main
           IRAS.F25  1.100e+01+/- 5.500e-01          Jy  6.195e-15  Jy 279.23  38.78  9.400             II/125/main
           IRAS.F60  9.510e+00+/- 7.608e-01          Jy  8.420e-16  Jy 279.23  38.78  9.400             II/125/main
          IRAS.F100  7.760e+00+/- 6.984e-01          Jy  2.349e-16  Jy 279.23  38.78  9.400             II/125/main
           TD1.1565  5.689e-09+/- 1.700e-11  10mW/m2/nm  5.689e-09  Jy 279.23  38.78 18.510          II/59B/catalog
           TD1.1965  4.928e-09+/- 1.300e-11  10mW/m2/nm  4.928e-09  Jy 279.23  38.78 18.510          II/59B/catalog
           TD1.2365  3.700e-09+/- 1.000e-11  10mW/m2/nm  3.700e-09  Jy 279.23  38.78 18.510          II/59B/catalog
           TD1.2740  3.123e-09+/- 9.000e-12  10mW/m2/nm  3.123e-09  Jy 279.23  38.78 18.510          II/59B/catalog
            ANS.15N -4.910e-01+/- 1.000e-03         mag  5.707e-09  Jy 279.23  38.78 14.000               II/97/ans
            ANS.15W -4.410e-01+/- 1.200e-02         mag  5.450e-09  Jy 279.23  38.78 14.000               II/97/ans
             ANS.18 -4.620e-01+/- 3.000e-03         mag  5.556e-09  Jy 279.23  38.78 14.000               II/97/ans
             ANS.25  4.600e-02+/- 4.000e-03         mag  3.480e-09  Jy 279.23  38.78 14.000               II/97/ans
             ANS.33  1.910e-01+/- 3.000e-03         mag  3.045e-09  Jy 279.23  38.78 14.000               II/97/ans
       HIPPARCOS.HP  8.680e-02+/- 2.100e-03         mag  3.840e-09  Jy 279.23  38.78  3.060          I/239/hip_main
            MIPS.24  8.900e+03+/- 8.900e+01         mJy  4.628e-15  Jy 279.23  38.78  0.010    J/ApJ/653/675/table1
            MIPS.70  1.142e+04+/- 2.283e+03         mJy  7.075e-16  Jy 279.23  38.78  0.010    J/ApJ/653/675/table1
          JOHNSON.B  1.900e-02+/- 1.000e-02         mag  6.216e-09  Jy 279.23  38.78  3.060             I/280B/ascc
          JOHNSON.V  7.400e-02+/- 2.000e-03         mag  3.428e-09  Jy 279.23  38.78  3.060             I/280B/ascc
          JOHNSON.V  3.300e-02+/- 1.200e-02         mag  3.560e-09  Jy 279.23  38.78  0.010         II/168/ubvmeans
        JOHNSON.B-V -1.000e-03+/- 5.000e-03         mag        nan  Jy 279.23  38.78  0.010         II/168/ubvmeans
        JOHNSON.U-B -6.000e-03+/- 6.000e-03         mag        nan  Jy 279.23  38.78  0.010         II/168/ubvmeans
          JOHNSON.B  3.200e-02+/- 1.300e-02         mag  6.142e-09  Jy 279.23  38.78  0.010         II/168/ubvmeans
          JOHNSON.U  2.600e-02+/- 1.432e-02         mag  4.086e-09  Jy 279.23  38.78  0.010         II/168/ubvmeans
          JOHNSON.K  1.300e-01+/- 1.900e-01         mag  3.764e-11  Jy 279.23  38.78  0.010 J/PASP/120/1128/catalog
          AKARI.N60  6.582e+00+/- 2.090e-01          Jy  4.614e-16  Jy 279.23  38.78  3.400              II/298/fis
        AKARI.WIDES  6.201e+00+/- 1.650e-01          Jy  2.566e-16  Jy 279.23  38.78  3.400              II/298/fis
        AKARI.WIDEL  4.047e+00+/- 3.500e-01          Jy  5.658e-17  Jy 279.23  38.78  3.400              II/298/fis
         AKARI.N160  3.221e+00+/- 2.550e-01          Jy  3.695e-17  Jy 279.23  38.78  3.400              II/298/fis
          AKARI.S9W  5.670e+01+/- 4.010e-01          Jy  2.169e-13  Jy 279.24  38.78  2.550              II/297/irc
         AKARI.L18W  1.254e+01+/- 7.770e-02          Jy  1.050e-14  Jy 279.24  38.78  2.550              II/297/irc
      STROMGREN.B-Y  3.000e-03+/- 3.000e-03         mag        nan  Jy 279.23  38.78 22.000          II/215/catalog
       STROMGREN.M1  1.570e-01+/- 3.000e-03         mag        nan  Jy 279.23  38.78 22.000          II/215/catalog
       STROMGREN.C1  1.088e+00+/- 4.000e-03         mag        nan  Jy 279.23  38.78 22.000          II/215/catalog
        STROMGREN.B  4.300e-02+/- 3.000e-03         mag  5.604e-09  Jy 279.23  38.78 22.000          II/215/catalog
         DIRBE.F2_2  6.217e+02+/- 9.500e+00          Jy  3.791e-11  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
         DIRBE.F3_5  2.704e+02+/- 1.400e+01          Jy  6.534e-12  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
         DIRBE.F4_9  1.504e+02+/- 6.200e+00          Jy  1.895e-12  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
          DIRBE.F12  2.910e+01+/- 1.570e+01          Jy  6.014e-14  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
          DIRBE.F25  1.630e+01+/- 3.120e+01          Jy  1.153e-14  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
          DIRBE.F60  1.220e+01+/- 5.610e+01          Jy  1.195e-15  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
         DIRBE.F100 -7.000e-01+/- 7.790e+01          Jy -2.238e-17  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
         DIRBE.F140  2.557e+02+/- 5.223e+03          Jy  3.577e-15  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
         DIRBE.F240  8.290e+01+/- 2.881e+03          Jy  4.152e-16  Jy 279.23  38.78  0.110    J/ApJS/154/673/DIRBE
            WISE.W2  1.143e+00+/- 1.900e-02         mag  8.428e-13  Jy 279.24  38.78  3.276             II/311/wise
            WISE.W3 -6.700e-02+/- 8.000e-03         mag  6.930e-14  Jy 279.24  38.78  3.276             II/311/wise
            WISE.W4 -1.270e-01+/- 6.000e-03         mag  5.722e-15  Jy 279.24  38.78  3.276             II/311/wise
                
    Make a quick plot:
    
    >>> p = pylab.figure()
    >>> p = pylab.loglog(eff_waves,myvalue,'ko')
    >>> p = pylab.show()
    
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
    dtypes = [('meas','f8'),('e_meas','f8'),('flag','a20'),
                  ('unit','a30'),('photband','a30'),('source','a50')]
    
    #-- extra can be added:
    names = list(results.dtype.names)

    if extra_fields is not None:
        for e_dtype in extra_fields:
            if e_dtype in names:
                dtypes.append((e_dtype,results.dtype[names.index(e_dtype)].str))
            else:
                dtypes.append((e_dtype,'f8'))
    
    #-- create empty master if not given
    newmaster = False
    if master is None or len(master)==0:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
        newmaster = True
    
    #-- add fluxes and magnitudes to the record array
    cols_added = 0
    for key in cat_info.options(source):
        if key in ['e_flag','q_flag','mag','bibcode']:
            continue
        photband = cat_info.get(source,key)
        key = key.replace('_[','[').replace(']_',']')
        #-- contains measurement, error, quality, units, photometric band and
        #   source
        if not take_mean or len(results)<=1: #take first
            cols = [results[key][:1],
                (e_flag+key in results.dtype.names and results[e_flag+key][:1] or np.ones(len(results[:1]))*np.nan),
                (q_flag+key in results.dtype.names and results[q_flag+key][:1] or np.ones(len(results[:1]))*np.nan),
                np.array(len(results[:1])*[units[key]],str),
                np.array(len(results[:1])*[photband],str),
                np.array(len(results[:1])*[source],str)]
        else:
            logger.warning('taking mean: {0} --> {1}+/-{2}'.format(results[key][0],results[key].mean(),results[key].std()))
            cols = [[results[key].mean()],[results[key].std()],
                (q_flag+key in results.dtype.names and results[q_flag+key][:1] or np.ones(len(results[:1]))*np.nan),
                np.array(len(results[:1])*[units[key]],str),
                np.array(len(results[:1])*[photband],str),
                np.array(len(results[:1])*[source],str)]
        #-- correct errors given in percentage
        if e_flag+key in results.dtype.names and units[e_flag+key]=='%':
            cols[1] = cols[1]/100.*cols[0]
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
    master_ = _breakup_colours(master[N:])
    master_ = _breakup_colours(master_)
    #-- combine and return
    master = np.core.records.fromrecords(master.tolist()[:N]+master_.tolist(),dtype=dtypes)
    
    #-- skip first line from building 
    if newmaster: master = master[1:]
    return master
    


def vizier2fund(source,results,units,master=None,e_flag='e_',q_flag='q_',extra_fields=None):
    """
    Convert/combine VizieR record arrays to measurement record arrays.
    
    This is probably only useful if C{results} contains only information on
    one target (or you have to give 'ID' as an extra field, maybe).
    
    The standard columns are:
    
        1. C{meas}: containing the measurement of a fundamental parameter
        2. C{e_meas}: the error on the measurement of a fundamental parameter
        3. C{flag}: an optional quality flag
        4. C{unit}: the unit of the measurement
        5. C{source}: name of the source catalog
    
    If a target appears more than once in a catalog, only the first match will
    be added.
    
    The result is a record array with each row a measurement.
    
    Example usage:
    
    >>> master = None
    >>> for source in cat_info_fund.sections():
    ...     results,units,comms = search(source,ID='AzV 79',radius=60.)
    ...     if results is not None:
    ...         master = vizier2fund(source,results,units,master,extra_fields=['_r','_RAJ2000','_DEJ2000'])
    >>> master = master[-np.isnan(master['meas'])]
    >>> for i in range(len(master)):
    ...     print '%10s %10.3e+/-%10.3e %11s %6.2f %6.2f %6.3f %23s'%(master[i]['name'],master[i]['meas'],master[i]['e_meas'],master[i]['unit'],master[i]['_RAJ2000'],master[i]['_DEJ2000'],master[i]['_r'],master[i]['source'])
          Teff  7.304e+03+/-       nan           K  12.67 -72.83  0.002         B/pastel/pastel
          logg  2.000e+00+/-       nan     [cm/s2]  12.67 -72.83  0.002         B/pastel/pastel
        [Fe/H] -8.700e-01+/-       nan       [Sun]  12.67 -72.83  0.002         B/pastel/pastel
    
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
    if cat_info_fund.has_option(source,'q_flag'):
        q_flag = cat_info_fund.get(source,'q_flag')
    if cat_info_fund.has_option(source,'e_flag'):
        e_flag = cat_info_fund.get(source,'e_flag')
    
    #-- basic dtypes
    dtypes = [('meas','f8'),('e_meas','f8'),('q_meas','f8'),('unit','a30'),
              ('source','a50'),('name','a50')]
    
    #-- extra can be added:
    names = list(results.dtype.names)
    if extra_fields is not None:
        for e_dtype in extra_fields:
            dtypes.append((e_dtype,results.dtype[names.index(e_dtype)].str))
    
    #-- create empty master if not given
    newmaster = False
    if master is None or len(master)==0:
        master = np.rec.array([tuple([('f' in dt[1]) and np.nan or 'none' for dt in dtypes])],dtype=dtypes)
        newmaster = True
    
    #-- add fluxes and magnitudes to the record array
    for key in cat_info_fund.options(source):
        if key in ['e_flag','q_flag']:
            continue
        photband = cat_info_fund.get(source,key)
        #-- contains measurement, error, quality, units, photometric band and
        #   source
        #if len(results[e_flag+key])>1:
        key = key.replace('_[','[').replace(']_',']')
        cols = [results[key][:1],
                (e_flag+key in results.dtype.names and results[e_flag+key][:1] or np.ones(len(results[:1]))*np.nan),
                (q_flag+key in results.dtype.names and results[q_flag+key][:1] or np.ones(len(results[:1]))*np.nan),
                np.array(len(results[:1])*[units[key]],str),
                np.array(len(results[:1])*[source],str),
                np.array(len(results[:1])*[key],str)]
        #-- add any extra fields if desired.
        if extra_fields is not None:
            for e_dtype in extra_fields:
                cols.append(results[:1][e_dtype])
        #-- add to the master
        rows = []
        for i in range(len(cols[0])):
            rows.append(tuple([col[i] for col in cols]))
        #print master
        master = np.core.records.fromrecords(master.tolist()+rows,dtype=dtypes)
    
    #-- skip first line from building 
    if newmaster: master = master[1:]
    return master


def catalog2bibcode(catalog):
    """
    Retrieve the ADS bibcode of a ViZieR catalog.
    
    @param catalog: name of the catalog (e.g. II/306/sdss8)
    @type catalog: str
    @return: bibtex code
    @rtype: str
    """
    catalog = catalog.split("/")
    N = max(2,len(catalog)-1)
    catalog = "/".join(catalog[:N])
    base_url = "http://cdsarc.u-strasbg.fr/viz-bin/Cat?{0}".format(catalog)
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url)

    code = None
    ff = open(filen,'r')
    for line in ff.readlines():
        if 'Keywords' in line: break
        if '=<A HREF' in line:
            code = line.split('>')[1].split('<')[0]
    ff.close()
    url.close()
    return code

def bibcode2bibtex(bibcode):
    """
    Retrieve the bibtex entry of an ADS bibcode.
    
    @param bibcode: bibcode (e.g. C{2011yCat.2306....0A})
    @type bibcode: str
    @return: bibtex entry
    @rtype: str
    """
    base_url = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode={0}&data_type=BIBTEX&db_key=AST&nocookieset=1".format(bibcode)
    url = urllib.URLopener()
    filen,msg = url.retrieve(base_url)

    bibtex = []
    ff = open(filen,'r')
    for line in ff.readlines():
        if (not bibtex and '@' in line) or bibtex:
            bibtex.append(line)
    ff.close()
    url.close()
    return "".join(bibtex)

def catalog2bibtex(catalog):
    """
    Retrieve the bibtex entry of a catalog.
    
    @param catalog: name of the catalog (e.g. II/306/sdss8)
    @type catalog: str
    @return: bibtex entry
    @rtype: str
    """
    bibcode = catalog2bibcode(catalog)
    bibtex = bibcode2bibtex(bibcode)
    return bibtex
    

#}

#{ Internal helper functions

def _get_URI(name=None,ID=None,ra=None,dec=None,radius=20.,
                     oc='deg',oc_eq='J2000',
                     out_all=True,out_max=1000000,
                     filetype='tsv',sort='_r',constraints=None,**kwargs):
    """
    Build Vizier URI from available options.
    
    kwargs are to catch unused arguments.
    
    @param name: name of a ViZieR catalog (e.g. 'II/246/out')
    @type name: str
    @param filetype: type of the retrieved file 
    @type filetype: str (one of 'tsv','csv','ascii'... see ViZieR site)
    @param oc: coordinates
    @type oc: str (one of 'deg'...)
    @param out_all: retrieve all or basic information
    @type out_all: boolean (True = all, None = basic)
    @param out_max: maximum number of rows
    @type out_max: integer
    @param ID: target name
    @type ID: str
    @param ra: target's right ascension
    @type ra: float
    @param dec: target's declination
    @type dec: float
    @param radius: search radius (arcseconds)
    @type radius: float
    @return: url
    @rtype: str
    """
    base_url = 'http://{0}/viz-bin/asu-{1}/VizieR?&-oc={2},eq={3}'.format(mirrors['current'],filetype,oc,oc_eq)
    if sort:
        base_url += '&-sort=%s'%(sort)
    if constraints is not None:
        for constr in constraints:
            base_url += '&%s'%(constr)
    
    if ID is not None:
        #-- if the ID is given in the form 'J??????+??????', derive the
        #   coordinates of the target from the name.
        if ID[0]=='J':
            ra = int(ID[1:3]),int(ID[3:5]),float(ID[5:10])
            dec = int(ID[10:13]),int(ID[13:15]),float(ID[15:])
            ra = '%02d+%02d+%.2f'%ra
            dec = '+%+02d+%02d+%.2f'%dec
            ID = None
    
    
    if name:    base_url += '&-source=%s'%(name)
    if out_all: base_url += '&-out.all'
    if out_max: base_url += '&-out.max=%s'%(out_max)
    if radius:  base_url += '&-c.rs=%s'%(radius)
    if ID is not None and ra is None: base_url += '&-c=%s'%(urllib.quote(ID))
    if ra is not None: base_url += '&-c.ra=%s&-c.dec=%s'%(ra,dec)
    logger.debug(base_url)
    #print base_url
    return base_url

def _breakup_colours(master):
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
        
        ########################################################################
        #-- NORMAL COLORS
        ########################################################################
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
        
        ########################################################################
        #-- STROMGREN COLORS
        ########################################################################
        #-- stromgren index M1
        elif color.upper()=='M1':
            # m1 = v - 2*b + y
            #-- in which column are the measurements (and error) located?
            index_meas, index_emeas = names.index('meas'),names.index('e_meas')
            index_band = names.index('photband')
            #-- this is the m1 row
            row = list(master[i])
            meas,e_meas = row[index_meas],row[index_emeas]
            #-- retrieve the measurements we need to calculate band-magnitudes
            my_photbands = list(master['photband'])
            if 'STROMGREN.Y' in my_photbands and 'STROMGREN.B' in my_photbands and (not 'STROMGREN.V' in my_photbands):
                index_b = my_photbands.index('STROMGREN.B')
                index_y = my_photbands.index('STROMGREN.Y')
                rowb,rowy = list(master[index_b]),list(master[index_y])
                b,e_b = rowb[index_meas],rowb[index_emeas]
                y,e_y = rowy[index_meas],rowy[index_emeas]
                #-- add extra row
                row1 = list(master[index_band])
                row1[index_band] = 'STROMGREN.V'
                row1[index_meas] = meas + 2*b - y
                row1[index_emeas] = np.sqrt(e_meas**2 + 2*e_b**2 + e_y**2)
                master = np.core.records.fromrecords(master.tolist()+[tuple(row1)],dtype=master.dtype)
                logger.debug("Added band STROMGREN.Y (m1)")
        #-- stromgren index C1
        elif color.upper()=='C1':
            # c1 = u - 2*v + b
            #-- in which column are the measurements (and error) located?
            index_meas, index_emeas = names.index('meas'),names.index('e_meas')
            index_band = names.index('photband')
            #-- this is the m1 row
            row = list(master[i])
            meas,e_meas = row[index_meas],row[index_emeas]
            #-- retrieve the measurements we need to calculate band-magnitudes
            my_photbands = list(master['photband'])
            if 'STROMGREN.V' in my_photbands and 'STROMGREN.B' in my_photbands and (not 'STROMGREN.U' in my_photbands):
                index_v = my_photbands.index('STROMGREN.V')
                index_b = my_photbands.index('STROMGREN.B')
                rowb,rowv = list(master[index_b]),list(master[index_v])
                b,e_b = rowb[index_meas],rowb[index_emeas]
                v,e_v = rowv[index_meas],rowv[index_emeas]
                #-- add extra row
                row1 = list(master[index_band])
                row1[index_band] = 'STROMGREN.U'
                row1[index_meas] = meas + 2*v - b
                row1[index_emeas] = np.sqrt(e_meas**2 + 2*e_v**2 + e_b**2)
                master = np.core.records.fromrecords(master.tolist()+[tuple(row1)],dtype=master.dtype)
                logger.debug("Added band STROMGREN.U (c1)")
        
        
        
        
    return master
    

def test():
    """
    Execute all docstrings.
    
    >>> import pylab
    >>> p = pylab.show()
    """
    import doctest
    doctest.testmod()

#}

if __name__=="__main__":
    test()
