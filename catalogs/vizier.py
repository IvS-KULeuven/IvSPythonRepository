"""
Interface to the ViZieR website.

Download or retrieve ViZieR catalogs.
"""

import numpy as np
import urllib
import logging
from divers import io

logger = logging.getLogger("VIZIER")

cat_info = {'II/169/main':{'Vmag':'GENEVA.V',     # Rufener Geneva photometric
                           'U-B' :'GENEVA.U-B',
                           'V-B' :'GENEVA.V-B',
                           'V1-B':'GENEVA.V1-B',
                           'B1-B':'GENEVA.B1-B',
                           'B2-B':'GENEVA.B2-B',
                           'G-B' :'GENEVA.G-B'},
            'II/246/out':{'Jmag':'2MASS.J',       # 2MASS point source
                          'Hmag':'2MASS.H',
                          'Kmag':'2MASS.KS'},
            'II/125/main':{'Fnu_12':'IRAS.F12',   # IRAS point source
                           'Fnu_25':'IRAS.F25',
                           'Fnu_60':'IRAS.F60',
                           'Fnu_100':'IRAS.F100'},
            'II/59B/catalog':{'F1565':'TD1.1565', # TD1 catalog
                              'F1965':'TD1.1965',
                              'F2365':'TD1.2365',
                              'F2740':'TD1.2740'},
            'II/97/ans':{'15N':'ANS.15N',         # ANS satellite
                         '15W':'ANS.15W',
                         '18' :'ANS.18',
                         '22' :'ANS.22',
                         '25' :'ANS.25',
                         '33' :'ANS.33',
                         'q_flag':'l_'},
            'I/259/tyc2':{'BTmag':'TYCHO2.BT',   # TYCHO2
                          'VTmag':'TYCHO2.VT'},
            'I/284/out':{'B1mag':'USNOB1.B1',    # USNOB1
                         'R1mag':'USNOB1.R1',
                         'B2mag':'USNOB1.B2',
                         'R2mag':'USNOB1.R2'},
            'II/293/glimpse':{'F(3.6)':'IRAC.36', # Spitzer space telescope
                              'F(4.5)':'IRAC.45',
                              'F(5.8)':'IRAC.58',
                              'F(8.0)':'IRAC.80'},
            'J/PASP/120/1128/catalog':{'Jmag':'JOHNSON.J',
                                       'Hmag':'JOHNSON.H',
                                       'Kmag':'JOHNSON.K'}
             }


def search(source,**kwargs):
    """
    Search and retrieve information from a ViZieR catalog.
    
    If you're looking for information on B{one target}, then give the target's
    C{ID} or coordinates (C{ra} and C{dec}), and a search C{radius}.
    
    If you're looking for information of B{a whole field}, then give the
    field's coordinates (C{ra} and C{dec}), and C{radius}.
    
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
        >>> print "Vega: Vmag = %.3f %s, RA = %.2f %s, DEC = %.2f %s"\
                     %(results['Vmag'],units['Vmag'],\
                       results['_RAJ2000'],units['_RAJ2000'],\
                       results['_DEJ2000'],units['_DEJ2000'])
        Vega: Vmag = 0.061 mag, RA = 279.24 deg, DEC = 38.77 deg
    
        2. Search for all targets in the 2MASS catalog in a particular field.
        Download the results to a FITS file, read the file, and plot the results
        to the screen.
    
        >>> filename = search('II/246/out',ra=100.79,dec=0.70,radius=1000.,\
                filetype='fits',filename='2mass_test',out_all=None)
    
        Now read in the FITS-file and plot the contents
        >>> import pyfits,pylab,os
        >>> ff = pyfits.open('2mass_test.fits')
        >>> p = pylab.gcf().canvas.set_window_title('test of <search>')
        >>> p = pylab.scatter(ff[1].data.field('_RAJ2000'),\
                             ff[1].data.field('_DEJ2000'),\
                             c=ff[1].data.field('Jmag'),\
                             s=(20-ff[1].data.field('Jmag'))**2,\
                             cmap=pylab.cm.hot_r,edgecolors='none')
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
    results = {}
    units = {}
    if filetype=='tsv':                
        data,comms = io.read(filen,return_string=True,skip_blanks=True,splitchar='\t',return_comments=True)
            
        if len(data)>0:
            #-- remove return characters at the end of each row
            for col in range(len(data)):
                data[col][-1] = data[col][-1].strip()
            #-- now convert this thing into a nice dictionary
            data = np.array(data)
            for i,key in enumerate(data[0]):
                 col = data[3:,i]
                 #-- fill empty values with nan
                 col = np.array([np.nan and (row.isspace() or not row) or row for row in col])
                 #-- try to make the column a floating array
                 try:
                     col = np.array(col,float)
                 except ValueError:
                     pass
                 #-- put in the dictionary
                 results[key] = col
                 units[key] = data[1,i]
    url.close()
    logger.info('Results converted to dictionary')
    return results,units,comms
    

def vizier2dict(source,results,units,e_flag='e_',q_flag='q_'):
    """
    """
    q_flag = cat_info[source].get('q_flag',q_flag)
    e_flag = cat_info[source].get('e_flag',e_flag)
    
    data = {}
    for key in cat_info[source]:
        if key in ['e_flag','q_flag']:
            continue
        
        system,band = cat_info[source][key].split('.')
        if not system in data:
            data[system] = {}
        #-- contains measurement, error, quality and units
        data[system][key] = [results[key],
                            (e_flag+key in results and results[e_flag+key] or nan),
                            (q_flag+key in results and results[q_flag+key] or nan),
                            units[key]]
    return data
    

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
    
    
    #from divers import mylogging
    #logger = mylogging.get_basic_logger("")
    #for source in cat_info:
        #results,units,comms = search(source,ID='vega',radius=60.)
        #if results: print vizier2dict(source,results,units)

        #results,units,comms = search(source,ID='HD180642',radius=5.)
        #if results: print vizier2dict(source,results,units)


#filename = search_vizier('II/246/out',ra=100.79441,dec=0.7033,radius=1.347395*60*60,filetype='fits',filename='2mass_test')
