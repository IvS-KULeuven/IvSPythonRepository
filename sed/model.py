# -*- coding: utf-8 -*-
"""
Interface to the SED library.

The most basic usage of this module is:

>>> wave,flux = get_table(teff=10000,logg=4.0)

This will retrieve the model SED with the specified B{effective temperature} and
B{logg}, from the standard B{grid}, in standard B{units} and with zero
B{reddening}. All these things can be specified though (see below).

Section 1. Available model grids
================================

    Section 1.1 Available grids
    ---------------------------

    - kurucz: The Kurucz model grids, (default setting)  reference: Kurucz 1993, yCat, 6039, 0
        - metallicity (z): m01 is -0.1 log metal abundance relative to solar (solar abundances from Anders and Grevesse 1989)
        - metallicity (z): p01 is +0.1 log metal abundance relative to solar (solar abundances from Anders and Grevesse 1989)
        - alpha enhancement (alpha): True means alpha enhanced (+0.4)
        - turbulent velocity (vturb): vturb in km/s
        - nover= True means no overshoot
        - odfnew=True means no overshoot but with better opacities and abundances

    - tmap: NLTE grids computed for sdB stars with the Tubingen NLTE Model
      Atmosphere package. No further parameters are available. Reference:
      Werner et al. 2003,


    Section 1.2 Plotting the domains of all spectral grids
    ------------------------------------------------------

    We make a plot of the domains of all spectral grids. Therefore, we first collect
    all the grid names

    >>> grids = get_gridnames()

    Preparation of the plot: set the color cycle of the current axes to the spectral
    color cycle.

    >>> p = pl.figure(figsize=(10,8))
    >>> color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(grids))]
    >>> p = pl.gca().set_color_cycle(color_cycle)

    To plot all the grid points, we run over all grid names (which are strings), and
    retrieve their dimensions. The dimensions are just two arrays giving the teff-
    and logg-coordinate of each SED in the grid. They can thus be easily plot:

    >>> for grid in grids:
    ...    teffs,loggs = get_grid_dimensions(grid=grid)
    ...    p = pl.plot(np.log10(teffs),loggs,'o',ms=7,label=grid)

    And we need to set some of the plotting details to make it look nicer.

    >>> p = pl.xlim(pl.xlim()[::-1])
    >>> p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Effective temperature [K]')
    >>> p = pl.ylabel('log( Surface gravity [cm s$^{-1}$]) [dex]')
    >>> xticks = [3000,5000,7000,10000,15000,25000,35000,50000,65000]
    >>> p = pl.xticks([np.log10(i) for i in xticks],['%d'%(i) for i in xticks])
    >>> p = pl.legend(loc='upper left',prop=dict(size='small'))
    >>> p = pl.grid()

    ]include figure]]ivs_sed_model_grid.png]

Section 2. Retrieval of model SEDs
==================================

Subsection 2.1 Default settings
-------------------------------

To get information on the grid that is currently defined, you can type the
following. Note that not all parameters are relevant for all grids, e.g. the
convection theory parameter C{ct} has no influence when the Kurucz grid is
chosen.

>>> print defaults
{'use_scratch': False, 'Rv': 3.1, 'co': 1.05, 'c': 0.5, 'grid': 'kurucz', 'alpha': False, 'odfnew': True, 'ct': 'mlt', 'a': 0.0, 'vturb': 2, 'law': 'fitzpatrick2004', 'm': 1.0, 't': 1.0, 'z': 0.0, 'nover': False, 'He': 97}

or

>>> print os.path.basename(get_file())
kurucz93_z0.0_k2odfnew_sed.fits

You can change the defaults with the function L{set_defaults}:

>>> set_defaults(z=0.5)
>>> print defaults
{'use_scratch': False, 'Rv': 3.1, 'co': 1.05, 'c': 0.5, 'grid': 'kurucz', 'alpha': False, 'odfnew': True, 'ct': 'mlt', 'a': 0.0, 'vturb': 2, 'law': 'fitzpatrick2004', 'm': 1.0, 't': 1.0, 'z': 0.5, 'nover': False, 'He': 97}

And reset the 'default' default values by calling L{set_defaults} without arguments

>>> set_defaults()
>>> print defaults
{'use_scratch': False, 'Rv': 3.1, 'co': 1.05, 'c': 0.5, 'grid': 'kurucz', 'alpha': False, 'odfnew': True, 'ct': 'mlt', 'a': 0.0, 'vturb': 2, 'law': 'fitzpatrick2004', 'm': 1.0, 't': 1.0, 'z': 0.0, 'nover': False, 'He': 97}

Subsection 2.2 Speeding up
--------------------------

When fitting an sed using the builder class, or repeatedly reading model seds,
or integrated photometry, the main bottleneck on the speed will be the disk access
This can be circumvented by using the scratch disk. To do this, call the function
copy2scratch() after setting the default settings as explained above. f.x.:

>>> set_defaults(grid='kurucz', z=0.5)
>>> copy2scratch()

You have to do this every time you change a grid setting. This function creates a
directory named 'your_username' on the scratch disk and works from there. So you
won`t disturbed other users.

After the fitting process use the function

>>> clean_scratch()

to remove the models that you used from the scratch disk. Be carefull with this,
because it will remove the models without checking if there is another process
using them. So if you have multiple scripts running that are using the same models,
only clean the scratch disk after the last process is finished.

The gain in speed can be up to 70% in single sed fitting, and up to 40% in binary
and multiple sed fitting.

For the sake of the examples, we'll set the defaults back to z=0.0:

>>> set_defaults()

Subsection 2.3 Model SEDs
-------------------------

Be careful when you supply parameters: e.g., not all grids are calculated for
the same range of metallicities. In L{get_table}, only the effective temperature
and logg are 'interpolatable' quantities. You have to set the metallicity to a
grid value. The reddening value can take any value: it is not interpolated but
calculated. You can thus also specify the type of reddening law (see L{reddening}).

>>> wave,flux = get_table(teff=12345,logg=4.321,ebv=0.12345,z=0.5)

but

>>> try:
...     wave,flux = get_table(teff=12345,logg=4.321,ebv=0.12345,z=0.6)
... except IOError,msg:
...     print msg
File sedtables/modelgrids/kurucz93_z0.6_k2odfnew_sed.fits not found in any of the specified data directories /STER/pieterd/IVSDATA/, /STER/kristofs/IVSdata

Since the Kurucz model atmospheres have not been calculated for the value of
C{z=0.6}.

Instead of changing the defaults of this module with L{set_defaults}, you can
also give extra arguments to L{get_table} to specify the grid you want to use.
The default settings will not change in this case.

>>> wave,flux = get_table(teff=16321,logg=4.321,ebv=0.12345,z=0.3,grid='tlusty')

The default B{units} of the SEDs are angstrom and erg/s/cm2/AA/sr. To change them,
do:

>>> wave,flux = get_table(teff=16321,logg=4.321,wave_units='micron',flux_units='Jy/sr')

To B{remove the steradian} from the units when you know the angular diameter of
your star in milliarcseconds, you can do (we have to convert diameter to surface):

>>> ang_diam = 3.21 # mas
>>> scale = conversions.convert('mas','sr',ang_diam/2.)
>>> wave,flux = get_table(teff=9602,logg=4.1,ebv=0.0,z=0.0,grid='kurucz')
>>> flux *= scale

The example above is representative for the case of Vega. So, if we now calculate
the B{synthetic flux} in the GENEVA.V band, we should end up with the zeropoint
magnitude of this band, which is close to zero:

>>> flam = synthetic_flux(wave,flux,photbands=['GENEVA.V'])
>>> print '%.3f'%(conversions.convert('erg/s/cm2/AA','mag',flam,photband='GENEVA.V')[0])
0.063

Compare this with the calibrated value

>>> print filters.get_info(['GENEVA.V'])['vegamag'][0]
0.061

Section 3. Retrieval of integrated photometry
=============================================

Instead of retrieving a model SED, you can immediately retrieve pre-calculated
integrated photometry. The benefit of this approach is that it is B{much} faster
than retrieving the model SED and then calculating the synthetic flux. Also,
you can supply arbitrary metallicities within the grid boundaries, as interpolation
is done in effective temperature, surface gravity, reddening B{and} metallicity.

Note that also the B{reddening law is fixed} now, you need to recalculate the
tables for different parameters if you need them.

The B{massive speed-up} is accomplished the following way: it may take a few tens
of seconds to retrieve the first pre-integrated SED, because all available
files from the specified grid will be loaded into memory, and a `markerarray'
will be made allowing a binary search in the grid. This makes it easy to retrieve
all models around the speficied point in N-dimensional space. Next, a linear
interpolation method is applied to predict the photometric values of the
specified point.

All defaults set for the retrieval of model SEDs are applicable for the integrated
photometry tables as well.

When retrieving integrated photometry, you also get the B{absolute luminosity}
(integration of total SED) as a bonus. This is the absolute luminosity assuming
the star has a radius of 1Rsol. Multiply by Rstar**2 to get the true luminosity.

Because photometric filters cannot trivially be assigned a B{wavelength} to (see
L{filters.eff_wave}), by default, no wavelength information is retrieved. If you
want to retrieve the effective wavelengths of the filters themselves (not taking
into account the model atmospheres), you can give an extra keyword argument
C{wave_units}. If you want to take into account the model atmosphere, use
L{filters.eff_wave}.

>>> photbands = ['GENEVA.U','2MASS.J']
>>> fluxes,Labs = get_itable(teff=16321,logg=4.321,ebv=0.12345,z=0.123,photbands=photbands)
>>> waves,fluxes,Labs = get_itable(teff=16321,logg=4.321,ebv=0.12345,z=0.123,photbands=photbands,wave_units='AA')

Note that the integration only gives you fluxes, and is thus independent from
the zeropoints of the filters (but dependent on the transmission curves). To
get the synthetic magnitudes, you can do

>>> mymags = [conversions.convert('erg/s/cm2/AA','mag',fluxes[i],photband=photbands[i]) for i in range(len(photbands))]

The mags don't mean anything in this case because they have not been corrected
for the distance to the star.

The retrieval of integrated photometry can go much faster if you want to do
it for a whole set of parameters. The L{get_itable_pix} function has a much
more flexible, reliable and fast interpolation scheme. It is possible to
interpolate also over doppler shift and interstellar Rv, as long as the grids
have been computed before. See L{get_itable_pix} for more information.

Subsection 3. Full example
==========================

We build an SED of Vega and compute synthetic magnitudes in the GENEVA and
2MASS bands.

These are the relevant parameters of Vega and photometric passbands

>>> ang_diam = 3.21 # mas
>>> teff = 9602
>>> logg = 4.1
>>> ebv = 0.0
>>> z = 0.0
>>> photbands = ['GENEVA.U','GENEVA.G','2MASS.J','2MASS.H','2MASS.KS']

We can compute (R/d) to scale the synthetic flux as

>>> scale = conversions.convert('mas','sr',ang_diam/2.)

We retrieve the SED

>>> wave,flux = get_table(teff=teff,logg=logg,ebv=ebv,z=z,grid='kurucz')
>>> flux *= scale

Then compute the synthetic fluxes, and compare them with the synthetic fluxes as
retrieved from the pre-calculated tables

>>> fluxes_calc = synthetic_flux(wave,flux,photbands)
>>> wave_int,fluxes_int,Labs = get_itable(teff=teff,logg=logg,ebv=ebv,z=z,photbands=photbands,wave_units='AA')
>>> fluxes_int *= scale

Convert to magnitudes:

>>> m1 = [conversions.convert('erg/s/cm2/AA','mag',fluxes_calc[i],photband=photbands[i]) for i in range(len(photbands))]
>>> m2 = [conversions.convert('erg/s/cm2/AA','mag',fluxes_int[i],photband=photbands[i]) for i in range(len(photbands))]

And make a nice plot

>>> p = pl.figure()
>>> p = pl.loglog(wave,flux,'k-',label='Kurucz model')
>>> p = pl.plot(wave_int,fluxes_calc,'ro',label='Calculated')
>>> p = pl.plot(wave_int,fluxes_int,'bx',ms=10,mew=2,label='Pre-calculated')
>>> p = [pl.annotate('%s: %.3f'%(b,m),(w,f),color='r') for b,m,w,f in zip(photbands,m1,wave_int,fluxes_calc)]
>>> p = [pl.annotate('%s: %.3f'%(b,m),(w-1000,0.8*f),color='b') for b,m,w,f in zip(photbands,m2,wave_int,fluxes_int)]
>>> p = pl.xlabel('Wavelength [Angstrom]')
>>> p = pl.ylabel('Flux [erg/s/cm2/AA]')

]include figure]]ivs_sed_model_example.png]

"""
import re
import os
import glob
import logging
import astropy.io.fits as pf
import numpy as np
try:
    from scipy.interpolate import LinearNDInterpolator
    from scipy.interpolate import griddata
    new_scipy = True
except ImportError:
    from Scientific.Functions.Interpolation import InterpolatingFunction
    new_scipy = False
from scipy.interpolate import interp1d

from ivs import config
from ivs.units import conversions
from ivs.units import constants
from ivs.aux import loggers
from ivs.aux.decorators import memoized,clear_memoization
import itertools
import functools
from ivs.aux import numpy_ext
from ivs.sed import filters
from ivs.inout import fits
from ivs.sigproc import interpol
from . import reddening
import getpass
import shutil

logger = logging.getLogger("SED.MODEL")
logger.addHandler(loggers.NullHandler)

caldir = os.sep.join(['sedtables','calibrators'])

#-- default values for grids
__defaults__ = dict(grid='kurucz',odfnew=True,z=+0.0,vturb=2,
                alpha=False,nover=False,                  # KURUCZ
                He=97,                                    # WD
                ct='mlt',                                 # NEMO (convection theory)
                t=1.0,a=0.0,c=0.5,m=1.0,co=1.05,          # MARCS and COMARCS
                Rv=3.1,law='fitzpatrick2004',             # reddening info for integrated grids
                use_scratch=False)

defaults = __defaults__.copy()
defaults_multiple = [defaults.copy(),defaults.copy()]
#-- relative location of the grids
basedir = 'sedtables/modelgrids/'
scratchdir = None

#{ Interface to library

def set_defaults(*args,**kwargs):
    """
    Set defaults of this module

    If you give no keyword arguments, the default values will be reset.
    """
    clear_memoization(keys=['ivs.sed.model'])
    #-- these are the default defaults
    if not kwargs:
        kwargs = __defaults__.copy()

    for key in kwargs:
        if key in defaults:
            defaults[key] = kwargs[key]
            logger.info('Set %s to %s'%(key,kwargs[key]))


def set_defaults_multiple(*args):
    """
    Set defaults for multiple stars
    """
    if not args:
        args = [defaults for i in range(len(defaults_multiple))]
    for i,arg in enumerate(args):
        for key in arg:
            if key in defaults_multiple[i]:
                defaults_multiple[i][key] = arg[key]
                logger.info('Set %s to %s (star %d)'%(key,arg[key],i))

def copy2scratch(**kwargs):
    """
    Copy the grids to the scratch directory to speed up the fitting process.
    Files are placed in the directory: /scratch/uname/ where uname is your username.

    This function checks the grids that are set with the functions set_defaults()
    and set_defaults_multiple(). Every time a grid setting is changed, this
    function needs to be called again.

    Don`t forget to remove the files from the scratch directory after the fitting
    process is completed with clean_scratch()

    It is possible to give z='*' and Rv='*' as an option; when you do that, the grids
    with all z, Rv values are copied. Don't forget to add that option to clean_scratch too!
    """
    global scratchdir
    uname = getpass.getuser()
    if not os.path.isdir('/scratch/%s/'%(uname)):
        os.makedirs('/scratch/%s/'%(uname))
    scratchdir = '/scratch/%s/'%(uname)

    #-- we have defaults for the single and multiple grid
    defaults_ = []
    defaults_.append(defaults)
    defaults_.extend(defaults_multiple)

    #-- now run over the defaults for the single and multiple grid, and
    #   copy the necessary files to the scratch disk
    for default in defaults_:
        default['use_scratch'] = False
        #-- set the z with the starred version '*' if asked for, but remember
        #   the original value to reset it after the loop is done.
        originalDefaults = {}
        for key in kwargs:
            if key in default:
                originalDefaults[key] = default[key]
                default[key] = kwargs[key]
                logger.debug('Using provided value for {0:s}={1:s} when copying to scratch'.format(key,str(kwargs[key])))
        #grid
        fname = get_file(integrated=False,**default)
        #-- we could have received a list (multiple files) or a string (single file)
        if isinstance(fname,str):
            fname = [fname]
        for ifname in fname:
            if not os.path.isfile(scratchdir + os.path.basename(ifname)):
                shutil.copy(ifname,scratchdir)
                logger.info('Copied grid: %s to scratch'%(ifname))
            else:
                logger.info('Using existing grid: %s from scratch'%(os.path.basename(ifname)))
        #integrated grid
        fname = get_file(integrated=True,**default)
        if isinstance(fname,str):
            fname = [fname]
        for ifname in fname:
            if not os.path.isfile(scratchdir + os.path.basename(ifname)):
                shutil.copy(ifname,scratchdir)
                logger.info('Copied grid: %s to scratch'%(ifname))
            else:
                logger.info('Using existing grid: %s from scratch'%(os.path.basename(ifname)))
        default['use_scratch'] = True
        for key in kwargs:
            if key in default:
                default[key] = originalDefaults[key]

def clean_scratch(**kwargs):
    """
    Remove the grids that were copied to the scratch directory by using the
    function copy2scratch(). Be carefull with this function, as it doesn't check
    if the models are still in use. If you are running multiple scripts that
    use the same models, only clean the scratch disk after the last script is
    finnished.
    """
    defaults_ = []
    defaults_.append(defaults)
    defaults_.extend(defaults_multiple)

    for default in defaults_:
        if default['use_scratch']:
            originalDefaults = {}
            for key in kwargs:
                if key in default:
                    originalDefaults[key] = default[key]
                    default[key] = kwargs[key]
                    logger.debug('Using provided value for {0:s}={1:s} when deleting from scratch'.format(key,str(kwargs[key])))

            #if z is not None:
                #previous_z = default['z']
                #default['z']
            fname = get_file(integrated=False,**default)
            if isinstance(fname,str):
                fname = [fname]
            for ifname in fname:
                if os.path.isfile(ifname):
                    logger.info('Removed file: %s'%(ifname))
                    os.remove(ifname)

            fname = get_file(integrated=True,**default)
            if isinstance(fname,str):
                fname = [fname]
            for ifname in fname:
                if os.path.isfile(ifname):
                    logger.info('Removed file: %s'%(ifname))
                    os.remove(ifname)
            default['use_scratch'] = False
            for key in kwargs:
                if key in default:
                    default[key] = originalDefaults[key]
            #if z is not None:
                #default['z'] = previous_z

def defaults2str():
    """
    Convert the defaults to a string, e.g. for saving files.
    """
    return '_'.join([str(i)+str(defaults[i]) for i in sorted(defaults.keys())])

def defaults_multiple2str():
    """
    Convert the defaults to a string, e.g. for saving files.
    """
    return '_'.join([str(i)+str(defaults[i]) for defaults in defaults_multiple for i in sorted(sorted(defaults.keys()))])


def get_gridnames(grid=None):
    """
    Return a list of available grid names.

    If you specificy the grid's name, you get two lists: one with all available
    original, non-integrated grids, and one with the pre-calculated photometry.

    @parameter grid: name of the type of grid (optional)
    @type grid: string
    @return: list of grid names
    @rtype: list of str
    """
    if grid is None:
        return ['kurucz','fastwind','cmfgen','sdb_uli','wd_boris','wd_da','wd_db',
                'tlusty','uvblue','atlas12','nemo','tkachenko','marcs','marcs2','tmap',
                ]
                #'marcs','marcs2','comarcs','tlusty','uvblue','atlas12']
    else:
        files = config.glob(basedir,'*%s*.fits'%(grid))
        integrated = [os.path.basename(ff) for ff in files if os.path.basename(ff)[0]=='i']
        original = [os.path.basename(ff) for ff in files if os.path.basename(ff)[0]!='i']
        return original,integrated



def get_file(integrated=False,**kwargs):
    """
    Retrieve the filename containing the specified SED grid.

    The keyword arguments are specific to the kind of grid you're using.

    Basic keywords are 'grid' for the name of the grid, and 'z' for metallicity.
    For other keywords, see the source code.

    Available grids and example keywords:
        - grid='kurucz93':
                    - metallicity (z): m01 is -0.1 log metal abundance relative to solar (solar abundances from Anders and Grevesse 1989)
                    - metallicity (z): p01 is +0.1 log metal abundance relative to solar (solar abundances from Anders and Grevesse 1989)
                    - alpha enhancement (alpha): True means alpha enhanced (+0.4)
                    - turbulent velocity (vturb): vturb in km/s
                    - nover= True means no overshoot
                    - odfnew=True means no overshoot but with better opacities and abundances
        - grid='tlusty':
                    - z: log10(Z/Z0)
        - grid='sdb_uli': metallicity and helium fraction (z, he=98)
        - grid='fastwind': no options
        - grid='wd_boris': no options
        - grid='stars': precomputed stars (vega, betelgeuse...)
        - grid='uvblue'
        - grid='marcs'
        - grid='marcs2'
        - grid='atlas12'
        - grid='tkachenko': metallicity z
        - grid='nemo': convection theory and metallicity (CM=Canuto and Mazzitelli 1991),
        (CGM=Canuto,Goldman,Mazzitelli 1996), (MLT=mixinglengththeory a=0.5)
        - grid='marcsjorissensp': high resolution spectra from 4000 to 25000 A of (online available) MARCS grid computed by A. Jorissen
        with turbospectrum v12.1.1 in late 2012, then converted to the Kurucz wavelength grid (by S. Bloemen and M. Hillen).

    @param integrated: choose integrated version of the gridcopy2scratch
    @type integrated: boolean
    @keyword grid: gridname (default Kurucz)
    @type grid: str
    @return: gridfile
    @rtype: str
    """

    #-- possibly you give a filename
    grid = kwargs.get('grid',defaults['grid'])
    use_scratch = kwargs.get('use_scratch',defaults['use_scratch'])
    if os.path.isfile(grid):
        logger.debug('Selected %s'%(grid))
        if integrated:
            basename = os.path.basename(grid)
            return os.path.join(os.path.dirname(grid),basename[0]=='i' and basename or 'i'+basename)
        logging.debug('Returning grid path: '+grid)
        return grid

    grid = grid.lower()

    #-- general
    z = kwargs.get('z',defaults['z'])
    Rv = kwargs.get('Rv',defaults['Rv'])
    #-- only for Kurucz
    vturb = int(kwargs.get('vturb',defaults['vturb']))
    odfnew = kwargs.get('odfnew',defaults['odfnew'])
    alpha = kwargs.get('alpha',defaults['alpha'])
    nover = kwargs.get('nover',defaults['nover'])
    #-- only for WD
    He = int(kwargs.get('He',defaults['He']))
    #-- only for Marcs and COMarcs
    t = kwargs.get('t',defaults['t'])
    a = kwargs.get('a',defaults['a'])
    c = kwargs.get('c',defaults['c'])
    m = kwargs.get('m',defaults['m'])
    co= kwargs.get('co',defaults['co'])
    #-- only for Nemo
    ct = kwargs.get('ct','mlt')

    #-- figure out what grid to use
    if grid=='fastwind':
        basename = 'fastwind_sed.fits'
    elif grid in ['kurucz','kurucz2']:
        if not isinstance(z,str): z = '%.1f'%(z)
        if not isinstance(vturb,str): vturb = '%d'%(vturb)
        if grid=='kurucz2' and integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        if not alpha and not nover and not odfnew:
            basename = 'kurucz93_z%s_k%s_sed%sls.fits'%(z,vturb,postfix)
        elif alpha and odfnew:
            basename = 'kurucz93_z%s_ak%sodfnew_sed%s.fits'%(z,vturb,postfix)
        elif odfnew:
            basename = 'kurucz93_z%s_k%sodfnew_sed%s.fits'%(z,vturb,postfix)
        elif nover:
            basename = 'kurucz93_z%s_k%snover_sed%s.fits'%(z,vturb,postfix)
    elif grid=='cmfgen':
        basename = 'cmfgen_sed.fits'
    elif grid=='sdb_uli':
        if not isinstance(z,str): z = '%.1f'%(z)
        if not isinstance(He,str): He = '%d'%(He)
        basename = 'SED_int_h%s_z%s.fits'%(He,z)
    elif grid=='wd_boris':
        basename = 'SED_WD_Gaensicke.fits'
    elif grid=='wd_da':
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'SED_WD_Koester_DA%s.fits'%(postfix)
    elif grid=='wd_db':
        basename = 'SED_WD_Koester_DB.fits'
    elif grid=='marcs':
        if not isinstance(z,str): z = '%.1f'%(z)
        if not isinstance(t,str): t = '%.1f'%(t)
        if not isinstance(a,str): a = '%.2f'%(a)
        if not isinstance(c,str): c = '%.2f'%(c)
        basename = 'marcsp_z%st%s_a%s_c%s_sed.fits'%(z,t,a,c)
    elif grid=='marcs2':
        basename = 'marcsp2_z0.00t2.0_m.1.0c0.00_sed.fits'
    elif grid=='comarcs':
        if not isinstance(z,str): z = '%.2f'%(z)
        if not isinstance(co,str): co = '%.2f'%(co)
        if not isinstance(m,str): m = '%.1f'%(m)
        basename = 'comarcsp_z%sco%sm%sxi2.50_sed.fits'%(z,co,m)
    elif grid=='stars':
        basename = 'kurucz_stars_sed.fits'
    elif grid=='tlusty':
        if not isinstance(z,str): z = '%.2f'%(z)
        basename = 'tlusty_z%s_sed.fits'%(z)
    elif grid=='uvblue':
        if not isinstance(z,str): z = '%.1f'%(z)
        basename = 'uvblue_z%s_k2_sed.fits'%(z)
    elif grid=='atlas12':
        if not isinstance(z,str): z = '%.1f'%(z)
        basename = 'atlas12_z%s_sed.fits'%(z)
    elif grid=='tkachenko':
        if not isinstance(z,str): z = '%.2f'%(z)
        basename = 'tkachenko_z%s.fits'%(z)
    elif grid=='nemo':
        ct = ct.lower()
        if ct=='mlt': ct = ct+'072'
        else: ct = ct+'288'
        basename = 'nemo_%s_z%.2f_v%d.fits'%(ct,z,vturb)
    elif grid=='tmap':
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'TMAP2012_lowres%s.fits'%(postfix) #only available for 1 metalicity
    elif grid=='heberb':
        basename = 'Heber2000_B_h909_extended.fits' #only 1 metalicity
    elif grid=='hebersdb':
        basename = 'Heber2000_sdB_h909_extended.fits' #only 1 metalicity

    elif grid=='tmapsdb':
        # grids for sdB star fitting (JorisV)
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'TMAP2012_sdB_extended%s.fits'%(postfix)
    elif grid=='kuruczsdb':
        # grids for sdB star fitting (JorisV)
        if not isinstance(z,str): z = '%.1f'%(z)
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'kurucz_z%s_sdB%s.fits'%(z,postfix)
    elif grid=='kuruczpagb':
        # grids for post-AGB star fitting (MichelH)
        if not isinstance(z,str): z = '%.1f'%(z)
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'kurucz_pAGB_z%s_sed%s.fits'%(z,postfix)

    elif grid=='tmaptest':
        """ Grids exclusively for testing purposes"""
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'TMAP2012_SEDtest%s.fits'%(postfix) #only available for 1 metalicity

    elif grid=='kurucztest':
        """ Grids exclusively for testing purposes"""
        if not isinstance(z,str): z = '%.1f'%(z)
        if integrated:
            postfix = '_lawfitzpatrick2004_Rv'
            if not isinstance(Rv,str): Rv = '{:.2f}'.format(Rv)
            postfix+= Rv
        else:
            postfix = ''
        basename = 'kurucz_SEDtest_z%s%s.fits'%(z,postfix) #only available for 1 metalicity
    elif grid=='marcsjorissensp':
        if not isinstance(z,str): z = '%.2f'%(z)
        if not isinstance(a,str): a = '%.2f'%(a)
        if not isinstance(m,str): m = '%.1f'%(m)
        basename = 'MARCSJorissenSp_m%s_z%s_a%s_sed.fits'%(m,z,a)
    else:
        raise ValueError("Grid {} is not recognized: either give valid descriptive arguments, or give an absolute filepath".format(grid))
    #-- retrieve the absolute path of the file and check if it exists:
    if not '*' in basename:
        if use_scratch:
            if integrated:
                grid = scratchdir+'i'+basename
            else:
                grid = scratchdir+basename
        else:
            if integrated:
                grid = config.get_datafile(basedir,'i'+basename)
            else:
                grid = config.get_datafile(basedir,basename)
    #-- we could also ask for a list of files, when wildcards are given:
    else:
        grid = config.glob(basedir,'i'+basename)
        if use_scratch:
            if integrated:
                grid = glob.glob(scratchdir+'i'+basename)
            else:
                grid = glob.glob(scratchdir+basename)
        else:
            if integrated:
                grid = config.glob(basedir,'i'+basename)
            else:
                grid = config.glob(basedir,basename)

    #grid.sort()
    logger.debug('Returning grid path(s): %s'%(grid))
    return grid


def _blackbody_input(fctn):
    """
    Prepare input and output for blackbody-like functions.

    If the user gives wavelength units and Flambda units, we only need to convert
    everything to SI (and back to the desired units in the end).

    If the user gives frequency units and Fnu units, we only need to convert
    everything to SI ( and back to the desired units in the end).

    If the user gives wavelength units and Fnu units, we need to convert
    the wavelengths first to frequency.
    """
    @functools.wraps(fctn)
    def dobb(x,T,**kwargs):
        wave_units = kwargs.get('wave_units','AA')
        flux_units = kwargs.get('flux_units','erg/s/cm2/AA')
        to_kwargs = {}
        #-- prepare input
        #-- what kind of units did we receive?
        curr_conv = constants._current_convention
        # X: wavelength/frequency
        x_unit_type = conversions.get_type(wave_units)
        x = conversions.convert(wave_units,curr_conv,x)
        # T: temperature
        if isinstance(T,tuple):
            T = conversions.convert(T[1],'K',T[0])
        # Y: flux
        y_unit_type = conversions.change_convention('SI',flux_units)
        #-- if you give Jy vs micron, we need to first convert wavelength to frequency
        if y_unit_type=='kg1 rad-1 s-2' and x_unit_type=='length':
            x = conversions.convert(conversions._conventions[curr_conv]['length'],'rad/s',x)
            x_unit_type = 'frequency'
        elif y_unit_type=='kg1 m-1 s-3' and x_unit_type=='frequency':
            x = conversions.convert('rad/s',conversions._conventions[curr_conv]['length'],x)
            x_unit_type = 'length'
        elif not y_unit_type in ['kg1 rad-1 s-2','kg1 m-1 s-3']:
            raise NotImplementedError(flux_units,y_unit_type)
        #-- correct for rad
        if x_unit_type=='frequency':
            x /= (2*np.pi)
            to_kwargs['freq'] = (x,'Hz')
        else:
            to_kwargs['wave'] = (x,conversions._conventions[curr_conv]['length'])
        #-- run function
        I = fctn((x,x_unit_type),T)

        #-- prepare output
        disc_integrated = kwargs.get('disc_integrated',True)
        ang_diam = kwargs.get('ang_diam',None)
        if disc_integrated:
            I *= np.sqrt(2*np.pi)
            if ang_diam is not None:
                scale = conversions.convert(ang_diam[1],'sr',ang_diam[0]/2.)
                I *= scale
        I = conversions.convert(curr_conv,flux_units,I,**to_kwargs)
        return I

    return dobb






@_blackbody_input
def blackbody(x,T,wave_units='AA',flux_units='erg/s/cm2/AA',disc_integrated=True,ang_diam=None):
    """
    Definition of black body curve.

    To get them into the same units as the Kurucz disc-integrated SEDs, they are
    multiplied by sqrt(2*pi) (set C{disc_integrated=True}).

    You can only give an angular diameter if disc_integrated is True.

    To convert the scale parameter back to mas, simply do:

    ang_diam = 2*conversions.convert('sr','mas',scale)

    See decorator L{blackbody_input} for details on how the input parameters
    are handled: the user is free to choose wavelength or frequency units, choose
    *which* wavelength or frequency units, and can even mix them. To be sure that
    everything is handled correctly, we need to do some preprocessing and unit
    conversions.

    Be careful when, e.g. during fitting, scale contains an error: be sure to set
    the option C{unpack=True} in the L{conversions.convert} function!

    >>> x = np.linspace(2.3595,193.872,500)
    >>> F1 = blackbody(x,280.,wave_units='AA',flux_units='Jy',ang_diam=(1.,'mas'))
    >>> F2 = rayleigh_jeans(x,280.,wave_units='micron',flux_units='Jy',ang_diam=(1.,'mas'))
    >>> F3 = wien(x,280.,wave_units='micron',flux_units='Jy',ang_diam=(1.,'mas'))


    >>> p = pl.figure()
    >>> p = pl.subplot(121)
    >>> p = pl.plot(x,F1)
    >>> p = pl.plot(x,F2)
    >>> p = pl.plot(x,F3)


    >>> F1 = blackbody(x,280.,wave_units='AA',flux_units='erg/s/cm2/AA',ang_diam=(1.,'mas'))
    >>> F2 = rayleigh_jeans(x,280.,wave_units='micron',flux_units='erg/s/cm2/AA',ang_diam=(1.,'mas'))
    >>> F3 = wien(x,280.,wave_units='micron',flux_units='erg/s/cm2/AA',ang_diam=(1.,'mas'))


    >>> p = pl.subplot(122)
    >>> p = pl.plot(x,F1)
    >>> p = pl.plot(x,F2)
    >>> p = pl.plot(x,F3)


    @param: wavelength
    @type: ndarray
    @param T: temperature, unit
    @type: tuple (float,str)
    @param wave_units: wavelength units (frequency or length)
    @type wave_units: str (units)
    @param flux_units: flux units (could be in Fnu-units or Flambda-units)
    @type flux_units: str (units)
    @param disc_integrated: if True, they are in the same units as Kurucz-disc-integrated SEDs
    @type disc_integrated: bool
    @param ang_diam: angular diameter (in mas or rad or something similar)
    @type ang_diam: (value, unit)
    @return: intensity
    @rtype: array
    """
    x,x_unit_type = x
    #-- make the appropriate black body
    if x_unit_type=='frequency': # frequency units
        factor = 2.0 * constants.hh / constants.cc**2
        expont = constants.hh / (constants.kB*T)
        I = factor * x**3 * 1. / (np.exp(expont*x) - 1.)
    elif x_unit_type=='length': # wavelength units
        factor = 2.0 * constants.hh * constants.cc**2
        expont = constants.hh*constants.cc / (constants.kB*T)
        I = factor / x**5. * 1. / (np.exp(expont/x) - 1.)
    else:
        raise ValueError(x_unit_type)
    return I


@_blackbody_input
def rayleigh_jeans(x,T,wave_units='AA',flux_units='erg/s/cm2/AA',disc_integrated=True,ang_diam=None):
    """
    Rayleigh-Jeans approximation of a black body.

    Valid at long wavelengths.

    For input details, see L{blackbody}.

    @return: intensity
    @rtype: array
    """
    x,x_unit_type = x
    #-- now make the appropriate model
    if x_unit_type=='frequency': # frequency units
        factor = 2.0 * constants.kB*T / constants.cc**2
        I = factor * x**2
    elif x_unit_type=='length': # wavelength units
        factor = 2.0 * constants.cc * constants.kB*T
        I = factor / x**4.
    else:
        raise ValueError(unit_type)
    return I


@_blackbody_input
def wien(x,T,wave_units='AA',flux_units='erg/s/cm2/AA',disc_integrated=True,ang_diam=None):
    """
    Wien approximation of a black body.

    Valid at short wavelengths.

    For input details, see L{blackbody}.

    @return: intensity
    @rtype: array
    """
    x,x_unit_type = x
    #-- now make the appropriate model
    if x_unit_type=='frequency': # frequency units
        factor = 2.0 * constants.hh / constants.cc**2
        expont = constants.hh / (constants.kB*T)
        I = factor * x**3 * 1. * np.exp(-expont*x)
    elif x_unit_type=='length': # wavelength units
        factor = 2.0 * constants.hh * constants.cc**2
        expont = constants.hh*constants.cc / (constants.kB*T)
        I = factor / x**5. * np.exp(-expont/x)
    else:
        raise ValueError(unit_type)
    return I


def get_table_single(teff=None,logg=None,ebv=None,rad=None,star=None,
              wave_units='AA',flux_units='erg/s/cm2/AA/sr',**kwargs):
    """
    Retrieve the spectral energy distribution of a model atmosphere.

    Wavelengths in A (angstrom)
    Fluxes in Ilambda = ergs/cm2/s/AA/sr, except specified via 'units',

    If you give 'units', and /sr is not included, you are responsible yourself
    for giving an extra keyword with the angular diameter C{ang_diam}, or other
    possibilities offered by the C{units.conversions.convert} function.

    Possibility to redden the fluxes according to the reddening parameter EB_V.

    Extra kwargs can specify the grid type.
    Extra kwargs can specify constraints on the size of the grid to interpolate.
    Extra kwargs can specify reddening law types.
    Extra kwargs can specify information for conversions.

    Example usage:

    >>> from pylab import figure,gca,subplot,title,gcf,loglog
    >>> p = figure(figsize=(10,6))
    >>> p=gcf().canvas.set_window_title('Test of <get_table>')
    >>> p = subplot(131)
    >>> p = loglog(*get_table(grid='FASTWIND',teff=35000,logg=4.0),**dict(label='Fastwind'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=35000,logg=4.0),**dict(label='Kurucz'))
    >>> p = loglog(*get_table(grid='TLUSTY',teff=35000,logg=4.0),**dict(label='Tlusty'))
    >>> p = loglog(*get_table(grid='MARCS',teff=5000,logg=2.0),**dict(label='Marcs'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=5000,logg=2.0),**dict(label='Kurucz'))
    >>> p = pl.xlabel('Wavelength [angstrom]');p = pl.ylabel('Flux [erg/s/cm2/AA/sr]')
    >>> p = pl.legend(loc='upper right',prop=dict(size='small'))
    >>> p = subplot(132)
    >>> p = loglog(*get_table(grid='FASTWIND',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='Fastwind'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='Kurucz'))
    >>> p = loglog(*get_table(grid='TLUSTY',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='Tlusty'))
    >>> p = loglog(*get_table(grid='MARCS',teff=5000,logg=2.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='Marcs'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=5000,logg=2.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='Kurucz'))
    >>> p = pl.xlabel('Wavelength [micron]');p = pl.ylabel('Flux [Jy/sr]')
    >>> p = pl.legend(loc='upper right',prop=dict(size='small'))
    >>> p = subplot(133);p = title('Kurucz')
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10000,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='10000'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10250,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='10250'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10500,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='10500'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10750,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='10750'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=11000,logg=4.0,wave_units='micron',flux_units='Jy/sr'),**dict(label='11000'))
    >>> p = pl.xlabel('Wavelength [micron]');p = pl.ylabel('Flux [Jy/sr]')
    >>> p = pl.legend(loc='upper right',prop=dict(size='small'))

    ]]include figure]]ivs_sed_model_comparison.png]

    @param teff: effective temperature
    @type teff: float
    @param logg: logarithmic gravity (cgs)
    @type logg: float
    @param ebv: reddening coefficient
    @type ebv: float
    @param wave_units: units to convert the wavelengths to (if not given, A)
    @type wave_units: str
    @param flux_units: units to convert the fluxes to (if not given, erg/s/cm2/AA/sr)
    @type flux_units: str
    @return: wavelength,flux
    @rtype: (ndarray,ndarray)
    """

    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)

    #-- read the file:
    ff = pf.open(gridfile)
    #-- a possible grid is the one where only selected stellar models are
    #   present. In that case, we there is no need for interpolation or
    #   other stuff.
    if star is not None:
        wave = ff[star.upper()].data.field('wavelength')
        flux = ff[star.upper()].data.field('flux')
    else:
        teff = float(teff)
        logg = float(logg)
        #-- if we have a grid model, no need for interpolation
        try:
            #-- extenstion name as in fits files prepared by Steven
            mod_name = "T%05d_logg%01.02f" %(teff,logg)
            mod = ff[mod_name]
            wave = mod.data.field('wavelength')
            flux = mod.data.field('flux')
            logger.debug('Model SED taken directly from file (%s)'%(os.path.basename(gridfile)))
        #-- if the teff/logg is not present, use the interpolation thing
        except KeyError:
            #-- it is possible we first have to set the interpolation function.
            #   This function is memoized, so if it will not be calculated
            #   twice.
            wave,teffs,loggs,flux,flux_grid = get_grid_mesh(**kwargs)
            logger.debug('Model SED interpolated from grid %s (%s)'%(os.path.basename(gridfile),kwargs))
            wave = wave + 0.
            flux = flux_grid(np.log10(teff),logg) + 0.

    #-- convert to arrays
    wave = np.array(wave,float)
    flux = np.array(flux,float)
    #-- redden if necessary
    if ebv is not None and ebv>0:
        if 'wave' in list(kwargs.keys()):
            removed = kwargs.pop('wave')
        flux = reddening.redden(flux,wave=wave,ebv=ebv,rtype='flux',**kwargs)
    if flux_units!='erg/s/cm2/AA/sr':
        flux = conversions.convert('erg/s/cm2/AA/sr',flux_units,flux,wave=(wave,'AA'),**kwargs)
    if wave_units!='AA':
        wave = conversions.convert('AA',wave_units,wave,**kwargs)

    ff.close()

    if rad != None:
        flux = rad**2 * flux

    return wave,flux


def get_itable_single(teff=None,logg=None,ebv=0,z=0,rad=None,photbands=None,
               wave_units=None,flux_units='erg/s/cm2/AA/sr',**kwargs):
    """
    Retrieve the spectral energy distribution of a model atmosphere in
    photometric passbands.

    Wavelengths in A (angstrom). If you set 'wavelengths' to None, no effective
    wavelengths will be calculated. Otherwise, the effective wavelength is
    calculated taking the model flux into account.
    Fluxes in Ilambda = ergs/cm2/s/AA/sr, except specified via 'units',

    If you give 'units', and /sr is not included, you are responsible yourself
    for giving an extra keyword with the angular diameter C{ang_diam}, or other
    possibilities offered by the C{units.conversions.convert} function.

    Possibility to redden the fluxes according to the reddening parameter EB_V.

    Extra kwargs can specify the grid type.
    Extra kwargs can specify constraints on the size of the grid to interpolate.
    Extra kwargs can specify reddening law types.
    Extra kwargs can specify information for conversions.

    @param teff: effective temperature
    @type teff: float
    @param logg: logarithmic gravity (cgs)
    @type logg: float
    @param ebv: reddening coefficient
    @type ebv: float
    @param photbands: photometric passbands
    @type photbands: list of photometric passbands
    @param wave_units: units to convert the wavelengths to (if not given, A)
    @type wave_units: str
    @param flux_units: units to convert the fluxes to (if not given, erg/s/cm2/AA/sr)
    @type flux_units: str
    @keyword clear_memory: flag to clear memory from previously loaded SED tables.
    If you set it to False, you can easily get an overloaded memory!
    @type clear_memory: boolean
    @return: (wave,) flux, absolute luminosity
    @rtype: (ndarray,)ndarray,float
    """
    if 'vrad' in kwargs:
        logger.debug('vrad is NOT taken into account when interpolating in get_itable()')
    if 'rv' in kwargs:
        logger.debug('Rv is NOT taken into account when interpolating in get_itable()')

    if photbands is None:
        raise ValueError('no photometric passbands given')
    ebvrange = kwargs.pop('ebvrange',(-np.inf,np.inf))
    zrange = kwargs.pop('zrange',(-np.inf,np.inf))
    clear_memory = kwargs.pop('clear_memory',True)
    #-- get the FITS-file containing the tables
    #c0 = time.time()
    #c1 = time.time() - c0
    #-- retrieve structured information on the grid (memoized)
    markers,(g_teff,g_logg,g_ebv,g_z),gpnts,ext = _get_itable_markers(photbands,ebvrange=ebvrange,zrange=zrange,
                            include_Labs=True,clear_memory=clear_memory,**kwargs)
    #c2 = time.time() - c0 - c1
    #-- if we have a grid model, no need for interpolation
    try:
        input_code = float('%3d%05d%03d%03d'%(int(round((z+5)*100)),int(round(teff)),int(round(logg*100)),int(round(ebv*100))))
        index = markers.searchsorted(input_code)
        output_code = markers[index]
        #-- if not available, go on and interpolate!
        #   we raise a KeyError for symmetry with C{get_table}.
        if not input_code==output_code:
            raise KeyError
        #c0_ = time.time()
        flux = ext[index]
        #c1_ = time.time()-c0_
        #flux = np.array([ext.data.field(photband)[index] for photband in photbands])
        #logger.debug('Model iSED taken directly from file (%s)'%(os.path.basename(gridfile)))
    #-- if the teff/logg is not present, use the interpolation thing
    except KeyError:
        #c1_ = 0
        #-- cheat edges in interpolating function
        if not new_scipy:
            teff = teff+1e-2
            logg = logg-1e-6
            ebv = ebv+1e-6
            #-- it is possible we have to interpolate: identify the grid points in
            #   the immediate vicinity of the given fundamental parameters
            i_teff = max(1,g_teff.searchsorted(teff))
            i_logg = max(1,g_logg.searchsorted(logg))
            i_ebv  = max(1,g_ebv.searchsorted(ebv))
            i_z  = max(1,g_z.searchsorted(z))
            if i_teff==len(g_teff): i_teff -= 1
            if i_logg==len(g_logg): i_logg -= 1
            if i_ebv==len(g_ebv): i_ebv -= 1
            if i_z==len(g_z): i_z -= 1
            #-- prepare fluxes matrix for interpolation, and x,y an z axis
            teffs_subgrid = g_teff[i_teff-1:i_teff+1]
            loggs_subgrid = g_logg[i_logg-1:i_logg+1]
            ebvs_subgrid = g_ebv[i_ebv-1:i_ebv+1]
            zs_subgrid = g_z[i_z-1:i_z+1]
            #-- iterates over df-1 values (df=degrees of freedom): we know that the
            #   grid is ordered via z in the last part (about twice as fast).
            #   Reducing the grid size to 2 increases the speed again with a factor 2.

            #-- if metallicity needs to be interpolated
            if not (z in g_z):
                fluxes = np.zeros((2,2,2,2,len(photbands)+1))
                for i,j,k in itertools.product(range(2),range(2),range(2)):
                    input_code = float('%3d%05d%03d%03d'%(int(round((zs_subgrid[i]+5)*100)),\
                                                    int(round(teffs_subgrid[j])),\
                                                    int(round(loggs_subgrid[k]*100)),\
                                                    int(round(ebvs_subgrid[1]*100))))
                    index = markers.searchsorted(input_code)
                    fluxes[i,j,k] = ext[index-1:index+1]
                myf = InterpolatingFunction([zs_subgrid,np.log10(teffs_subgrid),
                                        loggs_subgrid,ebvs_subgrid],np.log10(fluxes),default=-100*np.ones_like(fluxes.shape[1]))
                flux = 10**myf(z,np.log10(teff),logg,ebv) + 0.

            #-- if only teff,logg and ebv need to be interpolated (faster)
            else:
                fluxes = np.zeros((2,2,2,len(photbands)+1))
                for i,j in itertools.product(range(2),range(2)):
                    input_code = float('%3d%05d%03d%03d'%(int(round((z+5)*100)),\
                                                    int(round(teffs_subgrid[i])),\
                                                    int(round(loggs_subgrid[j]*100)),\
                                                    int(round(ebvs_subgrid[1]*100))))
                    index = markers.searchsorted(input_code)
                    fluxes[i,j] = ext[index-1:index+1]
                myf = InterpolatingFunction([np.log10(teffs_subgrid),
                                        loggs_subgrid,ebvs_subgrid],np.log10(fluxes),default=-100*np.ones_like(fluxes.shape[1]))
                flux = 10**myf(np.log10(teff),logg,ebv) + 0.
        #-- new scipy version
        else:
            #-- take care of inner edge of grid
            i_teff = max(1,g_teff.searchsorted(teff))
            i_logg = max(1,g_logg.searchsorted(logg))
            i_ebv  = max(1,g_ebv.searchsorted(ebv))
            i_z  = max(1,g_z.searchsorted(z))
            #-- take care of outer edge of grid
            if i_teff==len(g_teff): i_teff -= 1
            if i_logg==len(g_logg): i_logg -= 1
            if i_ebv==len(g_ebv): i_ebv -= 1
            if i_z==len(g_z): i_z -= 1
            if not (z in g_z):
                #-- prepare fluxes matrix for interpolation, and x,y an z axis
                myflux = np.zeros((16,4+len(photbands)+1))
                mygrid = itertools.product(g_teff[i_teff-1:i_teff+1],g_logg[i_logg-1:i_logg+1],g_z[i_z-1:i_z+1])
                for i,(t,g,zz) in enumerate(mygrid):
                    myflux[2*i,:4] = t,g,g_ebv[i_ebv-1],zz
                    myflux[2*i+1,:4] = t,g,g_ebv[i_ebv],zz
                    input_code = float('%3d%05d%03d%03d'%(int(round((zz+5)*100)),\
                                                        int(round(t)),int(round(g*100)),\
                                                        int(round(g_ebv[i_ebv]*100))))
                    index = markers.searchsorted(input_code)
                    myflux[2*i,4:] = ext[index-1]
                    myflux[2*i+1,4:] = ext[index]
                #-- interpolate in log10 of temperature
                myflux[:,0] = np.log10(myflux[:,0])
                flux = 10**griddata(myflux[:,:4],np.log10(myflux[:,4:]),(np.log10(teff),logg,ebv,z))
            else:
                #-- prepare fluxes matrix for interpolation, and x,y axis
                myflux = np.zeros((8,3+len(photbands)+1))
                mygrid = itertools.product(g_teff[i_teff-1:i_teff+1],g_logg[i_logg-1:i_logg+1])
                for i,(t,g) in enumerate(mygrid):
                    myflux[2*i,:3] = t,g,g_ebv[i_ebv-1]
                    myflux[2*i+1,:3] = t,g,g_ebv[i_ebv]
                    input_code = float('%3d%05d%03d%03d'%(int(round((z+5)*100)),\
                                                        int(round(t)),int(round(g*100)),\
                                                        int(round(g_ebv[i_ebv]*100))))
                    index = markers.searchsorted(input_code)
                    myflux[2*i,3:] = ext[index-1]
                    myflux[2*i+1,3:] = ext[index]
                #-- interpolate in log10 of temperature
                myflux[:,0] = np.log10(myflux[:,0])
                flux = 10**griddata(myflux[:,:3],np.log10(myflux[:,3:]),(np.log10(teff),logg,ebv))
    except IndexError:
        #-- probably metallicity outside of grid
        raise ValueError('point outside of grid (teff={teff}, logg={logg}, ebv={ebv}, z={z}'.format(**locals()))
    except ValueError:
        #-- you tried to make a code of a negative number
        raise ValueError('point outside of grid (teff={teff}, logg={logg}, ebv={ebv}, z={z}'.format(**locals()))
    if np.any(np.isnan(flux)):
        #-- you tried to make a code of a negative number
        raise ValueError('point outside of grid (teff={teff}, logg={logg}, ebv={ebv}, z={z}'.format(**locals()))
    if np.any(np.isinf(flux)):
        flux = np.zeros(fluxes.shape[-1])
    #return flux[:-1],flux[-1]#,np.array([c1_,c2,c3])

    #-- convert to arrays: remember that last column of the fluxes is actually
    #   absolute luminosity
    flux,Labs = np.array(flux[:-1],float),flux[-1]

    #-- Take radius into account when provided
    if rad != None:
        flux,Labs = flux*rad**2, Labs*rad**2

    if flux_units!='erg/s/cm2/AA/sr':
        flux = conversions.nconvert('erg/s/cm2/AA/sr',flux_units,flux,photband=photbands,**kwargs)

    if wave_units is not None:
        model = get_table(teff=teff,logg=logg,ebv=ebv,**kwargs)
        wave = filters.eff_wave(photbands,model=model)
        if wave_units !='AA':
            wave = wave = conversions.convert('AA',wave_units,wave,**kwargs)

        return wave,flux,Labs
    else:
        return flux,Labs

def get_itable(photbands=None, wave_units=None, flux_units='erg/s/cm2/AA/sr',
                                                        grids=None, **kwargs):
    """
    Retrieve the integrated spectral energy distribution of a combined model
    atmosphere.

    >>> teff1,teff2 = 20200,5100
    >>> logg1,logg2 = 4.35,2.00
    >>> ebv = 0.2,0.2
    >>> photbands = ['JOHNSON.U','JOHNSON.V','2MASS.J','2MASS.H','2MASS.KS']

    >>> wave1,flux1 = get_table(teff=teff1,logg=logg1,ebv=ebv[0])
    >>> wave2,flux2 = get_table(teff=teff2,logg=logg2,ebv=ebv[1])
    >>> wave3,flux3 = get_table_multiple(teff=(teff1,teff2),logg=(logg1,logg2),ebv=ebv,radius=[1,20])

    >>> iwave1,iflux1,iLabs1 = get_itable(teff=teff1,logg=logg1,ebv=ebv[0],photbands=photbands,wave_units='AA')
    >>> iflux2,iLabs2 = get_itable(teff=teff2,logg=logg2,ebv=ebv[1],photbands=photbands)
    >>> iflux3,iLabs3 = get_itable_multiple(teff=(teff1,teff2),logg=(logg1,logg2),z=(0,0),ebv=ebv,radius=[1,20.],photbands=photbands)

    >>> p = pl.figure()
    >>> p = pl.gcf().canvas.set_window_title('Test of <get_itable_multiple>')
    >>> p = pl.loglog(wave1,flux1,'r-')
    >>> p = pl.loglog(iwave1,iflux1,'ro',ms=10)
    >>> p = pl.loglog(wave2,flux2*20**2,'b-')
    >>> p = pl.loglog(iwave1,iflux2*20**2,'bo',ms=10)
    >>> p = pl.loglog(wave3,flux3,'k-',lw=2)
    >>> p = pl.loglog(iwave1,iflux3,'kx',ms=10,mew=2)

    @param teff: effective temperature
    @type teff: tuple floats
    @param logg: logarithmic gravity (cgs)
    @type logg: tuple floats
    @param ebv: reddening coefficient
    @type ebv: tuple floats
    @param z: metallicity
    @type z: tuple floats
    @param radius: ratio of R_i/(R_{i-1})
    @type radius: tuple of floats
    @param photbands: photometric passbands
    @type photbands: list
    @param flux_units: units to convert the fluxes to (if not given, erg/s/cm2/AA/sr)
    @type flux_units: str
    @param grids: specifications for grid1
    @type grids: list of dict
    @param full_output: return all individual SEDs
    @type full_output: boolean
    @return: wavelength,flux
    @rtype: (ndarray,ndarray)
    """
    #-- Find the parameters provided and store them separately.
    values, parameters, components = {}, set(), set()
    for key in list(kwargs.keys()):
        if re.search("^(teff|logg|ebv|z|rad)\d?$", key):
            par, comp = re.findall("^(teff|logg|ebv|z|rad)(\d?)$", key)[0]
            values[key] = kwargs.pop(key)
            parameters.add(par)
            components.add(comp)

    #-- If there is only one component, we can directly return the result
    if len(components) == 1:
        kwargs.update(values)
        return get_itable_single(photbands=photbands,wave_units=wave_units,
                                     flux_units=flux_units,**kwargs)

    #-- run over all fluxes and sum them, we do not need to multiply with the radius
    #   as the radius is provided as an argument to itable_single_pix.
    fluxes, Labs = [],[]
    for i, (comp, grid) in enumerate(zip(components,defaults_multiple)):
        trash = grid.pop('z',0.0), grid.pop('Rv',0.0)
        kwargs_ = kwargs
        kwargs_.update(grid)
        for par in parameters:
            kwargs_[par] = values[par+comp] if par+comp in values else values[par]

        f,L = get_itable_single(photbands=photbands,wave_units=None,**kwargs_)

        fluxes.append(f)
        Labs.append(L)

    fluxes = np.sum(fluxes,axis=0)
    Labs = np.sum(Labs,axis=0)

    if flux_units!='erg/s/cm2/AA/sr':
        fluxes = np.array([conversions.convert('erg/s/cm2/AA/sr',flux_units,fluxes[i],photband=photbands[i]) for i in range(len(fluxes))])

    if wave_units is not None:
        model = get_table_multiple(teff=teff,logg=logg,ebv=ebv, grids=grids,**kwargs)
        wave = filters.eff_wave(photbands,model=model)
        if wave_units !='AA':
            wave = wave = conversions.convert('AA',wave_units,wave)
        return wave,fluxes,Labs

    return fluxes,Labs


    ##-- set default parameters
    #if grids is None:
        #grids = [defaults_multiple[i] for i in range(len(teff))]
    #if radius is None:
        #radius = tuple([1. for i in teff])
    ##-- gather all the SEDs from the individual components
    #fluxes,Labs = [],[]
    #for i in range(len(teff)):
        #iteff,ilogg,iz,irrad,iebv = teff[i],logg[i],z[i],radius[i],ebv[0]
        #mykwargs = dict(list(grids[i].items()) + list(kwargs.items()))
        #if 'z' in mykwargs:
            #thrash = mykwargs.pop('z')
        ##mykwargs = dict(list(kwargs.items()))
        #iflux,iLabs = get_itable(teff=iteff,logg=ilogg,ebv=iebv,z=iz,photbands=photbands,clear_memory=False,**mykwargs)
        #fluxes.append(iflux*irrad**2)
        #Labs.append(iLabs*irrad**2)
    #fluxes = np.sum(fluxes,axis=0)
    #Labs = np.sum(Labs)
    #if flux_units!='erg/s/cm2/AA/sr':
        #fluxes = np.array([conversions.convert('erg/s/cm2/AA/sr',flux_units,fluxes[i],photband=photbands[i]) for i in range(len(fluxes))])

    #if wave_units is not None:
        #model = get_table_multiple(teff=teff,logg=logg,ebv=ebv, grids=grids,**kwargs)
        #wave = filters.eff_wave(photbands,model=model)
        #if wave_units !='AA':
            #wave = wave = conversions.convert('AA',wave_units,wave)
        #return wave,fluxes,Labs
    #return fluxes,Labs

def get_itable_single_pix(teff=None,logg=None,ebv=None,z=0,rv=3.1,vrad=0,photbands=None,
               wave_units=None,flux_units='erg/s/cm2/AA/sr',**kwargs):
    """
    Super fast grid interpolator.

    Possible kwargs are teffrange,loggrange etc.... that are past on to
    L{_get_pix_grid}. You should probably use these options when you want to
    interpolate in many variables; supplying these ranges will make the grid
    smaller and thus decrease memory usage.

    It is possible to fix C{teff}, C{logg}, C{ebv}, C{z}, C{rv} and/or C{vrad}
    to one value, in which case it B{has} to be a point in the grid. If you want
    to retrieve a list of fluxes with the same ebv value that is not in the grid,
    you need to give an array with all equal values. The reason is that the
    script can try to minimize the number of interpolations, by fixing a
    variable on a grid point. The fluxes on the other gridpoints will then not
    be interpolated over. These parameter also have to be listed with the
    additional C{exc_interpolpar} keyword.

    >>> teffs = np.linspace(5000,7000,100)
    >>> loggs = np.linspace(4.0,4.5,100)
    >>> ebvs = np.linspace(0,1,100)
    >>> zs = np.linspace(-0.5,0.5,100)
    >>> rvs = np.linspace(2.2,5.0,100)

    >>> set_defaults(grid='kurucz2')
    >>> flux,labs = get_itable_pix(teffs,loggs,ebvs,zs,rvs,photbands=['JOHNSON.V'])

    >>> names = ['teffs','loggs','ebvs','zs','rvs']
    >>> p = pl.figure()
    >>> for i in range(len(names)):
    ...     p = pl.subplot(2,3,i+1)
    ...     p = pl.plot(locals()[names[i]],flux[0],'k-')
    ...     p = pl.xlabel(names[i])


    Thanks to Steven Bloemen for the core implementation of the interpolation
    algorithm.

    The addition of the exc_interpolpar keyword was done by Michel Hillen (Jan 2016).
    """

    #-- setup some standard values when they are not provided
    ebv = np.array([0 for i in teff]) if ebv is None else ebv
    z = np.array([0.for i in teff]) if z is None else z
    rv = np.array([3.1 for i in teff]) if rv is None else rv
    vrad = np.array([0 for i in teff]) if vrad is None else vrad

    #for var in ['teff','logg','ebv','z','rv','vrad']:
        #if not hasattr(locals()[var],'__iter__'):
            #print var, locals()[var]
            #locals()[var] = np.array([ locals()[var] ])
    #print locals()
    vrad = 0
    N = 1
    clear_memory = kwargs.pop('clear_memory',False)
    #variables = kwargs.pop('variables',['teff','logg','ebv','z','rv','vrad'])                             # !!!
    for var in ['teff','logg','ebv','z','rv','vrad']:                                       # !!!
        if not hasattr(locals()[var],'__iter__'):
            kwargs.setdefault(var+'range',(locals()[var],locals()[var]))
        else:
            N = len(locals()[var])

    #-- retrieve structured information on the grid (memoized)
    axis_values,gridpnts,pixelgrid,cols = _get_pix_grid(photbands,
                            include_Labs=True,clear_memory=clear_memory,**kwargs)     # !!!

    #-- Remove parameters from the grid if it is requested that these should not be interpolated
    #-- (with the exc_interpolpar keyword). This can only work if the requested values of
    #-- these parameters all correspond to a single point in the original grid!
    #-- we check whether this condition is fulfilled
    #-- if not, then the parameter is not excluded from the interpolation
    #-- and a warning is raised to the log
    for var in kwargs.get('exc_interpolpar',[]): # e.g. for Kurucz, var can be anything in ['teff','logg','ebv','z']
        # retrieve the unique values in var
        var_uniquevalue = np.unique(np.array(locals()[var]))
        # if there is more than one unique value in var, then our condition is not fulfilled
        if len(var_uniquevalue) > 1:
            logger.warning('{} is requested to be excluded from interpolation, although fluxes for more than one value are requested!?'.format(var))
        else:
            # retrieve the index of var in the 'pixelgrid' and 'cols' arrays of the original grid
            var_index = np.where(cols == var)[0]
            # retrieve the index of the unique value in the original grid
            var_uniquevalue_index = np.where(axis_values[var_index] == var_uniquevalue[0])[0]
            # if the unique value does not correspond to a grid point of the original grid, then we only raise a warning
            if len(var_uniquevalue_index) == 0:
                logger.warning('{} can only be excluded from interpolation, as requested, if its values are all equal to an actual grid point!'.format(var))
            else:
                # remove var from the list of variables in the original grid
                trash = axis_values.pop(var_index)
                cols = np.delete(cols,[var_index])
                # since we do not know the axis of var in advance, we devise a clever way to
                # bring it to the first axis by transposing the array
                indices = [x for x in range(pixelgrid.ndim)]
                indices.remove(var_index)
                indices.insert(0,var_index)
                pixelgrid = np.transpose(pixelgrid, indices)
                # now we select the subgrid corresponding to the requested value of var
                pixelgrid = pixelgrid[var_uniquevalue_index[0]]

    #-- prepare input:
    values = np.zeros((len(cols),N))
    for i,col in enumerate(cols):
        values[i] = locals()[col]
    pars = 10**interpol.interpolate(values,axis_values,pixelgrid)
    flux,Labs = pars[:-1],pars[-1]

    #-- Take radius into account when provided
    if 'rad' in kwargs:
        flux,Labs = flux*kwargs['rad']**2, Labs*kwargs['rad']**2

    #-- change flux and wavelength units if needed
    if flux_units!='erg/s/cm2/AA/sr':
        flux = conversions.nconvert('erg/s/cm2/AA/sr',flux_units,flux,photband=photbands,**kwargs)
    if wave_units is not None:
        model = get_table(teff=teff,logg=logg,ebv=ebv,**kwargs)
        wave = filters.eff_wave(photbands,model=model)
        if wave_units !='AA':
            wave = conversions.convert('AA',wave_units,wave,**kwargs)
        return wave,flux,Labs
    else:
        return flux,Labs

def get_itable_pix(photbands=None, wave_units=None, flux_units='erg/s/cm2/AA/sr',
                grids=None, **kwargs):
    """
    Super fast grid interpolator for multiple tables, completely based on get_itable_pix.
    """
    #-- Find the parameters provided and store them separately.
    values, parameters, components = {}, set(), set()
    for key in list(kwargs.keys()):
        if re.search("^(teff|logg|ebv|z|rv|vrad|rad)\d?$", key):
            par, comp = re.findall("^(teff|logg|ebv|z|rv|vrad|rad)(\d?)$", key)[0]
            values[key] = kwargs.pop(key)
            parameters.add(par)
            components.add(comp)
    #-- If there is only one component, we can directly return the result
    if len(components) == 1:
        kwargs.update(values)
        return get_itable_single_pix(photbands=photbands,wave_units=wave_units,
                                     flux_units=flux_units,**kwargs)
    #-- run over all fluxes and sum them, we do not need to multiply with the radius
    #   as the radius is provided as an argument to itable_single_pix.
    fluxes, Labs = [],[]
    for i, (comp, grid) in enumerate(zip(components,defaults_multiple)):
        trash = grid.pop('z',0.0), grid.pop('Rv',0.0)
        kwargs_ = kwargs
        kwargs_.update(grid)
        for par in parameters:
            kwargs_[par] = values[par+comp] if par+comp in values else values[par]

        f,L = get_itable_single_pix(photbands=photbands,wave_units=None,**kwargs_)

        fluxes.append(f)
        Labs.append(L)
    fluxes = np.sum(fluxes,axis=0)
    Labs = np.sum(Labs,axis=0)

    if flux_units!='erg/s/cm2/AA/sr':
        fluxes = np.array([conversions.convert('erg/s/cm2/AA/sr',flux_units,fluxes[i],photband=photbands[i]) for i in range(len(fluxes))])

    if wave_units is not None:
        model = get_table_multiple(teff=teff,logg=logg,ebv=ebv, grids=grids,**kwargs)
        wave = filters.eff_wave(photbands,model=model)
        if wave_units !='AA':
            wave = conversions.convert('AA',wave_units,wave)
        return wave,fluxes,Labs
    return fluxes,Labs


#def get_table_multiple(teff=None,logg=None,ebv=None,radius=None,
              #wave_units='AA',flux_units='erg/cm2/s/AA/sr',grids=None,full_output=False,**kwargs):
def get_table(wave_units='AA',flux_units='erg/cm2/s/AA/sr',grids=None,full_output=False,**kwargs):
    """
    Retrieve the spectral energy distribution of a combined model atmosphere.

    Example usage:

    >>> teff1,teff2 = 20200,5100
    >>> logg1,logg2 = 4.35,2.00
    >>> wave1,flux1 = get_table(teff=teff1,logg=logg1,ebv=0.2)
    >>> wave2,flux2 = get_table(teff=teff2,logg=logg2,ebv=0.2)
    >>> wave3,flux3 = get_table_multiple(teff=(teff1,teff2),logg=(logg1,logg2),ebv=(0.2,0.2),radius=[1,20])

    >>> p = pl.figure()
    >>> p = pl.gcf().canvas.set_window_title('Test of <get_table_multiple>')
    >>> p = pl.loglog(wave1,flux1,'r-')
    >>> p = pl.loglog(wave2,flux2,'b-')
    >>> p = pl.loglog(wave2,flux2*20**2,'b--')
    >>> p = pl.loglog(wave3,flux3,'k-',lw=2)

    @param teff: effective temperature
    @type teff: tuple floats
    @param logg: logarithmic gravity (cgs)
    @type logg: tuple floats
    @param ebv: tuple reddening coefficients
    @type ebv: tuple floats
    @param radius: radii of the stars
    @type radius: tuple of floats
    @param wave_units: units to convert the wavelengths to (if not given, A)
    @type wave_units: str
    @param flux_units: units to convert the fluxes to (if not given, erg/s/cm2/AA/sr)
    @type flux_units: str
    @param grids: specifications for grid1
    @type grids: list of dict
    @param full_output: return all individual SEDs
    @type full_output: boolean
    @return: wavelength,flux
    @rtype: (ndarray,ndarray)
    """
    values, parameters, components = {}, set(), set()
    for key in list(kwargs.keys()):
        if re.search("^(teff|logg|ebv|z|rad)\d?$", key):
            par, comp = re.findall("^(teff|logg|ebv|z|rad)(\d?)$", key)[0]
            values[key] = kwargs.pop(key)
            parameters.add(par)
            components.add(comp)

    #-- If there is only one components we can directly return the result
    if len(components) == 1:
        kwargs.update(values)
        return get_table_single(wave_units=wave_units, flux_units=flux_units,
                                     full_output=full_output,**kwargs)

    #-- Run over all fluxes and sum them, we do not need to multiply with the radius
    #   as the radius is provided as an argument to get_table_single.
    waves, fluxes = [],[]
    for i, (comp, grid) in enumerate(zip(components,defaults_multiple)):
        trash = grid.pop('z',0.0), grid.pop('Rv',0.0)
        kwargs_ = kwargs.copy()
        kwargs_.update(grid)
        for par in parameters:
            kwargs_[par] = values[par+comp] if par+comp in values else values[par]

        w,f = get_table_single(**kwargs_)

        waves.append(w)
        fluxes.append(f)

    ##-- set default parameters
    #if grids is None:
        #grids = [defaults_multiple[i] for i in range(len(teff))]
    #if radius is None:
        #radius = tuple([1. for i in teff])
    ##-- gather all the SEDs from the individual components
    #waves,fluxes = [],[]
    #for i in range(len(teff)):
        #iteff,ilogg,iebv = teff[i],logg[i],ebv[i]
        #mykwargs = dict(list(grids[i].items()) + list(kwargs.items()))
        #iwave,iflux = get_table(teff=iteff,logg=ilogg,ebv=iebv,**mykwargs)
        #waves.append(iwave)
        #fluxes.append(iflux)

    #-- what's the total wavelength range? Merge all wavelength arrays and
    #   remove double points
    waves_ = np.sort(np.hstack(waves))
    waves_ = np.hstack([waves_[0],waves_[1:][np.diff(waves_)>0]])
    #   cut out the part which is common to every wavelength range
    wstart = max([w[0] for w in waves])
    wend = min([w[-1] for w in waves])
    waves_ = waves_[( (wstart<=waves_) & (waves_<=wend))]
    if full_output:
        fluxes_ = []
    else:
        fluxes_ = 0.

    ##-- interpolate onto common grid in log!
    #for i,(wave,flux) in enumerate(zip(waves,fluxes)):
        #intf = interp1d(np.log10(wave),np.log10(flux),kind='linear')
        #if full_output:
            #fluxes_.append(radius[i]**2*10**intf(np.log10(waves_)))
        #else:
            #fluxes_ += radius[i]**2*10**intf(np.log10(waves_))
    #-- interpolate onto common grid in log!
    for i,(wave,flux) in enumerate(zip(waves,fluxes)):
        intf = interp1d(np.log10(wave),np.log10(flux),kind='linear')
        if full_output:
            fluxes_.append(10**intf(np.log10(waves_)))
        else:
            fluxes_ += 10**intf(np.log10(waves_))

    if flux_units!='erg/cm2/s/AA/sr':
        fluxes_ = conversions.convert('erg/s/cm2/AA/sr',flux_units,fluxes_,wave=(waves_,'AA'),**kwargs)
    if wave_units!='AA':
        waves_ = conversions.convert('AA',wave_units,waves_,**kwargs)
    #-- where the fluxes are zero, log can do weird
    if full_output:
        fluxes_ = np.vstack(fluxes_)
        keep = -np.isnan(np.sum(fluxes_,axis=0))
        waves_ = waves_[keep]
        fluxes_ = fluxes_[:,keep]
    else:
        keep = -np.isnan(fluxes_)
        waves_ = waves_[keep]
        fluxes_ = fluxes_[keep]
    return waves_,fluxes_


def get_grid_dimensions(**kwargs):
    """
    Retrieve possible effective temperatures and gravities from a grid.

    E.g. kurucz, sdB, fastwind...

    @rtype: (ndarray,ndarray)
    @return: effective temperatures, gravities
    """
    gridfile = get_file(**kwargs)
    ff = pf.open(gridfile)
    teffs = []
    loggs = []
    for mod in ff[1:]:
        teffs.append(float(mod.header['TEFF']))
        loggs.append(float(mod.header['LOGG']))
    ff.close()

    #-- maybe the fits extensions are not in right order...
    matrix = np.vstack([np.array(teffs),np.array(loggs)]).T
    matrix = numpy_ext.sort_order(matrix,order=[0,1])
    teffs,loggs = matrix.T

    return teffs,loggs




def get_igrid_dimensions(**kwargs):
    """
    Retrieve possible effective temperatures, surface gravities and reddenings
    from an integrated grid.

    E.g. kurucz, sdB, fastwind...

    @rtype: (ndarray,ndarray,ndarray)
    @return: effective temperatures, surface, gravities, E(B-V)s
    """
    gridfile = get_file(integrated=True,**kwargs)
    ff = pf.open(gridfile)
    teffs = ff[1].data.field('TEFF')
    loggs = ff[1].data.field('LOGG')
    ebvs  = ff[1].data.field('EBV')
    ff.close()

    #correct = (teffs==14000) & (loggs==2.0)
    #teffs[correct] = 12000

    return teffs,loggs,ebvs







@memoized
def get_grid_mesh(wave=None,teffrange=None,loggrange=None,**kwargs):
    """
    Return InterpolatingFunction spanning the available grid of atmosphere models.

    WARNING: the grid must be entirely defined on a mesh grid, but it does not
    need to be equidistant.

    It is thus the user's responsibility to know whether the grid is evenly
    spaced in logg and teff (e.g. this is not so for the CMFGEN models).

    You can supply your own wavelength range, since the grid models'
    resolution are not necessarily homogeneous. If not, the first wavelength
    array found in the grid will be used as a template.

    It might take a long a time and cost a lot of memory if you load the entire
    grid. Therefor, you can also set range of temperature and gravity.

    WARNING: 30000,50000 did not work out for FASTWIND, since we miss a model!

    Example usage:

    @param wave: wavelength to define the grid on
    @type wave: ndarray
    @param teffrange: starting and ending of the grid in teff
    @type teffrange: tuple of floats
    @param loggrange: starting and ending of the grid in logg
    @type loggrange: tuple of floats
    @return: wavelengths, teffs, loggs and fluxes of grid, and the interpolating
    function
    @rtype: (3x1Darray,3Darray,interp_func)
    """
    #-- get the dimensions of the grid
    teffs,loggs = get_grid_dimensions(**kwargs)
    #-- build flux grid, assuming a perfectly sampled grid (needs not to be
    #   equidistant)
    if teffrange is not None:
        sa = (teffrange[0]<=teffs) & (teffs<=teffrange[1])
        teffs = teffs[sa]
    if loggrange is not None:
        sa = (loggrange[0]<=loggs) & (loggs<=loggrange[1])
        loggs = loggs[sa]

    #-- ScientificPython interface
    if not new_scipy:
        logger.warning('SCIENTIFIC PYTHON')
        #-- clip if necessary
        teffs = list(set(list(teffs)))
        loggs = list(set(list(loggs)))
        teffs = np.sort(teffs)
        loggs = np.sort(loggs)
        if wave is not None:
            flux = np.ones((len(teffs),len(loggs),len(wave)))
        #-- run over teff and logg, and interpolate the models onto the supplied
        #   wavelength range
        gridfile = get_file(**kwargs)
        ff = pf.open(gridfile)
        for i,teff in enumerate(teffs):
            for j,logg in enumerate(loggs):
                try:
                    mod_name = "T%05d_logg%01.02f" %(teff,logg)
                    mod = ff[mod_name]
                    wave_ = mod.data.field('wavelength')#array(mod.data.tolist())[:,0]
                    flux_ = mod.data.field('flux')#array(mod.data.tolist())[:,1]
                    #-- if there is no wavelength range given, we assume that
                    #   the whole grid has the same resolution, and the first
                    #   wave-array will be used as a template
                    if wave is None:
                        wave = wave_
                        flux = np.ones((len(teffs),len(loggs),len(wave)))
                except KeyError:
                    continue
                #-- it could be that we're lucky and the grid is completely
                #   homogeneous. In that case, there is no need for interpolation
                try:
                    flux[i,j,:] = flux_
                except:
                    flux[i,j,:] = np.interp(wave,wave_,flux_)
        ff.close()
        flux_grid = InterpolatingFunction([np.log10(teffs),loggs],flux)
        logger.info('Constructed SED interpolation grid')

    #-- Scipy interface
    else:
        logger.warning('SCIPY')
        #-- run over teff and logg, and interpolate the models onto the supplied
        #   wavelength range
        gridfile = get_file(**kwargs)
        ff = pf.open(gridfile)
        if wave is not None:
            fluxes = np.zeros((len(teffs),len(wave)))
        for i,(teff,logg) in enumerate(zip(teffs,loggs)):
            mod_name = "T%05d_logg%01.02f" %(teff,logg)
            mod = ff[mod_name]
            wave_ = mod.data.field('wavelength')
            flux_ = mod.data.field('flux')
            #-- if there is no wavelength range given, we assume that
            #   the whole grid has the same resolution, and the first
            #   wave-array will be used as a template
            if wave is None:
                wave = wave_
                flux = np.ones((len(teffs),len(wave)))
            try:
                flux[i] = flux_
            except:
                flux[i] = np.interp(wave,wave_,flux_)
        ff.close()
        flux_grid = LinearNDInterpolator(np.array([np.log10(teffs),loggs]).T,flux)
    return wave,teffs,loggs,flux,flux_grid

#}

#{ Calibration

def list_calibrators(library='calspec'):
    """
    Print and return the list of calibrators

    @parameter library: name of the library (calspec, ngsl, stelib)
    @type library: str
    @return: list of calibrator names
    @rtype: list of str
    """
    files = config.glob(os.path.join(caldir,library),'*.fits')
    targname = dict(calspec='targetid',ngsl='targname',stelib='object')[library]

    names = []
    for ff in files:
        name = pf.getheader(ff)[targname]
        #star_info = sesame.search(name)
        #if not star_info:
        #    continue
        #else:
        names.append(name)
    return names



def get_calibrator(name='alpha_lyr',version=None,wave_units=None,flux_units=None,library='calspec'):
    """
    Retrieve a calibration SED

    If C{version} is None, get the last version.

    Example usage:

    >>> wave,flux = get_calibrator(name='alpha_lyr')
    >>> wave,flux = get_calibrator(name='alpha_lyr',version='003')

    @param name: calibrator name
    @type name: str
    @param version: version of the calibration file
    @type version: str
    @param wave_units: units of wavelength arrays (default: AA)
    @type wave_units: str (interpretable by C{units.conversions.convert})
    @param flux_units: units of flux arrays (default: erg/s/cm2/AA)
    @type flux_units: str (interpretable by C{units.conversions.convert})
    @return: wavelength and flux arrays of calibrator
    @rtype: (ndarray,ndarray)
    """
    #-- collect calibration files
    files = config.glob(os.path.join(caldir,library),'*.fits')
    targname = dict(calspec='targetid',ngsl='targname',stelib='object')[library]

    calfile = None
    for ff in files:
        #-- check if the name matches with the given one
        fits_file = pf.open(ff)
        header = fits_file[0].header
        if name in ff or name in header[targname]:
            #-- maybe the target is correct, but the 'model version' is not
            if version is not None and version not in ff:
                fits_file.close()
                continue
            #-- extract the wavelengths and flux
            calfile = ff
            if library in ['calspec','ngsl']:
                wave = fits_file[1].data.field('wavelength')
                flux = fits_file[1].data.field('flux')
            elif library in ['stelib']:
                wave,flux = fits.read_spectrum(ff)
            else:
                raise ValueError("Don't know what to do with files from library {}".format(library))
        fits_file.close()

    if calfile is None:
        raise ValueError('Calibrator %s (version=%s) not found'%(name,version))

    if flux_units is not None:
        flux = conversions.convert('erg/s/cm2/AA',flux_units,flux,wave=(wave,'AA'))
    if wave_units is not None:
        wave = conversions.convert('AA',wave_units,wave)


    logger.info('Calibrator %s selected'%(calfile))

    return wave,flux


@memoized
def read_calibrator_info(library='ngsl'):
    filename = config.get_datafile('sedtables/calibrators','{}.ident'.format(library))
    names = []
    fits_files = []
    phot_files = []
    with open(filename,'r') as ff:
        for line in ff.readlines():
            line = line.strip().split(',')
            try:
                fits_file = config.get_datafile('sedtables/calibrators',line[1])
                phot_file = config.get_datafile('sedtables/calibrators',line[2])
            #-- it can happen that there is no photfile for a target
            except IOError:
                continue

            names.append(line[0])
            fits_files.append(fits_file)
            phot_files.append(phot_file)

    return names,fits_files,phot_files


def calibrate():
    """
    Calibrate photometry.

    Not finished!

    ABmag = -2.5 Log F_nu - 48.6 with F_nu in erg/s/cm2/Hz
    Flux computed as 10**(-(meas-mag0)/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    F0 = 3.6307805477010029e-20 erg/s/cm2/Hz

    STmag = -2.5 Log F_lam - 21.10 with F_lam in erg/s/cm2/AA
    Flux computed as 10**(-(meas-mag0)/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    F0 = 3.6307805477010028e-09 erg/s/cm2/AA

    Vegamag = -2.5 Log F_lam - C with F_lam in erg/s/cm2/AA
    Flux computed as 10**(-meas/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    """
    F0ST = 3.6307805477010028e-09
    F0AB = 3.6307805477010029e-20
    #-- get calibrator
    wave,flux = get_calibrator(name='alpha_lyr')
    zp = filters.get_info()

    #-- calculate synthetic fluxes
    syn_flux = synthetic_flux(wave,flux,zp['photband'])
    syn_flux_fnu = synthetic_flux(wave,flux,zp['photband'],units='Fnu')
    Flam0_lit = conversions.nconvert(zp['Flam0_units'],'erg/s/cm2/AA',zp['Flam0'],photband=zp['photband'])
    Fnu0_lit = conversions.nconvert(zp['Fnu0_units'],'erg/s/cm2/Hz',zp['Fnu0'],photband=zp['photband'])

    #-- we have Flam0 but not Fnu0: compute Fnu0
    keep = (zp['Flam0_lit']==1) & (zp['Fnu0_lit']==0)
    Fnu0 = conversions.nconvert(zp['Flam0_units'],'erg/s/cm2/Hz',zp['Flam0'],photband=zp['photband'])
    zp['Fnu0'][keep] = Fnu0[keep]
    zp['Fnu0_units'][keep] = 'erg/s/cm2/Hz'

    #-- we have Fnu0 but not Flam0: compute Flam0
    keep = (zp['Flam0_lit']==0) & (zp['Fnu0_lit']==1)
    Flam0 = conversions.nconvert(zp['Fnu0_units'],'erg/s/cm2/AA',zp['Fnu0'],photband=zp['photband'])

    #   set everything in correct units for convenience:
    Flam0 = conversions.nconvert(zp['Flam0_units'],'erg/s/cm2/AA',zp['Flam0'])
    Fnu0 = conversions.nconvert(zp['Fnu0_units'],'erg/s/cm2/Hz',zp['Fnu0'])

    #-- as a matter of fact, set Flam0 and Fnu for all the stuff for which we
    #   have no literature values
    keep = (zp['Flam0_lit']==0) & (zp['Fnu0_lit']==0)
    zp['Flam0'][keep] = syn_flux[keep]
    zp['Flam0_units'][keep] = 'erg/s/cm2/AA'
    zp['Fnu0'][keep] = syn_flux_fnu[keep]
    zp['Fnu0_units'][keep] = 'erg/s/cm2/Hz'

    keep = np.array(['DENIS' in photb and True or False for photb in zp['photband']])

    #-- we have no Flam0, only ZP vegamags
    keep = (zp['vegamag_lit']==1) & (zp['Flam0_lit']==0)
    zp['Flam0'][keep] = syn_flux[keep]
    zp['Flam0_units'][keep] = 'erg/s/cm2/AA'

    #-- we have no Flam0, no ZP vegamas but STmags
    keep = (zp['STmag_lit']==1) & (zp['Flam0_lit']==0)
    m_vega = 2.5*np.log10(F0ST/syn_flux) + zp['STmag']
    zp['vegamag'][keep] = m_vega[keep]

    #-- we have no Fnu0, no ZP vegamas but ABmags
    keep = (zp['ABmag_lit']==1) & (zp['Flam0_lit']==0)
    F0AB_lam = conversions.convert('erg/s/cm2/Hz','erg/s/cm2/AA',F0AB,photband=zp['photband'])
    m_vega = 2.5*np.log10(F0AB_lam/syn_flux) + zp['ABmag']
    zp['vegamag'][keep] = m_vega[keep]

    #-- set the central wavelengths of the bands
    set_wave = np.isnan(zp['eff_wave'])
    zp['eff_wave'][set_wave] = filters.eff_wave(zp['photband'][set_wave])

    return zp


#}

#{ Synthetic photometry

def synthetic_flux(wave,flux,photbands,units=None):
    """
    Extract flux measurements from a synthetic SED (Fnu or Flambda).

    The fluxes below 4micron are calculated assuming PHOTON-counting detectors
    (e.g. CCDs).

    Flam = int(P_lam * f_lam * lam, dlam) / int(P_lam * lam, dlam)

    When otherwise specified, we assume ENERGY-counting detectors (e.g. bolometers)

    Flam = int(P_lam * f_lam, dlam) / int(P_lam, dlam)

    Where P_lam is the total system dimensionless sensitivity function, which
    is normalised so that the maximum equals 1. Also, f_lam is the SED of the
    object, in units of energy per time per unit area per wavelength.

    The PHOTON-counting part of this routine has been thoroughly checked with
    respect to johnson UBV, geneva and stromgren filters, and only gives offsets
    with respect to the Kurucz integrated files (.geneva and stuff on his websites). These could be
    due to different normalisation.

    You can also readily integrate in Fnu instead of Flambda by suppling a list
    of strings to 'units'. This should have equal length of photbands, and
    should contain the strings 'flambda' and 'fnu' corresponding to each filter.
    In that case, the above formulas reduce to

    Fnu = int(P_nu * f_nu / nu, dnu) / int(P_nu / nu, dnu)

    and

    Fnu = int(P_nu * f_nu, dnu) / int(P_nu, dnu)

    Small note of caution: P_nu is not equal to P_lam according to
    Maiz-Apellaniz, he states that P_lam = P_nu/lambda. But in the definition
    we use above here, it *is* the same!

    The model fluxes should B{always} be given in Flambda (erg/s/cm2/AA). The
    program will convert them to Fnu where needed.

    The output is a list of numbers, equal in length to the 'photband' inputs.
    The units of the output are erg/s/cm2/AA where Flambda was given, and
    erg/s/cm2/Hz where Fnu was given.

    The difference is only marginal for 'blue' bands. For example, integrating
    2MASS in Flambda or Fnu is only different below the 1.1% level:

    >>> wave,flux = get_table(teff=10000,logg=4.0)
    >>> energys = synthetic_flux(wave,flux,['2MASS.J','2MASS.J'],units=['flambda','fnu'])
    >>> e0_conv = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',energys[0],photband='2MASS.J')
    >>> np.abs(energys[1]-e0_conv)/energys[1]<0.012
    True

    But this is not the case for IRAS.F12:

    >>> energys = synthetic_flux(wave,flux,['IRAS.F12','IRAS.F12'],units=['flambda','fnu'])
    >>> e0_conv = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',energys[0],photband='IRAS.F12')
    >>> np.abs(energys[1]-e0_conv)/energys[1]>0.1
    True

    If you have a spectrum in micron vs Jy and want to calculate the synthetic
    fluxes in Jy, a little bit more work is needed to get everything in the
    right units. In the following example, we first generate a constant flux
    spectrum in micron and Jy. Then, we convert flux to erg/s/cm2/AA using the
    wavelengths (this is no approximation) and convert wavelength to angstrom.
    Next, we compute the synthetic fluxes in the IRAS band in Fnu, and finally
    convert the outcome (in erg/s/cm2/Hz) to Jansky.

    >>> wave,flux = np.linspace(0.1,200,10000),np.ones(10000)
    >>> flam = conversions.convert('Jy','erg/s/cm2/AA',flux,wave=(wave,'micron'))
    >>> lam = conversions.convert('micron','AA',wave)
    >>> energys = synthetic_flux(lam,flam,['IRAS.F12','IRAS.F25','IRAS.F60','IRAS.F100'],units=['Fnu','Fnu','Fnu','Fnu'])
    >>> energys = conversions.convert('erg/s/cm2/Hz','Jy',energys)

    You are responsible yourself for having a response curve covering the
    model fluxes!

    Now. let's put this all in practice in a more elaborate example: we want
    to check if the effective wavelength is well defined. To do that we will:

        1. construct a model (black body)
        2. make our own weird, double-shaped filter (CCD-type and BOL-type detector)
        3. compute fluxes in Flambda, and convert to Fnu via the effective wavelength
        4. compute fluxes in Fnu, and compare with step 3.

    In an ideal world, the outcome of step (3) and (4) must be equal:

    Step (1): We construct a black body model.


    WARNING: OPEN.BOL only works in Flambda for now.

    See e.g. Maiz-Apellaniz, 2006.

    @param wave: model wavelengths (angstrom)
    @type wave: ndarray
    @param flux: model fluxes (erg/s/cm2/AA)
    @type flux: ndarray
    @param photbands: list of photometric passbands
    @type photbands: list of str
    @param units: list containing Flambda or Fnu flag (defaults to all Flambda)
    @type units: list of strings or str
    @return: model fluxes (erg/s/cm2/AA or erg/s/cm2/Hz)
    @rtype: ndarray
    """
    if isinstance(units,str):
        units = [units]*len(photbands)
    energys = np.zeros(len(photbands))

    #-- only keep relevant information on filters:
    filter_info = filters.get_info()
    keep = np.searchsorted(filter_info['photband'],photbands)
    filter_info = filter_info[keep]

    for i,photband in enumerate(photbands):
        #if filters.is_color
        waver,transr = filters.get_response(photband)
        #-- make wavelength range a bit bigger, otherwise F25 from IRAS has only
        #   one Kurucz model point in its wavelength range... this is a bit
        #   'ad hoc' but seems to work.
        region = ((waver[0]-0.4*waver[0])<=wave) & (wave<=(2*waver[-1]))
        #-- if we're working in infrared (>4e4A) and the model is not of high
        #   enough resolution (100000 points over wavelength range), interpolate
        #   the model in logscale on to a denser grid (in logscale!)
        if filter_info['eff_wave'][i]>=4e4 and sum(region)<1e5 and sum(region)>1:
            logger.debug('%10s: Interpolating model to integrate over response curve'%(photband))
            wave_ = np.logspace(np.log10(wave[region][0]),np.log10(wave[region][-1]),1e5)
            flux_ = 10**np.interp(np.log10(wave_),np.log10(wave[region]),np.log10(flux[region]),)
        else:
            wave_ = wave[region]
            flux_ = flux[region]
        if not len(wave_):
            energys[i] = np.nan
            continue
        #-- perhaps the entire response curve falls in between model points (happends with
        #   narrowband UV filters), or there's very few model points covering it
        if (np.searchsorted(wave_,waver[-1])-np.searchsorted(wave_,waver[0]))<5:
            wave__ = np.sort(np.hstack([wave_,waver]))
            flux_ = np.interp(wave__,wave_,flux_)
            wave_ = wave__
        #-- interpolate response curve onto model grid
        transr = np.interp(wave_,waver,transr,left=0,right=0)

        #-- integrated flux: different for bolometers and CCDs
        #-- WE WORK IN FLAMBDA
        if units is None or ((units is not None) and (units[i].upper()=='FLAMBDA')):
            if photband=='OPEN.BOL':
                energys[i] = np.trapz(flux_,x=wave_)
            elif filter_info['type'][i]=='BOL':
                energys[i] = np.trapz(flux_*transr,x=wave_)/np.trapz(transr,x=wave_)
            elif filter_info['type'][i]=='CCD':
                energys[i] = np.trapz(flux_*transr*wave_,x=wave_)/np.trapz(transr*wave_,x=wave_)

        #-- we work in FNU
        elif units[i].upper()=='FNU':
            #-- convert wavelengths to frequency, Flambda to Fnu
            freq_ = conversions.convert('AA','Hz',wave_)
            flux_f = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',flux_,wave=(wave_,'AA'))
            #-- sort again!
            sa = np.argsort(freq_)
            transr = transr[sa]
            freq_ = freq_[sa]
            flux_f = flux_f[sa]
            if filter_info['type'][i]=='BOL':
                energys[i] = np.trapz(flux_f*transr,x=freq_)/np.trapz(transr,x=freq_)
            elif filter_info['type'][i]=='CCD':
                energys[i] = np.trapz(flux_f*transr/freq_,x=wave_)/np.trapz(transr/freq_,x=wave_)
        else:
            raise ValueError('units %s not understood'%(units))

    #-- that's it!
    return energys


def synthetic_color(wave,flux,colors,units=None):
    """
    Construct colors from a synthetic SED.

    @param wave: model wavelengths (angstrom)
    @type wave: ndarray
    @param flux: model fluxes (erg/s/cm2/AA)
    @type flux: ndarray
    @param colors: list of photometric passbands
    @type colors: list of str
    @param units: list containing Flambda or Fnu flag (defaults to all Flambda)
    @type units: list of strings or str
    @return: flux ratios or colors
    @rtype: ndarray
    """
    if units is None:
        units = [None for color in colors]

    syn_colors = np.zeros(len(colors))
    for i,(color,unit) in enumerate(zip(colors,units)):
        #-- retrieve the passbands necessary to construct the color, and the
        #   function that defines the color
        photbands,color_func = filters.make_color(color)
        #-- compute the synthetic fluxes to construct the color
        fluxes = synthetic_flux(wave,flux,photbands,units=unit)
        #-- construct the color
        syn_colors[i] = color_func(*list(fluxes))

    return syn_colors



def luminosity(wave,flux,radius=1.):
    """
    Calculate the bolometric luminosity of a model SED.

    Flux should be in cgs per unit wavelength (same unit as wave).
    The latter is integrated out, so it is of no importance. After integration,
    flux, should have units erg/s/cm2.

    Returned luminosity is in solar units.

    If you give radius=1 and want to correct afterwards, multiply the obtained
    Labs with radius**2.

    @param wave: model wavelengths
    @type wave: ndarray
    @param flux: model fluxes (Flam)
    @type flux: ndarray
    @param radius: stellar radius in solar units
    @type radius: float
    @return: total bolometric luminosity
    @rtype: float
    """
    Lint = np.trapz(flux,x=wave)
    Labs = Lint*4*np.pi/constants.Lsol_cgs*(radius*constants.Rsol_cgs)**2
    return Labs

#}

@memoized
def _get_itable_markers(photbands,
                    teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf),
                    ebvrange=(-np.inf,np.inf),zrange=(-np.inf,np.inf),
                    include_Labs=True,clear_memory=True,**kwargs):
    """
    Get a list of markers to more easily retrieve integrated fluxes.
    """
    if clear_memory:
        clear_memoization(keys=['ivs.sed.model'])
    gridfiles = get_file(z='*',integrated=True,**kwargs)
    if isinstance(gridfiles,str):
        gridfiles = [gridfiles]
    #-- sort gridfiles per metallicity
    metals_sa = np.argsort([pf.getheader(ff,1)['z'] for ff in gridfiles])
    gridfiles = np.array(gridfiles)[metals_sa]
    flux = []
    gridpnts = []
    grid_z = []
    markers = []

    #-- collect information
    for gridfile in gridfiles:
        ff = pf.open(gridfile)
        ext = ff[1]
        z = ff[1].header['z']
        if z<zrange[0] or zrange[1]<z:
            continue

        teffs = ext.data.field('teff')
        loggs = ext.data.field('logg')
        ebvs = ext.data.field('ebv')
        keep = (ebvrange[0]<=ebvs) & (ebvs<=ebvrange[1])

        #-- for some reason, the Kurucz grid has a lonely point at Teff=14000,logg=2
        #   which messes up our interpolation
        #correct = (teffs==14000) & (loggs==2.0)
        #teffs[correct] = 12000

        teffs,loggs,ebvs = teffs[keep],loggs[keep],ebvs[keep]
        grid_teffs = np.sort(list(set(teffs)))
        grid_loggs = np.sort(list(set(loggs)))
        grid_ebvs = np.sort(list(set(ebvs)))
        grid_z.append(z)

        #-- we construct an array representing the teff-logg-ebv-z content, but
        #   in one number: 5000040031500 means:
        #   T=50000,logg=4.0,E(B-V)=0.31 and Z = 0.00
        # Note that Z is Z+5 so that we avoid minus signs...
        markers.append(np.zeros(len(teffs)))
        gridpnts.append(np.zeros((len(teffs),4)))

        for i,(it,il,ie) in enumerate(zip(teffs,loggs,ebvs)):
            markers[-1][i] = float('%3d%05d%03d%03d'%(int(round((z+5)*100)),int(round(it)),int(round(il*100)),int(round(ie*100))))
            gridpnts[-1][i]= it,il,ie,z

        flux.append(_get_flux_from_table(ext,photbands,include_Labs=include_Labs))
        ff.close()

    flux = np.vstack(flux)
    markers = np.hstack(markers)
    gridpnts = np.vstack(gridpnts)
    grid_z = np.sort(grid_z)

    return np.array(markers),(grid_teffs,grid_loggs,grid_ebvs,grid_z),gridpnts,flux


@memoized
def _get_pix_grid(photbands,
                    teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf),
                    ebvrange=(-np.inf,np.inf),zrange=(-np.inf,np.inf),
                    rvrange=(-np.inf,np.inf),vradrange=(-np.inf,np.inf),
                    include_Labs=True,clear_memory=True,
                    variables=['teff','logg','ebv','z','rv','vrad'],**kwargs):
    """
    Prepare the pixalted grid.

    In principle, it should be possible to return any number of free parameters
    here. I'm thinking about:

        teff, logg, ebv, z, Rv, vrad.
    """
    if clear_memory:
        clear_memoization(keys=['ivs.sed.model'])

    #-- remove Rv and z from the grid keywords
    trash = kwargs.pop('Rv', '*')
    trash = kwargs.pop('z', '*')
    gridfiles = get_file(z='*',Rv='*',integrated=True,**kwargs)
    if isinstance(gridfiles,str):
        gridfiles = [gridfiles]
    flux = []
    grid_pars = []
    grid_names = np.array(variables)
    #-- collect information from all the grid files
    for gridfile in gridfiles:
        with pf.open(gridfile) as ff:
            #-- make an alias for further reference
            ext = ff[1]
            #-- we already cut the grid here, in order not to take too much memory
            keep = np.ones(len(ext.data),bool)
            for name in variables:
                #-- we need to be carefull for rounding errors
                low,high = locals()[name+'range']
                in_range = (low<=ext.data.field(name)) & (ext.data.field(name)<=high)
                on_edge  = np.allclose(ext.data.field(name),low) | np.allclose(ext.data.field(name),high)
                #on_edge_low = np.less_equal(np.abs(ext.data.field(name)-low),1e-8 + 1e-5*np.abs(low))
                #on_edge_high = np.less_equal(np.abs(ext.data.field(name)-high),1e-8 + 1e-5*np.abs(high))
                #keep_this = (in_range | on_edge_low | on_edge_high)
                #if not sum(keep_this):
                    #logger.warning("_get_pix_grid: No selection done in axis {}".format(name))
                    #continue
                #keep = keep & keep_this
                keep = keep & (in_range | on_edge)
            partial_grid = np.vstack([ext.data.field(name)[keep] for name in variables])
            if sum(keep):
                grid_pars.append(partial_grid)
                #-- the flux grid:
                flux.append(_get_flux_from_table(ext,photbands,include_Labs=include_Labs)[keep])
    #-- make the entire grid: it consists of fluxes and grid parameters
    flux = np.vstack(flux)
    grid_pars = np.hstack(grid_pars)
    #-- this is also the place to put some stuff in logarithmic scale if
    #   this is needed
    #grid_pars[0] = np.log10(grid_pars[0])
    flux = np.log10(flux)

    #-- don't take axes into account if it has only one value
    keep = np.ones(len(grid_names),bool)
    for i in range(len(grid_names)):
        if np.all(grid_pars[i]==grid_pars[i][0]):
            keep[i] = False
    grid_pars = grid_pars[keep]

    #-- we need to know what variable parameters we have in the grid
    grid_names = grid_names[keep]

    #-- create the pixeltype grid
    axis_values, pixelgrid = interpol.create_pixeltypegrid(grid_pars,flux.T)
    return axis_values,grid_pars.T,pixelgrid,grid_names




def _get_flux_from_table(fits_ext,photbands,index=None,include_Labs=True):
    """
    Retrieve flux and flux ratios from an integrated SED table.

    @param fits_ext: fits extension containing integrated flux
    @type fits_ext: FITS extension
    @param photbands: list of photometric passbands
    @type photbands: list of str
    @param index: slice or index of rows to retrieve
    @type index: slice or integer
    @return: fluxes or flux ratios
    #@rtype: list
    """
    if index is None:
        index = slice(None) #-- full range
    fluxes = []
    for photband in photbands:
        try:
            if not filters.is_color(photband):
                fluxes.append(fits_ext.data.field(photband)[index])
            else:
                system,color = photband.split('.')
                if '-' in color:
                    band0,band1 = color.split('-')
                    fluxes.append(fits_ext.data.field('%s.%s'%(system,band0))[index]/fits_ext.data.field('%s.%s'%(system,band1))[index])
                elif color=='M1':
                    fv = fits_ext.data.field('STROMGREN.V')[index]
                    fy = fits_ext.data.field('STROMGREN.Y')[index]
                    fb = fits_ext.data.field('STROMGREN.B')[index]
                    fluxes.append(fv*fy/fb**2)
                elif color=='C1':
                    fu = fits_ext.data.field('STROMGREN.U')[index]
                    fv = fits_ext.data.field('STROMGREN.V')[index]
                    fb = fits_ext.data.field('STROMGREN.B')[index]
                    fluxes.append(fu*fb/fv**2)
        except KeyError:
            logger.warning('Passband %s missing from table'%(photband))
            fluxes.append(np.nan*np.ones(len(fits_ext.data)))
    #-- possibly include absolute luminosity
    if include_Labs:
        fluxes.append(fits_ext.data.field("Labs")[index])
    fluxes = np.array(fluxes).T
    if index is not None:
        fluxes = fluxes
    return fluxes



if __name__=="__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()
