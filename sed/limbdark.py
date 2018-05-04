# -*- coding: utf-8 -*-
"""
Interface to the limb-darkening library.

Section 1. Basic interface
==========================

Retrieve limb darkening coeffcients and goodness of fit parameters:

>>> coefs,fac,ssr,idiff = get_coeffs(teff=10000,logg=4.0,'JOHNSON.V',law='claret')

When you need this for a large grid of parameters, consider using L{get_itable}.
The available laws are C{claret}, C{linear}, C{logarithmic}, C{quadratic} and
C{power}.

Retrieve a normalised passband-integrated limb darkening profile via:

>>> mu,intensities = get_limbdarkening(teff=10000,logg=4.0,photbands=['JOHNSON.V'],normalised=True)

and fit a particular law via:

>>> coeffs,ssr,idiff = fit_law(mu,intensities[:,0],law='claret',fitmethod='equidist_mu_leastsq')
>>> print(coeffs)
[ 2.00943266 -3.71261279  4.25469623 -1.70466934]

You can evaluate this law on a particular angle-sampling grid:

>>> mus = np.linspace(0,1,1000)
>>> mu_law = ld_claret(mus,coeffs)

Or on a disk-sampling grid:

>>> rs = np.linspace(0,1,1000)
>>> r_law = ld_claret(_mu(rs),coeffs)

Make a simple plot:

>>> p = pl.figure()
>>> p,q = pl.subplot(121),pl.title('Angle-sampled')
>>> p = pl.plot(mu,intensities,'ko-')
>>> p = pl.plot(mus,mu_law,'r-',lw=2)
>>> p,q = pl.subplot(122),pl.title('Disk-sampled')
>>> p = pl.plot(_r(mu),intensities,'ko-')
>>> p = pl.plot(rs,r_law,'r-',lw=2)

]include figure]]ivs_limbdark_basic.png]

Section 2. Fit comparison
=========================

You can fit a limb darkening law in various ways: using different laws, different
optimizers and different x-coordinates (limb angle mu or radius). In the figure
below, you can see the influence of some of these choices:

>>> p = pl.figure()

>>> laws = ['claret','linear','logarithmic','quadratic','power']
>>> fitmethods = ['equidist_mu_leastsq','equidist_r_leastsq','leastsq']

>>> r,rs = _r(mu),_r(mus)
>>> integ_mu = np.trapz(intensities[:,0],x=mu)
>>> integ_r = np.trapz(intensities[:,0],x=r)
>>> for i,fitmethod in enumerate(fitmethods):
...     p = pl.subplot(3,2,2*i+1)
...     p = pl.title(fitmethod)
...     p = pl.plot(mu,intensities[:,0],'ko-')
...     p = pl.xlabel('$\mu$ = cos(limb angle)')
...     p = pl.ylabel('Normalised intensity')
...     p = pl.subplot(3,2,2*i+2)
...     p = pl.title(fitmethod)
...     p = pl.plot(r,intensities[:,0],'ko-')
...     p = pl.xlabel('r')
...     p = pl.ylabel('Normalised intensity')
...     for law in laws:
...         cl,ssr,idiff = fit_law(mu,intensities[:,0],law=law,fitmethod=fitmethod)
...         print("{:12s} ({:19s}): {}".format(law,fitmethod,cl))
...         Imus = globals()['ld_{}'.format(law)](mus,cl)
...         p = pl.subplot(3,2,2*i+1)
...         idiff = (np.trapz(Imus,x=mus)-integ_mu)/integ_mu
...         p = pl.plot(mus,Imus,'-',lw=2,label="{} ({:.3f},{:.3f}%)".format(law,ssr,idiff*100))
...         p = pl.subplot(3,2,2*i+2)
...         idiff = (np.trapz(Imus,x=rs)-integ_r)/integ_r
...         p = pl.plot(rs,Imus,'-',lw=2,label="{} ({:.3f},{:.3}%)".format(law,ssr,idiff*100))
...     p = pl.subplot(3,2,2*i+1)
...     leg = pl.legend(loc='best',fancybox=True,prop=dict(size='small'))
...     leg.get_frame().set_alpha(0.5)
...     p = pl.subplot(3,2,2*i+2)
...     leg = pl.legend(loc='best',fancybox=True,prop=dict(size='small'))
...     leg.get_frame().set_alpha(0.5)
claret       (equidist_mu_leastsq): [ 2.00943266 -3.71261279  4.25469623 -1.70466934]
linear       (equidist_mu_leastsq): [ 0.48655636]
divide by zero encountered in log
invalid value encountered in multiply
logarithmic  (equidist_mu_leastsq): [ 0.6  0. ]
quadratic    (equidist_mu_leastsq): [ 0.21209044  0.36591797]
power        (equidist_mu_leastsq): [ 0.29025864]
claret       (equidist_r_leastsq ): [ 1.99805089 -3.04504952  3.02780048 -1.09105603]
linear       (equidist_r_leastsq ): [ 0.42641628]
logarithmic  (equidist_r_leastsq ): [ 0.6  0. ]
quadratic    (equidist_r_leastsq ): [ 0.28650391  0.24476527]
power        (equidist_r_leastsq ): [ 0.31195553]
claret       (leastsq            ): [  3.49169665  -8.96767096  11.16873314  -4.72869922]
linear       (leastsq            ): [ 0.57180245]
logarithmic  (leastsq            ): [ 0.6  0. ]
quadratic    (leastsq            ): [ 0.00717934  0.65416204]
power        (leastsq            ): [ 0.27032565]

]include figure]]ivs_limbdark_fitcomp.png]

Author: Pieter Degroote, with thanks to Steven Bloemen.
"""
import logging
import os
import itertools
import astropy.io.fits as pf
import numpy as np
from scipy.optimize import leastsq,fmin
from scipy.interpolate import splrep, splev
from scipy.interpolate import LinearNDInterpolator
# from Scientific.Functions.Interpolation import InterpolatingFunction
from ivs.aux import loggers
from ivs.sed import reddening
from ivs.sed import model
from ivs.spectra import tools
from ivs.units import constants
from ivs.aux.decorators import memoized,clear_memoization
from ivs import config

logger = logging.getLogger("SED.LIMBDARK")
logger.addHandler(loggers.NullHandler())

#-- default values for grids
defaults = dict(grid='kurucz',odfnew=False,z=+0.0,vturb=2,
                alpha=False,nover=False)
#-- relative location of the grids
basedir = 'ldtables'

#{ Interface to the library

def set_defaults(**kwargs):
    """
    Set defaults of this module

    If you give no keyword arguments, the default values will be reset.
    """
    clear_memoization(keys=['ivs.sed.ld'])
    if not kwargs:
        kwargs = dict(grid='kurucz',odfnew=True,z=+0.0,vturb=2,
                alpha=False,nover=False,                  # KURUCZ
                He=97,                                    # WD
                t=1.0,a=0.0,c=0.5,m=1.0,co=1.05)          # MARCS and COMARCS

    for key in kwargs:
        if key in defaults:
            defaults[key] = kwargs[key]
            logger.info('Set %s to %s'%(key,kwargs[key]))



def get_gridnames():
    """
    Return a list of available grid names.

    @return: list of grid names
    @rtype: list of str
    """
    return ['kurucz']


def get_grid_dimensions(**kwargs):
    """
    Retrieve the possible effective temperatures and gravities from a grid.

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
    return np.array(teffs),np.array(loggs)


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

    #-- run over teff and logg, and interpolate the models onto the supplied
    #   wavelength range
    gridfile = get_file(**kwargs)

    if wave is not None:
        table = np.zeros((len(teffs),len(wave),18))
    ff = pf.open(gridfile)
    for i,(teff,logg) in enumerate(zip(teffs,loggs)):
        mod_name = "T%05d_logg%01.02f" %(teff,logg)
        mod = ff[mod_name]
        imu = np.array(mod.columns.names[1:], float)
        itable = np.array(mod.data.tolist())[:,1:]
        iwave = mod.data.field('wavelength')
        #-- if there is no wavelength range given, we assume that
        #   the whole grid has the same resolution, and the first
        #   wave-array will be used as a template
        if wave is None:
            wave = iwave
            table = np.zeros((len(teffs),len(wave),len(imu)))
        try:
            table[i] = itable
        except:
            for j,jimu in enumerate(imu):
                table[i,:,j] = np.interp(wave,iwave,itable[:,j])
    ff.close()
    table_grid = LinearNDInterpolator(np.array([np.log10(teffs),loggs]).T,table)
    return imu,wave,teffs,loggs,table,table_grid



def get_file(integrated=False,**kwargs):
    """
    Retrieve the filename containing the specified SED grid.

    The keyword arguments are specific to the kind of grid you're using.

    Basic keywords are 'grid' for the name of the grid, and 'z' for metallicity.
    For other keywords, see the source code.

    Available grids and example keywords:
        - grid='kurucz93':
                    * metallicity (z): m01 is -0.1 log metal abundance relative to solar (solar abundances from Anders and Grevesse 1989)
                    * metallicity (z): p01 is +0.1 log metal abundance relative to solar (solar abundances from Anders and Grevesse 1989)
                    * alpha enhancement (alpha): True means alpha enhanced (+0.4)
                    * turbulent velocity (vturb): vturb in km/s
                    * nover= True means no overshoot
                    * odfnew=True means no overshoot but with better opacities and abundances

    @param integrated: choose integrated version of the grid
    @type integrated: boolean
    @keyword grid: gridname (default Kurucz)
    @type grid: str
    @return: gridfile
    @rtype: str
    """
    #-- possibly you give a filename
    grid = kwargs.get('grid',defaults['grid'])
    if os.path.isfile(grid):
        return grid

    #-- general
    z = kwargs.get('z',defaults['z'])
    #-- only for Kurucz
    vturb = int(kwargs.get('vturb',defaults['vturb']))
    odfnew = kwargs.get('odfnew',defaults['odfnew'])
    alpha = kwargs.get('alpha',defaults['alpha'])
    nover = kwargs.get('nover',defaults['nover'])

    #-- figure out what grid to use
    if grid=='kurucz':
        if not isinstance(z,str): z = '%.1f'%(z)
        if not isinstance(vturb,str): vturb = '%d'%(vturb)
        if not alpha and not nover and not odfnew:
            basename = 'kurucz93_z%s_k%s_ld.fits'%(z,vturb)
        elif alpha and odfnew:
            basename = 'kurucz93_z%s_ak%sodfnew_ld.fits'%(z,vturb)
        elif odfnew:
            basename = 'kurucz93_z%s_k%sodfnew_ld.fits'%(z,vturb)
        elif nover:
            basename = 'kurucz93_z%s_k%snover_ld.fits'%(z,vturb)
    else:
        basename = grid

    #-- retrieve the absolute path of the file and check if it exists:
    if not '*' in basename:
        if integrated:
            grid = config.get_datafile(basedir,'i'+basename)
        else:
            grid = config.get_datafile(basedir,basename)
    #-- we could also ask for a list of files, when wildcards are given:
    else:
        grid = config.glob(basedir,'i'+basename)
        if integrated:
            grid = config.glob(basedir,'i'+basename)
        else:
            grid = config.glob(basedir,basename)
    logger.debug('Selected %s'%(grid))

    return grid



def get_table(teff=None,logg=None,ebv=None,vrad=None,star=None,
              wave_units='AA',flux_units='erg/cm2/s/AA/sr',**kwargs):
    """
    Retrieve the specific intensity of a model atmosphere.

    ebv is reddening
    vrad is radial velocity: positive is redshift, negative is blueshift (km/s!)

    extra kwargs are for reddening

    You get limb angles, wavelength and a table. The shape of the table is
    (N_wave,N_mu).

    WARNING: wave and flux units cannot be specificed for the moment.

    >>> mu,wave,table = get_table(10000,4.0)

    >>> p = pl.figure()
    >>> ax1 = pl.subplot(221)
    >>> p = pl.title('E(B-V)=0, vrad=0')
    >>> ax = pl.gca()
    >>> ax.set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0,1,len(mu))])
    >>> p = pl.loglog(wave,table)
    >>> p = pl.xlim(700,15000)
    >>> p = pl.ylim(1e3,1e8)

    >>> mu,wave,table = get_table(10000,4.0,ebv=0.5)

    >>> p = pl.subplot(222,sharex=ax1,sharey=ax1)
    >>> p = pl.title('E(B-V)=0.5, vrad=0')
    >>> ax = pl.gca()
    >>> ax.set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0,1,len(mu))])
    >>> p = pl.loglog(wave,table)
    >>> p = pl.xlim(700,15000)
    >>> p = pl.ylim(1e3,1e8)

    >>> mu,wave,table = get_table(10000,4.0,vrad=-1000.)

    >>> p = pl.subplot(223,sharex=ax1,sharey=ax1)
    >>> p = pl.title('E(B-V)=0., vrad=-10000.')
    >>> ax = pl.gca()
    >>> ax.set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0,1,len(mu))])
    >>> p = pl.loglog(wave,table)
    >>> p = pl.xlim(700,15000)
    >>> p = pl.ylim(1e3,1e8)

    >>> mu,wave,table = get_table(10000,4.0,vrad=-1000.,ebv=0.5)

    >>> p = pl.subplot(224,sharex=ax1,sharey=ax1)
    >>> p = pl.title('E(B-V)=0.5, vrad=-10000.')
    >>> ax = pl.gca()
    >>> ax.set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0,1,len(mu))])
    >>> p = pl.loglog(wave,table)
    >>> p = pl.xlim(700,15000)
    >>> p = pl.ylim(1e3,1e8)

    >>> mu,wave,table = get_table(10050,4.12)

    >>> p = pl.figure()
    >>> pl.gca().set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0,1,len(mu))])
    >>> p = pl.loglog(wave,table)


    @param teff: effective temperature (K)
    @type teff: float
    @param logg: log surface gravity (cgs, dex)
    @type logg: float
    @param ebv: reddening (mag)
    @type ebv: float
    @param vrad: radial velocity (for doppler shifting) (km/s)
    @type vrad: float
    @return: mu angles, wavelengths, table (Nwave x Nmu)
    @rtype: array, array, array
    """
    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)

    #-- read the file:
    ff = pf.open(gridfile)

    teff = float(teff)
    logg = float(logg)

    #-- if we have a grid model, no need for interpolation
    try:
        #-- extenstion name as in fits files prepared by Steven
        mod_name = "T%05d_logg%01.02f" %(teff,logg)
        mod = ff[mod_name]
        mu = np.array(mod.columns.names[1:], float)
        table = np.array(mod.data.tolist())[:,1:]
        wave = mod.data.field('wavelength')
        logger.debug('Model LD taken directly from file (%s)'%(os.path.basename(gridfile)))
    except KeyError:
        mu,wave,teffs,loggs,flux,flux_grid = get_grid_mesh(**kwargs)
        logger.debug('Model LD interpolated from grid %s (%s)'%(os.path.basename(gridfile),kwargs))
        wave = wave + 0.
        table = flux_grid(np.log10(teff),logg) + 0.

    ff.close()

    #-- velocity shift if necessary
    if vrad is not None and vrad!=0:
        cc = constants.cc/1000. #speed of light in km/s
        for i in range(len(mu)):
            flux_shift = tools.doppler_shift(wave,vrad,flux=table[:,i])
            table[:,i] = flux_shift - 5.*vrad/cc*table[:,i]

    #-- redden if necessary
    if ebv is not None and ebv>0:
        for i in range(len(mu)):
            table[:,i] = reddening.redden(table[:,i],wave=wave,ebv=ebv,rtype='flux',**kwargs)

    #-- that's it!
    return mu,wave,table








def get_itable(teff=None,logg=None,theta=None,mu=1,photbands=None,absolute=False,**kwargs):
    """
    mu=1 is center of disk
    """
    if theta is not None:
        mu = np.cos(theta)

    try:
        out = get_ld_grid(photbands,integrated=True,**kwargs)(teff,logg)
    except ValueError:
        print('Used teff and logg',teff,logg)
        raise
    a1x_,a2x_,a3x_,a4x_, I_x1 = out.reshape((len(photbands),5)).T
    Imu = ld_eval(mu,[a1x_,a2x_,a3x_,a4x_])
    if absolute:
        return Imu*I_x1
    else:
        return Imu

@memoized
def _get_itable_markers(photband,gridfile,
                    teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf)):
    """
    Get a list of markers to more easily retrieve integrated fluxes.
    """
    clear_memoization()

    ff = pf.open(gridfile)
    ext = ff[photband]
    columns = ext.columns.names
    names = ['teff','logg']
    grid = [ext.data.field('teff'),
            ext.data.field('logg')]
    if 'ebv' in columns:
        names.append('ebv')
        grid.append(ext.data.field('ebv'))
    if 'vrad' in columns:
        names.append('vrad')
        grid.append(ext.data.field('vrad'))

    grid_axes = [np.sort(list(set(i))) for i in grid]

    #-- we construct an array representing the teff-logg-ebv content, but
    #   in one number: 50000400 means:
    #   T=50000,logg=4.0
    N = len(grid[0])
    markers = np.zeros(N,float)
    gridpnts = np.zeros((N,len(grid)))
    pars = np.zeros((N,5))

    for i,entries in enumerate(zip(*grid)):
        markers[i] = encode(**{key:entry for key,entry in zip(names,entries)})
        gridpnts[i]= entries
        pars[i] = list(ext.data[i][-5:])
    ff.close()
    sa = np.argsort(markers)
    print('read in gridfile',gridfile)
    pars[:,-1] = np.log10(pars[:,-1])
    return markers[sa],grid_axes,gridpnts[sa],pars[sa]




def get_limbdarkening(teff=None,logg=None,ebv=None,vrad=None,z=None,photbands=None,normalised=False,**kwargs):
    """
    Retrieve a limb darkening law for a specific star and specific bandpass.

    Possibility to include reddening via EB-V parameter. If not given,
    reddening will not be performed...

    You choose your own reddening function.

    See e.g. Heyrovsky et al., 2007

    If you specify one angle (mu in radians), it will take the closest match
    from the grid.

    Mu = cos(theta) where theta is the angle between the surface and the line
    of sight. mu=1 means theta=0 means center of the star.

    Example usage:

    >>> teff,logg = 5000,4.5
    >>> mu,intensities = get_limbdarkening(teff=teff,logg=logg,photbands=['JOHNSON.V'])

    @keyword teff: effective temperature
    @type teff: float
    @keyword logg: logarithmic gravity (cgs)
    @type logg: float
    @keyword system: bandpass system
    @type system: string
    @keyword photbands: bandpass filters
    @type photbands: list of strings
    @keyword ebv: reddening coefficient
    @type ebv: float
    @keyword vrad: radial velocity (+ is redshift, - is blueshift)
    @type vrad: float
    @keyword mu: specificy specific angle
    @type mu: float
    """
    if z is None:
        z = defaults['z']
    #-- retrieve model atmosphere for a given teff and logg
    mus,wave,table = get_table(teff=teff,logg=logg,ebv=ebv,vrad=vrad,z=z,**kwargs)
    #-- compute intensity over the stellar disk, and normalise
    intensities = np.zeros((len(mus),len(photbands)))
    for i in range(len(mus)):
        intensities[i] = model.synthetic_flux(wave,table[:,i],photbands)
    if normalised:
        intensities/= intensities.max(axis=0)
    #-- or compute the intensity only for one angle:
    logger.info('Calculated LD')
    return mus,intensities

def get_coeffs(teff,logg,photband,ebv=None,vrad=None,law='claret',fitmethod='equidist_r_leastsq'):
    """
    Retrieve limb darkening coefficients on the fly.

    This is particularly useful if you only need a few and speed is not really
    a problem (although this function is nearly instantaneous, looping over it
    will expose that it's actually pretty slow).

    @keyword teff: effective temperature
    @type teff: float
    @keyword logg: logarithmic gravity (cgs)
    @type logg: float
    @keyword photbands: bandpass filters
    @type photbands: list of strings
    @keyword ebv: reddening coefficient
    @type ebv: float
    @keyword vrad: radial velocity (+ is redshift, - is blueshift)
    @type vrad: float
    """
    mu,intensities = get_limbdarkening(teff=teff,logg=logg,ebv=ebv,vrad=vrad,photbands=[photband])
    norm_factor = intensities.max(axis=0)
    coeffs,ssr,idiff = fit_law(mu,intensities[:,0]/norm_factor,law=law,fitmethod=fitmethod)
    return coeffs,norm_factor,ssr,idiff




def get_ld_grid_dimensions(**kwargs):
    """
    Returns the gridpoints of the limbdarkening grid (not unique values).
    """
    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)
    #-- the teff and logg range is the same for every extension, so just
    #   take the first one
    ff = pf.open(gridfile)
    teff_,logg_ = ff[1].data.field('Teff'),ff[1].data.field('logg')
    ff.close()
    return teff_,logg_

def intensity_moment(coeffs,ll=0,**kwargs):
    """
    Calculate the intensity moment (see Townsend 2003, eq. 39).

    Test analytical versus numerical implementation:

    #>>> photband = 'JOHNSON.V'
    #>>> l,logg = 2.,4.0
    #>>> gridfile = 'tables/ikurucz93_z0.0_k2_ld.fits'
    #>>> mu = np.linspace(0,1,100000)
    #>>> P_l = legendre(l)(mu)
    #>>> check = []
    #>>> for teff in np.linspace(9000,12000,19):
    #...    a1x,a2x,a3x,a4x,I_x1 = get_itable(gridfile,teff=teff,logg=logg,mu=1,photband=photband,evaluate=False)
    #...    coeffs = [a1x,a2x,a3x,a4x,I_x1]
    #...    Ix_mu = ld_claret(mu,coeffs)
    #...    Ilx1 = I_x1*np.trapz(P_l*mu*Ix_mu,x=mu)
    #...    a0x = 1 - a1x - a2x - a3x - a4x
    #...    limb_coeffs = [a0x,a1x,a2x,a3x,a4x]
    #...    Ilx2 = 0
    #...    for r,ar in enumerate(limb_coeffs):
    #...        s = 1 + r/2.
    #...        myIls = _I_ls(l,s)
    #...        Ilx2 += ar*myIls
    #...    Ilx2 = I_x1*Ilx2
    #...    Ilx3 = intensity_moment(teff,logg,photband,coeffs)
    #...    check.append(abs(Ilx1-Ilx2)<0.1)
    #...    check.append(abs(Ilx1-Ilx3)<0.1)
    #>>> np.all(np.array(check))
    #True


    @parameter teff: effecitve temperature (K)
    @type teff: float
    @parameter logg: log of surface gravity (cgs)
    @type logg: float
    @parameter ll: degree of the mode
    @type ll: float
    @parameter photband: photometric passbands
    @type photband: list of strings ( or iterable container)
    @keyword full_output: if True, returns intensity, coefficients and integrals
                         separately
    @type full_output: boolean
    @return: intensity moment
    """
    full_output = kwargs.get('full_output',False)
    #-- notation of Townsend 2002: coeffs in code are hat coeffs in the paper
    #   (for i=1..4 they are the same)
    #-- get the LD coefficients at the given temperature and logg
    a1x_,a2x_,a3x_,a4x_, I_x1 = coeffs
    a0x_ = 1 - a1x_ - a2x_ - a3x_ - a4x_
    limb_coeffs = np.array([a0x_,a1x_,a2x_,a3x_,a4x_])

    #-- compute intensity moment via helper coefficients
    int_moms = np.array([_I_ls(ll,1 + r/2.) for r in range(0,5,1)])
    #int_moms = np.outer(int_moms,np.ones(1))
    I_lx = I_x1 * (limb_coeffs * int_moms).sum(axis=0)
    #-- return value depends on keyword
    if full_output:
        return I_x1,limb_coeffs,int_moms
    else:
        return I_lx



#}
#{ Limbdarkening laws

def ld_eval(mu,coeffs):
    """
    Evaluate Claret's limb darkening law.
    """
    return ld_claret(mu,coeffs)

def ld_claret(mu,coeffs):
    """
    Claret's limb darkening law.
    """
    a1,a2,a3,a4 = coeffs
    Imu = 1-a1*(1-mu**0.5)-a2*(1-mu)-a3*(1-mu**1.5)-a4*(1-mu**2.)
    return Imu

def ld_linear(mu,coeffs):
    """
    Linear or linear cosine law
    """
    return 1-coeffs[0]*(1-mu)

def ld_nonlinear(mu,coeffs):
    """
    Nonlinear or logarithmic law
    """
    return 1-coeffs[0]*(1-mu)-coeffs[1]*mu*np.log(mu)

def ld_logarithmic(mu,coeffs):
    """
    Nonlinear or logarithmic law
    """
    return ld_nonlinear(mu,coeffs)

def ld_quadratic(mu,coeffs):
    """
    Quadratic law
    """
    return 1-coeffs[0]*(1-mu)-coeffs[1]*(1-mu)**2.0

def ld_uniform(mu,coeffs):
    """
    Uniform law.
    """
    return 1.

def ld_power(mu,coeffs):
    """
    Power law.
    """
    return mu**coeffs[0]

#}

#{ Fitting routines (thanks to Steven Bloemen)

def fit_law(mu,Imu,law='claret',fitmethod='equidist_r_leastsq'):
    """
    Fit an LD law to a sampled set of limb angles/intensities.

    In my (Pieter) experience, C{fitmethod='equidist_r_leastsq' seems
    appropriate for the Kurucz models.

    Make sure the intensities are normalised!

    >>> mu,intensities = get_limbdarkening(teff=10000,logg=4.0,photbands=['JOHNSON.V'],normalised=True)
    >>> p = pl.figure()
    >>> p = pl.plot(mu,intensities[:,0],'ko-')
    >>> cl,ssr,rdiff = fit_law(mu,intensities[:,0],law='claret')

    >>> mus = np.linspace(0,1,1000)
    >>> Imus = ld_claret(mus,cl)
    >>> p = pl.plot(mus,Imus,'r-',lw=2)



    @return: coefficients, sum of squared residuals,relative flux difference between prediction and model integrated intensity
    @rtype: array, float, float
    """
    #-- prepare array for coefficients and set the initial guess
    Ncoeffs = dict(claret=4,linear=1,nonlinear=2,logarithmic=2,quadratic=2,
                   power=1)
    c0 = np.zeros(Ncoeffs[law])
    c0[0] = 0.6
    #-- do the fitting
    if fitmethod=='leastsq':
        (csol, ierr)  = leastsq(ldres_leastsq, c0, args=(mu,Imu,law))
    elif fitmethod=='fmin':
        csol  = fmin(ldres_fmin, c0, maxiter=1000, maxfun=2000,args=(mu,Imu,law),disp=0)
    elif fitmethod=='equidist_mu_leastsq':
        mu_order = np.argsort(mu)
        tck = splrep(mu[mu_order],Imu[mu_order],s=0.0, k=2)
        mu_spl = np.linspace(mu[mu_order][0],1,5000)
        Imu_spl = splev(mu_spl,tck,der=0)
        (csol, ierr)  = leastsq(ldres_leastsq, c0, args=(mu_spl,Imu_spl,law))
    elif fitmethod=='equidist_r_leastsq':
        mu_order = np.argsort(mu)
        tck = splrep(mu[mu_order],Imu[mu_order],s=0., k=2)
        r_spl = np.linspace(mu[mu_order][0],1,5000)
        mu_spl = np.sqrt(1-r_spl**2)
        Imu_spl = splev(mu_spl,tck,der=0)
        (csol,ierr)  = leastsq(ldres_leastsq, c0, args=(mu_spl,Imu_spl,law))
    elif fitmethod=='equidist_mu_fmin':
        mu_order = np.argsort(mu)
        tck = splrep(mu[mu_order],Imu[mu_order],k=2, s=0.0)
        mu_spl = np.linspace(mu[mu_order][0],1,5000)
        Imu_spl = splev(mu_spl,tck,der=0)
        csol  = fmin(ldres_fmin, c0, maxiter=1000, maxfun=2000,args=(mu_spl,Imu_spl,law),disp=0)
    elif fitmethod=='equidist_r_fmin':
        mu_order = np.argsort(mu)
        tck = splrep(mu[mu_order],Imu[mu_order],k=2, s=0.0)
        r_spl = np.linspace(mu[mu_order][0],1,5000)
        mu_spl = np.sqrt(1-r_spl**2)
        Imu_spl = splev(mu_spl,tck,der=0)
        csol  = fmin(ldres_fmin, c0, maxiter=1000, maxfun=2000,args=(mu_spl,Imu_spl,law),disp=0)
    else:
        raise ValueError("Fitmethod {} not recognised".format(fitmethod))
    myfit = globals()['ld_%s'%(law)](mu,csol)
    res =  np.sum(Imu - myfit)**2
    int1,int2 = np.trapz(Imu,x=mu),np.trapz(myfit,x=mu)
    dflux = (int1-int2)/int1
    return csol,res,dflux


def fit_law_to_grid(photband,vrads=[0],ebvs=[0],zs=[0],
             law='claret',fitmethod='equidist_r_leastsq',**kwargs):
    """
    Gets the grid and fits LD law to all the models.

    Does not use mu=0 point in the fitting process.
    """
    teffs, loggs = get_grid_dimensions(**kwargs)
    teffs=teffs
    loggs=loggs
    grid_pars = []
    grid_coeffs = []
    Imu1s = []
    logger.info('Fitting photband {}'.format(photband))
    for teff_, logg_ in zip(teffs, loggs):
        for ebv_,vrad_,z_ in itertools.product(ebvs,vrads,zs):
            logger.info('teff={}, logg={}, ebv={}, vrad={}, z={}'.format(teff_, logg_,ebv_,vrad_,z_))
            mu, Imu = get_limbdarkening(teff=teff_, logg=logg_, ebv=ebv_,vrad=vrad_,z=z_,photbands=[photband],**kwargs)
            Imu1 = Imu.max()
            Imu = Imu[:,0]/Imu1
            coeffs,res,dflux = fit_law(mu[mu>0], Imu[mu>0],law=law,fitmethod=fitmethod)
            grid_coeffs.append(coeffs)
            grid_pars.append([teff_,logg_,ebv_,vrad_,z_])
            Imu1s.append([Imu1,res,dflux])
    #- wrap up results in nice arrays
    grid_pars = np.array(grid_pars)
    grid_coeffs = np.array(grid_coeffs)
    Imu1s = np.array(Imu1s)

    return grid_pars, grid_coeffs, Imu1s

def generate_grid(photbands,vrads=[0],ebvs=[0],zs=[0],
             law='claret',fitmethod='equidist_r_leastsq',outfile='mygrid.fits',**kwargs):

    if os.path.isfile(outfile):
        hdulist = pf.open(outfile,mode='update')
        existing_bands = [ext.header['extname'] for ext in hdulist[1:]]
    else:
        hdulist = pf.HDUList([])
        hdulist.append(pf.PrimaryHDU(np.array([[0,0]])))
        existing_bands = []

    hd = hdulist[0].header
    hd.update('FIT', fitmethod, 'FIT ROUTINE')
    hd.update('LAW', law, 'FITTED LD LAW')
    hd.update('GRID', kwargs.get('grid',defaults['grid']), 'GRID')

    for photband in photbands:
        if photband in existing_bands:
            logger.info('BAND {} already exists: skipping'.format(photband))
            continue
        pars,coeffs,Imu1s = fit_law_to_grid(photband,vrads=vrads,ebvs=ebvs,zs=zs,**kwargs)
        cols = []

        cols.append(pf.Column(name='Teff', format='E', array=pars[:,0]))
        cols.append(pf.Column(name="logg", format='E', array=pars[:,1]))
        cols.append(pf.Column(name="ebv" , format='E', array=pars[:,2]))
        cols.append(pf.Column(name="vrad", format='E', array=pars[:,3]))
        cols.append(pf.Column(name="z"   , format='E', array=pars[:,4]))
        for col in range(coeffs.shape[1]):
            cols.append(pf.Column(name='a{:d}'.format(col+1), format='E', array=coeffs[:,col]))
        cols.append(pf.Column(name='Imu1', format='E', array=Imu1s[:,0]))
        cols.append(pf.Column(name='SRS', format='E', array=Imu1s[:,1]))
        cols.append(pf.Column(name='dint', format='E', array=Imu1s[:,2]))

        newtable = pf.new_table(pf.ColDefs(cols))
        newtable.header.update('EXTNAME', photband, "SYSTEM.FILTER")
        newtable.header.update('SYSTEM', photband.split('.')[0], 'PASSBAND SYSTEM')
        newtable.header.update('FILTER', photband.split('.')[1], 'PASSBAND FILTER')
        hdulist.append(newtable)

    if os.path.isfile(outfile):
        hdulist.close()
    else:
        hdulist.writeto(outfile)
#}


#{ Internal functions

def _r(mu):
    """
    Convert mu to r coordinates
    """
    return np.sqrt(1.-mu**2.)

def _mu(r_):
    """
    Convert r to mu coordinates
    """
    return np.sqrt(1.-r_**2.)

def ldres_fmin(coeffs, mu, Imu, law):
    """
    Residual function for Nelder-Mead optimizer.
    """
    return sum((Imu - globals()['ld_%s'%(law)](mu,coeffs))**2)

def ldres_leastsq(coeffs, mu, Imu, law):
    """
    Residual function for Levenberg-Marquardt optimizer.
    """
    return Imu - globals()['ld_%s'%(law)](mu,coeffs)

def _I_ls(ll,ss):
    """
    Limb darkening moments (Dziembowski 1977, Townsend 2002)

    Recursive implementation.

    >>> _I_ls(0,0.5)
    0.6666666666666666
    >>> _I_ls(1,0.5)
    0.4
    >>> _I_ls(2,0.5)
    0.09523809523809523
    """
    if ll==0:
        return 1./(1.+ss)
    elif ll==1:
        return 1./(2.+ss)
    else:
        return (ss-ll+2)/(ss+ll-2+3.)*_I_ls(ll-2,ss)

@memoized
def get_ld_grid(photband,**kwargs):
    """
    Retrieve an interpolating grid for the LD coefficients

    Check outcome:

    #>>> bands = ['GENEVA.U', 'GENEVA.B', 'GENEVA.G', 'GENEVA.V']
    #>>> f_ld_grid = get_ld_grid(bands)
    #>>> ff = pf.open(_atmos['file'])
    #>>> all(ff['GENEVA.U'].data[257][2:]==f_ld_grid(ff['GENEVA.U'].data[257][0],ff['GENEVA.U'].data[257][1])[0:5])
    #True
    #>>> all(ff['GENEVA.G'].data[257][2:]==f_ld_grid(ff['GENEVA.G'].data[257][0],ff['GENEVA.G'].data[257][1])[10:15])
    #True
    #>>> ff.close()

    #Make some plots:

    #>>> photband = ['GENEVA.V']
    #>>> f_ld = get_ld_grid(photband)
    #>>> logg = 4.0
    #>>> mu = linspace(0,1,100)
    #>>> p = figure()
    #>>> p = gcf().canvas.set_window_title('test of function <get_ld_grid>')
    #>>> for teff in linspace(9000,12000,19):
    #...    out = f_ld(teff,logg)
    #...    a1x,a2x,a3x,a4x, I_x1 = out.reshape((len(photband),5)).T
    #...    p = subplot(221);p = title('Interpolation of absolute intensities')
    #...    p = plot(teff,I_x1,'ko')
    #...    p = subplot(222);p = title('Interpolation of LD coefficients')
    #...    p = scatter(4*[teff],[a1x,a2x,a3x,a4x],c=range(4),vmin=0,vmax=3,cmap=cm.spectral,edgecolors='none')
    #...    p = subplot(223);p = title('Without absolute intensity')
    #...    p = plot(mu,ld_eval(mu,[a1x,a2x,a3x,a4x]),'-')
    #...    p = subplot(224);p = title('With absolute intensity')
    #...    p = plot(mu,I_x1*ld_eval(mu,[a1x,a2x,a3x,a4x]),'-')

    """
    #-- retrieve the grid points (unique values)
    teffs,loggs = get_ld_grid_dimensions(**kwargs)
    teffs_grid = np.sort(np.unique1d(teffs))
    loggs_grid = np.sort(np.unique1d(loggs))
    coeff_grid = np.zeros((len(teffs_grid),len(loggs_grid),5*len(photband)))

    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)
    #-- fill the grid
    ff = pf.open(gridfile)
    for pp,iband in enumerate(photband):
        teffs = ff[iband].data.field('Teff')
        loggs = ff[iband].data.field('logg')
        for ii,(iteff,ilogg) in enumerate(zip(teffs,loggs)):
            indext = np.searchsorted(teffs_grid,iteff)
            indexg = np.searchsorted(loggs_grid,ilogg)
            #-- array and list are added for backwards compatibility with some
            #   pyfits versions
            coeff_grid[indext,indexg,5*pp:5*(pp+1)] = np.array(list(ff[iband].data[ii]))[2:]
    ff.close()
    #-- make an interpolating function
    f_ld_grid = InterpolatingFunction([teffs_grid,loggs_grid],coeff_grid)
    return f_ld_grid


@memoized
def _get_itable_markers(photband,gridfile,
                    teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf)):
    """
    Get a list of markers to more easily retrieve integrated fluxes.
    """
    clear_memoization()

    ff = pf.open(gridfile)
    ext = ff[photband]
    columns = ext.columns.names
    names = ['teff','logg']
    grid = [ext.data.field('teff'),
            ext.data.field('logg')]
    if 'ebv' in columns:
        names.append('ebv')
        grid.append(ext.data.field('ebv'))
    if 'vrad' in columns:
        names.append('vrad')
        grid.append(ext.data.field('vrad'))

    grid_axes = [np.sort(list(set(i))) for i in grid]

    #-- we construct an array representing the teff-logg-ebv content, but
    #   in one number: 50000400 means:
    #   T=50000,logg=4.0
    N = len(grid[0])
    markers = np.zeros(N,float)
    gridpnts = np.zeros((N,len(grid)))
    pars = np.zeros((N,5))

    for i,entries in enumerate(zip(*grid)):
        markers[i] = encode(**{key:entry for key,entry in zip(names,entries)})
        gridpnts[i]= entries
        pars[i] = list(ext.data[i][-5:])
    ff.close()
    sa = np.argsort(markers)
    print('read in gridfile',gridfile)
    pars[:,-1] = np.log10(pars[:,-1])
    return markers[sa],grid_axes,gridpnts[sa],pars[sa]

#}


if __name__=="__main__":
    import sys
    if not sys.argv[1:]:
        import pylab as pl
        import doctest
        doctest.testmod()
        pl.show()
    else:
        #from ivs.aux import argkwargparser
        #method,args,kwargs = argkwargparser.parse()
        #print("Calling {} with args=({}) and kwargs=({})".format(method,args,kwargs)
        photbands = ['JOHNSON.U','JOHNSON.B','JOHNSON.V','KEPLER.V','COROT.SIS','COROT.EXO']
        photbands+= ['2MASS.J','2MASS.H','2MASS.KS']
        #generate_grid(photbands,vrads=np.arange(-500,501,50),ebvs=np.arange(0,2.005,0.01),law='claret',outfile='claret.fits')
        #generate_grid(photbands,vrads=np.arange(-500,501,50),ebvs=np.arange(0,2.005,0.01),law='linear',outfile='linear.fits')
        generate_grid(['MOST.V','IRAC.36','COROT.EXO'],vrads=[0],ebvs=[0.06],zs=[0,0],
             law='claret',fitmethod='equidist_r_leastsq',outfile='HD261903.fits')
