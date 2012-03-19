# -*- coding: utf-8 -*-
"""
Interface to the limb-darkening library.
"""
import logging
import os
import pyfits
import numpy as np
from Scientific.Functions.Interpolation import InterpolatingFunction
from ivs.aux import loggers
from ivs.sed import reddening
from ivs.sed import model
from ivs.spectra import tools
from ivs.units import constants
from ivs.aux.decorators import memoized,clear_memoization
from ivs import config

logger = logging.getLogger("SED.LIMBDARK")
logger.addHandler(loggers.NullHandler)

#-- default values for grids
defaults = dict(grid='kurucz',odfnew=False,z=+0.0,vturb=2,
                alpha=False,nover=False)
#-- relative location of the grids
basedir = 'ldtables'

#{ Interface to library

def set_defaults(**kwargs):
    """
    Set defaults of this module
    
    If you give no keyword arguments, the default values will be reset.
    """
    clear_memoization(keys=['ivs.sed.ld'])
    if not kwargs:
        kwargs = dict(grid='kurucz',odfnew=False,z=+0.0,vturb=2,
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
    ff = pyfits.open(gridfile)
    teffs = []
    loggs = []
    for mod in ff[1:]:
        teffs.append(float(mod.header['TEFF']))
        loggs.append(float(mod.header['LOGG']))
    ff.close()
    return np.array(teffs),np.array(loggs)


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
              wave_units='A',flux_units='erg/cm2/s/A/sr',**kwargs):
    """
    Retrieve the specific intensity of a model atmosphere.
    
    ebv is reddening
    vrad is radial velocity: positive is redshift, negative is blueshift (km/s!)
    
    extra kwargs are for reddening
    
    wave and flux units cannot be specificed
    """
    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)
    
    #-- read the file:
    ff = pyfits.open(gridfile)
    
    teff = float(teff)
    logg = float(logg)
    
    #-- if we have a grid model, no need for interpolation
    #try:
    if True:
        #-- extenstion name as in fits files prepared by Steven
        mod_name = "T%05d_logg%01.02f" %(teff,logg)
        mod = ff[mod_name]
        mu = np.array(mod.columns.names[1:], float)
        table = np.array(mod.data.tolist())[:,1:]
        wave = mod.data.field('wavelength')
        logger.debug('Model LD taken directly from file (%s)'%(os.path.basename(gridfile)))
    
    #except KeyError:
    #    raise NotImplementedError
    
    ff.close()
    
    #-- velocity shift if necessary
    if vrad is not None and vrad>0:
        cc = constants.cc/1000. #speed of light in cc
        for i in range(len(mu)):
            wave_shift,flux_shift = tools.doppler_shift(wave,vrad,flux=table[:,i])
            table[:,i] = flux_shift - 5.*vrad/cc*table[:,i]
    
    #-- redden if necessary
    if ebv is not None and ebv>0:
        for i in range(len(mu)):
            table[:,i] = reddening.redden(table[:,i],wave=wave,ebv=ebv,rtype='flux',**kwargs)
    
    
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
        print 'Used teff and logg',teff,logg
        raise
    a1x_,a2x_,a3x_,a4x_, I_x1 = out.reshape((len(photbands),5)).T
    Imu = ld_eval(mu,[a1x_,a2x_,a3x_,a4x_])
    if absolute:
        return Imu*I_x1
    else:
        return Imu


def get_limbdarkening(teff=None,logg=None,ebv=None,vrad=None,photbands=None,**kwargs):
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
        >>> mu,intensities = get_limbdarkening(teff=teff,logg=logg,system='JOHNSON',band='V')

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
    #-- retrieve model atmosphere for a given teff and logg
    mus,wave,table = get_table(teff=teff,logg=logg,ebv=ebv,vrad=vrad,**kwargs)
    #-- compute intensity over the stellar disk, and normalise
    intensities = np.zeros((len(mus),len(photbands)))
    for i in range(len(mus)):
        intensities[i] = model.synthetic_flux(wave,table[:,i],photbands)
    #-- or compute the intensity only for one angle:
    logger.info('Calculated LD')
    return mus,intensities



def get_ld_grid_dimensions(**kwargs):
    """
    Returns the gridpoints of the limbdarkening grid (not unique values).
    """
    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)
    #-- the teff and logg range is the same for every extension, so just
    #   take the first one
    ff = pyfits.open(gridfile)
    teff_,logg_ = ff[1].data.field('Teff'),ff[1].data.field('logg')
    ff.close()
    return teff_,logg_

#}
#{ Laws

def ld_eval(mu,coeffs):
    """
    Evaluate Claret's limb darkening law.
    """
    a1,a2,a3,a4 = coeffs
    Imu = 1-a1*(1-mu**0.5)-a2*(1-mu)-a3*(1-mu**1.5)-a4*(1-mu**2.)    
    return Imu
    
def ld_claret(mu,coeffs):
    return ld_eval(mu,coeffs)
    
def ld_nonlinear(mu,coeffs):
    return 1-coeffs[0]*(1-mu)-coeffs[1]*mu*np.log(mu)

def ld_quadratic(mu,coeffs):
    return 1-coeffs[0]*(1-mu)-coeffs[1]*(1-mu)**2.0

def ld_uniform(mu,coeffs):
    return 1. 
    
#}


@memoized
def get_ld_grid(photband,**kwargs):
    """
    Retrieve an interpolating grid for the LD coefficients
    
    Check outcome:
    
    >>> bands = ['GENEVA.U', 'GENEVA.B', 'GENEVA.G', 'GENEVA.V']
    >>> f_ld_grid = get_ld_grid(bands)
    >>> ff = pyfits.open(_atmos['file'])
    >>> all(ff['GENEVA.U'].data[257][2:]==f_ld_grid(ff['GENEVA.U'].data[257][0],ff['GENEVA.U'].data[257][1])[0:5])
    True
    >>> all(ff['GENEVA.G'].data[257][2:]==f_ld_grid(ff['GENEVA.G'].data[257][0],ff['GENEVA.G'].data[257][1])[10:15])
    True
    >>> ff.close()
    
    Make some plots:
    
    >>> photband = ['GENEVA.V']
    >>> f_ld = get_ld_grid(photband)
    >>> logg = 4.0
    >>> mu = linspace(0,1,100)
    >>> p = figure()
    >>> p = gcf().canvas.set_window_title('test of function <get_ld_grid>')
    >>> for teff in linspace(9000,12000,19):
    ...    out = f_ld(teff,logg)
    ...    a1x,a2x,a3x,a4x, I_x1 = out.reshape((len(photband),5)).T
    ...    p = subplot(221);p = title('Interpolation of absolute intensities')
    ...    p = plot(teff,I_x1,'ko')
    ...    p = subplot(222);p = title('Interpolation of LD coefficients')
    ...    p = scatter(4*[teff],[a1x,a2x,a3x,a4x],c=range(4),vmin=0,vmax=3,cmap=cm.spectral,edgecolors='none')
    ...    p = subplot(223);p = title('Without absolute intensity')
    ...    p = plot(mu,ld_eval(mu,[a1x,a2x,a3x,a4x]),'-')
    ...    p = subplot(224);p = title('With absolute intensity')
    ...    p = plot(mu,I_x1*ld_eval(mu,[a1x,a2x,a3x,a4x]),'-')    
    
    """
    #-- retrieve the grid points (unique values)
    teffs,loggs = get_ld_grid_dimensions(**kwargs)
    teffs_grid = np.sort(np.unique1d(teffs))
    loggs_grid = np.sort(np.unique1d(loggs))
    coeff_grid = np.zeros((len(teffs_grid),len(loggs_grid),5*len(photband)))
    
    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)
    #-- fill the grid
    ff = pyfits.open(gridfile)
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

if __name__=="__main__":
    import doctest
    doctest.testmod()
    