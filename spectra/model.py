# -*- coding: utf-8 -*-
"""
Interface to stellar spectra library and functions to manipulate them.

Make a plot of the domains of all spectral grids. First collect all the grid
names

>>> grids = get_gridnames()

Prepare the plot

>>> p = pl.figure()
>>> color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(grids))]
>>> p = pl.gca().set_color_cycle(color_cycle)

And plot all the grid points. We have to set some custom default values for
some grids.

>>> for grid in grids:
...    vturb = 'ostar' in grid and 10 or 2
...    t = 0.
...    if 'marcs' in grid: t = 1.
...    teffs,loggs = get_grid_dimensions(grid=grid,vturb=vturb,t=t)
...    p = pl.plot(np.log10(teffs),loggs,'o',ms=7,label=grid)

Now take care of the plot details

>>> p = pl.xlim(pl.xlim()[::-1])
>>> p = pl.ylim(pl.ylim()[::-1])
>>> p = pl.xlabel('Effective temperature [K]')
>>> p = pl.ylabel('log( Surface gravity [cm s$^{-1}$]) [dex]')
>>> xticks = [3000,5000,7000,10000,15000,25000,35000,50000,65000]
>>> p = pl.xticks([np.log10(i) for i in xticks],['%d'%(i) for i in xticks])
>>> p = pl.legend(loc='upper left',prop=dict(size='small'))
>>> p = pl.grid()

]include figure]]ivs_spectra_model_grid.png]

"""
import os
import logging
import copy
import astropy.io.fits as pf
import inspect
from ivs import config
from ivs.spectra import tools
from ivs.aux import loggers

import numpy as np
from Scientific.Functions.Interpolation import InterpolatingFunction

logger = logging.getLogger("SPEC.MODEL")
logger.addHandler(loggers.NullHandler)

defaults = dict(grid='atlas',z=+0.0,vturb=2,band='vis',
                t=0.0,a=0.0,c=0.5,atm='p')          # MARCS and COMARCS
basedir = 'spectables'

#{ Interface to library

def set_defaults(**kwargs):
    """
    Set defaults of this module.

    If you give no keyword arguments, the default values will be reset.
    """
    if not kwargs:
        kwargs = dict(grid='atlas',z=+0.0,vturb=2,band='vis',
                t=0.0,a=0.0,c=0.5,atm='p')          # MARCS and COMARCS

    for key in kwargs:
        if key in defaults:
            defaults[key] = kwargs[key]


def get_gridnames():
    """
    Return a list of available grid names.

    @return: list of grid names
    @rtype: list of str
    """
    return ['cmfgen','ostar2002','bstar2006','atlas','marcs', 'heberb', 'hebersdb','tmapsdb']



def get_file(**kwargs):
    """
    Retrieve the filename containing the specified SED grid.

    Available grids and example keywords:
        - grid='fastwind': no options
        - grid='cmfgen': no options
        - grid='marcs': options c, atm (p or s)
        - grid='ostar2002': options z,v
        - grid='bstar2006a': options z,v
        - grid='bstar2006b': options z,v
        - grid='bstar2006': options z,v,a
        - grid='atlas': options z
        - grid='heberb': no options
        - grid='hebersdb': no options
        - grid='tmapsdb': no options

    Details for grid 'bstar2006':
        - metallicity in Z/Z0 with Z0 solar. z=0,0.001,0.01,0.033,0.1,0.2,0.5,1.,2
        - microturbulent velocity: v=2 or v=10 km/s
        - abundance: a='' or a='CN'. In the latter, the Helium abundance is
        increased to He/H=0.2, the nitrogen abundance is increased
        by a factor of 5, and the carbon abundance is halved (CNO cycle processed
        material brought to the stellar surface)

    Details for grid 'heberb': LTE Grid computed for B-type stars by Uli Heber,
    reff: Heber et al. 2000

    Details for grid 'hebersdb': LTE Grid computed for sdB stars by Uli Heber,
    reff: Heber et al. 2000

    Details for grid 'tmapsdb': NLTE Grid computed for sdB stars using the TMAP
    (TUEBINGEN NLTE MODEL ATMOSPHERE PACKAGE) code. reff: Werner K., et al. 2003
    and Rauch T., Deetjen J.L. 2003

    @param grid: gridname, or path to grid.
    @type grid: string
    """
    #-- possibly you give a filename
    grid = kwargs.get('grid',defaults['grid'])#.lower()
    if os.path.isfile(grid):
        return grid

    #-- general
    z = kwargs.get('z',defaults['z'])
    #-- only for Kurucz
    vturb = int(kwargs.get('vturb',defaults['vturb']))
    #-- only for Marcs
    t = kwargs.get('t',defaults['t'])
    a = kwargs.get('a',defaults['a'])
    c = kwargs.get('c',defaults['c'])
    atm = kwargs.get('atm',defaults['atm'])
    #-- only for TLUSTY
    band = kwargs.get('band',defaults['band'])

    #-- figure out what grid to use
    if grid=='cmfgen':
        basename = 'cmfgen_spectra.fits'
    elif grid=='marcs':
        basename = 'marcsp_%sz%.1ft%.1f_a%.2f_c%.2f_spectra.fits'%(atm,z,t,a,c)
    elif grid=='ostar2002':
        basename = 'OSTAR2002_z%.3fv%d_%s_spectra.fits'%(z,vturb,band)
    elif grid=='bstar2006':
        basename = 'BSTAR2006_z%.3fv%d_%s_spectra.fits'%(z,vturb,band)
    elif grid=='atlas':
        basename = 'ATLASp_z%.1ft%.1f_a%.2f_spectra.fits'%(z,t,a)
    elif grid=='heberb':
        basename = 'Heber2000_B_h909.fits'
    elif grid=='hebersdb':
        basename = 'Heber2000_sdB_h909.fits'
    elif grid=='tmapsdb':
        basename = 'TMAP2011_sdB.fits'
    else:
        raise ValueError("grid %s does not exist"%(grid))

    #-- retrieve the absolute path of the file and check if it exists:
    grid = config.get_datafile(basedir,basename)

    return grid


def get_table(teff=None,logg=None,vrad=0,vrot=0,**kwargs):
    """
    Retrieve synthetic spectrum.

    Wavelengths in angstrom, fluxes in eddington flux: erg/s/cm2/A.

    It is possible to rotationally broaden the spectrum by supplying at least
    'vrot' and optionally also other arguments for the rotation.rotin function.
    Supply vrot in km/s

    It is possibility to include a radial velocity shift in the spectrum. Supply
    'vrad' in km/s. (+vrad means redshift). Uses spectra.tools.doppler_shift
    """
    gridfile = get_file(**kwargs)
    template_wave = kwargs.pop('wave',None)

    ff = pf.open(gridfile)

    try:
        #-- extenstion name as in fits files prepared by Steven
        mod_name = "T%05d_logg%01.02f" %(teff,logg)
        mod = ff[mod_name]
        wave = mod.data.field('wavelength')
        flux = mod.data.field('flux')
        cont = mod.data.field('cont')
        logger.debug('Model spectrum taken directly from file (%s)'%(os.path.basename(gridfile)))
        if template_wave is not None:
            flux = np.interp(template_wave,wave,flux)
            cont = np.interp(template_wave,wave,cont)
            wave = template_wave
        #-- if the teff/logg is not present, use the interpolation thing
    except KeyError:
        #-- it is possible we first have to set the interpolation function.
        #   This function is memoized, so if it will not be calculated
        #   twice.
        meshkwargs = copy.copy(kwargs)
        meshkwargs['wave'] = kwargs.get('wave',None)
        meshkwargs['teffrange'] = kwargs.get('teffrange',None)
        meshkwargs['loggrange'] = kwargs.get('loggrange',None)
        wave,teffs,loggs,flux,flux_grid,cont_grid = get_grid_mesh(wave=template_wave,**kwargs)
        logger.debug('Model spectrum interpolated from grid %s (%s)'%(os.path.basename(gridfile),meshkwargs))
        wave = wave + 0.
        try:
            flux = flux_grid(np.log10(teff),logg) + 0.
            cont = cont_grid(np.log10(teff),logg) + 0.
        except ValueError:
            teffs,loggs = get_grid_dimensions(**kwargs)
            index = np.argmin(np.abs(  (teffs-teff)**2 + (loggs-logg)**2 ))
            #logger.error('teff=%f-->%f, logg=%f-->%f'%(teff,teffs[index],logg,loggs[index]))
            flux = flux_grid(np.log10(teffs[index]),loggs[index]) + 0.
            cont = cont_grid(np.log10(teffs[index]),loggs[index]) + 0.

    #-- convert to arrays
    wave = np.array(wave,float)
    flux = np.array(flux,float)

    if vrot>0:
        #-- calculate rotational broadening but reinterpolate on original
        #   wavelength grid. First we need to check which arguments we can pass
        #   through
        argspec = inspect.getargspec(tools.rotational_broadening)
        mykwargs = dict(list(zip(argspec.args[-len(argspec.defaults):],argspec.defaults)))
        for key in kwargs:
            if key in mykwargs:
                mykwargs[key] = kwargs[key]
        wave_,flux_ = tools.rotational_broadening(wave,flux,vrot,**mykwargs)
        flux = np.interp(wave,wave_,flux_)
    if vrad!=0:
        wave_ = tools.doppler_shift(wave,vrad)
        flux = np.interp(wave,wave_,flux)

    ff.close()
    return wave,flux,cont









def get_grid_dimensions(**kwargs):
    """
    Retrieve possible effective temperatures and gravities from a grid.

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






#@memoized
def get_grid_mesh(wave=None,teffrange=None,loggrange=None,**kwargs):
    """
    Return InterpolatingFunction spanning the available grid of spectrum models.

    WARNING: the grid must be entirely defined on a mesh grid, but it does not
    need to be equidistant.

    It is thus the user's responsibility to know whether the grid is evenly
    spaced in logg and teff

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
    @rtype: (1Darray,1Darray,1Darray,3Darray,InterpolatingFunction)
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
    #-- clip if necessary
    teffs = list(set(list(teffs)))
    loggs = list(set(list(loggs)))
    teffs = np.sort(teffs)
    loggs = np.sort(loggs)
    if wave is not None:
        flux = np.ones((len(teffs),len(loggs),len(wave)))
        cont = np.ones((len(teffs),len(loggs),len(wave)))
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
                cont_ = mod.data.field('cont')#array(mod.data.tolist())[:,1]
                #-- if there is no wavelength range given, we assume that
                #   the whole grid has the same resolution, and the first
                #   wave-array will be used as a template
                if wave is None:
                    wave = wave_
                    flux = np.ones((len(teffs),len(loggs),len(wave)))
                    cont = np.ones((len(teffs),len(loggs),len(wave)))
            except KeyError:
                continue
            #-- it could be that we're lucky and the grid is completely
            #   homogeneous. In that case, there is no need for interpolation
            try:
                flux[i,j,:] = flux_
                cont[i,j,:] = cont_
            except:
                flux[i,j,:] = np.interp(wave,wave_,flux_)
                cont[i,j,:] = np.interp(wave,wave_,cont_)
    flux_grid = InterpolatingFunction([np.log10(teffs),loggs],flux)
    cont_grid = InterpolatingFunction([np.log10(teffs),loggs],cont)
    #logger.info('Constructed spectrum interpolation grid')
    return wave,teffs,loggs,flux,flux_grid,cont_grid



#}

if __name__=="__main__":
    from ivs.aux import loggers
    logger = loggers.get_basic_logger()
    logger.setLevel(10)
    #import doctest
    import pylab as pl
    import numpy as np
    #doctest.testmod()

    grids = get_gridnames()

    p = pl.figure()
    color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(grids))]
    p = pl.gca().set_color_cycle(color_cycle)

    for grid in grids:
        vturb = 'ostar' in grid and 10 or 2
        t = 0.
        if 'marcs' in grid: t = 1.
        teffs,loggs = get_grid_dimensions(grid=grid,vturb=vturb,t=t)
        p = pl.plot(np.log10(teffs),loggs,'o',ms=7,label=grid)

    p = pl.xlim(pl.xlim()[::-1])
    p = pl.ylim(pl.ylim()[::-1])
    p = pl.xlabel('Effective temperature [K]')
    p = pl.ylabel('log( Surface gravity [cm s$^{-1}$]) [dex]')
    xticks = [3000,5000,7000,10000,15000,25000,35000,50000,65000]
    p = pl.xticks([np.log10(i) for i in xticks],['%d'%(i) for i in xticks])
    p = pl.legend(loc='upper left',prop=dict(size='small'))
    p = pl.grid()

    pl.show()
