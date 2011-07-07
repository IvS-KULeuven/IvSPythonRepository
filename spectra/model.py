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
import pyfits
import pyrotin4
from ivs import config
from ivs.units import constants
from ivs.units import conversions
from ivs.aux.decorators import memoized

import numpy as np
from Scientific.Functions.Interpolation import InterpolatingFunction

logger = logging.getLogger("SPEC.MODEL")

defaults = dict(grid='atlas',z=+0.0,vturb=2,band='vis',
                t=0.0,a=0.0,c=0.5,atm='p')          # MARCS and COMARCS
basedir = 'spectables'

#{ Interface to library

def set_defaults(**kwargs):
    """
    Set defaults of this module
    
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
    return ['cmfgen','ostar2002','bstar2006','atlas','marcs']



def get_file(**kwargs):
    """
    Retrieve the filename containing the specified SED grid.
    
    Details for grid 'bstar2006':
        - metallicity in Z/Z0 with Z0 solar. z=0,0.001,0.01,0.033,0.1,0.2,0.5,1.,2
        - microturbulent velocity: v=2 or v=10 km/s
        - abundance: a='' or a='CN'. In the latter, the Helium abundance is
                    increased to He/H=0.2, the nitrogen abundance is increased
                    by a factor of 5, and the carbon abundance is halved (CNO cycle processed
                    material brought to the stellar surface)
    
    Available grids and example keywords:
        - grid='fastwind': no options
        - grid='cmfgen': no options
        - grid='marcs': options c, atm (p or s)
        - grid='ostar2002': options z,v
        - grid='bstar2006a': options z,v
        - grid='bstar2006b': options z,v
        - grid='bstar2006: options z,v,a
                          
        - grid='atlas': options z
    
    @keyword grid: gridname
    """
    #-- possibly you give a filename
    grid = kwargs.get('grid',defaults['grid']).lower()
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
    else:
        raise ValueError, "grid %s does not exist"%(grid)

    #-- retrieve the absolute path of the file and check if it exists:
    grid = config.get_datafile(basedir,basename)
    
    return grid







def get_table(teff=None,logg=None,vrad=0,vrot=0,**kwargs):
    """
    Retrieve synthetic spectrum.
        
    Wavelengths in angstrom.
    
    Possibility to rotationally broaden the spectrum: supply at least 'vrot' and
    optionally also other keyword arguments for the rotation.rotin function.
    
    Possibility to include radial velocity shift (+vrad means redshift)
    
    vrot and vrad in km/s
    
    """
    gridfile = get_file(**kwargs)
    template_wave = kwargs.pop('wave',None)
    
    ff = pyfits.open(gridfile)
    
    try:
        #-- extenstion name as in fits files prepared by Steven
        mod_name = "T%05d_logg%01.02f" %(teff,logg)
        mod = ff[mod_name]
        wave = mod.data.field('wavelength')
        flux = mod.data.field('flux')
        cont = mod.data.field('cont')
        logger.debug('Model spectrum taken directly from file (%s)'%(os.path.basename(gridfile)))
        #-- if the teff/logg is not present, use the interpolation thing
    except KeyError:
        #-- it is possible we first have to set the interpolation function.
        #   This function is memoized, so if it will not be calculated
        #   twice.
        meshkwargs = copy.copy(kwargs)
        meshkwargs['wave'] = kwargs.get('wave',None)
        meshkwargs['teffrange'] = kwargs.get('teffrange',None)
        meshkwargs['loggrange'] = kwargs.get('loggrange',None)
        wave,teffs,loggs,flux,flux_grid,cont_grid = get_grid_mesh(wave=template_wave)
        logger.debug('Model spectrum interpolated from grid %s (%s)'%(os.path.basename(gridfile),meshkwargs))
        wave = wave + 0.
        try:
            flux = flux_grid(np.log10(teff),logg) + 0.
            cont = cont_grid(np.log10(teff),logg) + 0.
        except ValueError:
            teffs,loggs = get_grid_dimensions(**kwargs)
            index = np.argmin(np.abs(  (teffs-teff)**2 + (loggs-logg)**2 ))
            logger.error('teff=%f-->%f, logg=%f-->%f'%(teff,teffs[index],logg,loggs[index]))
            flux = flux_grid(np.log10(teffs[index]),loggs[index]) + 0.
            cont = cont_grid(np.log10(teffs[index]),loggs[index]) + 0.
    
    #-- convert to arrays
    wave = np.array(wave,float)
    flux = np.array(flux,float)
    
    if vrot>0:
        #-- calculate rotational broadening but reinterpolate on original
        #   wavelength grid
        wave_,flux_ = rotational_broadening(wave,flux,vrot,**kwargs)
        flux = np.interp(wave,wave_,flux_)
    if vrad!=0:
        wave_ = doppler_shift(wave,vrad)
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
    ff = pyfits.open(gridfile)
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
    ff = pyfits.open(gridfile)
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
    logger.info('Constructed spectrum interpolation grid')
    return wave,teffs,loggs,flux,flux_grid,cont_grid

#}

#{ Analysis

def doppler_shift(wave,vrad,vrad_units='km/s'):
    """
    Shift a spectrum with towards the red or blue side with some radial velocity.
    
    You can give units with the extra keywords C{vrad_units} (units of
    wavelengths are not important). The shifted wavelengths will be in the same
    units as the input wave array.
    
    If units are not supplied, the radial velocity is assumed to be in km/s.
    
    If you want to apply a barycentric correction, you'd probably want to
    reverse the sign!
    
    Example usage: shift a spectrum to the blue ('left') with 20 km/s.
    
    >>> wave = np.linspace(3000,8000,1000)
    >>> wave_shift1 = doppler_shift(wave,20.)
    >>> wave_shift2 = doppler_shift(wave,20000.,vrad_units='m/s')
    >>> print(wave_shift1[0],wave_shift1[-1])
    (3000.200138457119, 8000.5337025523177)
    >>> print(wave_shift2[0],wave_shift2[-1])
    (3000.200138457119, 8000.5337025523177)
    
    @param wave: wavelength array
    @type wave: ndarray
    @param vrad: radial velocity (negative shifts towards blue, positive towards red)
    @type vrad: float (units: km/s) or tuple (float,'units')
    @param vrad_units: units of radial velocity (default: km/s)
    @type vrad_units: str (interpretable for C{units.conversions.convert})
    @return: shifted wavelength array
    @rtype: ndarray
    """ 
    cc = constants.cc
    cc = conversions.convert('m/s',vrad_units,cc)
    wave_out = wave * (1+vrad/cc)
    return wave_out


def rotational_broadening(wave_spec,flux_spec,vrot,fwhm=0.25,epsilon=0.6,
                         chard=None,stepr=0,stepi=0,alam0=None,alam1=None,
                         irel=0,cont=None):
    """
    Apply rotational broadening to a spectrum assuming a linear limb darkening
    law.
    
    This function is based on the ROTIN program. See Fortran file for
    explanations of parameters.
    
    Limb darkening law is linear, default value is epsilon=0.6
    
    Possibility to normalize as well by giving continuum in 'cont' parameter.
    2. parameters for rotational convolution 

    VROT  - v sin i (in km/s):
             - if VROT=0 - rotational convolution is 
                 a) either not calculated,
                 b) or, if simultaneously FWHM is rather large
                    (vrot/c*lambda < FWHM/20.),
                    vrot is set to  FWHM/20*c/lambda;
             - if VROT >0 but the previous condition b) applies, the
                     value of VROT is changed as  in the previous case
             - if VROT<0 - the value of abs(VROT) is used regardless of
                     how small compared to FWHM it is
     CHARD - characteristic scale of the variations of unconvolved:
           - stellar spectrum (basically, characteristic distance
             between two neighbouring wavelength points) - in A
           - if =0 - program sets up default (0.01 A)
     STEPR - wavelength step for evaluation rotational convolution;
           - if =0, the program sets up default (the wavelength
                    interval corresponding to the rotational velocity
                    devided by 3.)                           
           - if <0, convolved spectrum calculated on the original
             (detailed) SYNSPEC wavelength mesh


     3. parameters for instrumental convolution

     FWHM  - full width at half maximum for Gaussian instrumental 
             profile
     STEPI - wavelength step for evaluating instrumental convolution
           - if =0, the program sets up default (FWHM/10.)
           - if <0, convolved spectrum calculated with the previous
                    wavelength mesh:
                    either the original (SYNSPEC) one if vrot=0,
                    or the one used in rotational convolution (vrot > 0)


     4. wavelength interval and normalization of spectra

     ALAM0 - initial wavelength
     ALAM1 - final wavelength
     IREL  - for =1 relative spectrum
                 =0 absolute spectrum
    
    @return: wavelength,flux
    @rtype: array, array
    """
    #-- set arguments
    if alam0 is None: alam0 = wave_spec[0]
    if alam1 is None: alam1 = wave_spec[-1]
    if cont is None: cont = (np.ones(1),np.ones(1))
    contw,contf = cont
    if chard is None:
        chard = np.diff(wave_spec).mean()
    
    #-- apply broadening
    w3,f3,ind = pyrotin4.pyrotin(wave_spec,flux_spec,contw,contf,
                  vrot,chard,stepr,fwhm,stepi,alam0,alam1,irel,epsilon)
    logger.info('ROTIN rot.broad. with vrot=%.3f (epsilon=%.2f)'%(vrot,epsilon))
    
    return w3[:ind],f3[:ind]
    
#}
if __name__=="__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()