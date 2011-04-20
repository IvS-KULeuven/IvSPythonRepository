# -*- coding: utf-8 -*-
"""
Interface to the SED library.

The most basic usage of this module is:

>>> wave,flux = get_table(teff=10000,logg=4.0)

This will retrieve the model SED with the specified effective temperature and
logg, from the standard grid, in standard units and with zero reddening. All
these things can be specified though.
"""
import os
import logging
import copy
import pyfits
import numpy as np
import pylab as pl
from Scientific.Functions.Interpolation import InterpolatingFunction

from ivs import config
from ivs.units import conversions
from ivs.misc import loggers
from ivs.misc.decorators import memoized
from ivs.sed import filters
from ivs.sed import reddening
from ivs.io import ascii

logger = logging.getLogger("SED.MODEL")
logger.addHandler(loggers.NullHandler)

caldir = os.sep.join(['sedtables','calibrators'])

#-- default values for grids
defaults = dict(grid='kurucz',odfnew=True,z=+0.0,vturb=2,
                alpha=False,nover=False,                  # KURUCZ
                He=97,                                    # WD
                t=1.0,a=0.0,c=0.5,m=1.0,co=1.05)          # MARCS and COMARCS
#-- relative location of the grids
basedir = 'sedtables/modelgrids/'

#{ Interface to library

def set_defaults(**kwargs):
    """
    Set defaults of this module
    """
    for key in kwargs:
        if key in defaults:
            defaults[key] = kwargs[key]

def reset_defaults():
    """
    Reset defaults to original defaults
    """
    mydefaults = dict(grid='kurucz',odfnew=True,z=+0.0,vturb=2,
                alpha=False,nover=False,                  # KURUCZ
                He=97,                                    # WD
                t=1.0,a=0.0,c=0.5,m=1.0,co=1.05)          # MARCS and COMARCS
    for key in mydefaults:
        defaults[key] = mydefaults[key]

def get_gridnames():
    """
    Return a list of available grid names.
    """
    return ['kurucz','fastwind','cmfgen','sdb_uli','wd_boris','wd_da','wd_db',
            'tlusty','uvblue','atlas12']
            #'marcs','marcs2','comarcs','tlusty','uvblue','atlas12']

def get_file(**kwargs):
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
        - grid='tlusty':
                    * z: log10(Z/Z0)
        - grid='sdb_uli': metallicity and helium fraction (z, he=98)
        - grid='fastwind': no options
        - grid='wd_boris': no options
        - grid='stars': precomputed stars (vega, betelgeuse...)
        - grid='uvblue'
        - grid='marcs'
        - grid='marcs2'
        - grid='atlas12'
    
    @keyword grid: gridname (default Kurucz)
    @type grid: str
    """
    #-- possibly you give a filename
    grid = kwargs.get('grid',defaults['grid'])
    if os.path.isfile(grid):
        return grid
    
    grid = grid.lower()
    
    #-- general
    z = kwargs.get('z',defaults['z'])
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
    
    #-- figure out what grid to use
    if grid=='fastwind':
        basename = 'fastwind_sed.fits'
        
    elif grid=='kurucz':
        if not alpha and not nover and not odfnew:
            basename = 'kurucz93_z%.1f_k%d_sed.fits'%(z,vturb)
        elif alpha and odfnew:
            basename = 'kurucz93_z%.1f_ak%dodfnew_sed.fits'%(z,vturb)
        elif odfnew:
            basename = 'kurucz93_z%.1f_k%dodfnew_sed.fits'%(z,vturb)
        elif nover:
            basename = 'kurucz93_z%.1f_k%dnover_sed.fits'%(z,vturb)            
    elif grid=='cmfgen':
        basename = 'cmfgen_sed.fits'        
    elif grid=='sdb_uli':
        basename = 'SED_int_h%d_z%.1f.fits'%(He,z)
    elif grid=='wd_boris':
        basename = 'SED_WD_Gaensicke.fits'
    elif grid=='wd_da':
        basename = 'SED_WD_Koester_DA.fits'
    elif grid=='wd_db':
        basename = 'SED_WD_Koester_DB.fits'        
    elif grid=='marcs':
        basename = 'marcsp_z%.1ft%.1f_a%.2f_c%.2f_sed.fits'%(z,t,a,c)
    elif grid=='marcs2':
        basename = 'marcsp2_z0.00t2.0_m.1.0c0.00_sed.fits'
    elif grid=='comarcs':
        basename = 'comarcsp_z%.2fco%.2fm%.1fxi2.50_sed.fits'%(z,co,m)        
    elif grid=='stars':
        basename = 'kurucz_stars_sed.fits'        
    elif grid=='tlusty':
        basename = 'tlusty_z%.2f_sed.fits'%(z)        
    elif grid=='uvblue':
        basename = 'uvblue_z%.1f_k2_sed.fits'%(z)        
    elif grid=='atlas12':
        basename = 'atlas12_z%.1f_sed.fits'%(z)
    
    #-- retrieve the absolute path of the file and check if it exists:
    grid = config.get_datafile(basedir,basename)
    
    return grid







def get_table(teff=None,logg=None,ebv=None,star=None,
              wave_units='A',flux_units='erg/cm2/s/A/sr',**kwargs):
    """
    Retrieve the spectral energy distribution of a model atmosphere.
        
    Wavelengths in A (angstrom)
    Fluxes in Ilambda = ergs/cm2/s/A/sr, except specified via 'units',
    
    If you give 'units', and /sr is not included, you are responsible yourself
    for giving an extra keyword with the angular diameter C{ang_diam}, or other
    possibilities offered by the C{units.conversions.convert} function.
        
    Possibility to redden the fluxes according to the reddening parameter EB_V.
    
    Extra kwargs can specify the grid type.
    Extra kwargs can specify constraints on the size of the grid to interpolate.
    Extra kwargs can specify reddening law types.
    Extra kwargs can specify information for conversions.
    
    Example usage:
    
    >>> from pylab import figure,plot,legend,gca,subplot,title,show,gcf,loglog
    >>> p = figure()
    >>> p=gcf().canvas.set_window_title('Test of <get_table>')
    >>> p = subplot(131)
    >>> p = loglog(*get_table(grid='FASTWIND',teff=35000,logg=4.0))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=35000,logg=4.0))
    >>> p = loglog(*get_table(grid='TLUSTY',teff=35000,logg=4.0))
    >>> p = loglog(*get_table(grid='MARCS',teff=5000,logg=2.0))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=5000,logg=2.0))
    >>> p = subplot(132)
    >>> p = loglog(*get_table(grid='FASTWIND',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='TLUSTY',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='MARCS',teff=5000,logg=2.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=5000,logg=2.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = subplot(133)
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10250,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10500,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=10750,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = loglog(*get_table(grid='KURUCZ',teff=11000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
    >>> p = show()
        
    @param teff: effective temperature
    @type teff: float
    @param logg: logarithmic gravity (cgs)
    @type logg: float
    @param ebv: reddening coefficient
    @type ebv: float
    @param wave_units: units to convert the wavelengths to (if not given, A)
    @type wave_units: str
    @param flux_units: units to convert the fluxes to (if not given, erg/s/cm2/A/sr)
    @type flux_units: str
    @return: wavelength,flux
    @rtype: (ndarray,ndarray)
    """
    
    #-- get the FITS-file containing the tables
    gridfile = get_file(**kwargs)
    
    #-- read the file:
    ff = pyfits.open(gridfile)
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
            meshkwargs = copy.copy(kwargs)
            meshkwargs['wave'] = kwargs.get('wave',None)
            meshkwargs['teffrange'] = kwargs.get('teffrange',None)
            meshkwargs['loggrange'] = kwargs.get('loggrange',None)
            wave,teffs,loggs,flux,flux_grid = get_grid_mesh(**meshkwargs)
            logger.debug('Model SED interpolated from grid %s (%s)'%(os.path.basename(gridfile),meshkwargs))
            wave = wave + 0.
            flux = flux_grid(np.log10(teff),logg) + 0.
    
    #-- convert to arrays
    wave = np.array(wave,float)
    flux = np.array(flux,float)
    #-- redden if necessary
    if ebv is not None:
        if 'wave' in kwargs.keys():
            removed = kwargs.pop('wave')
        flux = reddening.redden(flux,wave=wave,ebv=ebv,rtype='flux',**kwargs)
    
    if flux_units!='erg/s/cm2/A/sr':
        flux = conversions.convert('erg/s/cm2/A/sr',flux_units,flux,wave=(wave,'A'),**kwargs)
    if wave_units!='A':
        wave = conversions.convert('A',wave_units,wave,**kwargs)
    
    ff.close()
    
    return wave,flux

def get_grid_dimensions(**kwargs):
    """
    Retrieve possible effective temperatures and gravities from a grid.
    
    E.g. kurucz, sdB, fastwind...
    
    Example usage:
    
    >>> from pylab import plot,legend,figure,fill,getp,gca,xlabel,ylabel,gcf,show
    >>> from pylab import title,xlim,ylim,cm
    >>> from sigproc.base import bplot
    >>> import matplotlib
    >>> matplotlib.axes.set_default_color_cycle(bplot.get_color_cycle(cm.spectral,res=2))
    >>> p = figure();p=title('Grid of SEDs')
    >>> p=gcf().canvas.set_window_title('Test of <get_grid_dimensions>')
    >>> sizes = range(15,1,-1)
    >>> gridnames = get_gridnames()
    >>> for s,gridname in zip(sizes,np.sort(gridnames)):
    ...    teffs,loggs = get_grid_dimensions(grid=gridname)
    ...    teffs = np.log10(teffs)
    ...    p=plot(teffs,loggs,'o',label=gridname,mew=0,ms=s)
    ...    hull_pts = bplot.convex_hull(np.array([teffs,loggs]))
    ...    p=fill(hull_pts[:,0],hull_pts[:,1],color=getp(p[0],'color'),alpha=0.2)
    >>> p=legend(loc='upper left')
    >>> p=xlabel('log(Teff)')
    >>> p=ylabel('logg')
    >>> p=xlim(xlim()[::-1])
    >>> p=ylim(ylim()[::-1])
    >>> reset_defaults()
    
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
    flux_grid = InterpolatingFunction([np.log10(teffs),loggs],flux)
    logger.info('Constructed SED interpolation grid')
    return wave,teffs,loggs,flux,flux_grid




#}

#{ Calibration

def list_calibrators():
    """
    Print and return the list of calibrators
    """
    files = config.glob(caldir,'*.fits')
    names = []
    for ff in files:
        fits_file = pyfits.open(ff)
        names.append(fits_file[0].header['targetid'])
        star_info = sesame.search(names[-1])
        spType = ('spType' in star_info) and star_info['spType'] or ''
        logger.info('%-15s: (%8s) %s'%(names[-1],spType,fits_file[0].header['descrip']))
        fits_file.close()
    return names
    

def get_calibrator(name='alpha_lyr',version=None,wave_units=None,flux_units=None):
    """
    Retrieve a calibration SED
    
    If C{version} is None, get the last version.
    
    Example usage:
    
    >>> wave,flux = get_calibrator(name='alpha_lyr')
    >>> wave,flux = get_calibrator(name='alpha_lyr',version='003')
    
    """
    #-- collect calibration files
    files = config.glob(caldir,'*.fits')
    calfile = None
    for ff in files:
        #-- check if the name matches with the given one
        fits_file = pyfits.open(ff)
        header = fits_file[0].header
        if name in ff or name in header['targetid']:
            #-- maybe the target is correct, but the 'model version' is not
            if version is not None and version not in ff:
                fits_file.close()
                continue
            #-- extract the wavelengths and flux
            calfile = ff
            wave = fits_file[1].data.field('wavelength')
            flux = fits_file[1].data.field('flux')
        fits_file.close()
    
    if calfile is None:
        raise ValueError, 'Calibrator %s (version=%s) not found'%(name,version)
    
    if flux_units is not None:
        flux = conversions.convert('erg/s/cm2/A',flux_units,flux,wave=(wave,'A'))
    if wave_units is not None:
        wave = conversions.convert('A',wave_units,wave)
    
    
    logger.info('Calibrator %s selected'%(calfile))
    
    return wave,flux



def calibrate():
    """
    Calibrate photometry.
    
    Not finished!
    
    ABmag = -2.5 Log F_nu - 48.6 with F_nu in erg/s/cm2/Hz
    Flux computed as 10**(-(meas-mag0)/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    F0 = 3.6307805477010029e-20 erg/s/cm2/Hz
    
    STmag = -2.5 Log F_lam - 21.10 with F_lam in erg/s/cm2/A
    Flux computed as 10**(-(meas-mag0)/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    F0 = 3.6307805477010028e-09 erg/s/cm2/A
    
    Vegamag = -2.5 Log F_lam - C with F_lam in erg/s/cm2/A
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
    Flam0_lit = conversions.nconvert(zp['Flam0_units'],'erg/s/cm2/A',zp['Flam0'],photband=zp['photband'])
    
    #-- we have Flam0 but not Fnu0: compute Fnu0
    keep = (zp['Flam0_lit']==1) & (zp['Fnu0_lit']==0)
    Fnu0 = conversions.nconvert(zp['Flam0_units'],'erg/s/cm2/Hz',zp['Flam0'],photband=zp['photband'])
    zp['Fnu0'][keep] = Fnu0[keep]
    zp['Fnu0_units'][keep] = 'erg/s/cm2/Hz'
    
    #-- we have Fnu0 but not Flam0: compute Flam0
    keep = (zp['Flam0_lit']==0) & (zp['Fnu0_lit']==1)
    Flam0 = conversions.nconvert(zp['Fnu0_units'],'erg/s/cm2/A',zp['Fnu0'],photband=zp['photband'])
    
    #   set everything in correct units for convenience:
    Flam0 = conversions.nconvert(zp['Flam0_units'],'erg/s/cm2/A',zp['Flam0'])
    Fnu0 = conversions.nconvert(zp['Fnu0_units'],'erg/s/cm2/Hz',zp['Fnu0'])
    
    #-- as a matter of fact, set Flam0 for all the stuff
    keep = (zp['Flam0_lit']==0) & (zp['Fnu0_lit']==0)
    zp['Flam0'][keep] = syn_flux[keep]
    zp['Flam0_units'][keep] = 'erg/s/cm2/A'
    
    keep = np.array(['DENIS' in photb and True or False for photb in zp['photband']])
    print zip(syn_flux[keep],Flam0[keep])
    
    #-- we have no Flam0, only ZP vegamags
    keep = (zp['vegamag_lit']==1) & (zp['Flam0_lit']==0)
    zp['Flam0'][keep] = syn_flux[keep]
    zp['Flam0_units'][keep] = 'erg/s/cm2/A'
    
    #-- we have no Flam0, no ZP vegamas but STmags
    keep = (zp['STmag_lit']==1) & (zp['Flam0_lit']==0)
    m_vega = 2.5*np.log10(F0ST/syn_flux) + zp['STmag']
    zp['vegamag'][keep] = m_vega[keep]
    
    #-- we have no Fnu0, no ZP vegamas but ABmags
    keep = (zp['ABmag_lit']==1) & (zp['Flam0_lit']==0)
    F0AB_lam = conversions.convert('erg/s/cm2/Hz','erg/s/cm2/A',F0AB,photband=zp['photband'])
    m_vega = 2.5*np.log10(F0AB_lam/syn_flux) + zp['ABmag']
    zp['vegamag'][keep] = m_vega[keep]
    
    return zp
    

#}

#{ Synthetic photometry

def synthetic_flux(wave,flux,photbands):
    """
    Extract flux measurements from a synthetic SED.
    
    The fluxes below 4micron are calculated assuming PHOTON-counting detectors
    (e.g. CCDs).
    
    F = int(P_lam * f_lam * lam, dlam) / int(P_lam * lam, dlam)
    
    When otherwise specified, we assume ENEGY-counting detectors (e.g. bolometers)
    
    F = int(P_lam * f_lam, dlam) / int(P_lam, dlam)
    
    Where P_lam is the total system dimensionless sensitivity function, which
    is normalised so that the maximum equals 1. Also, f_lam is the SED of the
    object, in units of energy per time per unit area per wavelength.
    
    The PHOTON-counting part of this routine has been thoroughly checked with
    respect to johnson UBV, geneva and stromgren filters, and only gives offsets
    with respect to the Kurucz integrated files (.geneva and stuff on his websites). These could be
    due to different normalisation.
    
    You are responsible yourself for having a response curve covering the
    model fluxes!
    
    See e.g. Maiz-Apellaniz, 2006.
    
    @param wave: model wavelengths (angstrom)
    @type wave: ndarray
    @param flux: model fluxes (Flam)
    @type flux: ndarray
    @param photbands: list of photometric passbands
    @type photbands: list of str
    @return: model fluxes (Flam)
    @rtype: ndarray
    """    
    energys = np.zeros(len(photbands))
    
    #-- only keep relevant information on filters:
    filter_info = filters.get_info()
    keep = np.searchsorted(filter_info['photband'],photbands)
    filter_info = filter_info[keep]
    
    for i,photband in enumerate(photbands):
        waver,transr = filters.get_response(photband)
        #-- make wavelength range a bit bigger, otherwise F25 from IRAS has only
        #   one Kurucz model point in its wavelength range... this is a bit
        #   'ad hoc' but seems to work.
        region = ((waver[0]-0.4*waver[0])<=wave) & (wave<=(waver[-1]+0.4*waver[-1]))
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
        #-- interpolate response curve onto model grid
        transr = np.interp(wave_,waver,transr,left=0,right=0)
        #-- integrated flux: different for bolometers and CCDs
        if filter_info['type'][i]=='BOL':
            energys[i] = np.trapz(flux_*transr,x=wave_)/np.trapz(transr,x=wave_)
        elif filter_info['type'][i]=='CCD':
            energys[i] = np.trapz(flux_*transr*wave_,x=wave_)/np.trapz(transr*wave_,x=wave_)
    
    #-- that's it!
    return energys


#}

if __name__=="__main__":
    import doctest
    doctest.testmod()
    