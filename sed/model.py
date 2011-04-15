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

from ivs import config
from ivs.units import conversions
from ivs.misc import loggers
from ivs.reduction.photometry import calibration

logger = logging.getLogger("SED.MODEL")
logger.addHandler(loggers.NullHandler)

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
    

def get_sed_file(**kwargs):
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







def get_table(**kwargs):
    """
    Retrieve the spectral energy distribution of a model atmosphere.
        
    Wavelengths in A (angstrom)
    Fluxes in Ilambda = ergs/cm2/s/A/sr, except specified via 'units',
    
    If you give 'units', and /sr is not included, you are responsible yourself
    for giving an extra keyword with the angular diameter C{ang_diam}, or other
    possibilities offered by the C{units.conversions.convert} function.
        
    Possibility to redden the fluxes according to the reddening parameter EB_V.
    
    Example usage:
        >>> from pylab import figure,plot,legend,gca,subplot,title,show
        >>> p = figure(); p = subplot(121)
        >>> p = plot(*get_table(grid='FASTWIND',teff=35000,logg=4.0))
        >>> p = plot(*get_table(grid='KURUCZ',teff=35000,logg=4.0))
        >>> p = plot(*get_table(grid='TLUSTY',teff=35000,logg=4.0))
        >>> p = plot(*get_table(grid='MARCS',teff=5000,logg=2.0))
        >>> p = plot(*get_table(grid='KURUCZ',teff=5000,logg=2.0))
        >>> p = gca().set_xscale('log')
        >>> p = gca().set_yscale('log')
        >>> p = subplot(122)
        >>> p = plot(*get_table(grid='FASTWIND',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
        >>> p = plot(*get_table(grid='KURUCZ',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
        >>> p = plot(*get_table(grid='TLUSTY',teff=35000,logg=4.0,wave_units='micron',flux_units='Jy/sr'))
        >>> p = plot(*get_table(grid='MARCS',teff=5000,logg=2.0,wave_units='micron',flux_units='Jy/sr'))
        >>> p = plot(*get_table(grid='KURUCZ',teff=5000,logg=2.0,wave_units='micron',flux_units='Jy/sr'))
        >>> p = gca().set_xscale('log')
        >>> p = gca().set_yscale('log')
        >>> p = show()
        
    @keyword teff: effective temperature
    @type teff: float
    @keyword logg: logarithmic gravity (cgs)
    @type logg: float
    @keyword EB_V: reddening coefficient
    @type EB_V: float
    @keyword wave_units: units to convert the wavelengths to (if not given, A)
    @type wave_units: str
    @keyword flux_units: units to convert the fluxes to (if not given, erg/s/cm2/A/sr)
    @type flux_units: str
    @return: wavelength,flux
    @rtype: (ndarray,ndarray)
    """
    teff = kwargs.pop('teff',None)
    logg = kwargs.pop('logg',None)
    EB_V = kwargs.pop('EB_V',None)
    star = kwargs.pop('star',None)
    wave_units = kwargs.get('wave_units',None)
    flux_units = kwargs.get('flux_units',None)
    
    
    gridfile = get_sed_file(**kwargs)
    
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
            #Table name as in fits files prepared by Steven
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
            meshkwargs = copy.copy(gridfile)
            meshkwargs['wave'] = kwargs.get('wave',None)
            meshkwargs['teffrange'] = kwargs.get('teffrange',None)
            meshkwargs['loggrange'] = kwargs.get('loggrange',None)
            wave,teffs,loggs,flux,flux_grid = get_sed_grid_mesh(**meshkwargs)
            wave = wave + 0.
            logger.debug('Model SED interpolated from grid %s (%s)'%(os.path.basename(gridfile),meshkwargs))
            flux = flux_grid(log10(teff),logg) + 0.
            
    wave = np.array(wave,float)
    flux = np.array(flux,float)
    #-- redden if necessary
    if EB_V is not None:
        if 'wave' in kwargs.keys():
            removed=kwargs.pop('wave')
        flux = photometry.deredden(wave,flux,EB_V=EB_V,**kwargs)
    
    if flux_units:
        flux = conversions.convert('erg/s/cm2/A/sr',flux_units,flux,wave=(wave,'A'),**kwargs)
    if wave_units:
        wave = conversions.convert('A',wave_units,wave,**kwargs)
    
    return wave,flux


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
    
    logger.info('Calibrator %s selected'%(calfile))
    
    return wave,flux



def calibrate():
    """
    Calibrate photometry.
    
    Not finished!
    
    ABmag = -2.5 Log F_nu - 48.6 with F_nu in erg/s/cm2/Hz
    Flux computed as 10**(-(meas-mag0)/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    
    STmag = -2.5 Log F_lam - 21.10 with F_lam in erg/s/cm2/A
    Flux computed as 10**(-(meas-mag0)/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    
    Vegamag = -2.5 Log F_lam - C with F_lam in erg/s/cm2/A
    Flux computed as 10**(-meas/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    """
    #-- get calibrator
    wave,flux = get_calibrator(name='alpha_lyr')
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp = ascii.read2recarray(zp_file,dtype=[('photband','a50'),('eff_wave','>f4'),('type','a3'),
                                            ('vegamag','>f4'),('vegamag_lit','i'),
                                            ('abmag','>f4'),('abmag_lit','i'),
                                            ('stmag','>f4'),('stmag_lit','i'),
                                            ('Flam0','>f4'),('Flam0_units','a50'),('Flam0_lit','i'),
                                            ('Fnu0','>f4'),('Fnu0_units','a50'),('Fnu0_lit','i'),
                                            ('source','a100')])
    #-- calculate zeropoint fluxes
    #   Flam0 from VEGA magnitudes:
    
    #----------------- INTEGRATE THE MODEL ------------------
    keep = zp['vegamag_lit']==1
    energy = model.synthetic_flux(wave,flux,zp['photband'][keep])
    print energy
    #-- calculate VEGA magnitudes:
    zp['vegamag'][zp['vegamag_lit']==0] = 0.
    return wave,flux,zp['eff_wave'][keep],energy
    

#}

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
    """    
    energys = np.zeros(len(photbands))
    
    #-- only keep relevant information on filters:
    filter_info = calibration.get_filter_info()
    keep = np.searchsorted(filter_info['photband'],photbands)
    filter_info = filter_info[keep]
    
    for i,photband in enumerate(photbands):
        waver,transr = calibration.get_response(photband)
        
        #-- make wavelength range a bit bigger, otherwise F25 from IRAS has only
        #   one Kurucz model point in its wavelength range... this is a bit
        #   'ad hoc' but seems to work.
        region = ((waver[0]-0.4*waver[0])<=wave) & (wave<=(waver[-1]+0.4*waver[-1]))
        #-- if we're working in infrared (>4e4A) and the model is not of high
        #   enough resolution (100000 points over wavelength range), interpolate
        #   the model in logscale on to a denser grid (in logscale!)
        if filter_info['wave_eff'][i]>=4e4 and sum(region)<1e5 and sum(region)>1:
            wave_ = np.logspace(np.log10(wave[region][0]),np.log10(wave[region][-1]),1e5)
            flux = 10**np.interp(np.log10(wave_),np.log10(wave[region]),np.log10(flux[region]),)
            wave = wave_
        else:
            wave = wave[region]
            flux = flux[region]
        #-- interpolate response curve onto model grid
        trans_res = np.interp1d(wave,waver,transr,left=0,right=0)
        #-- integrated flux: different for bolometers and CCDs
        if filter_info['type'][i]=='BOL':
            energys[i] = trapz(flux*transr,x=wave)/trapz(transr,x=wave)
        elif filter_info['type'][i]=='CCD':
            energys[i] = trapz(flux*transr*wave,x=wave)/trapz(transr*wave,x=wave)
    
    #-- that's it!
    return energys


#}

if __name__=="__main__":
    import doctest
    doctest.testmod()
    