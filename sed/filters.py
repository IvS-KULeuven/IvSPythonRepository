# -*- coding: utf-8 -*-
"""
Functions relevant for photometric calibration

Table of contents:

    1. Available response functions
    2. Adding filters on the fly
        - Defining a new filter
        - Temporarily modifying an existing filter
    3. Adding filters permanently

Section 1. Available response functions
=======================================

Short list of available systems:

>>> responses = list_response()
>>> systems = [response.split('.')[0] for response in responses]
>>> set_responses = sorted(set([response.split('.')[0] for response in systems]))
>>> for i,resp in enumerate(set_responses):
...    print '%10s (%3d filters)'%(resp,systems.count(resp))
     2MASS (  3 filters)
    ACSHRC ( 17 filters)
    ACSSBC (  6 filters)
    ACSWFC ( 12 filters)
     AKARI ( 13 filters)
       ANS (  6 filters)
      APEX (  1 filters)
     ARGUE (  3 filters)
    BESSEL (  6 filters)
   BESSELL (  6 filters)
     COROT (  2 filters)
   COUSINS (  3 filters)
       DDO (  7 filters)
     DENIS (  3 filters)
     DIRBE ( 10 filters)
   EEV4280 (  1 filters)
     ESOIR ( 10 filters)
      GAIA (  4 filters)
     GALEX (  2 filters)
    GENEVA (  7 filters)
 HIPPARCOS (  1 filters)
     IPHAS (  3 filters)
      IRAC (  4 filters)
      IRAS (  4 filters)
    ISOCAM ( 21 filters)
   JOHNSON ( 25 filters)
    KEPLER ( 43 filters)
      KRON (  2 filters)
   LANDOLT (  6 filters)
      MIPS (  3 filters)
      MOST (  1 filters)
       MSX (  6 filters)
    NARROW (  1 filters)
    NICMOS (  6 filters)
      OAO2 ( 12 filters)
      PACS (  3 filters)
      SAAO ( 13 filters)
     SCUBA (  6 filters)
      SDSS ( 10 filters)
     SLOAN (  2 filters)
     SPIRE (  3 filters)
  STEBBINS (  6 filters)
   STISCCD (  2 filters)
   STISFUV (  4 filters)
   STISNUV (  7 filters)
 STROMGREN (  6 filters)
       TD1 (  4 filters)
     TYCHO (  2 filters)
    TYCHO2 (  2 filters)
  ULTRACAM (  5 filters)
    USNOB1 (  2 filters)
      UVEX (  5 filters)
   VILNIUS (  7 filters)
     VISIR ( 13 filters)
  WALRAVEN (  5 filters)
     WFCAM (  5 filters)
     WFPC2 ( 21 filters)
      WISE (  4 filters)
      WOOD ( 12 filters)

Plots of all passbands of all systems:

]include figure]]ivs_sed_filters_2MASS.png]

]include figure]]ivs_sed_filters_ACSHRC.png]

]include figure]]ivs_sed_filters_ACSSBC.png]

]include figure]]ivs_sed_filters_ACSWFC.png]

]include figure]]ivs_sed_filters_AKARI.png]

]include figure]]ivs_sed_filters_ANS.png]

]include figure]]ivs_sed_filters_APEX.png]

]include figure]]ivs_sed_filters_ARGUE.png]

]include figure]]ivs_sed_filters_BESSEL.png]

]include figure]]ivs_sed_filters_BESSELL.png]

]include figure]]ivs_sed_filters_COROT.png]

]include figure]]ivs_sed_filters_COUSINS.png]

]include figure]]ivs_sed_filters_DDO.png]

]include figure]]ivs_sed_filters_DENIS.png]

]include figure]]ivs_sed_filters_DIRBE.png]

]include figure]]ivs_sed_filters_ESOIR.png]

]include figure]]ivs_sed_filters_EEV4280.png]

]include figure]]ivs_sed_filters_GAIA.png]

]include figure]]ivs_sed_filters_GALEX.png]

]include figure]]ivs_sed_filters_GENEVA.png]

]include figure]]ivs_sed_filters_HIPPARCOS.png]

]include figure]]ivs_sed_filters_IPHAS.png]

]include figure]]ivs_sed_filters_IRAC.png]

]include figure]]ivs_sed_filters_IRAS.png]

]include figure]]ivs_sed_filters_ISOCAM.png]

]include figure]]ivs_sed_filters_JOHNSON.png]

]include figure]]ivs_sed_filters_KEPLER.png]

]include figure]]ivs_sed_filters_KRON.png]

]include figure]]ivs_sed_filters_LANDOLT.png]

]include figure]]ivs_sed_filters_MIPS.png]

]include figure]]ivs_sed_filters_MOST.png]

]include figure]]ivs_sed_filters_MSX.png]

]include figure]]ivs_sed_filters_NARROW.png]

]include figure]]ivs_sed_filters_NICMOS.png]

]include figure]]ivs_sed_filters_OAO2.png]

]include figure]]ivs_sed_filters_PACS.png]

]include figure]]ivs_sed_filters_SAAO.png]

]include figure]]ivs_sed_filters_SCUBA.png]

]include figure]]ivs_sed_filters_SDSS.png]

]include figure]]ivs_sed_filters_SLOAN.png]

]include figure]]ivs_sed_filters_SPIRE.png]

]include figure]]ivs_sed_filters_STEBBINS.png]

]include figure]]ivs_sed_filters_STISCCD.png]

]include figure]]ivs_sed_filters_STISFUV.png]

]include figure]]ivs_sed_filters_STISNUV.png]

]include figure]]ivs_sed_filters_STROMGREN.png]

]include figure]]ivs_sed_filters_TD1.png]

]include figure]]ivs_sed_filters_TYCHO.png]

]include figure]]ivs_sed_filters_TYCHO2.png]

]include figure]]ivs_sed_filters_ULTRACAM.png]

]include figure]]ivs_sed_filters_USNOB1.png]

]include figure]]ivs_sed_filters_UVEX.png]

]include figure]]ivs_sed_filters_VILNIUS.png]

]include figure]]ivs_sed_filters_VISIR.png]

]include figure]]ivs_sed_filters_WALRAVEN.png]

]include figure]]ivs_sed_filters_WFPC2.png]

]include figure]]ivs_sed_filters_WISE.png]

]include figure]]ivs_sed_filters_WOOD.png]

Section 2: Adding filters on the fly
====================================

Section 2.1: Defining a new filter
----------------------------------

You can add custom filters on the fly using L{add_custom_filter}. In this
example we add a weird-looking filter and check the definition of Flambda and
Fnu and its relation to the effective wavelength of a passband:

Prerequisites: some modules that come in handy:

>>> from ivs.sigproc import funclib
>>> from ivs.sed import model
>>> from ivs.units import conversions

First, we'll define a double peakd Gaussian profile on the wavelength grid of
the WISE.W3 response curve:

>>> wave = get_response('WISE.W3')[0]
>>> trans = funclib.evaluate('gauss',wave,[1.5,76000.,10000.,0.])
>>> trans+= funclib.evaluate('gauss',wave,[1.0,160000.,25000.,0.])

This is what it looks like:

>>> p = pl.figure()
>>> p = pl.plot(wave/1e4,trans,'k-')
>>> p = pl.xlabel("Wavelength [micron]")
>>> p = pl.ylabel("Transmission [arbitrary units]")

]include figure]]ivs_sed_filters_weird_trans.png]

We can add this filter to the list of predefined filters in the following way
(for the doctests to work, we have to do a little work around and call
filters via that module, this is not needed in a normal workflow):

>>> model.filters.add_custom_filter(wave,trans,photband='LAMBDA.CCD',type='CCD')
>>> model.filters.add_custom_filter(wave,trans,photband='LAMBDA.BOL',type='BOL')

Note that we add the filter twice, once assuming that it is mounted on a
bolometer, and once on a CCD device. We'll call the filter C{LAMBDA.CCD} and
C{LAMBDA.BOL}. From now on, they are available within functions as L{get_info}
and L{get_response}. For example, what is the effective (actually pivot)
wavelength?

>>> effwave_ccd = model.filters.eff_wave('LAMBDA.CCD')
>>> effwave_bol = model.filters.eff_wave('LAMBDA.BOL')
>>> print(effwave_ccd,effwave_bol)
(119263.54911400242, 102544.27931275869)

Let's do some synthetic photometry now. Suppose we have a black body atmosphere:

>>> bb = model.blackbody(wave,5777.)

We now calculate the synthetic flux, assuming the CCD and BOL device. We
compute the synthetic flux both in Flambda and Fnu:

>>> flam_ccd,flam_bol = model.synthetic_flux(wave,bb,['LAMBDA.CCD','LAMBDA.BOL'])
>>> fnu_ccd,fnu_bol = model.synthetic_flux(wave,bb,['LAMBDA.CCD','LAMBDA.BOL'],units=['FNU','FNU'])

You can see that the fluxes can be quite different when you assume photon or
energy counting devices!

>>> flam_ccd,flam_bol
(897.68536911320564, 1495.248213834755)
>>> fnu_ccd,fnu_bol
(4.2591095543803019e-06, 5.2446332430111098e-06)

Can we now readily convert Flambda to Fnu with assuming the pivot wavelength?

>>> fnu_fromflam_ccd = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',flam_ccd,wave=(effwave_ccd,'A'))
>>> fnu_fromflam_bol = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',flam_bol,wave=(effwave_bol,'A'))

Which is equivalent with:

>>> fnu_fromflam_ccd = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',flam_ccd,photband='LAMBDA.CCD')
>>> fnu_fromflam_bol = conversions.convert('erg/s/cm2/AA','erg/s/cm2/Hz',flam_bol,photband='LAMBDA.BOL')

Apparently, with the definition of pivot wavelength, you can safely convert from
Fnu to Flambda:

>>> print(fnu_ccd,fnu_fromflam_ccd)
(4.2591095543803019e-06, 4.259110447428463e-06)
>>> print(fnu_bol,fnu_fromflam_bol)
(5.2446332430111098e-06, 5.2446373530017525e-06)

The slight difference you see is numerical.

Section 2.2: Temporarily modifying an existing filter
-----------------------------------------------------

Under usual conditions, you are prohibited from overwriting an existing predefined
response curve. That is, if you try to L{add_custom_filter} with a C{photband}
that already exists as a file, a C{ValueError} will be raised (this is not the
case for a custom defined filter, which you can overwrite without problems!).
If, for testing purposes, you want to use another definition of a predefined
response curve, you need to set C{force=True} in L{add_custom_filter}, and then
call

>>> set_prefer_file(False)

To reset and use the original definitions again, do

>>> set_prefer_file(True)

Section 3.: Adding filters permanently
--------------------------------------

Add a new response curve file to the ivs/sed/filters directory. The file should
contain two columns, the first column is the wavelength in angstrom, the second
column is the transmission curve. The units of the later are not important.

Then, call L{update_info}. The contents of C{zeropoints.dat} will automatically
be updated. Make sure to add any additional information on the new filters
manually in that file (e.g. is t a CCD or bolometer, what are the zeropoint
magnitudes etc).

"""
import os
import glob
import pyfits
import logging
import numpy as np

from ivs import config
from ivs.aux.decorators import memoized
from ivs.aux import decorators
from ivs.aux import loggers
from ivs.io import ascii

basedir = os.path.dirname(__file__)

logger = logging.getLogger("SED.FILT")
logger.addHandler(loggers.NullHandler())

custom_filters = {'_prefer_file':True}

#{ response curves
@memoized
def get_response(photband):
    """
    Retrieve the response curve of a photometric system 'SYSTEM.FILTER'
    
    OPEN.BOL represents a bolometric open filter.
    
    Example usage:
    
    >>> p = pl.figure()
    >>> for band in ['J','H','KS']:
    ...    p = pl.plot(*get_response('2MASS.%s'%(band)))
    
    If you defined a custom filter with the same name as an existing one and
    you want to use that one in the future, set C{prefer_file=False} in the
    C{custom_filters} module dictionary.
    
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @return: (wavelength [A], response)
    @rtype: (array, array)
    """
    photband = photband.upper()
    prefer_file = custom_filters['_prefer_file']
    if photband=='OPEN.BOL':
        return np.array([1,1e10]),np.array([1/(1e10-1),1/(1e10-1)])    
    #-- either get from file or get from dictionary
    photfile = os.path.join(basedir,'filters',photband)
    photfile_is_file = os.path.isfile(photfile)
    #-- if the file exists and files have preference
    if photfile_is_file and prefer_file:
        wave, response = ascii.read2array(photfile).T[:2]
    #-- if the custom_filter exist
    elif photband in custom_filters:
        wave, response = custom_filters[photband]['response']
    #-- if the file exists but custom filters have preference
    elif photfile_is_file:
        wave, response = ascii.read2array(photfile).T[:2]
    else:
        raise IOError,('{0} does not exist {1}'.format(photband,custom_filters.keys()))
    sa = np.argsort(wave)
    return wave[sa],response[sa]

def create_custom_filter(wave,peaks,range=(3000,4000),sigma=3.):
    """
    Create a custom filter as a sum of Gaussians.
    
    @param wave: wavelength to evaluate the profile on
    @type wave: ndarray
    @param peaks: heights of the peaks
    @type peaks: ndarray of length N, with N peaks
    @param range: wavelength range of the peaks
    @type range: tuple
    @param sigma: width of the peaks in units of (range/N)
    @type sigma: float
    @return: filter profile
    @rtype: ndarray
    """
    wpnts = np.linspace(range[0],range[1],len(peaks)+2)[1:-1]
    sigma = (range[1]-range[0])/(sigma*len(peaks))
    gauss = lambda x,p: p[0] * np.exp( -(x-p[1])**2 / (2.0*p[2]**2))
    els = [gauss(wave,[pk,mu,sigma]) for pk,mu in zip(peaks,wpnts)]
    profile = np.array(els).sum(axis=0)
    return profile


def add_custom_filter(wave,response,**kwargs):
    """
    Add a custom filter to the set of predefined filters.
    
    Extra keywords are:
        'eff_wave', 'type',
        'vegamag', 'vegamag_lit',
        'ABmag', 'ABmag_lit',
        'STmag', 'STmag_lit',
        'Flam0', 'Flam0_units', 'Flam0_lit',
        'Fnu0', 'Fnu0_units', 'Fnu0_lit',
        'source'
        
    default C{type} is 'CCD'.
    default C{photband} is 'CUSTOM.FILTER'
    
    @param wave: wavelength (angstrom)
    @type wave: ndarray
    @param response: response
    @type response: ndarray
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    """
    kwargs.setdefault('photband','CUSTOM.FILTER')
    kwargs.setdefault('copy_from','JOHNSON.V')
    kwargs.setdefault('force',False)
    photband = kwargs['photband']
    #-- check if the filter already exists:
    photfile = os.path.join(basedir,'filters',photband)
    if os.path.isfile(photfile) and not kwargs['force']:
        raise ValueError,'bandpass {0} already exists'.format(photfile)
    elif photband in custom_filters:
        logger.debug('Overwriting previous definition of {0}'.format(photband))
    custom_filters[photband] = dict(response=(wave,response))
    #-- set effective wavelength
    kwargs.setdefault('type','CCD')
    kwargs.setdefault('eff_wave',eff_wave(photband,det_type=kwargs['type']))
    #-- add info for zeropoints.dat file: make sure wherever "lit" is part of
    #   the name, we replace it with "0". Then, we overwrite any existing
    #   information with info given
    myrow = get_info([kwargs['copy_from']])
    for name in myrow.dtype.names:
        if 'lit' in name:
            myrow[name] = 0
        myrow[name] = kwargs.pop(name,myrow[name])
    del decorators.memory[__name__]
    #-- add info:
    custom_filters[photband]['zp'] = myrow
    logger.debug('Added photband {0} to the predefined set'.format(photband))

def set_prefer_file(prefer_file=True):
    """
    Set whether to prefer files or custom filters when both exist.
    
    @param prefer_file: boolean
    @type prefer_file: bool
    """
    custom_filters['_prefer_file'] = prefer_file
    logger.info("Prefering {}".format(prefer_file and 'files' or 'custom filters'))


def add_spectrophotometric_filters(R=200.,lambda0=950.,lambdan=3350.):
    #-- STEP 1: define wavelength bins
    Delta = np.log10(1.+1./R)
    x = np.arange(np.log10(lambda0),np.log10(lambdan)+Delta,Delta)
    x = 10**x
    photbands = []
    for i in range(len(x)-1):
        wave = np.linspace(x[i],x[i+1],100)
        resp = np.ones_like(wave)
        dw = wave[1]-wave[0]
        wave = np.hstack([wave[0]-dw,wave,wave[-1]+dw])
        resp = np.hstack([0,resp,0])
        photband = 'BOXCAR_R{0:d}.{1:d}'.format(int(R),int(x[i]))
        try:
            add_custom_filter(wave,resp,photband=photband)
        except ValueError:
            logger.info('{0} already exists, skipping'.format(photband))
        photbands.append(photband)
    logger.info('Added spectrophotometric filters')
    return photbands


def list_response(name='*',wave_range=(-np.inf,+np.inf)):
    """
    List available response curves.
    
    Specify a glob string C{name} and/or a wavelength range to make a selection
    of all available curves. If nothing is supplied, all curves will be returned.
    
    @param name: list all curves containing this string
    @type name: str 
    @param wave_range: list all curves within this wavelength range (A)
    @type wave_range: (float, float)
    @return: list of curve files
    @rtype: list of str
    """
    #-- collect all curve files but remove human eye responses
    if not '*' in name:
        name_ = '*' + name + '*'
    else:
        name_ = name
    curve_files = sorted(glob.glob(os.path.join(basedir,'filters',name_.upper())))
    curve_files = sorted(curve_files+[key for key in custom_filters.keys() if ((name in key) and not (key=='_prefer_file'))])
    curve_files = [cf for cf in curve_files if not ('HUMAN' in cf or 'EYE' in cf) ]
    #-- select in correct wavelength range
    curve_files = [os.path.basename(curve_file) for curve_file in curve_files if (wave_range[0]<=eff_wave(os.path.basename(curve_file))<=wave_range[1])]
    #-- log to the screen and return
    for curve_file in curve_files: logger.info(curve_file)
    return curve_files



def is_color(photband):
    """
    Return true if the photometric passband is actually a color.
    
    @param photband: name of the photometric passband
    @type photband: string
    @return: True or False
    @rtype: bool
    """
    if '-' in photband.split('.')[1]:
        return True
    elif photband.split('.')[1].upper() in ['M1','C1']:
        return True
    else:
        return False


def get_color_photband(photband):
    """
    Retrieve the photometric bands from color
    
    @param photband: name of the photometric passband
    @type photband: string
    @return: tuple of strings
    @rtype: tuple
    """
    system,band = photband.split('.')
    band = band.strip() # remove extra spaces
    if '-' in band:
        bands = tuple(['%s.%s'%(system,iband) for iband in band.split('-')])
    elif band.upper()=='M1':
        bands = tuple(['%s.%s'%(system,iband) for iband in ['V','B','Y']])
    elif band.upper()=='C1':
        bands = tuple(['%s.%s'%(system,iband) for iband in ['V','B','U']])
    else:
        raise ValueError('cannot recognize color {}'.format(photband))
    return bands

def make_color(photband):
    """
    Make a color from a color name and fluxes.
    
    You get two things: a list of photbands that are need to construct the color,
    and a function which you need to pass fluxes to compute the color.
    
    >>> bands, func = make_color('JOHNSON.B-V')
    >>> print(bands)
    ('JOHNSON.B', 'JOHNSON.V')
    >>> print(func(2,3.))
    0.666666666667
    
    @return: photbands, function to construct color
    @rtype: tuple,callable
    """
    system,band = photband.split('.')
    band = band.strip() # remove extra spaces
    photbands = get_color_photband(photband)
    if len(band.split('-'))==2:
        func = lambda f0,f1: f0/f1
    elif band=='M1':
        func = lambda fv,fb,fy: fv*fy/fb**2
    elif band=='C1':
        func = lambda fv,fb,fu: fu*fb/fv**2
    else:
        raise ValueError('cannot recognize color {}'.format(photband))
    return photbands,func


def eff_wave(photband,model=None,det_type=None):
    """
    Return the effective wavelength of a photometric passband.
    
    The effective wavelength is defined as the average wavelength weighed with
    the response curve.
    
    >>> eff_wave('2MASS.J')
    12393.093155655277
    
    If you give model fluxes as an extra argument, the wavelengths will take
    these into account to calculate the `true' effective wavelength (e.g., 
    Van Der Bliek, 1996), eq 2.
    
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER') or array/list of str
    @param model: model wavelength and fluxes
    @type model: tuple of 1D arrays (wave,flux)
    @return: effective wavelength [A]
    @rtype: float or numpy array
    """
    
    #-- if photband is a string, it's the name of a photband: put it in a container
    #   but unwrap afterwards
    if isinstance(photband,unicode):
        photband = str(photband)
    if isinstance(photband,str):
        single_band = True
        photband = [photband]
    #-- else, it is a container
    else:
        single_band = False
        
    my_eff_wave = []
    for iphotband in photband:
        try:
            wave,response = get_response(iphotband)
            #-- bolometric or ccd?
            if det_type is None and len(get_info([iphotband])):
                det_type = get_info([iphotband])['type'][0]
            elif det_type is None:
                det_type = 'CCD'
            if model is None:
                #this_eff_wave = np.average(wave,weights=response)
                if det_type=='BOL':
                    this_eff_wave = np.sqrt(np.trapz(response,x=wave)/np.trapz(response/wave**2,x=wave))
                else:
                    this_eff_wave = np.sqrt(np.trapz(wave*response,x=wave)/np.trapz(response/wave,x=wave))
            else:
                #-- interpolate response curve onto higher resolution model and
                #   take weighted average
                is_response = response>1e-10
                start_response,end_response = wave[is_response].min(),wave[is_response].max()
                fluxm = np.sqrt(10**np.interp(np.log10(wave),np.log10(model[0]),np.log10(model[1])))
                
                if det_type=='CCD':
                    this_eff_wave = np.sqrt(np.trapz(wave*fluxm*response,x=wave) / np.trapz(fluxm*response/wave,x=wave))
                elif det_type=='BOL':
                    this_eff_wave = np.sqrt(np.trapz(fluxm*response,x=wave) / np.trapz(fluxm*response/wave**2,x=wave))     
        #-- if the photband is not defined:
        except IOError:
            this_eff_wave = np.nan
        my_eff_wave.append(this_eff_wave)
    
    if single_band:
        my_eff_wave = my_eff_wave[0]
    else:
        my_eff_wave = np.array(my_eff_wave,float)
    
    return my_eff_wave

@memoized
def get_info(photbands=None):
    """
    Return a record array containing all filter information.
    
    The record arrays contains following columns:
        - photband
        - eff_wave
        - type
        - vegamag, vegamag_lit
        - ABmag, ABmag_lit
        - STmag, STmag_lit
        - Flam0, Flam0_units, Flam0_lit
        - Fnu0, Fnu0_units, Fnu0_lit,
        - source
    
    @param photbands: list of photbands to get the information from. The input
    order is equal to the output order. If C{None}, all filters are returned.
    @type photbands: iterable container (list, tuple, 1Darray)
    @return: record array containing all information on the requested photbands.
    @rtype: record array
    """
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp = ascii.read2recarray(zp_file)
    for iph in custom_filters:
        if iph=='_prefer_file': continue
        if 'zp' in custom_filters[iph]:
            zp = np.hstack([zp,custom_filters[iph]['zp']])
    zp = zp[np.argsort(zp['photband'])]
    
    #-- list photbands in order given, and remove those that do not have
    #   zeropoints etc.
    if photbands is not None:
        order = np.searchsorted(zp['photband'],photbands)
        zp = zp[order]
        keep = (zp['photband']==photbands)
        zp = zp[keep]
    
    return zp





def update_info(zp=None):
    """
    Update information in zeropoint file, e.g. after calibration.
    
    Call first L{ivs.sed.model.calibrate} without arguments, and pass the output
    to this function.
    
    @param zp: updated contents from C{zeropoints.dat}
    @type zp: recarray
    """
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp_,comms = ascii.read2recarray(zp_file,return_comments=True)
    existing = [str(i.strip()) for i in zp_['photband']]
    resp_files = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),'filters/*')))
    resp_files = [os.path.basename(ff) for ff in resp_files if not os.path.basename(ff) in existing]
    resp_files.remove('HUMAN.EYE')
    resp_files.remove('HUMAN.CONES')
    resp_files.remove('CONES.EYE')
    if zp is None:
        zp = zp_
        logger.info('No new calibrations; previous information on existing response curves is copied')
    else:
        logger.info('Received new calibrations contents of zeropoints.dat will be updated')
    
    #-- update info on previously non existing response curves
    new_zp = np.zeros(len(resp_files),dtype=zp.dtype)
    logger.info('Found {} new response curves, adding them with default information'.format(len(resp_files)))
    for i,respfile in enumerate(resp_files):
        new_zp[i]['photband'] = respfile
        new_zp[i]['eff_wave'] = float(eff_wave(respfile))
        new_zp[i]['type'] = 'CCD'
        new_zp[i]['vegamag'] = np.nan
        new_zp[i]['ABmag'] = np.nan
        new_zp[i]['STmag'] = np.nan
        new_zp[i]['Flam0_units'] = 'erg/s/cm2/AA'
        new_zp[i]['Fnu0_units'] = 'erg/s/cm2/AA'
        new_zp[i]['source'] = 'nan'
    zp = np.hstack([zp,new_zp])
    sa = np.argsort(zp['photband'])
    ascii.write_array(zp[sa],'zeropoints.dat',header=True,auto_width=True,comments=['#'+line for line in comms[:-2]],use_float='%g')
    


if __name__=="__main__":
    import sys
    import pylab as pl
    if not sys.argv[1:]:
        import doctest
        doctest.testmod()
        pl.show()
    
    else:
        import itertools
        responses = list_response()
        systems = [response.split('.')[0] for response in responses]
        set_responses = sorted(set([response.split('.')[0] for response in systems]))
        this_filter = 0
        for i,resp in enumerate(responses):
            # what system is this, and how many filters are in this system?
            this_system = resp.split('.')[0]
            nr_filters = systems.count(this_system)
            # call the plot containing the filters of the same system. If this is the
            # the first time the plot is called (first filter of system), then set
            # the title and color cycle
            p = pl.figure(set_responses.index(this_system),figsize=(10,4.5))
            if not hasattr(pl.gca(),'color_cycle'):
                color_cycle = itertools.cycle([pl.cm.spectral(j) for j in np.linspace(0, 1.0, nr_filters)])
                p = pl.gca().color_cycle = color_cycle
            color = pl.gca().color_cycle.next()
            p = pl.title(resp.split('.')[0])
            # get the response curve and plot it
            wave,trans = get_response(resp)
            p = pl.plot(wave/1e4,trans,label=resp,color=color)
            # and set labels
            p = pl.xlabel('Wavelength [micron]')
            p = pl.ylabel('Transmission')
            # if there are not more filters in this systems, save the plot to a file
            # and close it
            this_filter+=1
            if this_filter==nr_filters:
                this_filter = 0
                p = pl.legend(prop=dict(size='small'))
                p = pl.savefig('/home/pieterd/python/ivs/doc/images/ivs_sed_filters_%s'%(this_system));p = pl.close()
