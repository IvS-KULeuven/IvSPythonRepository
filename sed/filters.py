# -*- coding: utf-8 -*-
"""
Functions relevant for photometric calibration

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
   COUSINS (  2 filters)
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

logger = logging.getLogger("CAT.VIZIER")
logger.addHandler(loggers.NullHandler())

custom_filters = {}

#{ response curves
@memoized
def get_response(photband):
    """
    Retrieve the response curve of a photometric system 'SYSTEM.FILTER'
    
    OPEN.BOL represents a bolometric open filter.
    
    Example usage:
    
    >>> from pylab import plot,show
    >>> for band in ['J','H','KS']:
    ...    p = plot(*get_response('2MASS.%s'%(band)))
    >>> p = show()
    
    @param photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @return: (wavelength [A], response)
    @rtype: (array, array)
    """
    photband = photband.upper()
    if photband=='OPEN.BOL':
        return np.array([1,1e10]),np.array([1/(1e10-1),1/(1e10-1)])    
    #-- either get from file or get from dictionary
    photfile = os.path.join(basedir,'filters',photband)
    if os.path.isfile(photfile):
        wave, response = ascii.read2array(photfile).T[:2]
    elif photband in custom_filters:
        wave, response = custom_filters[photband]['response']
    else:
        raise IOError
    sa = np.argsort(wave)
    return wave[sa],response[sa]


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
    photband = kwargs['photband']
    #-- check if the filter already exists:
    photfile = os.path.join(basedir,'filters',photband)
    if os.path.isfile(photfile):
        raise ValueError,'bandpass {0} already exists'.format(photfile)
    elif photband in custom_filters:
        logger.warning('Overwriting previous definition of {0}'.format(photband))
    custom_filters[photband] = dict(response=(wave,response))
    #-- set effective wavelength
    kwargs.setdefault('eff_wave',eff_wave(photband))
    #-- add info for zeropoints.dat file: make sure wherever "lit" is part of
    #   the name, we replace it with "0". Then, we overwrite any existing
    #   information with info given
    myrow = get_info(['JOHNSON.V'])
    for name in myrow.dtype.names:
        if 'lit' in name:
            myrow[name] = 0
        myrow[name] = kwargs.pop(name,myrow[name])
    del decorators.memory['ivs.sed.filters']
    #-- add info:
    custom_filters[photband]['zp'] = myrow
    logger.info('Added photband {0} to the predefined set'.format(photband))


def add_spectrophotometric_filters(R=200.,lambda0=950.,lambdan=3350.):
    #-- STEP 1: define wavelength bins
    Delta = np.log10(1.+1./R)
    x = np.arange(np.log10(lambda0),np.log10(lambdan)+Delta,Delta)
    x = 10**x
    photbands = []
    for i in range(len(x)-1):
        wave = np.linspace(x[i],x[i+1],100)
        resp = np.ones_like(wave)
        photband = 'BOXCAR_R{0:d}.{1:d}'.format(int(R),int(x[i]))
        try:
            add_custom_filter(wave,resp,photband=photband)
        except ValueError:
            logger.info('{0} already exists, skipping'.format(photband))
        photbands.append(photband)
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
    curve_files = sorted(curve_files+[key for key in custom_filters.keys() if name in key])
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
    
    return bands


def eff_wave(photband,model=None):
    """
    Return the effective wavelength of a photometric passband.
    
    The effective wavelength is defined as the average wavelength weighed with
    the response curve.
    
    >>> eff_wave('2MASS.J')
    12412.136241640892
    
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
            if model is None:
                #this_eff_wave = np.average(wave,weights=response)
                this_eff_wave = np.sqrt(np.trapz(wave*response,x=wave)/np.trapz(response/wave,x=wave))
            else:
                #-- interpolate response curve onto higher resolution model and
                #   take weighted average
                is_response = response>1e-10
                start_response,end_response = wave[is_response].min(),wave[is_response].max()
                fluxm = np.sqrt(10**np.interp(np.log10(wave),np.log10(model[0]),np.log10(model[1])))
                #-- bolometric or ccd?
                det_type = get_info([iphotband])['type'][0]
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





def update_info(zp):
    """
    Update information in zeropoint file, e.g. after calibration.
    
    Call first L{ivs.sed.model.calibrate} without arguments, and pass the output
    to this function.
    
    @param zp: updated contents from C{zeropoints.dat}
    @type zp: recarray
    """
    zp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'zeropoints.dat')
    zp_,comms = ascii.read2recarray(zp_file,return_comments=True)
    ascii.write_array(zp,'zeropoints.dat',header=True,auto_width=True,comments=['#'+line for line in comms[:-2]],use_float='%g')
    


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
