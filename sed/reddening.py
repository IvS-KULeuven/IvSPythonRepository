# -*- coding: utf-8 -*-
"""
Definitions of interstellar reddening curves

Example usage:
    
>>> import pylab as pl

Use the general interface to get different curves:

>>> wave = np.r_[1e3:1e5:10]
>>> for name in ['chiar2006','fitzpatrick1999','fitzpatrick2004','cardelli1989','seaton1979']:
...   wave_,mag_ = get_law(name,wave=wave)
...   p = pl.plot(1e4/wave_,mag_,label=name)
>>> p = pl.xlim(0,10)
>>> p = pl.ylim(0,12)
>>> p = pl.legend()

Use the general interface to get the same curves but with different Rv:

>>> for Rv in [2.0,3.1,5.1]:
...   wave_,mag_ = get_law('cardelli1989',wave=wave,Rv=Rv)
...   p = pl.plot(1e4/wave_,mag_,'--',lw=2)
>>> p = pl.xlim(0,10)
>>> p = pl.ylim(0,12)

Get the curves seperately:

>>> wave1,mag1 = cardelli1989()
>>> wave2,mag2 = chiar2006()
>>> wave3,mag3 = seaton1979()
>>> wave4,mag4 = fitzpatrick1999()
>>> wave5,mag5 = fitzpatrick2004()

And plot them:

>>> p = pl.figure()
>>> p = pl.plot(1e4/wave1,mag1*3.1)
>>> p = pl.plot(1e4/wave2,mag2*3.1)
>>> p = pl.plot(1e4/wave3,mag3*3.1)
>>> p = pl.plot(1e4/wave4,mag4*3.1)
>>> p = pl.plot(1e4/wave5,mag5*3.1)
>>> p = pl.xlim(0,10)
>>> p = pl.ylim(0,12)
>>> p = pl.show()
"""

import os
import numpy as np
import logging

from ivs.io import ascii
from ivs.misc import loggers
from ivs.misc.decorators import memoized
from ivs.units import conversions
from ivs.sed import filters

logger = logging.getLogger("SED.RED")
logger.addHandler(loggers.NullHandler)

basename = os.path.join(os.path.dirname(__file__),'redlaws')

#{ Main interface

def get_law(name,norm='E(B-V)',wave_units='A',**kwargs):
    """
    Retrieve an interstellar reddening law.
    
    Parameter C{name} must be the function name of one of the laws defined in
    this module.
    
    By default, the law will be interpolated on a grid from 100 angstrom to
    10 micron in steps of 10 angstrom. This can be adjusted with the parameter
    C{wave} (array), which B{must} be in angstrom. You can change the units
    ouf the returned wavelength array via C{wave_units}.
    
    By default, the curve is normalised with respect to E(B-V) (you get
    A(l)/E(B-V)). You can set the C{norm} keyword to Av if you want A(l)/E(B-V).
    Remember that
    
    A(V) = Rv * E(B-V)
    
    The parameter C{Rv} is by default 3.1, other reasonable values lie between
    2.0 and 5.1
    
    Extra accepted keywords depend on the type of reddening law used.
    
    Example usage:
    
    >>> wave = np.r_[1e3:1e5:10]
    >>> wave,mag = get_law('cardelli1989',wave=wave,Rv=3.1)
    
    @param name: name of the interstellar law
    @type name: str, one of the functions defined here
    @param norm: type of normalisation of the curve
    @type norm: str (one of E(B-V), Av)
    @param wave_units: wavelength units
    @type wave_units: str (interpretable for units.conversions.convert)
    @keyword wave: wavelength array to interpolate the law on
    @type wave: ndarray
    @return: wavelength, reddening magnitude
    @rtype: (ndarray,ndarray)
    """
    #-- get the inputs
    Rv = kwargs.setdefault('Rv',3.1)
    
    #-- get the curve
    wave,mag = globals()[name.lower()](**kwargs)
    
    #-- interpolate on user defined grid
    wave_ = kwargs.get('wave',None)
    if wave_ is not None:
        mag = np.interp(wave_,wave,mag,right=0)
        wave = wave_
    
    #-- convert to A(lambda)/Av if needed
    if norm.lower()=='e(b-v)':
        mag *= Rv
    elif norm.lower()!='av':
        raise ValueError, "do not understand normalization %s for reddening"%(norm)
    
    #-- set the units of the wavelengths
    if wave_units != 'A':
        wave = conversions.convert('A',wave_units,wave)
    
    return wave,mag


def redden(flux,wave=None,photbands=None,ebv=0.,rtype='flux',law='cardelli1989',**kwargs):
    """
    Redden flux or magnitudes
    
    The reddening parameters C{ebv} means E(B-V).
    
    If it is negative, we B{deredden}.
    
    If you give the keyword C{wave}, it is assumed that you want to (de)redden
    a B{model}, i.e. a spectral energy distribution.
    
    If you give the keyword C{photbands}, it is assumed that you want to (de)redden
    B{photometry}, i.e. integrated fluxes.
    
    @param flux: fluxes to (de)redden (magnitudes if C{rtype='mag'})
    @type flux: ndarray (floats)
    @param wave: wavelengths matching the fluxes (or give C{photbands})
    @type wave: ndarray (floats)
    @param photbands: photometry bands matching the fluxes (or give C{wave})
    @type photbands: ndarray of str
    @param ebv: reddening parameter E(B-V)
    @type ebv: float
    @param rtype: type of dereddening (magnituds or fluxes)
    @type rtype: str ('flux' or 'mag')
    @return: (de)reddened flux/magnitude
    @rtype: ndarray (floats)
    """
    if photbands is not None:
        wave = filters.get_info(photbands)['eff_wave']
        
    wave,mag = get_law(law,wave=wave,**kwargs)
    if rtype=='flux':
        flux_dered = flux / 10**(mag*ebv/2.5)
    elif rtype=='mag':
        flux_dered = flux - mag*ebv
    return flux_dered

def deredden(flux,wave=None,photbands=None,ebv=0.,rtype='flux',**kwargs):
    """
    Deredden flux or magnitudes.
    
    @param flux: fluxes to (de)redden (NOT magnitudes)
    @type flux: ndarray (floats)
    @param wave: wavelengths matching the fluxes (or give C{photbands})
    @type wave: ndarray (floats)
    @param photbands: photometry bands matching the fluxes (or give C{wave})
    @type photbands: ndarray of str
    @param ebv: reddening parameter E(B-V)
    @type ebv: float
    @param rtype: type of dereddening (magnituds or fluxes)
    @type rtype: str ('flux' or 'mag')
    @return: (de)reddened flux
    @rtype: ndarray (floats)
    """
    return redden(flux,wave=wave,photbands=photbands,ebv=-ebv,rtype=rtype,**kwargs)
    

#}

#{ Curve definitions
@memoized
def chiar2006(Rv=3.1,curve='ism',**kwargs):
    """
    Extinction curve at infrared wavelengths from Chiar and Tielens (2006)
    
    We return A(lambda)/E(B-V), by multiplying A(lambda)/Av with Rv.
    
    This is only defined for Rv=3.1. If it is different, this will raise an
    AssertionError
    
    Extra kwags are to catch unwanted keyword arguments.
    
    UNCERTAIN NORMALISATION
    
    @param Rv: Rv
    @type Rv: float
    @param curve: extinction curve
    @type curve: string (one of 'gc' or 'ism', galactic centre or local ISM)
    @return: wavelengths (A), A(lambda)/Av
    @rtype: (ndarray,ndarray)
    """
    source = os.path.join(basename,'Chiar2006.red')
    
    #-- check Rv
    assert(Rv==3.1)
    wavelengths,gc,ism = ascii.read2array(source).T
    if curve=='gc':
        alam_ak = gc
    elif curve=='ism':
        keep = ism>0
        alam_ak = ism[keep]
        wavelengths = wavelengths[keep]
    else:
        raise ValueError,'no curve %s'%(curve)
    alam_aV = alam_ak * 0.09
    #plot(1/wavelengths,alam_aV,'o-')
    return wavelengths*1e4,alam_aV



@memoized
def fitzpatrick1999(Rv=3.1,**kwargs):
    """
    From Fitzpatrick 1999 (downloaded from ASAGIO database)
    
    This function returns A(lambda)/A(V).
    
    To get A(lambda)/E(B-V), multiply the return value with Rv (A(V)=Rv*E(B-V))
    
    Extra kwags are to catch unwanted keyword arguments.
    
    @param Rv: Rv (2.1, 3.1 or 5.0)
    @type Rv: float
    @return: wavelengths (A), A(lambda)/Av
    @rtype: (ndarray,ndarray)
    """
    filename = 'Fitzpatrick1999_Rv_%.1f'%(Rv)
    filename = filename.replace('.','_') + '.red'
    myfile = os.path.join(basename,filename)
    wave,alam_ebv = ascii.read2array(myfile).T
    alam_av = alam_ebv/Rv
    
    logger.info('Fitzpatrick1999 curve with Rv=%.2f'%(Rv))
    
    return wave,alam_av

@memoized
def fitzpatrick2004(Rv=3.1,**kwargs):
    """
    From Fitzpatrick 2004 (downloaded from FTP)
    
    This function returns A(lambda)/A(V).
    
    To get A(lambda)/E(B-V), multiply the return value with Rv (A(V)=Rv*E(B-V))
    
    Extra kwags are to catch unwanted keyword arguments.
    
    @param Rv: Rv (2.1, 3.1 or 5.0)
    @type Rv: float
    @return: wavelengths (A), A(lambda)/Av
    @rtype: (ndarray,ndarray)
    """
    filename = 'Fitzpatrick2004_Rv_%.1f.red'%(Rv)
    myfile = os.path.join(basename,filename)
    wave_inv,elamv_ebv = ascii.read2array(myfile,skip_lines=15).T
    
    logger.info('Fitzpatrick2004 curve with Rv=%.2f'%(Rv))
    
    return 1e4/wave_inv[::-1],((elamv_ebv+Rv)/Rv)[::-1]


@memoized
def donnell1994(**kwargs):
    """
    Small improvement on Cardelli 1989 by James E. O'Donnell (1994).
    
    Extra kwags are to catch unwanted keyword arguments.
    
    @keyword Rv: Rv
    @type Rv: float
    @keyword wave: wavelengths to compute the curve on
    @type wave: ndarray
    @return: wavelengths (A), A(lambda)/Av
    @rtype: (ndarray,ndarray)
    """
    return cardelli1989(curve='donnell',**kwargs)


@memoized
def cardelli1989(Rv=3.1,curve='cardelli',wave=None,**kwargs):
    """
    Construct extinction laws from Cardelli (1989).
    
    Improvement in optical by James E. O'Donnell (1994)
    
    wavelengths in Angstrom!
    
    This function returns A(lambda)/A(V).
    
    To get A(lambda)/E(B-V), multiply the return value with Rv (A(V)=Rv*E(B-V))
    
    Extra kwags are to catch unwanted keyword arguments.
    
    @param Rv: Rv
    @type Rv: float
    @param curve: extinction curve
    @type curve: string (one of 'cardelli' or 'donnell')
    @param wave: wavelengths to compute the curve on
    @type wave: ndarray
    @return: wavelengths (A), A(lambda)/Av
    @rtype: (ndarray,ndarray)
    """
    if wave is None:
        wave = np.r_[100.:100000.:10]
    
    all_x = 1./(wave/1.0e4)
    alam_aV = np.zeros_like(all_x)
    
    #-- infrared
    infrared = all_x<1.1
    x = all_x[infrared]
    ax = +0.574*x**1.61
    bx = -0.527*x**1.61
    alam_aV[infrared] = ax + bx/Rv
    
    #-- optical
    optical = (1.1<=all_x) & (all_x<3.3)
    x = all_x[optical]
    y = x-1.82
    if curve=='cardelli':
        ax = 1 + 0.17699*y    - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 \
               + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        bx =     1.41338*y    + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 \
               - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    elif curve=='donnell':
        ax = 1 + 0.104*y    - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 \
               - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8
        bx =     1.952*y    + 2.908*y**2 - 3.989*y**3 - 7.985*y**4 \
              + 11.102*y**5 + 5.491*y**6 -10.805*y**7 + 3.347*y**8
    else:
        raise ValueError,'curve %s not found'%(curve)
    alam_aV[optical] = ax + bx/Rv
    
    #-- ultraviolet
    ultraviolet = (3.3<=all_x) & (all_x<8.0)
    x = all_x[ultraviolet]
    Fax = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
    Fbx = +0.21300*(x-5.9)**2 + 0.120700*(x-5.9)**3
    Fax[x<5.9] = 0
    Fbx[x<5.9] = 0
    ax = +1.752 - 0.316*x - 0.104 / ((x-4.67)**2 + 0.341) + Fax
    bx = -3.090 + 1.825*x + 1.206 / ((x-4.62)**2 + 0.263) + Fbx
    alam_aV[ultraviolet] = ax + bx/Rv
    
    #-- far UV
    fuv = 8.0<=all_x
    x = all_x[fuv]
    ax = -1.073 - 0.628*(x-8) + 0.137*(x-8)**2 - 0.070*(x-8)**3
    bx = 13.670 + 4.257*(x-8) - 0.420*(x-8)**2 + 0.374*(x-8)**3
    alam_aV[fuv] = ax + bx/Rv
    
    logger.info('%s curve with Rv=%.2f'%(curve.title(),Rv))
    
    return wave,alam_aV


@memoized
def seaton1979(Rv=3.1,wave=None,**kwargs):
    """
    Extinction curve from Seaton, 1979.
    
    This function returns A(lambda)/A(V).
    
    To get A(lambda)/E(B-V), multiply the return value with Rv (A(V)=Rv*E(B-V))
    
    Extra kwags are to catch unwanted keyword arguments.
    
    @param Rv: Rv
    @type Rv: float
    @param wave: wavelengths to compute the curve on
    @type wave: ndarray
    @return: wavelengths (A), A(lambda)/Av
    @rtype: (ndarray,ndarray)
    """
    if wave is None: wave = np.r_[1000.:10000.:10]
    
    all_x = 1e4/(wave)
    alam_aV = np.zeros_like(all_x)
    
    #-- far infrared
    x_ = np.r_[1.0:2.8:0.1]
    X_ = np.array([1.36,1.44,1.84,2.04,2.24,2.44,2.66,2.88,3.14,3.36,3.56,3.77,3.96,4.15,4.26,4.40,4.52,4.64])
    fir = all_x<=2.7
    alam_aV[fir] = np.interp(all_x[fir][::-1],x_,X_,left=0)[::-1]
    
    #-- infrared
    infrared = (2.70<=all_x) & (all_x<3.65)
    x = all_x[infrared]
    alam_aV[infrared] = 1.56 + 1.048*x + 1.01 / ( (x-4.60)**2 + 0.280)
    
    #-- optical
    optical = (3.65<=all_x) & (all_x<7.14)
    x = all_x[optical]
    alam_aV[optical] = 2.29 + 0.848*x + 1.01 / ( (x-4.60)**2 + 0.280)
    
    #-- ultraviolet
    ultraviolet = (7.14<=all_x) & (all_x<=10)
    x = all_x[ultraviolet]
    alam_aV[ultraviolet] = 16.17 - 3.20*x + 0.2975*x**2
    
    logger.info('Seaton curve with Rv=%.2f'%(Rv))
    
    return wave,alam_aV/Rv

#}

if __name__=="__main__":
    import doctest
    doctest.testmod()