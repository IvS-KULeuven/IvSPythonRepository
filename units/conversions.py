# -*- coding: utf-8 -*-
"""
Convert one unit (and uncertainty) to another.

Some of the many possibilities include:
    
    1. Conversions between equal-type units: meter to nano-lightyears, erg/s
    to W, cy/d to muHz, but also erg/s/cm2/A to W/m2/mum, sr to deg2, etc...
    2. Conversions between unequal-type units: angstrom to km/s via the speed
    of light, F(lambda) to F(nu), F(nu) to lambdaF(lambda)/sr, meter to
    cycles/arcsec (interferometry), etc...
    3. Nonlinear conversions: vegamag to erg/s/cm2/A or Jy, Celcius to
    Fahrenheit or Kelvin, calender date to (any kind of modified) Julian Day, ...
    4. Conversions of magnitude to flux amplitudes via 'Amag' and 'ppt' or 'ampl'
    5. Conversions of magnitude colors to flux ratios via 'mag_color' and 'flux_ratio'
    6. Coordinate transformations (equatorial to galactic etc.)
    7. Logarithmic conversions, e.g. from logTeff to Teff via '[K]' and [K]
    8. Inclusion of uncertainties, both in input values and/or reference values
    when converting between unequal-type units, automatically recognised when
    giving two positional argument (value, error) instead of one (value).
    9. Currency exchange rates. First, you need to call L{set_exchange_rates}
    for this to work, since the latest currency definitions and rates need to
    queried from the European Central Bank (automatically done when using the
    terminal tool).

B{Warning:} frequency units are technically given in cycles per time (cy).
This means that if you want to convert e.g. muHz to d-1 (or 1/d), you need to
ask for cy/d.
    
This module can be used as:

    - a Python module
    - a standalone terminal tool
    
Section 1. The Python module
============================
    
The main function L{convert} (see link for a full list of examples) does all the
work and is called via

>>> result = convert('km','m',1.)

or when the error is known

>>> result,e_result = convert('km','m',1.,0.1)

Be B{careful} when mixing nonlinear conversions (e.g. magnitude to flux) with
linear conversions (e.g. Jy to W/m2/m).

Note: when your favorite conversion is not implemented, there are five places
where you can add information:

    1. C{_scalings}: if your favorite prefix (e.g. Tera, nano...) is not
    available
    2. C{_aliases}: if your unit is available but not under the name you are
    used to.
    3. C{_factors}: if your unit is not available.
    4. C{_switch}: if your units are available, but conversions from one to
    another is not straightforward and extra infromation is needed (e.g. to go
    from angstrom to km/s in a spectrum, a reference wavelength is needed).
    5. C{_convention}: if your favorite base unit system is not defined.

If you need to add a linear factor, just give the factor in SI units, and the
SI base units it consists of. If you need to add a nonlinear factor, you have
to give a class definition (see the examples).

Section 2. The Terminal tool
============================

For help and a list of all defined units and abbreviations, do::

    $:> python conversions.py --help
    Usage: List of available units:
    =====================================           | =====================================           
    =   Units of angle                  =           | =   Units of area                   =           
    =====================================           | =====================================           
    am              = arcminute                     | a               = are                           
    as              = arcsecond                     | ac              = acre (international)          
    cy              = cycle                         | =====================================           
    deg             = degree                        | =   Units of distance               =           
    rad             = radian                        | =====================================           
    rpm             = revolutions per minute        | A               = angstrom                      
    sr              = sterradian                    | AU              = astronomical unit             
    =====================================           | Rearth          = Earth radius                  
    =   Units of coordinate             =           | Rsol            = Solar radius                  
    =====================================           | a0              = Bohr radius                   
    complex_coord   = <own unit>                    | ell             = ell                           
    deg_coord       = degrees                       | ft              = foot (international)          
    ecl             = ecliptic                      | in              = inch (international)          
    equ             = equatorial                    | ly              = light year                    
    gal             = galactic                      | m               = meter                         
    rad_coord       = radians                       | mi              = mile (international)          
    =====================================           | pc              = parsec                        
    =   Units of energy/power           =           | yd              = yard (international)          
    =====================================           | =====================================           
    J               = Joule                         | =   Units of flux                   =           
    Lsol            = Solar luminosity              | =====================================           
    W               = Watt                          | ABmag           = AB magnitude                  
    cal             = calorie (international table) | Amag            = amplitude in magnitude        
    eV              = electron volt                 | Jy              = Jansky                        
    erg             = ergon                         | STmag           = ST magnitude                  
    =====================================           | ampl            = fractional amplitude          
    =   Units of force                  =           | flux_ratio      = flux ratio                    
    =====================================           | mag             = magnitude                     
    N               = Newton                        | mag_color       = color                         
    dyn             = dyne                          | pph             = amplitude in parts per hundred
    =====================================           | ppm             = amplitude in parts per million
    =   Units of pressure               =           | ppt             = amplitude in parts per thousand
    =====================================           | vegamag         = Vega magnitude                
    Pa              = Pascal                        | =====================================           
    at              = atmosphere (technical)        | =   Units of mass                   =           
    atm             = atmosphere (standard)         | =====================================           
    bar             = baros                         | Mearth          = Earth mass                    
    psi             = pound per square inch         | Mjup            = Jupiter mass                  
    torr            = Torricelli                    | Mlun            = Lunar mass                    
    =====================================           | Msol            = Solar mass                    
    =   Units of time                   =           | g               = gram                          
    =====================================           | lbs             = pound                         
    CD              = calender day                  | st              = stone                         
    JD              = Julian day                    | =====================================           
    MJD             = modified Julian day           | =   Units of temperature            =           
    cr              = century                       | =====================================           
    d               = day                           | C               = Celcius                       
    h               = hour                          | F               = Fahrenheit                    
    hz              = Hertz                         | K               = Kelvin                        
    j               = jiffy                         | 
    min             = minute                        | 
    mo              = month                         | 
    s               = second                        | 
    sidereal        = sidereal day                  | 
    wk              = week                          | 
    yr              = year                          | 

    Usage: conversions.py --from=<unit> --to=<unit> [options] value [error]

    Options:
    -h, --help            show this help message and exit
    --from=_FROM          units to convert from
    --to=_TO              units to convert to

    Extra quantities when changing base units (e.g. Flambda to Fnu):
        -w WAVE, --wave=WAVE
                            wavelength with units (e.g. used to convert Flambda to
                            Fnu)
        -f FREQ, --freq=FREQ
                            frequency (e.g. used to convert Flambda to Fnu)
        -p PHOTBAND, --photband=PHOTBAND
                            photometric passband

To convert from one unit to another, do::

    $:> python conversions.py --from=nRsol/h --to cm/s=1.2345
    1.2345 nRsol/h    =    0.0238501 cm/s

In fact, the C{to} parameter is optional (despite it not being a positional
argument, C{from} is not optional). The script will assume you want to convert
to SI::
    
    $:> python conversions.py --from=nRsol/h 1.2345
    1.2345 nRsol/h    =    0.000238501 m1 s-1

It is also possible to compute with uncertainties, and you can give extra
keyword arguments if extra information is needed for conversions::

    $:> python conversions.py --from=mag --to=erg/s/cm2/A --photband=GENEVA.U 7.84 0.02
    7.84 +/- 0.02 mag    =    4.12191e-12 +/- 7.59283e-14 erg/s/cm2/A
    
If you want to do coordinate transformations, e.g. from fractional radians to 
degrees/arcminutes/arcseconds, you can do::
    
    $:> python conversions.py --from=rad_coord --to=equ 5.412303,0.123
      (5.412303, 0.123) rad_coord    =    20:40:24.51,7:02:50.6 equ

Section 3. Fundamental constants and base unit systems
======================================================

Although fundamental constants are supposed to be constants, there are two
reasons why one might to change their value:

    1. You want to work in a different base unit system (e.g. cgs instead of SI)
    2. You disagree with some of the literature values and want to use your own.
    
Section 3.1. Changing base unit system
--------------------------------------

A solution to the first option might be to define all constants in SI, and define
them again in cgs, adding a postfix C{_cgs} to the constants' name. This is done
in the module L{ivs.units.constants}, but it doesn't offer any flexibility, and
one has to redefine all the values manually for any given convention. Also, you
maybe want to work in cgs by default, which means that typing the postfix every
time is time consuming. For this purpose, you can redefine all constants in a
different base system with the simply command

>>> set_convention(units='cgs')
('SI', 'standard')

and resetting to the default SI can be done by calling L{set_convention} without
any arguments, or simply

>>> set_convention(units='SI')
('cgs', 'standard')

B{Warning:} The value of the constants are changed in B{all} modules, where
an import statement C{from ivs.units import constants} is made. It will B{not}
change the values in modules where the import is done via
C{from ivs.units.constants import *}. If you want to get the value of a
fundamental constant regardless of the preference of base system, you might want
to call

>>> Msol = get_constant('Msol','cgs')

Section 3.2. Changing the values of the fundamental constants
-------------------------------------------------------------

If you disagree with some of the literature values, you can also redefine the
values of the fundamental constants. For example, to use the values defined in
the stellar evolution code B{MESA}, simply do

>>> set_convention(units='cgs',values='mesa')
('SI', 'standard')

But for the sake of the examples, we'll switch back to the default SI...

>>> set_convention()
('cgs', 'mesa')


Section 4.  Examples
====================

B{Example 1:} The following is an exam question on the Biophysics exam (1st year
Bachelor) about error propagation.

Question: Suppose there is a party to celebrate the end of the Biophysics exam.
You want to invite 4 persons, and you want to know if 1 liter of champagne is
enough to fill 5 glasses. The glasses are cylinders with a circumference of
15.5+/-0.5cm, and a height of 10.0+/-0.5cm. Calculate the volume of one glass
and its uncertainty. Can you pour 5 glasses of champagne from the 1 liter
bottle?

Answer:

>>> r = Unit( (15.5/(2*np.pi), 0.5/(2*np.pi)), 'cm')
>>> h = Unit( (10.0,0.5), 'cm')
>>> V = np.pi*r**2*h
>>> print V
0.000191184875389+/-1.56051027314e-05 m3
>>> print (5*V).convert('dm3')
0.955924376946+/-0.0780255136569

It is not sufficient within about 1 sigma.

"""
#-- standard libraries
import itertools
import re
import os
import sys
import logging
import urllib
import numpy as np

#-- optional libraries: WARNING: when these modules are not installed, the
#   module's use is restricted
try: import ephem
except ImportError: print("Unable to load pyephem, stellar coordinate transformations unavailable")

#-- from IVS repository
from ivs.units import constants
#from ivs.units.constants import *
from ivs.units.uncertainties import unumpy,AffineScalarFunc,ufloat
from ivs.units.uncertainties.unumpy import log10,sqrt,sin,cos
from ivs.sed import filters
from ivs.io import ascii
from ivs.aux import loggers
from ivs.aux.decorators import memoized

logger = logging.getLogger("UNITS.CONV")
logger.addHandler(loggers.NullHandler)

#{ Main functions

def convert(_from,_to,*args,**kwargs):
    """
    Convert one unit to another.
    
    Basic explanation
    =================
    
    The unit strings C{_from} and C{_to} should by default be given in the form
    
    C{erg s-1 cm-2 A-1}
    
    Common alternatives are also accepted (see below).
    
    Square brackets '[]' denote a logarithmic value.
    
    If one positional argument is given, it can be either a scalar, numpy array
    or C{uncertainties} object. The function will also return one argument of
    the same type.
    
    If two positional arguments are given, the second argument is assumed to be
    the uncertainty on the first (scalar or numpy array). The function will also
    return two arguments.

    Basic examples:
    
    >>> convert('km','cm',1.)
    100000.0
    >>> convert('m/s','km/h',1,0.1)
    (3.5999999999999996, 0.36)
    
    Keyword arguments can give extra information, for example when converting
    from Flambda to Fnu, and should be tuples (float(,error),'unit'):
    
    >>> convert('A','km/s',4553,0.1,wave=(4552.,0.1,'A'))
    (65.85950307557613, 9.314963362464114)
    
    Extra
    =====
    
    The unit strings C{_from} and C{_to} should by default be given in the form
    
    C{erg s-1 cm-2 A-1}
    
    Common alternatives are also accepted, but don't drive this too far:
    
    C{erg/s/cm2/A}
    
    The crasiest you're allowed to go is something like
    
    >>> print(convert('10mW m-2/nm','erg s-1 cm-2 A-1',1.))
    1.0
    
    But there is a limit on the interpretation of this prefactor also. Floats
    will probably not work, as are exponentials require two digits.
    
    Parentheses are in no circumstances accepted. Some common aliases are also
    resolved (for a full list, see dictionary C{_aliases}):
    
    C{erg/s/cm2/angstrom}
    
    You don't really have to spell both units if one is consistently within one
    base unit. But of course you have to give at least one!:
    
    >>> convert('kg','cgs',1.)
    1000.0
    >>> convert('g','SI',1.)
    0.001
    >>> convert('SI','g',1.)
    1000.0
    
    B{WARNINGS}:
        1. the conversion involving sr and pixels is B{not tested}.
        2. the conversion involving magnitudes is calibrated but not fully tested
        3. non-integer powers are not functioning yet
    
    Examples:
    
    B{Spectra}:
    
    >>> convert('A','km/s',4553.,wave=(4552.,'A'))
    65.85950307557613
    >>> convert('A','km/s',4553.,wave=(4552.,0.1,'A'))
    (65.85950307557613, 6.587397133195861)
    >>> convert('nm','m/s',455.3,wave=(0.4552,'mum'))
    65859.50307564587
    >>> convert('km/s','A',65.859503075576129,wave=(4552.,'A'))
    4553.0
    >>> convert('nm','Ghz',1000.)
    299792.4579999999
    >>> convert('km h-1','nRsol s-1',1.)
    0.3993883287866966
    >>> convert('erg s-1 cm-2 A-1','SI',1.)
    10000000.0
    
    B{Fluxes}:
    
    >>> convert('erg/s/cm2/A','Jy',1e-10,wave=(10000.,'angstrom'))
    333.564095198152
    >>> convert('erg/s/cm2/A','Jy',1e-10,freq=(constants.cc/1e-6,'hz'))
    333.564095198152
    >>> convert('erg/s/cm2/A','Jy',1e-10,freq=(constants.cc,'Mhz'))
    333.564095198152
    >>> convert('Jy','erg/s/cm2/A',333.56409519815202,wave=(10000.,'A'))
    1e-10
    >>> convert('Jy','erg/s/cm2/A',333.56409519815202,freq=(constants.cc,'Mhz'))
    1e-10
    >>> convert('W/m2/mum','erg/s/cm2/A',1e-10,wave=(10000.,'A'))
    1.0000000000000003e-11
    >>> convert('Jy','W/m2/Hz',1.)
    1e-26
    >>> print convert('W/m2/Hz','Jy',1.)
    1e+26
    >>> print convert('Jy','erg/cm2/s/Hz',1.)
    1e-23
    >>> print convert('erg/cm2/s/Hz','Jy',1.)
    1e+23
    >>> convert('Jy','erg/s/cm2',1.,wave=(2.,'micron'))
    1.49896229e-09
    >>> convert('erg/s/cm2','Jy',1.,wave=(2.,'micron'))
    667128190.3963041
    >>> convert('Jy','erg/s/cm2/micron/sr',1.,wave=(2.,'micron'),ang_diam=(3.,'mas'))
    4511059.829810158
    >>> convert('Jy','erg/s/cm2/micron/sr',1.,wave=(2.,'micron'),pix=(3.,'mas'))
    3542978.105308904
    >>> convert('erg/s/cm2/micron/sr','Jy',1.,wave=(2.,'micron'),ang_diam=(3.,'mas'))
    2.2167739682629828e-07
    >>> convert('Jy','erg/s/cm2/micron',1.,wave=(2,'micron'))
    7.49481145e-10
    >>> print(convert('10mW m-2 nm-1','erg s-1 cm-2 A-1',1.))
    1.0
    >>> print convert('Jy','erg s-1 cm-2 micron-1 sr-1',1.,ang_diam=(2.,'mas'),wave=(1.,'micron'))
    40599538.4683
    
    B{Angles}:
    
    >>> convert('sr','deg2',1.)
    3282.806350011744
    
    B{Magnitudes and amplitudes}:
    
    >>> print(convert('ABmag','Jy',0.,photband='SDSS.U'))
    3767.03798984
    >>> print(convert('Jy','erg cm-2 s-1 A-1',3630.7805477,wave=(1.,'micron')))
    1.08848062485e-09
    >>> print(convert('ABmag','erg cm-2 s-1 A-1',0.,wave=(1.,'micron'),photband='SDSS.G'))
    4.93934836475e-09
    >>> print(convert('erg cm-2 s-1 A-1','ABmag',1e-8,wave=(1.,'micron'),photband='SDSS.G'))
    -0.765825856568
    >>> print(convert('ppm','muAmag',1.))
    1.0857356618
    >>> print(convert('mAmag','ppt',1.,0.1))
    (0.9214583192957981, 0.09218827316735488)
    >>> convert('mag_color','flux_ratio',0.599,0.004,photband='GENEVA.U-B')
    (1.1391327795013375, 0.004196720251233045)
    
    B{Frequency analysis}:
    
    >>> convert('cy/d','muHz',1.)
    11.574074074074074
    >>> convert('muhz','cy/d',11.574074074074074)
    1.0
    
    B{Interferometry}:
    
    >>> convert('m','cy/arcsec',85.,wave=(2.2,'micron'))
    187.3143767923207
    >>> convert('cm','cy/arcmin',8500.,wave=(2200,'nm'))/60.
    187.3143767923207
    >>> convert('cy/arcsec','m',187.,wave=(2.2,'mum'))
    84.85734129005546
    >>> convert('cyc/arcsec','m',187.,wave=(1,'mum'))
    38.57151876820703
    >>> convert('cycles/arcsec','m',187.,freq=(300000.,'Ghz'))
    38.54483473437972
    >>> convert('cycles/mas','m',0.187,freq=(300000.,'Ghz'))
    38.54483473437971
    
    B{Temperature}:
    
    >>> print(convert('F','K',123.))
    323.705555556
    >>> print(convert('kF','kK',0.123))
    0.323705555556
    >>> print(convert('K','F',323.7))
    122.99
    >>> print(convert('C','K',10.))
    283.15
    >>> print(convert('C','F',10.))
    50.0
    >>> print(convert('dC','kF',100.))
    0.05
    
    B{Time and Dates}:
    
    >>> convert('sidereal d','d',1.)
    1.0027379093
    >>> convert('JD','CD',2446257.81458)
    (1985.0, 7.0, 11.314580000005662)
    >>> convert('CD','JD',(1985,7,11.31))
    2446257.81
    >>> convert('CD','JD',(1985,7,11,7,31,59))
    2446257.813877315
    >>> convert('MJD','CD',0.,jtype='corot')
    (2000.0, 1.0, 1.5)
    >>> convert('JD','MJD',2400000.5,jtype='mjd')
    0.0
    >>> convert('MJD','CD',0,jtype='mjd')
    (1858.0, 11.0, 17.0)
    
    B{Coordinates}: When converting coordinates with pyephem, you get pyephem
    coordinates back. They have built-in string represenations and float
    conversions:
    
    >>> x,y = convert('equ','gal',('17:45:40.4','-29:00:28.1'),epoch='2000')
    >>> print (x,y)
    (6.282224277178722, -0.000825178833899176)
    >>> print x,y
    359:56:41.8 -0:02:50.2
    >>> x,y = convert('gal','equ',('00:00:00.00','00:00:00.0'),epoch='2000')
    >>> print (x,y)
    (4.64964430303663, -0.5050315085342665)
    >>> print x,y
    17:45:37.20 -28:56:10.2
    
    It is also possible to immediately convert to radians or degrees in floats
    (this can be useful for plotting):
    >>> convert('equ','deg_coord',('17:45:40.4','-29:00:28.1'))
    (266.41833333333335, -29.00780555555556)
    >>> convert('equ','rad_coord',('17:45:40.4','-29:00:28.1'))
    (4.649877104342426, -0.5062817157227474)
    
    
    @param _from: units to convert from
    @type _from: str
    @param _to: units to convert to
    @type _to: str
    @keyword unpack: set to True if you don't want 'uncertainty objects'. If True
    and uncertainties are given, they will be returned as a tuple (value, error)
    instead of uncertainty object. Set to False probably only for internal uses
    @type unpack: boolean, defaults to True
    @return: converted value
    @rtype: float
    """
    #-- remember if user wants to unpack the results to have no trace of
    #   uncertainties, or wants to get uncertainty objects back
    unpack = kwargs.pop('unpack',True)
    
    #-- get the input arguments: if only one is given, it is either an
    #   C{uncertainty} from the C{uncertainties} package, or it is just a float
    if len(args)==1:
        start_value = args[0]
    #   if two arguments are given, we assume the first is the actual value and
    #   the second is the error on the value
    elif len(args)==2:
        start_value = unumpy.uarray([args[0],args[1]])
    else:
        raise ValueError,'illegal input'
    
    #-- (un)logarithmicize (denoted by '[]')
    m_in = re.search(r'\[(.*)\]',_from)
    m_out = re.search(r'\[(.*)\]',_to)
    
    if m_in is not None:
        _from = m_in.group(1)
        start_value = 10**start_value
    if m_out is not None:
        _to = m_out.group(1)
    
    #-- It is possible the user gave a convention for either the from or to
    #   units (but not both!)
    #-- break down the from and to units to their basic elements
    if _from in _conventions:
        _from = change_convention(_from,_to)
    elif _to in _conventions:
        _to = change_convention(_to,_from)
    fac_from,uni_from = breakdown(_from)
    fac_to,uni_to = breakdown(_to)
    
        
    
    #-- convert the kwargs to SI units if they are tuples (make a distinction
    #   when uncertainties are given)
    kwargs_SI = {}
    for key in kwargs:
        if isinstance(kwargs[key],tuple) and len(kwargs[key])==2:
            fac_key,uni_key = breakdown(kwargs[key][1])
            kwargs_SI[key] = fac_key*kwargs[key][0]
        elif isinstance(kwargs[key],tuple) and len(kwargs[key])==3:
            fac_key,uni_key = breakdown(kwargs[key][2])
            kwarg_value = unumpy.uarray([kwargs[key][0],kwargs[key][1]])
            kwargs_SI[key] = fac_key*kwarg_value
        else:
            kwargs_SI[key] = kwargs[key]
    
    #-- add some default values if necessary
    if uni_from!=uni_to and uni_from=='m1' and not ('wave' in kwargs_SI):# or 'freq' in kwargs_SI or 'photband' in kwargs_SI):
        kwargs_SI['wave'] = fac_from*start_value
        logger.warning('Assumed input value to serve also for "wave" key')
    elif uni_from!=uni_to and uni_from=='cy1 s-1' and not ('wave' in kwargs_SI):# or 'freq' in kwargs_SI or 'photband' in kwargs_SI):
        kwargs_SI['freq'] = fac_from*start_value
        logger.warning('Assumed input value to serve also for "freq" key')
        
    logger.debug('Convert %s to %s'%(uni_from,uni_to))
    
    #-- conversion is easy if same units
    ret_value = 1.
    if uni_from==uni_to:
        #-- if nonlinear conversions from or to:
        if isinstance(fac_from,NonLinearConverter):
            ret_value *= fac_from(start_value,**kwargs_SI)
        else:
            ret_value *= fac_from*start_value
            
    #-- otherwise a little bit more complicated
    else:
        #-- first check where the unit differences are
        uni_from_ = uni_from.split()
        uni_to_ = uni_to.split()
        only_from_c,only_to_c = sorted(list(set(uni_from_) - set(uni_to_))),sorted(list(set(uni_to_) - set(uni_from_)))
        only_from_c,only_to_c = [list(components(i))[1:] for i in only_from_c],[list(components(i))[1:] for i in only_to_c]
        #-- push them all bach to the left side (change sign of right hand side components)
        left_over = " ".join(['%s%d'%(i,j) for i,j in only_from_c])
        left_over+= " "+" ".join(['%s%d'%(i,-j) for i,j in only_to_c])
        left_over = breakdown(left_over)[1]
        only_from = "".join(left_over.split())
        only_to = ''

        #-- first we remove any differences concerning (ster)radians
        #   we recently added fac_from* to all these things, maybe this needs to 
        #   change?
        if 'rad2' in only_from:
            start_value = fac_from*_switch['rad2_to_'](start_value,**kwargs_SI)
            only_from = only_from.replace('rad2','')
            logger.debug('Switching to /sr')
            fac_from = 1.
        elif 'rad-2' in only_from:
            start_value = fac_from*_switch['rad-2_to_'](start_value,**kwargs_SI)
            only_from = only_from.replace('rad-2','')
            logger.debug('Switching from /sr')
            fac_from = 1.
        elif 'rad1' in only_from:
            start_value = fac_from*_switch['rad1_to_'](start_value,**kwargs_SI)
            only_from = only_from.replace('rad1','')
            logger.debug('Switching to /rad')
            fac_from = 1.
        elif 'rad-1' in only_from:
            start_value = fac_from*_switch['rad-1_to_'](start_value,**kwargs_SI)
            only_from = only_from.replace('rad-1','')
            logger.debug('Switching from /rad')
            fac_from = 1.
        
        #-- then we do what is left over (if anything is left over)
        if only_from or only_to:
            logger.debug("Convert %s to %s"%(only_from,only_to))
            
            #-- nonlinear conversions need a little tweak
            try:
                key = '%s_to_%s'%(only_from,only_to)
                logger.debug('Switching from %s to %s'%(only_from,only_to))
                if isinstance(fac_from,NonLinearConverter):
                    ret_value *= _switch[key](fac_from(start_value,**kwargs_SI),**kwargs_SI)
                #-- linear conversions are easy
                else:
                    logger.debug('fac_from=%s, start_value=%s'%(fac_from,start_value))
                    ret_value *= _switch[key](fac_from*start_value,**kwargs_SI)
            except KeyError:
                logger.critical('cannot convert %s to %s: no %s definition in dict _switch'%(_from,_to,key))
                raise
        else:
            ret_value *= start_value
    #-- final step: convert to ... (again distinction between linear and
    #   nonlinear converters)
    if isinstance(fac_to,NonLinearConverter):
        ret_value = fac_to(ret_value,inv=True,**kwargs_SI)
    else:
        ret_value /= fac_to
    
    #-- logarithmicize
    if m_out is not None:
        ret_value = log10(ret_value)
    
    #-- unpack the uncertainties if: 
    #    1. the input was not given as an uncertainty
    #    2. the input was without uncertainties, but extra keywords had uncertainties
    if unpack and (len(args)==2 or (len(args)==1 and isinstance(ret_value,AffineScalarFunc))):
        ret_value = unumpy.nominal_values(ret_value),unumpy.std_devs(ret_value)
        #-- convert to real floats if real floats were given
        if not ret_value[0].shape:
            ret_value = np.asscalar(ret_value[0]),np.asscalar(ret_value[1])
    
    
    return ret_value


def nconvert(_froms,_tos,*args,**kwargs):
    """
    Convert a list/array/tuple of values with different units to other units.
    
    This silently catches some exceptions and replaces the value with nan!
    """
    if len(args)==1:
        ret_value = np.zeros((len(args[0])))
    elif len(args)==2:
        ret_value = np.zeros((len(args[0]),2))
    if isinstance(_tos,str):
        _tos = [_tos for i in _froms]
    elif isinstance(_froms,str):
        _froms = [_froms for i in _tos]
    
    for i,(_from,_to) in enumerate(zip(_froms,_tos)):
        myargs = [iarg[i] for iarg in args]
        mykwargs = {}
        for key in kwargs:
            if not isinstance(kwargs[key],str) and hasattr(kwargs[key],'__iter__'):
                mykwargs[key] = kwargs[key][i]
            else:
                mykwargs[key] = kwargs[key]
        try:
            ret_value[i] = convert(_from,_to,*myargs,**mykwargs)
        except ValueError: #no calibration
            ret_value[i] = np.nan
    
    if len(args)==2:
        ret_value = ret_value.T
    return ret_value


def change_convention(to_,units):
    """
    Change units from one convention to another.
        
    Example usage:
    
    >>> units1 = 'kg m2 s-2 K-1 mol-1'
    >>> print change_convention('cgs',units1)
    K-1 g1 cm2 mol-1 s-2
    >>> print change_convention('sol','g cm2 s-2 K-1 mol-1')
    Tsol-1 Msol1 Rsol2 mol-1 s-2
    
    @param to_: convention name to change to
    @type to_: str
    @param units: units to change
    @type units: str
    @return: units in new convention
    @rtype: str
    """
    #-- break down units in base units, and breakdown that string in
    #   whole units and non-alpha digits (we'll weave them back in after
    #   translation)
    factor,units = breakdown(units)
    new_units = re.findall(r'[a-z]+',units, re.I)
    powers = re.findall(r'[0-9\W]+',units, re.I)
    #-- make the translation dictionary
    translator = {}
    for key in sorted(_conventions['SI'].keys()):
        translator[_conventions['SI'][key]] = _conventions[to_][key]
    #-- translate
    new_units = [unit in translator and translator[unit] or unit for unit in new_units]
    #-- weave them back in
    new_units = "".join(["".join([i,j]) for i,j in zip(new_units,powers)])
    return new_units
    
def set_convention(units='SI',values='standard'):
    """
    Consistently change the values of the fundamental constants to the paradigm
    of other system units or programs.
    
    This can be important in pulsation codes, where e.g. the value of the
    gravitational constant and the value of the solar radius can have siginficant
    influences on the frequencies.
    
    @param units: name of the new units base convention
    @type units: str
    @param values: name of the value set for fundamental constants
    @type values: str
    @return: name of old convention
    @rtype: str
    """
    to_return = constants._current_convention,constants._current_values
    values = values.lower()
    #-- first reload the constants to their original values (i.e. SI based)
    reload(constants)
    #-- then, where possible, replace all the constants with the value from
    #   the other convention:
    cvars = dir(constants)
    cvars = [i for i in cvars if i+'_units' in cvars]
    for const in cvars:
        #-- break down the old units in their basic components, but replace
        #   the value of the constant from those of the C{values} set if
        #   possible
        old_units = getattr(constants,const+'_units')
        const_ = const+'_'+values
        if hasattr(constants,const_):
            logger.info('Override {0} value to {1}'.format(const,values))
            old_value = getattr(constants,const_)
        else:
            old_value = getattr(constants,const)
        factor,old_units = breakdown(old_units)
        #-- replace old base units with new conventions's base units
        new_units = change_convention(units,old_units)
        #-- convert the value from the old to the new convenction
        new_value = convert(old_units,new_units,old_value)
        #-- and attach to the constants module
        setattr(constants,const,new_value)
        setattr(constants,const+'_units',new_units)
    constants._current_convention = units
    constants._current_values = values
    logger.info('Changed convention to {0} with values from {1} set'.format(units,values))
    return to_return

def get_constant(constant_name,units='SI',value='standard'):
    """
    Convenience function to retrieve the value of a constant in a particular
    system.
    
    >>> Msol = get_constant('Msol','SI')
    >>> Msol = get_constant('Msol','kg')
    
    @param constant_name: name of the constant
    @type constant: str
    @param units_system: name of the unit base system
    @type units_system: str
    @param value: name of the parameter set the get the value from
    @type value: str
    @return: value of the constant in the unit base system
    @rtype: float
    """
    value = value.lower()
    #-- see if we need (and have) a different value than the standard one
    const_ = constant_name+'_'+value
    if hasattr(constants,const_):
        old_value = getattr(constants,const_)
    else:
        old_value = getattr(constants,constant_name)
    old_units = getattr(constants,constant_name+'_units')
    new_value = convert(old_units,units,old_value)
    return new_value

def get_constants(units='SI',values='standard'):
    """
    Convenience function to retrieve the value of all constants in a particular
    system.
    
    This can be helpful when you want to attach a copy of the fundamental constants
    to some module or class instances, and be B{sure} these values never change.
    
    Yes, I know, there's a lot that can go wrong with values that are supposed
    to be CONSTANT!
    """
    cvars = dir(constants)
    cvars = [i for i in cvars if i+'_units' in cvars]
    myconstants = {}
    for cvar in cvars:
        myconstants[cvar] = get_constant(cvar,units,values)
    return myconstants
    
    
    
    
def round_arbitrary(x, base=5):
    """
    Round to an arbitrary base.
    
    Example usage:
    
    >>> round_arbitrary(1.24,0.25)
    1.25
    >>> round_arbitrary(1.37,0.75)
    1.5
    
    
    @param x: number to round
    @type x: float
    @param base: base to round to
    @type base: integer or float
    @return: rounded number
    @rtype: float
    """
    return base * round(float(x)/base)
#}
#{ Conversions basics
def solve_aliases(unit):
    """
    Resolve simple aliases in a unit's name.
    
    Resolves aliases and replaces division signs with negative power.
    
    @param unit: unit (e.g. erg s-1 angstrom-1)
    @type unit: str
    @return: aliases-resolved unit (e.g. erg s-1 A-1)
    @rtype: str
    """
    #-- resolve aliases
    for alias in _aliases:
        unit = unit.replace(alias[0],alias[1])
    
    #-- replace slash-forward with negative powers
    if '/' in unit:
        unit_ = [uni.split('/') for uni in unit.split()]
        for i,uni in enumerate(unit_):
            for j,after_div in enumerate(uni[1:]):
                if not after_div[-1].isdigit(): after_div += '1'
                #m = re.search(r'(\d*)(.+?)(-{0,1}\d+)',after_div)
                m = re.search(r'(\d*)(.+?)(-{0,1}[\d\.]+)',after_div)
                if m is not None:
                    factor,basis,power = m.group(1),m.group(2),m.group(3)
                    if not '.' in power: power = int(power)
                    else:                power = float(power)
                    if factor: factor = float(factor)
                    else: factor = 1.
                else:
                    factor,basis,power = 1.,after_div,1
                uni[1+j] = '%s%g'%(basis,-power)
                if factor!=1: uni[1+j] = '%d%s'%(factor,uni[1+j])
        ravelled = []
        for uni in unit_:
            ravelled += uni
        unit = " ".join(ravelled)
    return unit
    

def components(unit):
    """
    Decompose a unit into: a factor, SI base unit, power.
    
    Examples:
    
    >>> print(components('m'))
    (1.0, 'm', 1)
    >>> print(components('g2'))
    (0.001, 'kg', 2)
    >>> print(components('hg3'))
    (0.1, 'kg', 3)
    >>> print(components('Mg4'))
    (1000.0, 'kg', 4)
    >>> print(components('mm'))
    (0.001, 'm', 1)
    >>> print(components('W3'))
    (1.0, 'kg m2 s-3', 3)
    >>> print(components('s-2'))
    (1.0, 's', -2)
    
    @param unit: unit name
    @type unit: str
    @return: 3-tuple with factor, SI base unit and power
    @rtype: (float,str,int)
    """
    factor = 1.
    #-- manually check if there is a prefactor of the form '10-14' or '10e-14'
    #   and expand it if necessary
    m = re.search('\d\d[-+]\d\d',unit[:5])
    if m is not None:
        factor *= float(m.group(0)[:2])**(float(m.group(0)[2]+'1')*float(m.group(0)[3:5]))
        unit = unit[5:]
    m = re.search('\d\d[[eE][-+]\d\d',unit[:6])
    if m is not None:
        factor *= float(m.group(0)[:2])**(float(m.group(0)[3]+'1')*float(m.group(0)[4:6]))
        unit = unit[6:]
    
    if not unit[-1].isdigit(): unit += '1'
    #-- decompose unit in base name and power
    #m = re.search(r'(\d*)(.+?)(-{0,1}\d+)',unit)
    m = re.search(r'(\d*)(.+?)(-{0,1}[\d\.]+)',unit)
    if m is not None:
        factor_,basis,power = m.group(1),m.group(2),m.group(3)
        #-- try to make power an integer, otherwise make it a float
        if not '.' in power: power = int(power)
        else:                power = float(power)
        if factor_:
            factor *= float(factor_)
    else:
        factor,basis,power = 1.,unit,1
    #-- decompose the base name (which can be a composition of a prefix
    #   (e.g., 'mu') and a unit name (e.g. 'm')) into prefix and unit name
    #-- check if basis is part of _factors dictionary. If not, find the
    #   combination of _scalings and basis which is inside the dictionary!
    for scale in _scalings:
        scale_unit,base_unit = basis[:len(scale)],basis[len(scale):]
        if scale_unit==scale and base_unit in _factors and not basis in _factors:
            #if basis in _factors:
            #    raise ValueError,'ambiguity between %s and %s-%s'%(basis,scale_unit,base_unit)
            factor *= _scalings[scale]
            basis = base_unit
            break
    #-- if we didn't find any scalings, check if the 'raw' unit is already
    #   a base unit
    else:
        if not basis in _factors:
            raise ValueError, 'Unknown unit %s'%(basis)
        
    #-- switch from base units to SI units
    if hasattr(_factors[basis][0],'__call__'):
        factor = factor*_factors[basis][0]()
    else:
        factor *= _factors[basis][0]
    basis = _factors[basis][1]
    
    return factor,basis,power

def breakdown(unit):
    """
    Decompose a unit into SI base units containing powers.
    
    Examples:
    
    >>> print(breakdown('erg s-1 W2 kg2 cm-2'))
    (0.001, 'kg5 m4 s-9')
    >>> print(breakdown('erg s-1 cm-2 A-1'))
    (10000000.0, 'kg1 m-1 s-3')
    >>> print(breakdown('W m-3'))
    (1.0, 'kg1 m-1 s-3')
    
    @param unit: unit's name
    @type unit: str
    @return: 2-tuple factor, unit's base name
    @rtype: (float,str)
    """
    #-- solve aliases
    unit = solve_aliases(unit)
    #-- break down in basic units
    units = unit.split()
    total_factor = 1.
    total_units = []
    total_power = []
    for unit in units:
        factor,basis,power = components(unit)
        total_factor = total_factor*factor**power
        basis = basis.split()
        for base in basis:
            factor_,basis_,power_ = components(base)
            if basis_ in total_units:
                index = total_units.index(basis_)
                total_power[index] += power_*power
            else:
                total_units.append(basis_)
                total_power.append(power_*power)
    
    #-- make sure to return a sorted version
    total_units = sorted(['%s%s'%(i,j) for i,j in zip(total_units,total_power) if j!=0])
    return total_factor," ".join(total_units)

def get_help():
    """
    Return a string with a list and explanation of all defined units
    """
    #set_exchange_rates()
    help_text = {}
    for fac in sorted(_factors.keys()):
        if _factors[fac][2] not in help_text:
            help_text[_factors[fac][2]] = []
        help_text[_factors[fac][2]].append('%-15s = %-30s'%(fac,_factors[fac][3]))
    text = [[],[]]
    bar = '%-48s'%('================='+20*'=')
    for i,key in enumerate(sorted(help_text.keys())):
        text[i%2] += [bar,'%-48s'%("=   Units of %-20s   ="%(key)),bar]
        text[i%2] += help_text[key]
    out = ''
    #for i,j in itertools.izip_longest(*text,fillvalue=''):
    #    out += '%s| %s\n'%(i,j)
    return out
        
#}
#{ Linear change-of-base conversions
        
def distance2velocity(arg,**kwargs):
    """
    Switch from distance to velocity via a reference wavelength.
    
    @param arg: distance (SI, m)
    @type arg: float
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @return: velocity (SI, m/s)
    @rtype: float
    """
    if 'wave' in kwargs:
        wave = kwargs['wave']
        velocity = (arg-wave) / wave * constants.cc
    else:
        raise ValueError,'reference wavelength (wave) not given'
    return velocity

def velocity2distance(arg,**kwargs):
    """
    Switch from velocity to distance via a reference wavelength.
    
    @param arg: velocity (SI, m/s)
    @type arg: float
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @return: distance (SI, m)
    @rtype: float
    """
    if 'wave' in kwargs:
        wave = kwargs['wave']
        distance = wave / constants.cc * arg + wave
    else:
        raise ValueError,'reference wavelength (wave) not given'
    return distance

def fnu2flambda(arg,**kwargs):
    """
    Switch from Fnu to Flambda via a reference wavelength.
    
    Flambda and Fnu are spectral irradiance in wavelength and frequency,
    respectively
    
    @param arg: spectral irradiance (SI,W/m2/Hz)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        wave = kwargs['wave']
        flambda = constants.cc/wave**2 * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        flambda = freq**2/constants.cc * arg
    else:
        raise ValueError,'reference wave/freq not given'
    return flambda

def flambda2fnu(arg,**kwargs):
    """
    Switch from Flambda to Fnu via a reference wavelength.
    
    Flambda and Fnu are spectral irradiance in wavelength and frequency,
    respectively
    
    @param arg: spectral irradiance (SI, W/m2/m)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI,W/m2/Hz)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        wave = kwargs['wave']
        fnu = wave**2/constants.cc * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        fnu = constants.cc/freq**2 * arg
    else:
        raise ValueError,'reference wave/freq not given'
    return fnu

def fnu2nufnu(arg,**kwargs):
    """
    Switch from Fnu to nuFnu via a reference wavelength.
    
    Flambda and Fnu are spectral irradiance in wavelength and frequency,
    respectively
    
    @param arg: spectral irradiance (SI,W/m2/Hz)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        wave = kwargs['wave']
        fnu = constants.cc/wave * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        fnu = freq * arg
    else:
        raise ValueError,'reference wave/freq not given'
    return fnu

def nufnu2fnu(arg,**kwargs):
    """
    Switch from nuFnu to Fnu via a reference wavelength.
    
    Flambda and Fnu are spectral irradiance in wavelength and frequency,
    respectively
    
    @param arg: spectral irradiance (SI,W/m2/Hz)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        wave = kwargs['wave']
        fnu = wave/constants.cc * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        fnu = arg / freq
    else:
        raise ValueError,'reference wave/freq not given'
    return fnu

def flam2lamflam(arg,**kwargs):
    """
    Switch from lamFlam to Flam via a reference wavelength.
    
    Flambda and Fnu are spectral irradiance in wavelength and frequency,
    respectively
    
    @param arg: spectral irradiance (SI,W/m2/Hz)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        wave = kwargs['wave']
        lamflam = wave * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        lamflam = cc/freq * arg
    else:
        raise ValueError,'reference wave/freq not given'
    return lamflam

def lamflam2flam(arg,**kwargs):
    """
    Switch from lamFlam to Flam via a reference wavelength.
    
    Flambda and Fnu are spectral irradiance in wavelength and frequency,
    respectively
    
    @param arg: spectral irradiance (SI,W/m2/Hz)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        wave = kwargs['wave']
        flam = arg / wave
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        flam = arg / (cc/freq)
    else:
        raise ValueError,'reference wave/freq not given'
    return flam

def distance2spatialfreq(arg,**kwargs):
    """
    Switch from distance to spatial frequency via a reference wavelength.
    
    @param arg: distance (SI, m)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spatial frequency (SI, cy/as)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        spatfreq = 2*np.pi*arg/kwargs['wave']
    elif 'freq' in kwargs:
        spatfreq = 2*np.pi*arg/(constants.cc/kwargs['freq'])
    else:
        raise ValueError,'reference wave/freq not given'
    return spatfreq

def spatialfreq2distance(arg,**kwargs):
    """
    Switch from spatial frequency to distance via a reference wavelength.
    
    @param arg: spatial frequency (SI, cy/as)
    @type arg: float
    @keyword photband: photometric passband
    @type photband: str ('SYSTEM.FILTER')
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: distance (SI, m)
    @rtype: float
    """
    if 'photband' in kwargs:
        lameff = filters.eff_wave(kwargs['photband'])
        lameff = convert('A','m',lameff)
        kwargs['wave'] = lameff
    if 'wave' in kwargs:
        distance = kwargs['wave']*arg/(2*np.pi)
    elif 'freq' in kwargs:
        distance = constants.cc/kwargs['freq']*arg/(2*np.pi)
    else:
        raise ValueError,'reference wave/freq not given'
    return distance

def per_sr(arg,**kwargs):
    """
    Switch from [Q] to [Q]/sr
    
    @param arg: some SI unit
    @type arg: float
    @return: some SI unit per steradian
    @rtype: float
    """
    if 'ang_diam' in kwargs:
        radius = kwargs['ang_diam']/2.
        surface = np.pi*radius**2
    elif 'radius' in kwargs:
        radius = kwargs['radius']
        surface = np.pi*radius**2
    elif 'pix' in kwargs:
        pix = kwargs['pix']
        surface = pix**2
    else:
        raise ValueError,'angular size (ang_diam/radius) not given'
    Qsr = arg/surface
    return Qsr

def times_sr(arg,**kwargs):
    """
    Switch from [Q]/sr to [Q]
    
    @param arg: some SI unit per steradian
    @type arg: float
    @return: some SI unit
    @rtype: float
    """
    if 'ang_diam' in kwargs:
        radius = kwargs['ang_diam']/2.
        surface = np.pi*radius**2
    elif 'radius' in kwargs:
        radius = kwargs['radius']
        surface = np.pi*radius**2
    elif 'pix' in kwargs:
        pix = kwargs['pix']
        surface = pix**2
    else:
        raise ValueError,'angular size (ang_diam/radius) not given'
    Q = arg*surface
    return Q

def per_cy(arg,**kwargs):
    return arg / (2*np.pi)

def times_cy(arg,**kwargs):
    return 2*np.pi*arg

def period2freq(arg,**kwargs):
    return 1./arg

def do_nothing(arg,**kwargs):
    logger.debug('Experimental: probably just dropped the "cy" unit, please check results')
    return arg


#}


#{ Stellar calibrations

def derive_luminosity(radius,temperature, unit='Lsol'):
    """
    Convert radius and effective temperature to stellar luminosity.
    
    Units given to radius and temperature must be understandable by C{convert}.
    
    Stellar luminosity is returned in solar units.
    
    
    @param luminosity: (Luminosity(, error), units)
    @type luminosity: 2 or 3 tuple
    @param temperature: (effective temperature(, error), units)
    @type temperature: 2 or 3 tuple
    @return: radius (and error) in SI units
    @rtype: 1- or 2-tuple
    """
    #-- take care of luminosity
    if len(radius)==3:
        radius = unumpy.uarray([radius[0],radius[1]]),radius[2]
    radius = convert(radius[-1],'SI',*radius[:-1],unpack=False)
    #-- take care of effective temperature
    if len(temperature)==3:
        temperature = unumpy.uarray([temperature[0],temperature[1]]),temperature[2]
    teff = convert(temperature[-1],'SI',*temperature[:-1],unpack=False)
    #-- calculate radius in SI units
    lumi = 4*np.pi*constants.sigma*radius**2*teff**4
    lumi = convert('W',unit,lumi)
    return lumi

def derive_radius(luminosity,temperature, unit='m'):
    """
    Convert luminosity and effective temperature to stellar radius.
    
    Units given to luminosity and temperature must be understandable by C{convert}.
    
    Stellar radius is returned in SI units.
    
    Example usage:
    
    #>>> calculate_radius((3.9,'[Lsol]'),(3.72,'[K]'))/Rsol
    107.994124114
    
    @param luminosity: (Luminosity(, error), units)
    @type luminosity: 2 or 3 tuple
    @param temperature: (effective temperature(, error), units)
    @type temperature: 2 or 3 tuple
    @return: radius (and error) in SI units
    @rtype: 1- or 2-tuple
    """
    #-- take care of luminosity
    if len(luminosity)==3:
        luminosity = unumpy.uarray([luminosity[0],luminosity[1]])
    lumi = convert(luminosity[-1],'SI',*luminosity[:-1],unpack=False)
    #-- take care of effective temperature
    if len(temperature)==3:
        temperature = unumpy.uarray([temperature[0],temperature[1]])
    teff = convert(temperature[-1],'SI',*temperature[:-1],unpack=False)
    #-- calculate radius in SI units
    R = sqrt(lumi / (teff**4)/(4*np.pi*constants.sigma))
    R = convert('m',unit,R)
    return R    

def derive_radius_slo(numax,Deltanu0,teff,unit='Rsol'):
    """
    Derive stellar radius from solar-like oscillations diagnostics.
    
    @param numax: (numax(, error), units)
    @type numax: 2 or 3 tuple
    @param Deltanu0: (large separation(, error), units)
    @type Deltanu0: 2 or 3 tuple
    @param teff: (effective temperature(, error), units)
    @type teff: 2 or 3 tuple
    @return: radius (and error) in whatever units
    @rtype: 1- or 2-tuple
    """
    numax_sol = convert(constants.numax_sol[-1],'mHz',*constants.numax_sol[:-1],unpack=False)
    Deltanu0_sol = convert(constants.Deltanu0_sol[-1],'muHz',*constants.Deltanu0_sol[:-1],unpack=False)
    #-- take care of temperature
    if len(teff)==3:
        teff = (unumpy.uarray([teff[0],teff[1]]),teff[2])
    teff = convert(teff[-1],'5777K',*teff[:-1],unpack=False)
    #-- take care of numax
    if len(numax)==3:
        numax = (unumpy.uarray([numax[0],numax[1]]),numax[2])
    numax = convert(numax[-1],'mHz',*numax[:-1],unpack=False)
    #-- take care of large separation
    if len(Deltanu0)==3:
        Deltanu0 = (unumpy.uarray([Deltanu0[0],Deltanu0[1]]),Deltanu0[2])
    Deltanu0 = convert(Deltanu0[-1],'muHz',*Deltanu0[:-1],unpack=False)
    R = sqrt(teff)/numax_sol * numax/Deltanu0**2 * Deltanu0_sol**2
    return convert('Rsol',unit,R)    
    

    
    
def derive_logg(mass,radius, unit='[cm/s2]'):
    """
    Convert mass and radius to stellar surface gravity.
    
    Units given to mass and radius must be understandable by C{convert}.
    
    Logarithm of surface gravity is returned in CGS units.
    
    @param mass: (mass(, error), units)
    @type mass: 2 or 3 tuple
    @param radius: (radius(, error), units)
    @type radius: 2 or 3 tuple
    @return: log g (and error) in CGS units
    @rtype: 1- or 2-tuple
    """
    #-- take care of mass
    if len(mass)==3:
        mass = (unumpy.uarray([mass[0],mass[1]]),mass[2])
    M = convert(mass[-1],'g',*mass[:-1],unpack=False)
    #-- take care of radius
    if len(radius)==3:
        radius = (unumpy.uarray([radius[0],radius[1]]),radius[2])
    R = convert(radius[-1],'cm',*radius[:-1],unpack=False)
    #-- calculate surface gravity in logarithmic CGS units
    logg = log10(constants.GG_cgs*M / (R**2))
    logg = convert('[cm/s2]',unit,logg)
    return logg

def derive_logg_slo(teff,numax, unit='[cm/s2]'):
    """
    Derive stellar surface gravity from solar-like oscillations diagnostics.
    
    Units given to teff and numax must be understandable by C{convert}.
    
    Logarithm of surface gravity is returned in CGS units.
    
    @param teff: (effective temperature(, error), units)
    @type teff: 2 or 3 tuple
    @param numax: (numax(, error), units)
    @type numax: 2 or 3 tuple
    @return: log g (and error) in CGS units
    @rtype: 1- or 2-tuple
    """
    numax_sol = convert(constants.numax_sol[-1],'mHz',*constants.numax_sol[:-1],unpack=False)
    #-- take care of temperature
    if len(teff)==3:
        teff = (unumpy.uarray([teff[0],teff[1]]),teff[2])
    teff = convert(teff[-1],'5777K',*teff[:-1],unpack=False)
    #-- take care of numax
    if len(numax)==3:
        numax = (unumpy.uarray([numax[0],numax[1]]),numax[2])
    numax = convert(numax[-1],'mHz',*numax[:-1],unpack=False)
    #-- calculate surface gravity in logarithmic CGS units
    GG = convert(constants.GG_units,'Rsol3 Msol-1 s-2',constants.GG)
    surf_grav = GG*sqrt(teff)*numax / numax_sol
    logg = convert('Rsol s-2',unit,surf_grav)
    return logg    



    

def derive_mass(surface_gravity,radius,unit='kg'):
    """
    Convert surface gravity and radius to stellar mass.
    """
    #-- take care of logg
    if len(surface_gravity)==3:
        surface_gravity = (unumpy.uarray([surface_gravity[0],surface_gravity[1]]),surface_gravity[2])
    grav = convert(surface_gravity[-1],'m/s2',*surface_gravity[:-1],unpack=False)
    #-- take care of radius
    if len(radius)==3:
        radius = (unumpy.uarray([radius[0],radius[1]]),radius[2])
    R = convert(radius[-1],'m',*radius[:-1],unpack=False)
    #-- calculate mass in SI
    M = grav*R**2/constants.GG
    return convert('kg',unit,M)

def derive_numax(mass,radius,temperature,unit='mHz'):
    """
    Derive the predicted nu_max according to Kjeldsen and Bedding (1995).
    
    Example: compute the predicted numax for the Sun in mHz
    >>> print derive_numax((1.,'Msol'),(1.,'Rsol'),(5777.,'K'))
    3.05
    """
    #-- take care of mass
    if len(mass)==3:
        mass = (unumpy.uarray([mass[0],mass[1]]),mass[2])
    M = convert(mass[-1],'Msol',*mass[:-1],unpack=False)
    #-- take care of radius
    if len(radius)==3:
        radius = (unumpy.uarray([radius[0],radius[1]]),radius[2])
    R = convert(radius[-1],'Rsol',*radius[:-1],unpack=False)
    #-- take care of effective temperature
    if len(temperature)==3:
        temperature = (unumpy.uarray([temperature[0],temperature[1]]),temperature[2])
    teff = convert(temperature[-1],'5777K',*temperature[:-1],unpack=False)
    #-- predict nu_max
    nu_max = M/R**2/sqrt(teff)*3.05
    return convert('mHz',unit,nu_max)

def derive_nmax(mass,radius,temperature):
    """
    Derive the predicted n_max according to Kjeldsen and Bedding (1995).
    
    Example: compute the predicted numax for the Sun in mHz
    >>> print derive_nmax((1.,'Msol'),(1.,'Rsol'),(5777.,'K'))
    21.0
    """
    #-- take care of mass
    if len(mass)==3:
        mass = (unumpy.uarray([mass[0],mass[1]]),mass[2])
    M = convert(mass[-1],'Msol',*mass[:-1])
    #-- take care of radius
    if len(radius)==3:
        radius = (unumpy.uarray([radius[0],radius[1]]),radius[2])
    R = convert(radius[-1],'Rsol',*radius[:-1])
    #-- take care of effective temperature
    if len(temperature)==3:
        temperature = unumpy.uarray([temperature[0],temperature[1]])
    teff = convert(temperature[-1],'5777K',*temperature[:-1])
    #-- predict n_max
    n_max = sqrt(M/teff/R)*22.6 - 1.6
    return n_max
    
def derive_Deltanu0(mass,radius,unit='mHz'):
    """
    Derive the predicted large spacing according to Kjeldsen and Bedding (1995).
    
    Example: compute the predicted large spacing for the Sun in mHz
    >>> print derive_Deltanu0((1.,'Msol'),(1.,'Rsol'),unit='muHz')
    134.9
    """
    #-- take care of mass
    if len(mass)==3:
        mass = (unumpy.uarray([mass[0],mass[1]]),mass[2])
    M = convert(mass[-1],'Msol',*mass[:-1])
    #-- take care of radius
    if len(radius)==3:
        radius = (unumpy.uarray([radius[0],radius[1]]),radius[2])
    R = convert(radius[-1],'Rsol',*radius[:-1])
    #-- predict large spacing
    Deltanu0 = sqrt(M)*R**(-1.5)*134.9
    return convert('muHz',unit,Deltanu0)

def derive_ampllum(luminosity,mass,temperature,wavelength,unit='ppm'):
    """
    Derive the luminosity amplitude around nu_max of solar-like oscillations.
    
    See Kjeldsen and Bedding (1995).
    
    >>> print derive_ampllum((1.,'Lsol'),(1,'Msol'),(5777.,'K'),(550,'nm'))
    (4.7, 0.3)
    """
    #-- take care of mass
    if len(mass)==3:
        mass = (unumpy.uarray([mass[0],mass[1]]),mass[2])
    M = convert(mass[-1],'Msol',*mass[:-1])
    #-- take care of effective temperature
    if len(temperature)==3:
        temperature = unumpy.uarray([temperature[0],temperature[1]])
    teff = convert(temperature[-1],'5777K',*temperature[:-1])
    #-- take care of luminosity
    if len(luminosity)==3:
        luminosity = unumpy.uarray([luminosity[0],luminosity[1]])
    lumi = convert(luminosity[-1],'Lsol',*luminosity[:-1])
    #-- take care of wavelength
    if len(wavelength)==3:
        wavelength = unumpy.uarray([wavelength[0],wavelength[1]])
    wave = convert(wavelength[-1],'550nm',*wavelength[:-1])
    ampllum = lumi / wave / teff**2 / M * ufloat((4.7,0.3))
    return convert('ppm',unit,ampllum)

def derive_amplvel(luminosity,mass,unit='cm/s'):
    """
    Derive the luminosity amplitude around nu_max of solar-like oscillations.
    
    See Kjeldsen and Bedding (1995).
    
    >>> print derive_amplvel((1.,'Lsol'),(1,'Msol'),unit='cm/s')
    (23.4, 1.4)
    """
    #-- take care of mass
    if len(mass)==3:
        mass = (unumpy.uarray([mass[0],mass[1]]),mass[2])
    M = convert(mass[-1],'Msol',*mass[:-1])
    #-- take care of luminosity
    if len(luminosity)==3:
        luminosity = unumpy.uarray([luminosity[0],luminosity[1]])
    lumi = convert(luminosity[-1],'Lsol',*luminosity[:-1])
    amplvel = lumi / M * ufloat((23.4,1.4))
    return convert('cm/s',unit,amplvel)
#}



#{ Nonlinear change-of-base functions

class NonLinearConverter():
    """
    Base class for nonlinear conversions
    
    This class keeps track of prefix-factors and powers.
    
    To have a real nonlinear converter, you need to define the C{__call__}
    attribute.
    """
    def __init__(self,prefix=1.,power=1.):
        self.prefix = prefix
        self.power = power
    def __rmul__(self,other):
        if type(other)==type(5) or type(other)==type(5.):
            return self.__class__(prefix=self.prefix*other)
    def __div__(self,other):
        if type(other)==type(5) or type(other)==type(5.):
            return self.__class__(prefix=self.prefix*other)
    def __pow__(self,other):
        if type(other)==type(5) or type(other)==type(5.):
            return self.__class__(prefix=self.prefix,power=self.power+other)

class Fahrenheit(NonLinearConverter):
    """
    Convert Fahrenheit to Kelvin and back
    """
    def __call__(self,a,inv=False):
        if not inv: return (a*self.prefix+459.67)*5./9.
        else:       return (a*9./5.-459.67)/self.prefix

class Celcius(NonLinearConverter):
    """
    Convert Celcius to Kelvin and back
    """
    def __call__(self,a,inv=False):
        if not inv: return a*self.prefix+273.15
        else:       return (a-273.15)/self.prefix


class AmplMag(NonLinearConverter):
    """
    Convert a Vega magnitude to W/m2/m (Flambda) and back
    """
    def __call__(self,meas,inv=False):
        #-- this part should include something where the zero-flux is retrieved
        if not inv: return 10**(meas*self.prefix/2.5) - 1.
        else:       return (2.5*log10(1.+meas))/self.prefix

class VegaMag(NonLinearConverter):
    """
    Convert a Vega magnitude to W/m2/m (Flambda) and back
    """
    def __call__(self,meas,photband=None,inv=False,**kwargs):
        #-- this part should include something where the zero-flux is retrieved
        zp = filters.get_info()
        match = zp['photband']==photband.upper()
        if sum(match)==0: raise ValueError, "No calibrations for %s"%(photband)
        F0 = convert(zp['Flam0_units'][match][0],'W/m3',zp['Flam0'][match][0])
        mag0 = float(zp['vegamag'][match][0])
        if not inv: return 10**(-(meas-mag0)/2.5)*F0
        else:       return -2.5*log10(meas/F0)+mag0

class ABMag(NonLinearConverter):
    """
    Convert an AB magnitude to W/m2/Hz (Fnu) and back
    """
    def __call__(self,meas,photband=None,inv=False,**kwargs):
        zp = filters.get_info()
        F0 = 3.6307805477010024e-23
        match = zp['photband']==photband.upper()
        if sum(match)==0: raise ValueError, "No calibrations for %s"%(photband)
        mag0 = float(zp['ABmag'][match][0])
        if np.isnan(mag0): mag0 = 0.
        if not inv:
            try:
                return 10**(-(meas-mag0)/2.5)*F0
            except OverflowError:
                return np.nan
        else:       return -2.5*log10(meas/F0)

class STMag(NonLinearConverter):
    """
    Convert an ST magnitude to W/m2/m (Flambda) and back
    
    mag = -2.5*log10(F) - 21.10
    
    F0 = 3.6307805477010028e-09 erg/s/cm2/A
    """
    def __call__(self,meas,photband=None,inv=False,**kwargs):
        zp = filters.get_info()
        F0 = 0.036307805477010027
        match = zp['photband']==photband.upper()
        if sum(match)==0: raise ValueError, "No calibrations for %s"%(photband)
        mag0 = float(zp['STmag'][match][0])
        if np.isnan(mag0): mag0 = 0.
        if not inv: return 10**(-(meas-mag0)/-2.5)*F0
        else:       return -2.5*log10(meas/F0)


class Color(NonLinearConverter):
    """
    Convert a color to a flux ratio and back
    
    B-V = -2.5log10(FB) + CB - (-2.5log10(FV) + CV)
    B-V = -2.5log10(FB) + CB + 2.5log10(FV) - CV
    B-V = -2.5log10(FB/FV) + (CB-CV)
    
    and thus
    
    FB/FV = 10 ** [((B-V) - (CB-CV)) / (-2.5)]
    
    where
    
    CB = 2.5log10[FB(m=0)]
    CV = 2.5log10[FV(m=0)]
    
    Stromgren colour indices:
    
    m1 = v - 2b + y
    c1 = u - 2v + b
    Hbeta = HBN - HBW
    """
    def __call__(self,meas,photband=None,inv=False,**kwargs):
        #-- we have two types of colours: the stromgren M1/C1 type, and the
        #   normal Band1 - Band2 type. We need to have conversions back and
        #   forth: this translates into four cases.
        system,band = photband.split('.')
        if '-' in band and not inv:
            band0,band1 = band.split('-')
            f0 = convert('mag','SI',meas,photband='.'.join([system,band0]),unpack=False)
            f1 = convert('mag','SI',0.00,photband='.'.join([system,band1]))
            return f0/f1
        elif '-' in band and inv:
            #-- the units don't really matter, we choose SI'
            #   the flux ratio is converted to color by assuming that the
            #   denominator flux is equal to one.
            band0,band1 = band.split('-')
            m0 = convert('W/m3','mag',meas,photband='.'.join([system,band0]),unpack=False)
            m1 = convert('W/m3','mag',1.00,photband='.'.join([system,band1]))
            return m0-m1
        elif photband=='STROMGREN.C1' and not inv:
            fu = convert('mag','SI',meas,photband='STROMGREN.U',unpack=False)
            fb = convert('mag','SI',0.00,photband='STROMGREN.B')
            fv = convert('mag','SI',0.00,photband='STROMGREN.V')
            return fu*fb/fv**2
        elif photband=='STROMGREN.C1' and inv:
            mu = convert('W/m3','mag',meas,photband='STROMGREN.U',unpack=False)
            mb = convert('W/m3','mag',1.00,photband='STROMGREN.B')
            mv = convert('W/m3','mag',1.00,photband='STROMGREN.V')
            return mu-2*mv+mb
        elif photband=='STROMGREN.M1' and not inv:
            fv = convert('mag','SI',meas,photband='STROMGREN.V',unpack=False)
            fy = convert('mag','SI',0.00,photband='STROMGREN.Y')
            fb = convert('mag','SI',0.00,photband='STROMGREN.B')
            return fv*fy/fb**2
        elif photband=='STROMGREN.M1' and inv:
            mu = convert('W/m3','mag',meas,photband='STROMGREN.V',unpack=False)
            mb = convert('W/m3','mag',1.00,photband='STROMGREN.Y')
            mv = convert('W/m3','mag',1.00,photband='STROMGREN.B')
            return mv-2*mb+my
        else:
            raise ValueError, "No color calibrations for %s"%(photband)

class JulianDay(NonLinearConverter):
    """
    Convert a calender date to Julian date and back
    """
    def __call__(self,meas,inv=False,**kwargs):
        if inv:
            L= meas+68569
            N= 4*L//146097
            L= L-(146097*N+3)//4
            I= 4000*(L+1)//1461001
            L= L-1461*I//4+31
            J= 80*L//2447
            day = L-2447*J//80+0.5
            L= J//11
            month = J+2-12*L
            year = 100*(N-49)+I+L
            
            return year,month,day
        else:
            year,month,day = meas[:3]
            hour = len(meas)>3 and meas[3] or 0.
            mint = len(meas)>4 and meas[4] or 0.
            secn = len(meas)>5 and meas[5] or 0.    
            a = (14 - month)//12
            y = year + 4800 - a
            m = month + 12*a - 3
            jd = day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045
            jd += hour/24.
            jd += mint/24./60.
            jd += secn/24./3600.
            jd -= 0.5
            return jd

class ModJulianDay(NonLinearConverter):
    """
    Convert a Modified Julian Day to Julian Day  and back
    
    The CoRoT conversion has been checked with the CoRoT data archive: it is
    correct at least to the second (the archive tool's precision).
    """
    ZP = {'COROT':2451545.,
            'HIP':2440000.,
            'MJD':2400000.5}
    def __call__(self,meas,inv=False,jtype='MJD'):
        if inv:
            return meas-self.ZP[jtype.upper()]
        else:
            return meas+self.ZP[jtype.upper()]

class GalCoords(NonLinearConverter):
    """
    Convert Galactic coords to complex coords and back
    """
    def __call__(self,mycoord,inv=False,epoch='2000'):
        if inv:
            x,y = mycoord.real,mycoord.imag
            equ = ephem.Equatorial(x,y,epoch=epoch)
            gal = ephem.Galactic(equ,epoch=epoch)
            return gal.lon,gal.lat
        else:
            x,y = mycoord
            gal = ephem.Galactic(x,y,epoch=epoch)
            equ = ephem.Equatorial(gal,epoch=epoch)
            return float(equ.ra) + 1j*float(equ.dec)

class EquCoords(NonLinearConverter):
    """
    Convert Equatorial coords to complex coords and back
    """
    def __call__(self,mycoord,inv=False,epoch='2000'):
        if inv:
            x,y = mycoord.real,mycoord.imag
            equ = ephem.Equatorial(x,y,epoch=epoch)
            return equ.ra,equ.dec
        else:
            x,y = mycoord
            equ = ephem.Equatorial(x,y,epoch=epoch)
            return float(equ.ra) + 1j*float(equ.dec)

class EclCoords(NonLinearConverter):
    """
    Convert Ecliptic coords to complex coords and back
    """
    def __call__(self,mycoord,inv=False,epoch='2000'):
        if inv:
            x,y = mycoord.real,mycoord.imag
            equ = ephem.Equatorial(x,y,epoch=epoch)
            ecl = ephem.Ecliptic(equ,epoch=epoch)
            return ecl.long,ecl.lat
        else:
            x,y = mycoord
            ecl = ephem.Ecliptic(x,y,epoch=epoch)
            equ = ephem.Equatorial(ecl,epoch=epoch)
            return float(equ.ra) + 1j*float(equ.dec)

class DegCoords(NonLinearConverter):
    """
    Convert Complex coords to degrees coordinates and back
    """
    def __call__(self,mycoord,inv=False):
        if inv:
            x,y = mycoord.real,mycoord.imag
            return x/np.pi*180,y/np.pi*180
        else:
            x,y = mycoord
            return x/180.*np.pi + 1j*y/180.*np.pi

class RadCoords(NonLinearConverter):
    """
    Convert Complex coords to degrees coordinates and back
    """
    def __call__(self,mycoord,inv=False):
        if inv:
            x,y = mycoord.real,mycoord.imag
            return x,y
        else:
            x,y = mycoord
            return x + 1j*y

#}
#{ Currencies
@memoized
def set_exchange_rates():
    """
    Download currency exchange rates from the European Central Bank.
    """
    myurl = 'http://www.ecb.europa.eu/stats/eurofxref/eurofxref-daily.xml'
    url = urllib.URLopener()
    logger.info('Downloading current exchanges rates from ecb.europa.eu')
    filen,msg = url.retrieve(myurl)
    ff = open(filen,'r')
    for line in ff.readlines():
        if '<Cube currency=' in line:
            prefix,curr,interfix,rate,postfix = line.split("'")
            _factors[curr] = (1/float(rate),'EUR','currency','<some currency>')
    ff.close()
    
            
            
#}

#{ Computations with units
class Unit(object):
    """
    Trial class to calculate with numbers containing units. Not yet possible
    to supply an error.
    
    Initalisation is done via:
    
    >>> a = Unit(2.,'m')
    >>> b = Unit(4.,'km')
    >>> c = Unit(3.,'cm2')
    
    And you can calculate via:
    
    >>> print a*b
    8000.0 m2
    >>> print a/c
    6666.66666667 m-1
    >>> print a+b
    4002.0 m1
    
    For example, when you want to calculated the equatorial velocity of the sun,
    you could do:
    
    >>> distance = Unit(2*np.pi,'Rsol')
    >>> time = Unit(22.,'d')
    >>> print (distance/time)
    2299.03495719 m1 s-1
    
    or directly to km/s:
    
    >>> print (distance/time).convert('km/s')
    2.29903495719
    
    >>> distance = Unit((2*np.pi,0.1),'Rsol')
    >>> print (distance/time).convert('km/s')
    2.29903495719+/-0.0365902777778
    
    To compute the surface gravity of the sun:
    
    >>> G = Unit(constants.GG,constants.GG_units)
    >>> M = Unit(constants.Msol,constants.Msol_units)
    >>> R = Unit(constants.Rsol,constants.Rsol_units)
    >>> print np.log10((G*M/R**2).convert('cgs'))
    4.43830739117
    
    or 
    
    >>> G = Unit(constants.GG,constants.GG_units)
    >>> M = Unit(1.,'Msol')
    >>> R = Unit(1.,'Rsol')
    >>> print np.log10((G*M/R**2).convert('cgs'))
    4.43830739117
    
    
    """
    def __init__(self,value,unit=None,**kwargs):
        #-- input values
        kwargs.setdefault('unpack',False)
        self.value = value
        self.unit = unit
        self.kwargs = kwargs
        
        
        if isinstance(self.value,tuple):
            self.value = ufloat(self.value)
        
        #-- values and units to work with
        if self.unit is not None:
            self._SI_value = convert(self.unit,'SI',self.value,**kwargs)
        else:
            self._SI_value = self.value
        self._basic_unit = breakdown(self.unit)[1]
    
    def convert(self,unit):
        return convert(self.unit,unit,self.value,**self.kwargs)
    
    def __lt__(self,other):
        return self._SI_value<other._SI_value
    
    def __add__(self,other):
        if self._basic_unit!=other._basic_unit:
            raise ValueError,'unequal units %s and %s'%(self._basic_unit,other._basic_unit)
        return Unit(self._SI_value+other._SI_value,self._basic_unit)
        
    def __sub__(self,other):
        if self._basic_unit!=other._basic_unit:
            raise ValueError,'unequal units %s and %s'%(self._basic_unit,other._basic_unit)
        return Unit(self._SI_value-other._SI_value,self._basic_unit)
    
    def __mul__(self,other):
        if hasattr(other,'_basic_unit'):
            new_unit = ' '.join([self._basic_unit,other._basic_unit])
            fac,new_unit = breakdown(new_unit)
            outcome = self._SI_value*other._SI_value
        else:
            outcome = other*self.value
            new_unit = self.unit
        return Unit(outcome,new_unit)
    
    def __rmul__(self,other):
        return self.__mul__(other)
    
    def __div__(self,other):
        if hasattr(other,'_basic_unit'):
            #-- reverse signs in second units
            uni_b_ = ''
            isalpha = True
            prev_char_min = False
            for i,char in enumerate(other._basic_unit):
                if char=='-':
                    prev_char_min = True
                    continue
                if isalpha and not char.isalpha() and not prev_char_min:
                    uni_b_ += '-'
                prev_char_min = False
                uni_b_ += char
                isalpha = char.isalpha()
            new_unit = ' '.join([self._basic_unit,uni_b_])
            fac,new_unit = breakdown(new_unit)
            outcome = self._SI_value/other._SI_value
        else:
            outcome = self.value/other
            new_unit = self.unit
        return Unit(outcome,new_unit)
    
    def __rdiv__(self,other):
        if hasattr(other,'_basic_unit'):
            #-- reverse signs in second units
            uni_b_ = ''
            isalpha = True
            for i,char in enumerate(other._basic_unit):
                if char=='-': continue
                if isalpha and not char.isalpha():
                    uni_b_ += '-'
                uni_b_ += char
                isalpha = char.isalpha()
            new_unit = ' '.join([self._basic_unit,uni_b_])
            fac,new_unit = breakdown(new_unit)
            outcome = other._SI_value/self._SI_value
        else:
            outcome = other/self.value
            new_unit = self.unit
        return Unit(outcome,new_unit)
    
    def __pow__(self,power):
        mycomps = [components(u) for u in self._basic_unit.split()]
        mycomps = [(u[0]**power,u[1],u[2]*power) for u in mycomps]
        factor = np.product([u[0] for u in mycomps])
        new_unit = ' '.join(['%s%g'%(u[1],u[2]) for u in mycomps])
        fac,new_unit = breakdown(new_unit)
        return Unit(self._SI_value**power,new_unit)
        #new_unit = ' '.join(power*[self._basic_unit])
        #fac,new_unit = breakdown(new_unit)
        #return Unit(self._SI_value**power,new_unit)
    
    def __str__(self):
        return '%s %s'%(self.value,self.unit)
    
def usin(unit):
    value = sin(unit.convert('rad'))
    return value
        
        

#}


            
_fluxcalib = os.path.join(os.path.abspath(os.path.dirname(__file__)),'fluxcalib.dat')
#-- basic units which the converter should know about
_factors = {
# DISTANCE
           'm':     (  1e+00,       'm','distance','meter'), # meter
           'A':     (  1e-10,       'm','distance','angstrom'), # Angstrom
           'AU':    (constants.au,    constants.au_units,'distance','astronomical unit'), # astronomical unit
           'pc':    (constants.pc,    constants.pc_units,'distance','parsec'), # parsec
           'ly':    (constants.ly,    constants.ly_units,'distance','light year'), # light year
           'Rsol':  (constants.Rsol,  constants.Rsol_units,'distance','Solar radius'), # Solar radius
           'Rearth':(constants.Rearth,constants.Rearth_units,'distance','Earth radius'), # Earth radius
           'ft':    (0.3048,        'm','distance','foot (international)'), # foot (international)
           'in':    (0.0254,        'm','distance','inch (international)'), # inch (international)
           'mi':    (1609.344,      'm','distance','mile (international)'), # mile (international)
           'a0':    (constants.a0,  constants.a0_units,'distance','Bohr radius'), # Bohr radius
           'ell':   (1.143,         'm','distance','ell'), # ell
           'yd':    (0.9144,        'm','distance','yard (international)'), # yard (international)
# MASS
           'g':     (  1e-03,       'kg','mass','gram'), # gram
           'amu':   (1.66053892173e-27,'kg','mass','atomic mass'), # atomic mass unit (wikipedia)
           'Msol':  (constants.Msol,   constants.Msol_units,'mass','Solar mass'), # Solar mass
           'Mearth':(constants.Mearth, constants.Mearth_units,'mass','Earth mass'), # Earth mass
           'Mjup':  (constants.Mjup,   constants.Mjup_units,'mass','Jupiter mass'), # Jupiter mass
           'Mlun':  (constants.Mlun,   constants.Mlun_units,'mass','Lunar mass'), # Lunar mass
           'lbs':   (0.45359237,    'kg','mass','pound'), # pound
           'st':    (6.35029318,    'kg','mass','stone'), # stone
           'ounce': (0.0283495231,  'kg','mass','ounce'), # ounce
           'mol':   (1./constants.NA,'mol','mass','molar mass'), # not really a mass...
# TIME
           's':     (  1e+00,       's','time','second'),     # second
           'min':   (  60.,         's','time','minute'),     # minute
           'h':     (3600.,         's','time','hour'),     # hour 
           'd':     (24*3600.,      's','time','day'),     # day
           'wk':    (7*24*3600.,    's','time','week'),     # week
           'mo':    (30*7*24*3600., 's','time','month'),     # month
           'sidereal': (1.0027379093,'','time','sidereal day'),     # sidereal
           'yr':    (3.1558149984e7,'s','time','year'),     # year
           'cr':    (100*365*24*3600,'s','time','century'),    # century
           'hz':    (1e+00,         'cy s-1','time','Hertz'),# Hertz
           'JD':    (1e+00,         'JD','time','Julian day'), # Julian Day
           'CD':    (JulianDay,     'JD','time','calender day'), # Calender Day
           'MJD':   (ModJulianDay,  'JD','time','modified Julian day'), # Modified Julian Day
           'j':     (1/60.,         's','time','jiffy'),  # jiffy
# ANGLES
           'rad':         (1e+00,               'rad','angle','radian'),  # radian
           'cy':          (1e+00,               'cy','angle','cycle'),   # cycle
           'deg':         (np.pi/180.,          'rad','angle','degree'),  # degree
           'am':          (np.pi/180./60.,      'rad','angle','arcminute'),  # arcminute
           'as':          (np.pi/180./3600.,    'rad','angle','arcsecond'),  # arcsecond
           'sr':          (1,                   'rad2','angle','sterradian'), # sterradian #1/39.4784176045
           'rpm':         (0.104719755,         'rad/s','angle','revolutions per minute'),# revolutions per minute
# COORDINATES
           'complex_coord':(1e+00+0*1j, 'complex_coord','coordinate','<own unit>'), # own unit
           'equ':          (EquCoords,  'complex_coord','coordinate','equatorial'), # Equatorial coordinates
           'gal':          (GalCoords,  'complex_coord','coordinate','galactic'), # Galactic coordinates
           'ecl':          (EclCoords,  'complex_coord','coordinate','ecliptic'), # Ecliptic coordinates
           'deg_coord':    (DegCoords,  'complex_coord','coordinate','degrees'), # Coordinates in degrees
           'rad_coord':    (RadCoords,  'complex_coord','coordinate','radians'), # Coordinates in radians
# FORCE
           'N':     (1e+00,         'kg m s-2','force','Newton'), # newton
           'dyn':   (1e-05,         'kg m s-2','force','dyne'), # dyne
# TEMPERATURE
           'K':      (1e+00,        'K','temperature','Kelvin'), # Kelvin
           'F':      (Fahrenheit,   'K','temperature','Fahrenheit'), # Fahrenheit
           'C':      (Celcius,      'K','temperature','Celcius'), # Celcius
           'Tsol':   (constants.Tsol,constants.Tsol_units,'temperature','Solar temperature'), # solar temperature
# ENERGY & POWER
           'J':     (  1e+00,       'kg m2 s-2','energy/power','Joule'), # Joule
           'W':     (  1e+00,       'kg m2 s-3','energy/power','Watt'), # Watt
           'erg':   (  1e-07,       'kg m2 s-2','energy/power','ergon'), # ergon
           'eV':    (1.60217646e-19,'kg m2 s-2','energy/power','electron volt'), # electron volt
           'cal':   (4.1868,        'kg m2 s-2','energy/power','calorie (international table)'),# calorie (International table)
           'Lsol':  (constants.Lsol, constants.Lsol_units,'energy/power','Solar luminosity'), # solar luminosity
           'hp':    (745.699872,    'kg m2 s-3','energy/power','Horsepower'), # horsepower
# PRESSURE
           'Pa':    (  1e+00,       'kg m-1 s-2','pressure','Pascal'), # Pascal
           'bar':   (  1e+05,       'kg m-1 s-2','pressure','baros'), # baros
           'at':    (  98066.5,     'kg m-1 s-2','pressure','atmosphere (technical)'), # atmosphere (technical)
           'atm':   ( 101325,       'kg m-1 s-2','pressure','atmosphere (standard)'), # atmosphere (standared)
           'torr':  (    133.322,   'kg m-1 s-2','pressure','Torricelli'), # Torricelli
           'psi':   (   6894.,      'kg m-1 s-2','pressure','pound per square inch'), # pound per square inch
# AREA
           'ac':    (4046.8564224,  'm2','area','acre (international)'), # acre (international)
           'a':     (100.,          'm2','area','are'), # are
# FLUX
# -- absolute magnitudes
           'Jy':      (1e-26,         'kg s-2 cy-1','flux','Jansky'), # W/m2/Hz
           'vegamag': (VegaMag,       'kg m-1 s-3','flux','Vega magnitude'),  # W/m2/m
           'mag':     (VegaMag,       'kg m-1 s-3','flux','magnitude'),  # W/m2/m
           'STmag':   (STMag,         'kg m-1 s-3','flux','ST magnitude'),  # W/m2/m
           'ABmag':   (ABMag,         'kg s-2 cy-1','flux','AB magnitude'), # W/m2/Hz
# -- magnitude differences (colors)
           'mag_color':(Color,         'flux_ratio','flux','color'),
           'flux_ratio':(1+00,         'flux_ratio','flux','flux ratio'),
# -- magnitude amplitudes
           'ampl':    (1e+00,         'ampl','flux','fractional amplitude'),
           'Amag':    (AmplMag,       'ampl','flux','amplitude in magnitude'),
           'pph':     (1e-02,         'ampl','flux','amplitude in parts per hundred'), # amplitude
           'ppt':     (1e-03,         'ampl','flux','amplitude in parts per thousand'), # amplitude
           'ppm':     (1e-06,         'ampl','flux','amplitude in parts per million'), # amplitude
# -- currency
           'EUR':     (1e+00,         'EUR','currency','EURO')
           }
#-- set of conventions:
_conventions = {'SI': dict(mass='kg',length='m', time='s',temperature='K',
                          electric_current='ampere',lum_intens='cd',amount='mol'), # International standard
               'cgs':dict(mass='g', length='cm',time='s',temperature='K',
                          electric_current='ampere',lum_intens='cd',amount='mol'), # Centi-gramme-second
               'sol':dict(mass='Msol',length='Rsol',time='s',temperature='Tsol',
                          electric_current='ampere',lum_intens='cd',amount='mol'), # solar
               'imperial':dict(mass='lbs',length='yd',time='s',temperature='K',
                          electric_current='ampere',lum_intens='cd',amount='mol'), # Imperial (UK/US) system
               }           
           
#-- scaling factors for prefixes            
_scalings ={'y':       1e-24, # yocto
            'z':       1e-21, # zepto
            'a':       1e-18, # atto
            'f':       1e-15, # femto
            'p':       1e-12, # pico
            'n':       1e-09, # nano
            'mu':      1e-06, # micro
            'u':       1e-06, # micro
            'm':       1e-03, # milli
            'c':       1e-02, # centi
            'd':       1e-01, # deci
            'da':      1e+01, # deca
            'h':       1e+02, # hecto
            'k':       1e+03, # kilo
            'M':       1e+06, # mega
            'G':       1e+09, # giga
            'T':       1e+12, # tera
            'P':       1e+15, # peta
            'E':       1e+18, # exa
            'Z':       1e+21, # zetta
            'Y':       1e+24  # yotta
            }
 
#-- some common aliases
_aliases = [('micron','mum'),('au','AU'),
            ('micro','mu'),('milli','m'),('kilo','k'),('mega','M'),('giga','G'),
            ('nano','n'),
            ('watt','W'),('Watt','W'),
            ('Hz','hz'),
            ('joule','J'),('Joule','J'),
            ('jansky','Jy'),('Jansky','Jy'),('jy','Jy'),
            ('arcsec','as'),('arcmin','am'),
            ('cycles','cy'),('cycle','cy'),('cyc','cy'),
            ('angstrom','A'),('Angstrom','A'),
            ('inch','in'),
            ('^',''),('**',''),
            ('celcius','C'),('fahrenheit','F'),('hr','h'),
            ('galactic','gal'),('equatorial','equ'),('ecliptic','ecl'),
            ('Vegamag','vegamag'),('mile','mi'),
            ('oz','ounce'),
            ('pk','hp'),('mph','mi/h')
            ]
 
#-- Change-of-base function definitions
_switch = {'s1_to_':       distance2velocity, # switch from wavelength to velocity
           's-1_to_':      velocity2distance, # switch from wavelength to velocity
           'cy-1m1_to_':   distance2spatialfreq,  # for interferometry
           'cy1m-1_to_':   spatialfreq2distance, # for interferometry
           'cy-1m1s1_to_': fnu2flambda, # switch from Fnu to Flam
           'cy1m-1s-1_to_':flambda2fnu, # switch from Flam to Fnu
           'cy-1s1_to_':   fnu2nufnu, # switch from Fnu to nuFnu
           'cy1s-1_to_':   nufnu2fnu, # switch from nuFnu to Fnu
           'm1_to_':       lamflam2flam, # switch from lamFlam to Flam
           'm-1_to_':      flam2lamflam, # switch from Flam to lamFlam
           'rad2_to_':     per_sr,
           'rad-2_to_':    times_sr,
           'rad1_to_':     per_cy,
           'rad-1_to_':    times_cy,
           'cy-1s2_to_':   period2freq, # same for both, since just inverse
           'cy1s-2_to_':   period2freq,
           'cy1_to_':      do_nothing,
           'cy-1_to_':     do_nothing,
           'cy2_to_':      do_nothing,
           'cy-2_to_':     do_nothing}
 
 
if __name__=="__main__":
    if not sys.argv[1:]:
        import doctest
        doctest.testmod()
        quit()
    from optparse import OptionParser, Option, OptionGroup
    import datetime
    import copy
    
    logger = loggers.get_basic_logger()
    
    #-- make sure we can parse strings as Python code
    def check_pythoncode(option, opt, value):
        try:
            return eval(value)
        except ValueError:
            raise OptionValueError(
                "option %s: invalid python code: %r" % (opt, value))

    #-- generate a custom help log
    class MyOption (Option):
        TYPES = Option.TYPES + ("pythoncode",)
        TYPE_CHECKER = copy.copy(Option.TYPE_CHECKER)
        TYPE_CHECKER["pythoncode"] = check_pythoncode   
    usage = "List of available units:\n" + get_help() + "\nUsage: %prog --from=<unit> --to=<unit> [options] value [error]"
    
    #-- define all the input parameters
    parser = OptionParser(option_class=MyOption,usage=usage)
    parser.add_option('--from',dest='_from',type='str',
                        help="units to convert from",default=None)
    parser.add_option('--to',dest='_to',type='str',
                        help="units to convert to",default='SI')
    
    group = OptionGroup(parser,'Extra quantities when changing base units (e.g. Flambda to Fnu)')
    group.add_option('--wave','-w',dest='wave',type='str',
                        help="wavelength with units (e.g. used to convert Flambda to Fnu)",default=None)
    group.add_option('--freq','-f',dest='freq',type='str',
                        help="frequency (e.g. used to convert Flambda to Fnu)",default=None)
    group.add_option('--photband','-p',dest='photband',type='str',
                        help="photometric passband",default=None)
    parser.add_option_group(group)
    
    #-- prepare inputs for functions
    (options, args) = parser.parse_args()
    options = vars(options)
    _from = options.pop('_from')
    _to = options.pop('_to')
    #-- in case of normal floats or floats with errors
    if not any([',' in i for i in args]):
        args = tuple([float(i) for i in args])
    #-- in case of tuples (like coordinates)
    else:
        args = tuple((eval(args[0]),))
    #-- remove None types
    for option in copy.copy(options):
        if options[option] is None:
            options.pop(option)
    #-- set type correctly
    for option in options:
        if isinstance(options[option],str) and ',' in options[option]:
            entry = options[option].split(',')
            options[option] = (float(entry[0]),entry[1])
    #-- check if currencies are asked. If so, download the latest exchange rates
    if (_from.isupper() and len(_from)==3) and (_to.isupper() and (len(_to)==3 or _to=='SI')):
        set_exchange_rates()
    #-- do the conversion
    output = convert(_from,_to,*args,**options)
    
    #-- and nicely print to the screen
    if _to=='SI':
        fac,_to = breakdown(_from)
    if isinstance(output,tuple) and len(output)==2 and len(args)==2:
        print "%g +/- %g %s    =    %g +/- %g %s"%(args[0],args[1],_from,output[0],output[1],_to)
    elif isinstance(output,tuple) and len(output)==2:
        print "%s %s    =    %s,%s %s"%(args[0],_from,output[0],output[1],_to)
    elif _to.lower()=='cd':
        year,month,day = output
        day,fraction = int(day),day-int(day)
        hour = fraction*24
        hour,fraction = int(hour),hour-int(hour)
        minute = fraction*60
        minute,fraction = int(minute),minute-int(minute)
        second = int(fraction*60)
        dt = datetime.datetime(year,month,day,hour,minute,second)
        print "%g %s    =    %s %s"%(args[0],_from,dt,_to)
    else:
        print "%g %s    =    %g %s"%(args[0],_from,output,_to)