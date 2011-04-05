# -*- coding: utf-8 -*-
"""
Convert one unit to another.

Some of the  many possibilities include:
    
    1. Conversions between equal-type units: meter to nano-lightyears, erg/s
    to W, cy/d to muHz, but also erg/s/cm2/A to W/m2/mum, sr to deg2, etc...
    2. Conversions between unequal-type units: angstrom to km/s via the speed
    of light, F(lambda) to F(nu), F(nu) to lambdaF(lambda)/sr, meter to
    cycles/arcsec (interferometry), etc...
    3. Nonlinear conversions: vegamag to erg/s/cm2/A or Jy, Celcius to
    Fahrenheit or Kelvin, calender date to (modified) Julian Day, Equatorial to
    Ecliptic coordinates, etc...

The main function C{convert} does all the work and is called via

C{result = convert('km','m',1.)}

Be B{careful} when mixing nonlinear conversions (e.g. magnitude to flux) with
linear conversions (e.g. Jy to W/m2/m).

Note: when your favorite conversion is not implemented, there are four places
where you can add information:

    1. C{_scalings}: if your favorite prefix (e.g. Tera, nano...) is not
    available
    2. C{_aliases}: if your unit is available but not under the name you are
    used to.
    3. C{_factors}: if your unit is not available.
    4. C{_switch}: if your units are available, but conversions from one to
    another is not straightforward and extra infromation is needed (e.g. to go
    from angstrom to km/s in a spectrum, a reference wavelength is needed).

If you need to add a linear factor, just give the factor in SI units, and the
SI base units it consists of. If you need to add a nonlinear factor, you have
to give a class definition (see the examples).
"""
#-- standard libraries
import re
import os
import numpy as np
try: import ephem
except ImportError: print("Unable to load pyephem, coordinate transfos unavailable")

#-- from IVS repository
from ivs.units.constants import *
from ivs.io import ascii

#{ Main functions

def convert(_from,_to,*args,**kwargs):
    """
    Convert one unit to another.
    
    Keyword arguments can give extra information, for example when converting
    from Flambda to Fnu, and should be tuples (float,'unit').
    
    The unit strings should by default be given in the form
    
    C{erg s-1 cm-2 A-1}
    
    Common alternatives are also accepted, but don't drive this to far:
    
    C{erg/s/cm2/A}
    
    The crasiest you're allowed to go is
    
    >>> print(convert('10mW m-2/nm','erg s-1 cm-2 A-1',1.))
    1.0
    
    Parentheses are in no circumstances accepted. Some common aliases are also
    resolved (for a full list, see dictionary C{_aliases}):
    
    C{erg/s/cm2/angstrom}
    
    B{WARNINGS}:
        1. the conversion involving sr and pixels is B{not tested}.
        2. the conversion involving magnitudes is not tested and not calibrated
    
    Examples:
    
    B{Spectra}:
    
    >>> convert('km','cm',1.)
    100000.0
    >>> convert('A','km/s',4553.,wave=(4552.,'A'))
    65.859503075576129
    >>> convert('nm','m/s',455.3,wave=(0.4552,'mum'))
    65859.503075645873
    >>> convert('km/s','A',65.859503075576129,wave=(4552.,'A'))
    4553.0
    >>> convert('nm','Ghz',1000.)
    299792.45799999993
    >>> convert('km h-1','nRsun s-1',1.)
    0.39939292275740873
    >>> convert('erg s-1 cm-2 A-1','SI',1.)
    10000000.0
    
    B{Fluxes}:
    
    >>> convert('erg/s/cm2/A','Jy',1e-10,wave=(10000.,'angstrom'))
    333.56409519815202
    >>> convert('erg/s/cm2/A','Jy',1e-10,freq=(cc/1e-6,'hz'))
    333.56409519815202
    >>> convert('erg/s/cm2/A','Jy',1e-10,freq=(cc,'Mhz'))
    333.56409519815202
    >>> convert('Jy','erg/s/cm2/A',333.56409519815202,wave=(10000.,'A'))
    1e-10
    >>> convert('Jy','erg/s/cm2/A',333.56409519815202,freq=(cc,'Mhz'))
    1e-10
    >>> convert('W/m2/mum','erg/s/cm2/A',1e-10,wave=(10000.,'A'))
    1.0000000000000001e-11
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
    667128190.39630413
    >>> convert('Jy','erg/s/cm2/micron/sr',1.,wave=(2.,'micron'),ang_diam=(3.,'mas'))
    4511059.8297938667
    >>> convert('Jy','erg/s/cm2/micron/sr',1.,wave=(2.,'micron'),pix=(3.,'mas'))
    3542978.1052961089
    >>> convert('erg/s/cm2/micron/sr','Jy',1.,wave=(2.,'micron'),ang_diam=(3.,'mas'))
    2.2167739682709884e-07
    >>> convert('Jy','erg/s/cm2/micron',1.,wave=(2,'micron'))
    7.4948114500000012e-10
    >>> print(convert('10mW m-2 nm-1','erg s-1 cm-2 A-1',1.))
    1.0
    >>> print convert('Jy','erg s-1 cm-2 micron-1 sr-1',1.,ang_diam=(2.,'mas'),wave=(1.,'micron'))
    40599538.4681
    
    B{Angles}:
    >>> convert('sr','deg2',1.)
    3282.8063499998884
    
    B{Magnitudes}:
    
    >>> print(convert('ABmag','Jy',0.,photband='SDSS.U'))
    3767.03798984
    >>> print(convert('Jy','erg cm-2 s-1 A-1',3630.7805477,wave=(1.,'micron')))
    1.08848062485e-09
    >>> print(convert('ABmag','erg cm-2 s-1 A-1',0.,wave=(1.,'micron'),photband='SDSS.G'))
    1.08848062485e-09
    >>> print(convert('erg cm-2 s-1 A-1','ABmag',1e-8,wave=(1.,'micron'),photband='SDSS.G'))
    -2.40794824268
    
    B{Frequency analysis}:
    
    >>> convert('cy/d','muHz',1.)
    11.574074074074074
    >>> convert('muhz','cy/d',11.574074074074074)
    1.0
    
    B{Interferometry}:
    
    >>> convert('m','cy/arcsec',85.,wave=(2.2,'micron'))
    187.3143767923207
    >>> convert('cm','cy/arcmin',8500.,wave=(2200,'nm'))/60.
    187.31437679232073
    >>> convert('cy/arcsec','m',187.,wave=(2.2,'mum'))
    84.857341290055444
    >>> convert('cyc/arcsec','m',187.,wave=(1,'mum'))
    38.571518768207014
    >>> convert('cycles/arcsec','m',187.,freq=(300000.,'Ghz'))
    38.54483473437972
    >>> convert('cycles/mas','m',0.187,freq=(300000.,'Ghz'))
    38.544834734379712
    
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
    2446257.8100000001
    >>> convert('CD','JD',(1985,7,11,7,31,59))
    2446257.8138773148
    >>> convert('MJD','CD',0.,jtype='corot')
    (2000.0, 1.0, 1.5)
    >>> convert('JD','MJD',2400000.5,jtype='mjd')
    0.0
    >>> convert('MJD','CD',0,jtype='mjd')
    (1858.0, 11.0, 17.0)
    
    B{Coordinates}:
    
    >>> convert('equ','gal',('17:45:40.4','-29:00:28.1'),epoch='2000')
    (6.282224277178722, -0.00082517883389919317)
    >>> convert('gal','equ',('00:00:00.00','00:00:00.0'),epoch='2000')
    (4.6496443030366299, -0.50503150853426648)
    
    @param _from: units to convert from
    @type _from: str
    @param _to: units to convert to
    @type _to: str
    @return: converted value
    @rtype: float
    """
    #-- break down the from and to units to their basic elements
    fac_from,uni_from = breakdown(_from)
    if _to!='SI':
        fac_to,uni_to = breakdown(_to)
    else:
        fac_to,uni_to = 1.,uni_from
    
    #-- convert the kwargs to SI units if they are tuples
    kwargs_SI = {}
    for key in kwargs:
        if isinstance(kwargs[key],tuple):
            fac_key,uni_key = breakdown(kwargs[key][1])
            kwargs_SI[key] = fac_key*kwargs[key][0]
        else:
            kwargs_SI[key] = kwargs[key]
    
    #-- easy if same units
    ret_value = 1.
    if uni_from==uni_to:
        #-- if nonlinear conversions from or to:
        if isinstance(fac_from,NonLinearConverter):
            ret_value *= fac_from(args[0],**kwargs_SI)
        else:
            ret_value *= fac_from*args[0]
    #-- otherwise a little bit more complicated
    else:
        uni_from_ = uni_from.split()
        uni_to_ = uni_to.split()
        only_from = "".join(sorted(list(set(uni_from_) - set(uni_to_))))
        only_to = "".join(sorted(list(set(uni_to_) - set(uni_from_))))
        
        #-- especially for conversion from and to sterradians
        if 'cy-2' in only_to:
            args = _switch['_to_cy-2'](args[0],**kwargs_SI),
            only_to = only_to.replace('cy-2','')
        if 'cy-2' in only_from:
            args = _switch['cy-2_to_'](args[0],**kwargs_SI),
            only_from = only_from.replace('cy-2','')
        
        #-- nonlinear conversions need a little tweak
        try:
            key = '%s_to_%s'%(only_from,only_to)
            if isinstance(fac_from,NonLinearConverter):
                ret_value *= _switch[key](fac_from(args[0],**kwargs_SI),**kwargs_SI)
            #-- linear conversions are easy
            else:
                ret_value *= _switch[key](fac_from*args[0],**kwargs_SI)
        except KeyError:
            raise KeyError,'cannot convert %s to %s: no %s definition in dict _switch'%(only_from,only_to and only_to or '[DimLess]',key)
    
    #-- final step: convert to ... (again distinction between linear and
    #   nonlinear converters)
    if isinstance(fac_to,NonLinearConverter):
        ret_value = fac_to(ret_value,inv=True,**kwargs_SI)
    else:
        ret_value /= fac_to
    
    return ret_value

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
                m = re.search(r'(\d*)(.+?)(-{0,1}\d+)',after_div)
                if m is not None:
                    factor,basis,power = m.group(1),m.group(2),int(m.group(3))
                    if factor: factor = float(factor)
                    else: factor = 1.
                else:
                    factor,basis,power = 1.,after_div,1
                uni[1+j] = '%s%d'%(basis,-power)
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
    (0.10000000000000001, 'kg', 3)
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
    if not unit[-1].isdigit(): unit += '1'
    #-- decompose unit in base name and power
    m = re.search(r'(\d*)(.+?)(-{0,1}\d+)',unit)
    if m is not None:
        factor,basis,power = m.group(1),m.group(2),int(m.group(3))
        if factor: factor = float(factor)
        else: factor = 1.
    else:
        factor,basis,power = 1.,unit,1
    #-- decompose the base name (which can be a composition of a prefix
    #   (e.g., 'mu') and a unit name (e.g. 'm')) into prefix and unit name
    #-- check if basis is part of _factors dictionary. If not, find the
    #   combination of _scalings and basis which is inside the dictionary!
    for scale in _scalings:
        scale_unit,base_unit = basis[:len(scale)],basis[len(scale):]
        if scale_unit==scale and base_unit in _factors and not basis in _factors:
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
        velocity = (arg-wave) / wave * cc
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
        distance = wave / cc * arg + wave
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
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'wave' in kwargs:
        wave = kwargs['wave']
        flambda = cc/wave**2 * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        flambda = freq**2/cc * arg
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
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI,W/m2/Hz)
    @rtype: float
    """
    if 'wave' in kwargs:
        wave = kwargs['wave']
        fnu = wave**2/cc * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        fnu = cc/freq**2 * arg
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
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'wave' in kwargs:
        wave = kwargs['wave']
        fnu = cc/wave * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        fnu = freq/cc * arg
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
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spectral irradiance (SI, W/m2/m)
    @rtype: float
    """
    if 'wave' in kwargs:
        wave = kwargs['wave']
        fnu = wave/cc * arg
    elif 'freq' in kwargs:
        freq = kwargs['freq']
        fnu = cc/freq * arg
    else:
        raise ValueError,'reference wave/freq not given'
    return fnu

def distance2frequency(arg,**kwargs):
    """
    Switch from distance to frequency via the speed of light, or vice versa.
    
    @param arg: distance (SI, m)
    @type arg: float
    @return: frequency (SI, Hz)
    @rtype: float
    """
    return cc/arg

def distance2spatialfreq(arg,**kwargs):
    """
    Switch from distance to spatial frequency via a reference wavelength.
    
    @param arg: distance (SI, m)
    @type arg: float
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: spatial frequency (SI, cy/as)
    @rtype: float
    """
    if 'wave' in kwargs:
        spatfreq = 2*np.pi*arg/kwargs['wave']
    elif 'freq' in kwargs:
        spatfreq = 2*np.pi*arg*cc*kwargs['freq']
    else:
        raise ValueError,'reference wave/freq not given'
    return spatfreq

def spatialfreq2distance(arg,**kwargs):
    """
    Switch from spatial frequency to distance via a reference wavelength.
    
    @param arg: spatial frequency (SI, cy/as)
    @type arg: float
    @keyword wave: reference wavelength (SI, m)
    @type wave: float
    @keyword freq: reference frequency (SI, Hz)
    @type freq: float
    @return: distance (SI, m)
    @rtype: float
    """
    if 'wave' in kwargs:
        distance = kwargs['wave']*arg/(2*np.pi)
    elif 'freq' in kwargs:
        distance = cc/kwargs['freq']*arg/(2*np.pi)
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
#}

#{ Nonlinear change-of-base functions
def read_fluxcalib():
    dtypes = [('PHOTBAND','a50'),
              ('VEGAMAG',np.float),
              ('ABMAG',np.float),
              ('STMAG',np.float),
              ('F0',np.float)]
    data = ascii.read2recarray(_fluxcalib,dtype=dtypes)
    return data

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

class VegaMag(NonLinearConverter):
    """
    Convert a Vega magnitude to W/m2/m (Flambda) and back
    """
    def __call__(self,meas,photband=None,inv=False):
        #-- this part should include something where the zero-flux is retrieved
        data = read_fluxcalib()
        match = data['PHOTBAND']==photband.upper()
        if sum(match)==0: raise ValueError, "No calibrations for %s"%(photband)
        F0 = data['F0'][match][0]
        if not inv: return 10**(-meas/2.5)*F0
        else:       return -2.5*np.log10(meas/F0)

class ABMag(NonLinearConverter):
    """
    Convert an AB magnitude to W/m2/Hz (Fnu) and back
    """
    def __call__(self,meas,photband=None,inv=False,**kwargs):
        data = read_fluxcalib()
        F0 = 3.6307805477010024e-23
        match = data['PHOTBAND']==photband.upper()
        if sum(match)==0: raise ValueError, "No calibrations for %s"%(photband)
        mag0 = data['ABMAG'][match][0]
        if np.isnan(mag0): mag0 = 0.
        if not inv: return 10**(-(meas-mag0)/2.5)*F0
        else:       return -2.5*np.log10(meas/F0)

class STMag(NonLinearConverter):
    """
    Convert an ST magnitude to W/m2/m (Flambda) and back
    """
    def __call__(self,meas,photband=None,inv=False):
        data = read_fluxcalib()
        F0 = 0.036307805477010027
        match = data['PHOTBAND']==photband.upper()
        if sum(match)==0: raise ValueError, "No calibrations for %s"%(photband)
        mag0 = data['STMAG'][match][0]
        if np.isnan(mag0): mag0 = 0.
        if not inv: return 10**(-(meas-mag0)/-2.5)*F0
        else:       return -2.5*np.log10(meas/F0)

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
    """
    ZP = {'COROT':2451545.,
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
            return gal.long,gal.lat
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

_fluxcalib = os.path.join(os.path.abspath(os.path.dirname(__file__)),'fluxcalib.dat')
#-- basic units which the converter should know about
_factors = {
# DISTANCE
           'm':     (  1e+00,       'm'),
           'A':     (  1e-10,       'm'),
           'AU':    (au,            'm'),
           'pc':    (pc,            'm'),
           'ly':    (ly,            'm'),
           'Rsun':  (Rsun,          'm'),
           'ft':    (0.3048,        'm'),
           'in':    (0.0254,        'm'),
           'mi':    (1609.344,      'm'),
# MASS
           'g':     (  1e-03,       'kg'),
           'Msun':  (Msun,          'kg'),
# TIME
           's':     (  1e+00,       's'),
           'min':   (  60.,         's'),
           'h':     (3600.,         's'),
           'd':     (24*3600.,      's'),
           'sidereal': (1.0027379093,''),
           'yr':    (365*24*3600.,  's'),
           'cr':    (100*365*24*3600,'s'),
           'hz':    (1e+00,         'cy s-1'),
           'JD':    (1e+00,         'JD'), # Julian Day
           'CD':    (JulianDay,     'JD'), # Calender Day
           'MJD':   (ModJulianDay,  'JD'), # Modified Julian Day
# ANGLES
           'rad':         (0.15915494309189535, 'cy'),
           'cy':          (1e+00,               'cy'),
           'deg':         (1./360.,             'cy'),
           'am':          (1./360./60.,         'cy'),
           'as':          (1./360./3600.,       'cy'),
           'sr':          (1/39.4784176045,     'cy2'),
# COORDINATES
           'complex_coord':(1e+00+0*1j, 'complex_coord'),
           'equ':          (EquCoords,  'complex_coord'),
           'gal':          (GalCoords,  'complex_coord'),
           'ecl':          (EclCoords,  'complex_coord'),
# FORCE
           'N':     (1e+00,         'kg m s-2'),
           'dy':    (1e-05,         'kg m s-2'),
# TEMPERATURE
           'K':      (1e+00,        'K'),
           'F':      (Fahrenheit,   'K'),
           'C':      (Celcius,      'K'),
# ENERGY & POWER
           'J':     (  1e+00,       'kg m2 s-2'),
           'W':     (  1e+00,       'kg m2 s-3'),
           'erg':   (  1e-07,       'kg m2 s-2'),
           'eV':    (1.60217646e-19,'kg m2 s-2'),
           'cal':   (4.184,         'kg m2 s-2'),
           'Lsun':  (Lsun,          'kg m2 s-3'),
# PRESSURE
           'Pa':    (  1e+00,       'kg m-1 s-2'),
           'bar':   (  1e+05,       'kg m-1 s-2'),
           'at':    (  98066.5,     'kg m-1 s-2'),
           'atm':   ( 101325,       'kg m-1 s-2'),
           'torr':  (    133.322,   'kg m-1 s-2'),
           'psi':   (   6894.,      'kg m-1 s-2'),
# FLUX
           'Jy':      (1e-26,         'kg s-2 cy-1'),
           'vegamag': (VegaMag,       'kg m-1 s-3'),  # in W/m2/m
           'mag':     (VegaMag,       'kg m-1 s-3'),  # in W/m2/m
           'STmag':   (STMag,         'kg m-1 s-3'),  # in W/m2/m
           'ABmag':   (ABMag,         'kg s-2 cy-1'), # in W/m2/Hz
           }
            
#-- scaling factors for prefixes            
_scalings ={'y':       1e-24, # yocto
            'z':       1e-21, # zepto
            'a':       1e-18, # atto
            'f':       1e-15, # femto
            'p':       1e-12, # pico
            'n':       1e-09, # nano
            'mu':      1e-06, # micro
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
_aliases = [('micron','mum'),
            ('micro','mu'),
            ('milli','m'),
            ('kilo','k'),
            ('mega','M'),
            ('giga','G'),
            ('nano','n'),
            ('watt','W'),
            ('Watt','W'),
            ('Hz','hz'),
            ('joule','J'),
            ('Joule','J'),
            ('jansky','Jy'),
            ('Jansky','Jy'),
            ('arcsec','as'),
            ('arcmin','am'),
            ('cycles','cy'),
            ('cycle','cy'),
            ('cyc','cy'),
            ('angstrom','A'),
            ('Angstrom','A'),
            ('inch','in'),
            ('^',''),
            ('**',''),
            ('galactic','gal'),
            ('equatorial','equ'),
            ('ecliptic','ecl'),
            ]
 
#-- Change-of-base function definitions
_switch = {'_to_s-1':distance2velocity, # switch from wavelength to velocity
           's-1_to_':velocity2distance, # switch from wavelength to velocity
           'm1_to_cy1s-1':distance2frequency,  # switch from wavelength to frequency
           'cy1s-1_to_m1':distance2frequency,  # switch from frequency to wavelength
           'm1_to_':distance2spatialfreq, # for interferometry
           '_to_m1':spatialfreq2distance, # for interferometry
           'cy-1s-2_to_m-1s-3':fnu2flambda,
           'm-1s-3_to_cy-1s-2':flambda2fnu,
           'cy-1s-2_to_s-3':fnu2nufnu,
           's-3_to_cy-1s-2':nufnu2fnu,
           '_to_cy-2':per_sr,
           'cy-2_to_':times_sr}
 
 
if __name__=="__main__":
    import doctest
    doctest.testmod()
    #-- timing
    #import time
    #c0 = time.time()
    #nr = 50000
    #for i in xrange(nr):
        #x = convert('ABmag','erg s-1 cm-2 micron-1 sr-1',0.,ang_diam=(2.,'mas'),wave=(1.,'micron'))
    #print "One conversion takes on average",(time.time()-c0)/nr*1000,'ms'
    #print convert('cy/d2','s/yr',1e-9)
    #print convert('cy/d2','Hz/yr',1e-9)
    #print(convert('ABmag','Jy',0.,photband='SDSS.U'))
    #3630.7805477
    #print(convert('Jy','erg cm-2 s-1 A-1',3630.7805477,wave=(1.,'micron')))
    #1.08848062485e-09