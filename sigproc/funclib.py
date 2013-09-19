"""
Database with model functions.

To be used with the L{ivs.sigproc.fit.minimizer} function or with the L{evaluate}
function in this module.

>>> p = plt.figure()
>>> x = np.linspace(-10,10,1000)
>>> p = plt.plot(x,evaluate('gauss',x,[5,1.,2.,0.5]),label='gauss')
>>> p = plt.plot(x,evaluate('voigt',x,[20.,1.,1.5,3.,0.5]),label='voigt')
>>> p = plt.plot(x,evaluate('lorentz',x,[5,1.,2.,0.5]),label='lorentz')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib01.png]

>>> p = plt.figure()
>>> x = np.linspace(0,10,1000)[1:]
>>> p = plt.plot(x,evaluate('power_law',x,[2.,3.,1.5,0,0.5]),label='power_law')
>>> p = plt.plot(x,evaluate('power_law',x,[2.,3.,1.5,0,0.5])+evaluate('gauss',x,[1.,5.,0.5,0,0]),label='power_law + gauss')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib02.png]

>>> p = plt.figure()
>>> x = np.linspace(0,10,1000)
>>> p = plt.plot(x,evaluate('sine',x,[1.,2.,0,0]),label='sine')
>>> p = plt.plot(x,evaluate('sine_linfreqshift',x,[1.,0.5,0,0,.5]),label='sine_linfreqshift')
>>> p = plt.plot(x,evaluate('sine_expfreqshift',x,[1.,0.5,0,0,1.2]),label='sine_expfreqshift')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib03.png]

>>> p = plt.figure()
>>> p = plt.plot(x,evaluate('sine',x,[1.,2.,0,0]),label='sine')
>>> p = plt.plot(x,evaluate('sine_orbit',x,[1.,2.,0,0,0.1,10.,0.1]),label='sine_orbit')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib03a.png]

>>> p = plt.figure()
>>> x_single = np.linspace(0,10,1000)
>>> x_double = np.vstack([x_single,x_single])
>>> p = plt.plot(x_single,evaluate('kepler_orbit',x_single,[2.5,0.,0.5,0,3,1.]),label='kepler_orbit (single)')
>>> y_double = evaluate('kepler_orbit',x_double,[2.5,0.,0.5,0,3,2.,-4,2.],type='double')
>>> p = plt.plot(x_double[0],y_double[0],label='kepler_orbit (double 1)')
>>> p = plt.plot(x_double[1],y_double[1],label='kepler_orbit (double 2)')
>>> p = plt.plot(x,evaluate('box_transit',x,[2.,0.4,0.1,0.3,0.5]),label='box_transit')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib04.png]

>>> p = plt.figure()
>>> x = np.linspace(-1,1,1000)
>>> gammas = [-0.25,0.1,0.25,0.5,1,2,4]
>>> y = np.array([evaluate('soft_parabola',x,[1.,0,1.,gamma]) for gamma in gammas])
divide by zero encountered in power
>>> for iy,gamma in zip(y,gammas): p = plt.plot(x,iy,label="soft_parabola $\gamma$={:.2f}".format(gamma))
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib05.png]

>>> p = plt.figure()
>>> x = np.logspace(-1,2,1000)
>>> blbo = evaluate('blackbody',x,[10000.,1.],wave_units='micron',flux_units='W/m3')
>>> raje = evaluate('rayleigh_jeans',x,[10000.,1.],wave_units='micron',flux_units='W/m3')
>>> wien = evaluate('wien',x,[10000.,1.],wave_units='micron',flux_units='W/m3')
>>> p = plt.subplot(221)
>>> p = plt.title(r'$\lambda$ vs $F_\lambda$')
>>> p = plt.loglog(x,blbo,label='Black Body')
>>> p = plt.loglog(x,raje,label='Rayleigh-Jeans')
>>> p = plt.loglog(x,wien,label='Wien')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

>>> blbo = evaluate('blackbody',x,[10000.,1.],wave_units='micron',flux_units='Jy')
>>> raje = evaluate('rayleigh_jeans',x,[10000.,1.],wave_units='micron',flux_units='Jy')
>>> wien = evaluate('wien',x,[10000.,1.],wave_units='micron',flux_units='Jy')
>>> p = plt.subplot(223)
>>> p = plt.title(r"$\lambda$ vs $F_\\nu$")
>>> p = plt.loglog(x,blbo,label='Black Body')
>>> p = plt.loglog(x,raje,label='Rayleigh-Jeans')
>>> p = plt.loglog(x,wien,label='Wien')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

>>> x = np.logspace(0.47,3.47,1000)
>>> blbo = evaluate('blackbody',x,[10000.,1.],wave_units='THz',flux_units='Jy')
>>> raje = evaluate('rayleigh_jeans',x,[10000.,1.],wave_units='THz',flux_units='Jy')
>>> wien = evaluate('wien',x,[10000.,1.],wave_units='THz',flux_units='Jy')
>>> p = plt.subplot(224)
>>> p = plt.title(r"$\\nu$ vs $F_\\nu$")
>>> p = plt.loglog(x,blbo,label='Black Body')
>>> p = plt.loglog(x,raje,label='Rayleigh-Jeans')
>>> p = plt.loglog(x,wien,label='Wien')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

>>> blbo = evaluate('blackbody',x,[10000.,1.],wave_units='THz',flux_units='W/m3')
>>> raje = evaluate('rayleigh_jeans',x,[10000.,1.],wave_units='THz',flux_units='W/m3')
>>> wien = evaluate('wien',x,[10000.,1.],wave_units='THz',flux_units='W/m3')
>>> p = plt.subplot(222)
>>> p = plt.title(r"$\\nu$ vs $F_\lambda$")
>>> p = plt.loglog(x,blbo,label='Black Body')
>>> p = plt.loglog(x,raje,label='Rayleigh-Jeans')
>>> p = plt.loglog(x,wien,label='Wien')
>>> leg = plt.legend(loc='best')
>>> leg.get_frame().set_alpha(0.5)

]include figure]]ivs_sigproc_fit_funclib06.png]

"""
import numpy as np
from numpy import pi,cos,sin,sqrt,tan,arctan
from scipy.special import erf,jn
from ivs.sigproc.fit import Model, Function
import ivs.timeseries.keplerorbit as kepler
from ivs.sed import model as sed_model
import logging

logger = logging.getLogger("SP.FUNCLIB")

#{ Function Library
def blackbody(wave_units='AA',flux_units='erg/s/cm2/AA',disc_integrated=True):
    """
    Blackbody (T, scale).
    
    @param wave_units: wavelength units
    @type wave_units: string
    @param flux_units: flux units
    @type flux_units: string
    @param disc_integrated: sets units equal to SED models
    @type disc_integrated: bool
    """
    pnames = ['T', 'scale']
    function = lambda p, x: p[1]*sed_model.blackbody(x, p[0], wave_units=wave_units,\
                                                flux_units=flux_units,\
                                                disc_integrated=disc_integrated)
    function.__name__ = 'blackbody'
    
    return Function(function=function, par_names=pnames)
    
def rayleigh_jeans(wave_units='AA',flux_units='erg/s/cm2/AA',disc_integrated=True):
    """
    Rayleigh-Jeans tail (T, scale).
    
    @param wave_units: wavelength units
    @type wave_units: string
    @param flux_units: flux units
    @type flux_units: string
    @param disc_integrated: sets units equal to SED models
    @type disc_integrated: bool
    """
    pnames = ['T', 'scale']
    function = lambda p, x: p[1]*sed_model.rayleigh_jeans(x, p[0], wave_units=wave_units,\
                                                flux_units=flux_units,\
                                                disc_integrated=disc_integrated)
    function.__name__ = 'rayleigh_jeans'
    
    return Function(function=function, par_names=pnames)
    
def wien(wave_units='AA',flux_units='erg/s/cm2/AA',disc_integrated=True):
    """
    Wien approximation (T, scale).
    
    @param wave_units: wavelength units
    @type wave_units: string
    @param flux_units: flux units
    @type flux_units: string
    @param disc_integrated: sets units equal to SED models
    @type disc_integrated: bool
    """
    pnames = ['T', 'scale']
    function = lambda p, x: p[1]*sed_model.wien(x, p[0], wave_units=wave_units,\
                                                flux_units=flux_units,\
                                                disc_integrated=disc_integrated)
    function.__name__ = 'wien'
    
    return Function(function=function, par_names=pnames)

def kepler_orbit(type='single'):
    """
    Kepler orbits ((p,t0,e,omega,K,v0) or (p,t0,e,omega,K1,v01,K2,v02))
    
    A single kepler orbit
    parameters are: [p, t0, e, omega, k, v0]
    
    A double kepler orbit
    parameters are: [p, t0, e, omega, k_1, v0_1, k_2, v0_2]
    Warning: This function uses 2d input and output!
    """
    if type == 'single':
        pnames = ['p','t0','e','omega','k','v0']
        function = lambda p, x: kepler.radial_velocity(p, times=x, itermax=8)
        function.__name__ = 'kepler_orbit_single'
    
        return Function(function=function, par_names=pnames)
        
    elif type == 'double':
        pnames = ['p','t0','e','omega','k1','v01','k2','v02' ]
        function = lambda p, x: [kepler.radial_velocity([p[0],p[1],p[2],p[3],p[4],p[5]], times=x[0], itermax=8),
                                 kepler.radial_velocity([p[0],p[1],p[2],p[3],p[6],p[7]], times=x[1], itermax=8)]
        def residuals(syn, data, weights=None, errors=None, **kwargs):
            return np.hstack( [( data[0] - syn[0] ) * weights[0],  ( data[1] - syn[1] ) * weights[1] ] )
        function.__name__ = 'kepler_orbit_double'
    
        return Function(function=function, par_names=pnames, resfunc=residuals)

def box_transit(t0=0.):
    """
    Box transit model (cont,freq,ingress,egress,depth)
    
    @param t0: reference time (defaults to 0)
    @type t0: float
    """
    pnames = 'cont','freq','ingress','egress','depth'
    def function(p,x):
        cont,freq,ingress,egress,depth = p
        model = np.ones(len(x))*cont
        phase = np.fmod((x - t0) * freq, 1.0)
        phase = np.where(phase<0,phase+1,phase)
        transit_place = (ingress<=phase) & (phase<=egress)
        model = np.where(transit_place,model-depth,model)
        return model
    function.__name__ = 'box_transit'
    return Function(function=function, par_names=pnames)


def polynomial(d=1):
    """
    Polynomial (a1,a0).
    
    y(x) = ai*x**i + a(i-1)*x**(i-1) + ... + a1*x + a0
    
    @param d: degree of the polynomial
    @type d: int
    """
    pnames = ['a{:d}'.format(i) for i in range(d,0,-1)]+['a0']
    function = lambda p, x: np.polyval(p,x)
    function.__name__ = 'polynomial'
    
    return Function(function=function, par_names=pnames)


def soft_parabola():
    """
    Soft parabola (ta,vlsr,vinf,gamma).
    
    See Olofsson 1993ApJS...87...267O.
    
    T_A(x) = T_A(0) * [ 1 - ((x- v_lsr)/v_inf)**2 ] ** (gamma/2)
    """
    pnames = ['ta','vlsr','vinf','gamma']
    def function(p,x):
        term = (x-p[1]) / p[2]
        y = p[0] * (1- term**2)**(p[3]/2.)
        if p[3]<=0: y[np.abs(term)>=1] = 0
        y[np.isnan(y)] = 0
        return y
    function.__name__ = 'soft_parabola'
    
    return Function(function=function, par_names=pnames)


def gauss(use_jacobian=True):
    """
    Gaussian (a,mu,sigma,c)
    
    f(x) = a * exp( - (x - mu)**2 / (2 * sigma**2) ) + c
    """
    pnames = ['a', 'mu', 'sigma', 'c']
    function = lambda p, x: p[0] * np.exp( -(x-p[1])**2 / (2.0*p[2]**2)) + p[3]
    function.__name__ = 'gauss'
    
    if not use_jacobian:
        return Function(function=function, par_names=pnames)
    else:
        def jacobian(p, x):
            ex = np.exp( -(x-p[1])**2 / (2.0*p[2]**2) )
            return np.array([-ex, -p[0] * (x-p[1]) * ex / p[2]**2, -p[0] * (x-p[1])**2 * ex / p[2]**3, [-1 for i in x] ]).T
        return Function(function=function, par_names=pnames, jacobian=jacobian)
    
def sine():
    """
    Sine (ampl,freq,phase,const)
    
    f(x) = ampl * sin(2pi*freq*x + 2pi*phase) + const
    """
    pnames = ['ampl', 'freq', 'phase', 'const']
    function = lambda p, x: p[0] * sin(2*pi*(p[1]*x + p[2])) + p[3]
    function.__name__ = 'sine'
    
    return Function(function=function, par_names=pnames)


def sine_linfreqshift(t0=0.):
    """
    Sine with linear frequency shift (ampl,freq,phase,const,D).
        
    Similar to C{sine}, but with extra parameter 'D', which is the linear
    frequency shift parameter.
    
    @param t0: reference time (defaults to 0)
    @type t0: float
    """
    pnames = ['ampl', 'freq', 'phase', 'const','D']
    def function(p,x):
        freq = (p[1] + p[4]/2.*(x-t0))*(x-t0)
        return p[0] * sin(2*pi*(freq + p[2])) + p[3]
    function.__name__ = 'sine_linfreqshift'
    return Function(function=function, par_names=pnames)
    
def sine_expfreqshift(t0=0.):
    """
    Sine with exponential frequency shift (ampl,freq,phase,const,K).
        
    Similar to C{sine}, but with extra parameter 'K', which is the exponential
    frequency shift parameter.
    
    frequency(x) = freq / log(K) * (K**(x-t0)-1)
    
    f(x) = ampl * sin( 2*pi * (frequency + phase))
    
    @param t0: reference time (defaults to 0)
    @type t0: float
    """
    pnames = ['ampl', 'freq', 'phase', 'const','K']
    def function(p,x):
        freq = p[1] / np.log(p[4]) * (p[4]**(x-t0)-1.)
        return p[0] * sin(2*pi*(freq + p[2])) + p[3]
    function.__name__ = 'sine_expfreqshift'
    return Function(function=function, par_names=pnames)

    
def sine_orbit(t0=0.,nmax=10):
    """
    Sine with a sinusoidal frequency shift (ampl,freq,phase,const,forb,asini,omega,(,ecc))
    
    Similar to C{sine}, but with extra parameter 'asini' and 'forb', which are
    the orbital parameters. forb in cycles/day or something similar, asini in au.
    
    For eccentric orbits, add longitude of periastron 'omega' (radians) and
    'ecc' (eccentricity).
        
    @param t0: reference time (defaults to 0)
    @type t0: float
    @param nmax: number of terms to include in series for eccentric orbit
    @type nmax: int
    """
    pnames = ['ampl', 'freq', 'phase', 'const','forb','asini','omega']
    def ane(n,e): return 2.*sqrt(1-e**2)/e/n*jn(n,n*e)    
    def bne(n,e): return 1./n*(jn(n-1,n*e)-jn(n+1,n*e))
    def function(p,x):
        ampl,freq,phase,const,forb,asini,omega = p[:7]
        ecc = None
        if len(p)==8:
            ecc = p[7]
        cc = 173.144632674 # speed of light in AU/d
        alpha = freq*asini/cc
        if ecc is None:
            frequency = freq*(x-t0) + alpha*(sin(2*pi*forb*x) - sin(2*pi*forb*t0))
        else:
            ns = np.arange(1,nmax+1,1)
            ans,bns = np.array([[ane(n,ecc),bne(n,ecc)] for n in ns]).T
            ksins = sqrt(ans**2*cos(omega)**2+bns**2*sin(omega)**2)
            thns = arctan(bns/ans*tan(omega))
            tau = -np.sum(bns*sin(omega))
            frequency = freq*(x-t0) + \
               alpha*(np.sum(np.array([ksins[i]*sin(2*pi*ns[i]*forb*(x-t0)+thns[i]) for i in range(nmax)]),axis=0)+tau)
        return ampl * sin(2*pi*(frequency + phase)) + const
    function.__name__ = 'sine_orbit'
    return Function(function=function, par_names=pnames)

#def generic(func_name):
    ##func = model.FUNCTION(function=getattr(evaluate,func_name), par_names)
    #raise NotImplementedError

def power_law():
    """
    Power law (A,B,C,f0,const)
    
    P(f) = A / (1+ B(f-f0))**C + const
    """
    pnames = ['ampl','b','c','f0','const']
    function = lambda p, x: p[0] / (1 + (p[1]*(x-p[3]))**p[2]) + p[4]
    function.__name__ = 'power_law'
    return Function(function=function, par_names=pnames)


def lorentz():
    """
    Lorentz profile (ampl,mu,gamma,const)
    
    P(f) = A / ( (x-mu)**2 + gamma**2) + const
    """
    pnames = ['ampl','mu','gamma','const']
    function = lambda p,x: p[0] / ((x-p[1])**2 + p[2]**2) + p[3]
    function.__name__ = 'lorentz'
    return Function(function=function, par_names=pnames)

def voigt():
    """
    Voigt profile (ampl,mu,sigma,gamma,const)
    
    z = (x + gamma*i) / (sigma*sqrt(2))
    V = A * Real[cerf(z)] / (sigma*sqrt(2*pi))
    """
    pnames = ['ampl','mu','sigma','gamma','const']
    def function(p,x):
        x = x-p[1]
        z = (x+1j*p[3])/(p[2]*sqrt(2))
        return p[0]*_complex_error_function(z).real/(p[2]*sqrt(2*pi))+p[4]
    function.__name__ = 'voigt'
    return Function(function=function, par_names=pnames)

#}

#{ Combination functions

def multi_sine(n=10):
    """
    Multiple sines.
    
    @param n: number of sines
    @type n: int
    """
    return Model(functions=[sine() for i in range(n)])

def multi_blackbody(n=3,**kwargs):
    """
    Multiple black bodies.
    
    @param n: number of blackbodies
    @type n: int
    """
    return Model(functions=[blackbody(**kwargs) for i in range(n)])

#}

#{ Internal Helper functions

def _complex_error_function(x):
    """
    Complex error function
    """
    cef_value = np.exp(-x**2)*(1-erf(-1j*x))
    if sum(np.isnan(cef_value))>0:
        logging.warning("Complex Error function: NAN encountered, results are biased")
        noisnans = np.compress(1-np.isnan(cef_value),cef_value)
        try:
            last_value = noisnans[-1]
        except:
            last_value = 0
            logging.warning("Complex Error function: all values are NAN, results are wrong")
        cef_value = np.where(np.isnan(cef_value),last_value,cef_value)

    return cef_value

#}


def evaluate(funcname, domain, parameters, **kwargs):
    """
    Evaluate a function on specified interval with specified parameters.
    
    Extra keywords are passed to the funcname.
    
    Example:
    
    >>> x = np.linspace(-5,5,1000)
    >>> y = evaluate('gauss',x,[1.,0.,2.,0.])
    
    @parameter funcname: name of the function
    @type funcname: str
    @parameter domain: domain to evaluate onto
    @type domain: array
    @parameter parameters: parameter of the function
    @type parameters: array
    """
    function = globals()[funcname](**kwargs)
    function.setup_parameters(parameters)
    return function.evaluate(domain)
    
if  __name__=="__main__":
    import doctest
    from matplotlib import pyplot as plt
    doctest.testmod()
    plt.show()