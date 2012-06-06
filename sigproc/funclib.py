"""
Database with model functions.

To be used with the Lsigproc.fit.minimizer function or with the L{evaluate}
function in this module.

>>> x = np.linspace(-10,10,1000)
>>> p = plt.plot(x,evaluate('gauss',x,[5,1.,2.,0.5]),label='gauss')
>>> p = plt.plot(x,evaluate('voigt',x,[20.,1.,1.5,3.,0.5]),label='voigt')
>>> p = plt.plot(x,evaluate('lorentz',x,[5,1.,2.,0.5]),label='lorentz')
"""
import numpy as np
from numpy import pi,cos,sin,sqrt,tan,arctan
from scipy.special import erf,jn
from ivs.sigproc.fit import Model, Function
import ivs.timeseries.keplerorbit as kepler

#{ Function Library

def kepler_orbit(type='single'):
    """
    A single kepler orbit
    parameters are: [P, T0, e, omega, K, v0]
    
    A double kepler orbit
    parameters are: [P, T0, e, omega, K_1, v0_1, K_2, v0_2]
    Warning: This function uses 2d input and output!
    """
    if type == 'single':
        function = lambda p, x: kepler.radial_velocity(p, times=x, itermax=8)
        pnames = ['p','t0','e','omega','k','v0']
    
        return Function(function=function, par_names=pnames)
        
    elif type == 'double':
        function = lambda p, x: np.vstack([kepler.radial_velocity([p[0],p[1],p[2],p[3],p[4],p[5]], times=x[0], itermax=8),
                                 kepler.radial_velocity([p[0],p[1],p[2],p[3],p[6],p[7]], times=x[1], itermax=8)])
        pnames = ['p','t0','e','omega','k1','v01','k2','v02' ]
    
        return Function(function=function, par_names=pnames)

def gauss():
    """
    Your standard gaussian
    parameters are: [A, mu, sigma, cte]
    f(x) = Amp * exp( - (x - mu)**2 / (2 * sigma**2) ) + cte
    """
    pnames = ['a', 'mu', 'sigma', 'c']
    function = lambda p, x: p[0] * np.exp( -(x-p[1])**2 / (2.0*p[2]**2)) + p[3]
    
    return Function(function=function, par_names=pnames)
    
def sine():
    """
    Your standard sine function
    parameters are: [Amp, omega, phi, cte]
    f(x) = Amp * sin(omega*x + phi) + cte
    """
    pnames = ['a', 'omega', 'phi', 'c']
    function = lambda p, x: p[0] * np.sin(p[1]*x + p[2]) + p[3]
    
    return Function(function=function, par_names=pnames)
    
#def generic(func_name):
    ##func = model.FUNCTION(function=getattr(evaluate,func_name), par_names)
    #raise NotImplementedError

def power_law():
    """
    Power law
    
    P(f) = A / (1+ Bf)**C + const
    """
    pnames = ['ampl','b','c','const']
    function = lambda p, x: p[0] / (1 + (par[1]*x)**par[2]) + par[3]
    return Function(function=function, par_names=pnames)


def lorentz():
    """
    Lorentz profile
    
    P(f) = A / ( (x-mu)**2 + gamma**2)
    """
    pnames = ['ampl','mu','gamma','const']
    function = lambda p,x: p[0] / ((x-p[1])**2 + p[2]**2) + p[3]
    return Function(function=function, par_names=pnames)

def voigt():
    """
    Voigt profile
    
    z = (x + gamma*i) / (sigma*sqrt(2))
    V = A * Real[cerf(z)] / (sigma*sqrt(2*pi))
    """
    pnames = ['ampl','mu','sigma','gamma','const']
    def function(p,x):
        x = x-p[1]
        z = (x+1j*p[3])/(p[2]*sqrt(2))
        return p[0]*_complex_error_function(z).real/(p[2]*sqrt(2*pi))+p[4]
    return Function(function=function, par_names=pnames)

#}

#{ Internal Helper functions

def _complex_error_function(x):
    """
    Complex error function
    """
    cef_value = np.exp(-x**2)*(1-erf(-1j*x))
    if sum(np.isnan(cef_value))>0:
        logger.warning("Complex Error function: NAN encountered, results are biased")
        noisnans = np.compress(1-np.isnan(cef_value),cef_value)
        try:
            last_value = noisnans[-1]
        except:
            last_value = 0
            logger.warning("Complex Error function: all values are NAN, results are wrong")
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