"""
Database with model functions to be used with the sigproc.fit.minimizer function. 
"""
import numpy as np
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
    Power law.
    
    Function:
    
    P(f) = A / (1+ Bf)**C + const
    """
    pnames = ['a','b','c','const']
    function = lambda p, x: p[0] / (1 + (par[1]*x)**par[2]) + par[3]
    return Function(function=function, par_names=pnames)



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
    
    