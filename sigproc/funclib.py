"""
Database with model functions
"""
import numpy as np
from ivs.sigproc.fit import Model, Function
import ivs.timeseries.keplerorbit as kepler
#from ivs.sigproc import evaluate

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
#}