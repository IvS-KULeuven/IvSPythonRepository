"""
    LMfit-py provides a Least-Squares Minimization routine and
    class with a simple, flexible approach to parameterizing a
    model for fitting to data.  Named Parameters can be held
    fixed or freely adjusted in the fit, or held between lower
    and upper bounds.  If the separate asteval module has been
    installed, parameters can be constrained as a simple
    mathematical expression of other Parameters.

    version: 0.7.2
    last update: 20-June-2013
    License: BSD
    Author:  Matthew Newville <newville@cars.uchicago.edu>
            Center for Advanced Radiation Sources,
            The University of Chicago

    Changes applied to lmfit and uncertainties to make it work with sigproc:

    fixed uncertainties import to accomodate its place in the ivs repository
    uncertainties.umath
        import uncertainties -> from ivs.sigproc.lmfit import uncertainties

    uncertainties.unumpy.__init__:
        from uncertainties.unumpy import core -> import core
        uncertainties.unumpy -> ivs.sigproc.lmfit.uncertainties.unumpy

    uncertainties.unumpy.core:
        import uncertainties -> from ivs.sigproc.lmfit import uncertainties

    uncertainties.unumpy.ulinalg
        from uncertainties import __author__ -> from ivs.sigproc.lmfit.uncertainties import __author__
        from uncertainties.unumpy import core -> import core

    Delete all tests as they are not nessesary at all.
    uncertainties.unumpy.test_unumpy, uncertainties.unumpy.test_ulinalg,
    uncertainties.test_umath, uncertainties.test_uncertainties

"""
__version__ = '0.7.2'
from .minimizer import minimize, Minimizer, MinimizerException, make_paras_and_func
from .parameter import Parameter, Parameters
from .confidence import conf_interval, conf_interval2d, ConfidenceInterval
from .printfuncs import (fit_report, ci_report,
                         report_fit, report_ci, report_errors)

import uncertainties
from uncertainties import ufloat, correlated_values


#======================================================
#{ Add nessessary methods to several classes from lmfit

import re
import numpy as np
from ivs.aux import decorators

@decorators.extend(Parameters)
def can_kick(self,pnames=None):
    """
    Checks if the given parameters can be kicked and returns the good ones in a list.
    """
    if pnames == None:
        pnames = list(self.keys())
    kick_pars = []
    for key in pnames:
        if self[key].can_kick():
            kick_pars.append(key)
    return kick_pars

@decorators.extend(Parameters)
def kick(self,pnames=None):
    """
    Kicks the given parameters to a new value chosen from the uniform
    distribution between max and min value.
    """
    if pnames == None:
        pnames = self.can_kick()
    for key in pnames:
        self[key].kick()

@decorators.extend(Parameters)
def __getattr__(self, name):
    """
    Allow easy access to the parameter attributes like value ect.
    fx. all values can be reached with Parameters.value
    """
    if name in ['name', 'value', 'min', 'max', 'vary', 'expr',
                'stderr', 'mcerr', 'cierr', 'correl']:
        return [getattr(p, name) for n, p in list(self.items())]
    else:
        raise AttributeError

@decorators.extend(Parameter)
def __init__(self, name=None, value=None, vary=True,
                 min=None, max=None, expr=None):
    self.name = name
    self._val = value
    self.user_value = value
    self.init_value = value
    self.min = min
    self.max = max
    self.vary = vary
    self.expr = expr
    self.deps   = None
    self.stderr = None
    self.correl = {}
    self.cierr = {}
    self.mcerr = None
    if self.max is not None and value > self.max:
        self._val = self.max
    if self.min is not None and value < self.min:
        self._val = self.min
    self.from_internal = lambda val: val

@decorators.extend(Parameter)
def can_kick(self):
    """
    Returns True is the parameter can be kicked (has defined min and max)
    """
    return self.min != None and self.max != None and \
           np.isfinite(self.min) and np.isfinite(self.max)

@decorators.extend(Parameter)
def kick(self):
    """
    Kick the starting value to a random value between min and max
    """
    value = np.random.uniform(low=self.min, high=self.max)
    self._val = value
    self.user_value = value

@decorators.extend(Parameter)
def __getattr__(self, name):
    """
    Allow to reach any cierr directly with coded attribute name (fx. ci95),
    Allow to get min and max with the bounds attribute
    Allow to obtain procentual error
    """
    ciexpr = r"^ci(\d\d+?)$"
    errexpr = r"(stderr|mcerr)pc"
    if name == 'bounds':
        return (self.min, self.max)
    if re.match(ciexpr, name):
        ci = '0.' + str(re.findall(ciexpr, name)[0])
        if ci in self.cierr:
            return self.cierr[ci]
        else:
            return (np.nan, np.nan)
    if re.match(errexpr, name):
        error = getattr( self, str(re.findall(errexpr, name)[0]) )
        if error == None or self.value == 0. or np.isnan(self.value):
            return np.nan
        else:
            return abs(error/self.value * 100.)
    else:
        raise AttributeError('%s not assigned'%(name))

@decorators.extend(Minimizer)
def start_minimize(self, engine, **kwargs):
    "Start the actual fitting"
    engine = engine.lower()

    # scalar minimize methods:
    scal_min = dict(powel='Powell', cg='CG', newton='Newton-CG',
                    cobyla='COBYLA', slsqp='SLSQP')
    if engine in ['powell', 'cg', 'newton', 'cobyla', 'slsqp']:
        engine = scal_min[engine]
        self.scalar_minimize(method=engine, **kwargs)

    # other methods
    elif engine == 'anneal':
        self.anneal(**kwargs)
    elif engine == 'lbfgsb':
        self.lbfgsb(**kwargs)
    elif engine == 'fmin' or engine == 'nelder':
        self.fmin(**kwargs)
    else:
        self.leastsq(**kwargs)

@decorators.extend(ConfidenceInterval)
def calc_all_ci(self):
    """
    Calculates all cis, and returns the cis in a structured dict as:
    ci = {pname = {'0.95'=(lower, upper), '0.66'=(lower, upper)} }
    """
    out = {}
    for p in self.p_names:

        lower = self.calc_ci(p, -1)
        upper = self.calc_ci(p, 1)

        o = {}
        for s, l, u in zip(self.sigmas, lower, upper):
            o[s] = (l[1], u[1])
        out[p] = o

    if self.trace:
        self.trace_dict = map_trace_to_names(self.trace_dict,
                                                self.minimizer.params)

    return out

#}
#======================================================

__xall__ = ['minimize', 'Minimizer', 'Parameter', 'Parameters',
           'conf_interval', 'conf_interval2d', 'make_paras_and_func',
           'fit_report', 'ci_report', 'report_errors',
           'report_fit', 'report_ci', 'ufloat',
           'correlated_values']
