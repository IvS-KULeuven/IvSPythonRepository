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

    lmfit.parameter: 
        added kick method to Parameter and Parameters class
    
    lmfit.confidence:
        change the way conf_interval returns the CIs to more usefull structure
        actual change is in ConfidenceInterval.calc_all_ci()
        ci = {pname = {'0.95'=(lower, upper), '0.66'=(lower, upper)} }
    
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
from .confidence import conf_interval, conf_interval2d
from .printfuncs import (fit_report, ci_report,
                         report_fit, report_ci, report_errors)

from . import uncertainties
from .uncertainties import ufloat, correlated_values

__xall__ = ['minimize', 'Minimizer', 'Parameter', 'Parameters',
           'conf_interval', 'conf_interval2d', 'make_paras_and_func',
           'fit_report', 'ci_report', 'report_errors',
           'report_fit', 'report_ci', 'ufloat',
           'correlated_values']
