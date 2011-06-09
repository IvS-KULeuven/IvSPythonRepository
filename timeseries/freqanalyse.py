# -*- coding: utf-8 -*-
"""
Frequency Analysis Routines.

Here are some examples:

Section 1. Pulsation frequency analysis
=======================================

Do a frequency analysis of the star HD129929, after Aerts 2004. We reproduce
their Fig. 8 here.

Import necessary modules:

>>> from ivs.catalogs import vizier

Read in the data and do an iterative prewhitening frequency analysis:

>>> data,units,comms = vizier.search('J/A+A/415/241/table1')
>>> params = iterative_prewhitening(data.HJD,data.Umag-data.Umag.mean(),f0=6.2,fn=7.2,maxiter=6)
>>> print pl.mlab.rec2txt(params,precision=6)
      const       ampl       freq       phase    e_const     e_ampl     e_freq    e_phase
   0.000317   0.014698   6.461699    0.315793   0.000401   0.001331   0.000006   0.090586
   0.000000   0.014891   6.978305    0.346205   0.000000   0.001223   0.000006   0.082163
   0.000000   0.011687   6.449590   -0.361744   0.000000   0.001087   0.000007   0.092995
   0.000000   0.011568   6.990430   -0.234580   0.000000   0.001000   0.000006   0.086491
   0.000000   0.009386   6.590939   -0.456697   0.000000   0.000941   0.000007   0.100237
   0.000000   0.007637   6.966174   -0.172689   0.000000   0.000868   0.000008   0.113673

Plot the results:

>>> p = pl.figure()
>>> p = pl.vlines(params['freq'],0,params['ampl'],color='k',linestyle='-')
>>> p = pl.xlabel('Frequency [c d$^{-1}$]')
>>> p = pl.ylabel('Amplitude [magnitude]')
>>> p = pl.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
>>> p = pl.annotate('$\ell=0$',(6.59,0.010),xycoords='data',va='bottom',ha='center',fontsize='large')
>>> p = pl.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
>>> p = pl.xlim(6.2,7.2)
>>> p = pl.ylim(0,0.018)

Author: Pieter Degroote
"""
import logging
import numpy as np
import pylab as pl
from ivs.timeseries import fit
from ivs.timeseries import evaluate
from ivs.timeseries import pergrams
from ivs.timeseries.decorators import defaults_pergram
from ivs.misc import numpy_ext

logger = logging.getLogger("TS.FREQANAL")

#{ Convenience functions

@defaults_pergram
def find_frequency(times,signal,
            f0=None,fn=None,df=None,optimize=0,
            max_loops=20, scale_region=0.1, scale_df=0.20,
            **kwargs):
    """
    Find one frequency, automatically going to maximum precision and return
    parameters & error estimates.
    
    This routine will make the frequency grid finer until it is finer than the
    estimated error on the frequency. After that, it will compute harmonic
    parameters and estimate errors.
    
    There is a possibility to escape this optimization by setting dfscale=0 or
    freqregscale=0.
    
    You can include a nonlinear least square update of the parameters, by
    setting the keyword C{optimize} to 1 (optimization outside loop) or
    2 (optimization after each iteration).
    
    Possible Extra keywords: see definition of the used periodogram function,
    evaluating and fitting functions, etc...
    
    Example keywords:
        - 'correlation_correction', default=True
        - 'freqregscale', default=0.5: factor for zooming in on frequency
        - 'dfscale', default = 0.25: factor for optimizing frequency resolution
    
    There is a possibility to enter frequency search range in units of the
    Nyquist frequency, by setting the keyword 'units' to 'relative'.
    
    Example usage:
    
    @rtype: record array
    @return: parameters and errors
    """
    
    #-- initial values
    e_f = 0
    freq_diff = np.inf
    prev_freq = -np.inf
    counter = 0
    
    f_max = fn + 0.
    f_min = f0 + 0.
    
    #-- calculate periodogram until frequency precision is
    #   under 1/10th of correlation corrected version
    while freq_diff>e_f/10.:
        #-- calculate periodogram
        freqs,ampls = pergrams.scargle(times,signal,f0=f0,fn=fn,df=df,**kwargs)
        frequency = freqs[np.argmax(ampls)]
        
        #-- estimate parameters and calculate a fit, errors and residuals
        params = fit.sine(times,signal,frequency)
        errors = fit.e_sine(times,signal,params)
        
        #-- optimize inside loop if necessary and if we gained prediction
        #   value:
        if optimize==2:
            params_,errors_,gain = fit.optimize(times,signal,params,'sine')
            if gain>0:
                params = params_
                logger.info('Accepted optimization (gained %g%%)'%gain)
        
        #-- improve precision
        freq_diff = abs(frequency-prev_freq)
        prev_freq = frequency
        freq_region = fn-f0
        f0 = max(f_min,frequency-freq_region*scale_region/2.)
        fn = min(f_max,frequency+freq_region*scale_region/2.)
        df *= scale_df
        
        #-- possibilities to escape iterative zoom in
        if scale_region==0 or scale_df==0: break
        if counter >= max_loops:
            logger.error("Frequency precision not reached in %d steps, breaking loop"%(max_loops))
            break
        counter += 1
        
    #-- optimize parameters outside of loop if necessary:
    if optimize==1:
        params_,errors_,gain = fit.optimize(times,signal,params,'sine')
        if gain>0:
            params = params_
            logger.info('Accepted optimization (gained %g%%)'%gain)
    
    params = numpy_ext.recarr_join(params,errors)
    
    return params



def iterative_prewhitening(times,signal,maxiter=1000,optimize=0,**kwargs):
    """
    Fit one or more functions to a timeseries via iterative prewhitening.
    
    This function will use C{find_frequency} to fit the function parameters to
    the original signal (including any previously found parameters),
    optimize the parameters if needed, and remove it from the data to search
    for a new frequency in the residuals.
    
    It is always the original signal that is used to fit all parameters again;
    only the (optimized) frequency is remembered from step to step.
    """
    
    residuals = signal.copy()
    frequencies = []
    
    while maxiter:
        #-- compute the next frequency from the residuals
        params = find_frequency(times,residuals,**kwargs)
        
        #-- do the fit including all frequencies
        frequencies.append(params['freq'][-1])
        allparams = fit.sine(times,signal,frequencies)
        
        #-- if there's a need to optimize, optimize the last n parameters
        if optimize>0:
            residuals_for_optimization = residuals
            if optimize<=len(params):
                model_fixed_params = evaluate.sine(times,allparams[:-optimize])
                residuals_for_optimization -= model_fixed_params
            uparams,e_uparams, gain = fit.optimize(times,residuals_for_optimization,allparams[-optimize:],'sine')
            #-- only accept the optimization if we gained prediction power
            if gain>0:
                allparams[-optimize:] = uparams
                logger.info('Accepted optimization (gained %g%%)'%gain)
        
        #-- compute the residuals to use in the next prewhitening step
        modelfunc = evaluate.sine(times,allparams)
        residuals = signal - modelfunc
        
        #-- exhaust the counter
        maxiter -= 1
        
    #-- calculate the errors
    e_allparams = fit.e_sine(times,signal,allparams)
    
    allparams = numpy_ext.recarr_join(allparams,e_allparams)
    
    return allparams
    
    


if __name__=="__main__":
    import doctest
    import pylab as pl
    import sys
    doctest.testmod()
    pl.show()
    sys.exit()












  


    