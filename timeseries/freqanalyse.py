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
>>> params = iterative_prewhitening(data.HJD,data.Umag-data.Umag.mean(),f0=6.2,fn=7.2,maxiter=6,\
          stopcrit=(stopcrit_scargle_snr,4.,6,))
>>> print pl.mlab.rec2txt(params,precision=6)
      const       ampl       freq       phase    e_const     e_ampl     e_freq    e_phase   stopcrit
   0.000310   0.014722   6.461700    0.443450   0.000401   0.001323   0.000006   0.089865   4.950496
   0.000000   0.014866   6.978306   -0.050189   0.000000   0.001224   0.000006   0.082305   5.677022
   0.000000   0.011687   6.449591    0.016647   0.000000   0.001090   0.000007   0.093280   5.747565
   0.000000   0.011546   6.990431   -0.482814   0.000000   0.001001   0.000006   0.086678   6.454955
   0.000000   0.009380   6.590938   -0.382048   0.000000   0.000938   0.000007   0.100016   6.510729
   0.000000   0.007645   6.966174    0.323627   0.000000   0.000863   0.000008   0.112876   5.703801

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


Section 2: Line profile variability
===================================

We generate some line profiles as Gaussian profiles, with varying average. This
as an analogue for purely radial variability.

>>> wave = np.linspace(4950.,5000.,200.)
>>> times = np.linspace(0,1.,100)
>>> A,sigma,mu = 0.5, 5., 4975.
>>> frequency = 1/0.2
>>> wshifts = 10.*np.sin(2*np.pi*frequency*times)
>>> spectra = np.zeros((len(times),len(wave)))
>>> for i,wshift in enumerate(wshifts):
...    spectra[i] = 1.-evaluate.gauss(wave,[A,mu+wshift,sigma])+np.random.normal(size=len(wave),scale=0.01)

Once the line profiles are constructed, we can compute the Fourier transform
of these lines:

>>> output = spectrum_2D(times,wave,spectra,f0=frequency/2.,fn=3*frequency,df=0.001,threads=2,method='scargle',subs_av=True,full_output=True)

With this output, we can find retrieve the frequency. We plot the periodograms
for each wavelength bin on the left, and on the right the average over all
wavelength bins:

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.imshow(output['pergram'][1][::-1],extent=[wave[0],wave[-1],frequency/2,3*frequency],aspect='auto')
>>> p = pl.xlabel('Wavelength (A)')
>>> p = pl.ylabel('Frequency (c/d)')
>>> cbar = pl.colorbar()
>>> cbar.set_label('Amplitude')
>>> p = pl.subplot(122)
>>> p = pl.plot(output['pergram'][1].mean(axis=1),output['pergram'][0],'k-')
>>> p = pl.ylim(frequency/2,3*frequency)
>>> p = pl.xlabel('Amplitude')
>>> p = pl.ylabel('Frequency (c/d)')

]]include figure]]timeseries_freqanalyse_02.png]

We can then fix the frequency and compute all the parameters. However, we now
choose not to subtract the average profile, but instead use the GLS periodogram
to include the constant

>>> output = spectrum_2D(times,wave,spectra,f0=frequency-0.05,fn=frequency+0.05,df=0.001,threads=2,method='gls',subs_av=False,full_output=True)

>>> p = pl.figure()
>>> p = pl.subplot(221)
>>> p = pl.errorbar(wave,output['pars']['const'],yerr=output['pars']['e_const'],fmt='ko-')
>>> p,q = pl.xlabel('Wavelength (A)'),pl.ylabel('Constant')
>>> p = pl.subplot(222)
>>> p = pl.errorbar(wave,output['pars']['ampl'],yerr=output['pars']['e_ampl'],fmt='ro-')
>>> p,q = pl.xlabel('Wavelength (A)'),pl.ylabel('Amplitude')
>>> p = pl.subplot(223)
>>> p = pl.errorbar(wave,output['pars']['freq'],yerr=output['pars']['e_freq'],fmt='ko-')
>>> p,q = pl.xlabel('Wavelength (A)'),pl.ylabel('Frequency (c/d)')
>>> p = pl.subplot(224)
>>> p = pl.errorbar(wave,output['pars']['phase'],yerr=output['pars']['e_phase'],fmt='ro-')
>>> p,q = pl.xlabel('Wavelength (A)'),pl.ylabel('Phase')

]]include figure]]timeseries_freqanalyse_03.png]

Section 3. Time-frequency analysis
==================================

Example usage: we generate some data where the frequency jumps to double
its value in the middle of the time series. The first half is thus given by

>>> params = np.rec.fromarrays([[10.],[1.],[0.],[0.]],names=['freq','ampl','const','phase'])
>>> times = np.linspace(0,15,1000)
>>> signal = evaluate.sine(times,params)

And the second half by

>>> params['freq'] = 20.
>>> signal[:len(times)/2] = evaluate.sine(times[:len(times)/2],params)

The time-frequency analysis is calculate via the command

>>> output = time_frequency(times,signal,fn=30.)

And the outcome can be plotted via

>>> p = pl.figure()
>>> p = pl.imshow(output['pergram'][1].T[::-1],aspect='auto',extent=[times[0],times[-1],output['pergram'][0][0],output['pergram'][0][-1]])
>>> p,q = pl.xlabel('Time (d)'),pl.ylabel('Frequency (c/d)')

]]include figure]]timeseries_freqanalyse_04.png]

or

>>> p = pl.figure()
>>> p = pl.subplot(221)
>>> p = pl.errorbar(output['times'],output['pars']['const'],yerr=output['pars']['e_const'],fmt='ko-')
>>> p,q = pl.xlabel('Time (d)'),pl.ylabel('Constant')
>>> p = pl.subplot(222)
>>> p = pl.errorbar(output['times'],output['pars']['ampl'],yerr=output['pars']['e_ampl'],fmt='ro-')
>>> p,q = pl.xlabel('Time (d)'),pl.ylabel('Amplitude')
>>> p = pl.subplot(223)
>>> p = pl.errorbar(output['times'],output['pars']['freq'],yerr=output['pars']['e_freq'],fmt='ko-')
>>> p,q = pl.xlabel('Time (d)'),pl.ylabel('Frequency (c/d)')
>>> p = pl.subplot(224)
>>> p = pl.errorbar(output['times'],output['pars']['phase'],yerr=output['pars']['e_phase'],fmt='ro-')
>>> p,q = pl.xlabel('Time (d)'),pl.ylabel('Phase')

]]include figure]]timeseries_freqanalyse_05.png]

Author: Pieter Degroote
"""
import logging
import numpy as np
import pylab as pl
from ivs.sigproc import fit
from ivs.sigproc import evaluate
from ivs.timeseries import pergrams
from ivs.timeseries.decorators import defaults_pergram
from ivs.aux import numpy_ext

logger = logging.getLogger("TS.FREQANAL")

#{ Convenience functions

def find_frequency(times,signal,method='scargle',model='sine',full_output=False,
            optimize=0,max_loops=20, scale_region=0.1, scale_df=0.20, model_kwargs=None,
            correlation_correction=True,prewhiteningorder_snr=False,
            prewhiteningorder_snr_window=1.,**kwargs):
    """
    Find one frequency, automatically going to maximum precision and return
    parameters & error estimates.

    This routine makes the frequency grid finer until it is finer than the
    estimated error on the frequency. After that, it will compute (harmonic)
    parameters and estimate errors.

    There is a possibility to escape this optimization by setting scale_df=0 or
    freqregscale=0.

    You can include a nonlinear least square update of the parameters, by
    setting the keyword C{optimize=1} (optimization outside loop) or
    C{optimize=2} (optimization after each iteration).

    The method with which to find the frequency can be set with the keyword
    C{method}, the model used to fit and optimize should be set with C{model}.
    Extra keywords for the model functions should go in C{model_kwargs}. If
    C{method} is a tuple, the first method will be used for the first frequency
    search only. This could be useful to take advantage of such methods as
    fasper which do not allow for iterative zoom-ins. By default, the function looks
    for the highest (or deepest in the case of the pdm method) peak, but instead it is
    possible to go for the peak having the highest SNR before prewhitening by setting
    C{prewhiteningorder_snr} to True. In this case, the noise spectrum is calculated
    using a convolution with a C{prewhiteningorder_snr_window} wide box.

    Possible extra keywords: see definition of the used periodogram function.

    B{Warning}: the timeseries must be B{sorted in time} and B{cannot contain
    the same timepoint twice}. Otherwise, a 'ValueError, concatenation problem'
    can occur.

    Example keywords:
        - 'correlation_correction', default=True
        - 'freqregscale', default=0.5: factor for zooming in on frequency
        - 'dfscale', default = 0.25: factor for optimizing frequency resolution

    Example usage: We generate a sine signal

    >>> times = np.linspace(0,100,1000)
    >>> signal = np.sin(2*np.pi*2.5*times) + np.random.normal(size=len(times))

    Compute the frequency

    >>> parameters, pgram, model = find_frequency(times,signal,full_output=True)

    Make a plot:

    >>> p = pl.figure()
    >>> p = pl.axes([0.1,0.3,0.85,0.65])
    >>> p = pl.plot(pgram[0],pgram[1],'k-')
    >>> p = pl.xlim(2.2,2.8)
    >>> p = pl.ylabel('Amplitude')
    >>> p = pl.axes([0.1,0.1,0.85,0.2])
    >>> p = pl.plot(pgram[0][:-1],np.diff(pgram[0]),'k-')
    >>> p = pl.xlim(2.2,2.8)
    >>> p,q = pl.xlabel('Frequency (c/d)'),pl.ylabel('Frequency resolution $\Delta F$')

    ]]include figure]]timeseries_freqanalyse_06.png]

    @rtype: record array(, 2x1Darray, 1Darray)
    @return: parameters and errors(, periodogram, model function)
    """
    if model_kwargs is None:
        model_kwargs = dict()
    #-- initial values
    e_f = 0
    freq_diff = np.inf
    prev_freq = -np.inf
    counter = 0

    f_max = np.inf
    f_min = 0.#-np.inf

    #-- calculate periodogram until frequency precision is
    #   under 1/10th of correlation corrected version of frequency error
    method_kwargs = kwargs.copy() # don't modify the dictionary the user gave

    while freq_diff>e_f/10.:
        #-- possibly, we might want to use different periodograms for the first
        #   calculation than for the zoom in, since some periodograms are faster
        #   than others but do not have the ability to 'zoom in' (e.g. the FFT)
        if freq_diff==np.inf and not isinstance(method,str):
            method_ = method[1]
            method = method[0]  # override method to be a string the next time
        #-- calculate periodogram
        freqs,ampls = getattr(pergrams,method)(times,signal,**method_kwargs)
        f0,fn,df = freqs[0],freqs[-1],freqs[1]-freqs[0]
        #-- now use the second method for the zoom-ins from now on
        if freq_diff==np.inf and not isinstance(method,str):
            method = method_
        #-- extract the frequency: this part should be generalized, but for now,
        #   it will do:
        if method in ['pdm']:
            frequency = freqs[np.argmin(ampls)]
        #-- instead of going for the highest peak, let's get the most significant one
        if prewhiteningorder_snr:
            if counter == 0: #we calculate a noise spectrum with a convolution in a 1 d-1 window
                windowlength = float(prewhiteningorder_snr_window)/(freqs[1]-freqs[0])
                window = np.ones(int(windowlength))/float(windowlength)
                ampls_ = np.concatenate((ampls[::-1],ampls,ampls[::-1])) #we mirror the amplitude spectrum on both ends so the convolution will be better near the edges
                noises_ = np.convolve(ampls_, window, 'same')
                noises = np.split(noises_,3)[1] #and we recut the resulted convolution to match the original frequency range
                freqs_old = np.copy(freqs)
                noises_old = np.copy(noises)
            else:
                noises = np.interp(freqs,freqs_old,noises_old) #we use the original noise spectrum in this narrower windows too, which should save some time, and avoid the problem of having a wider window for the SNR calculation than the width of the zoom-in window
            frequency = freqs[np.argmax(ampls/noises)]
        else:
            frequency = freqs[np.argmax(ampls)]
        if full_output and counter==0:
            freqs_,ampls_ = freqs,ampls
        #-- estimate parameters and calculate a fit, errors and residuals
        params = getattr(fit,model)(times,signal,frequency,**model_kwargs)
        if hasattr(fit,'e_'+model):
            errors = getattr(fit,'e_'+model)(times,signal,params,correlation_correction=correlation_correction)
            e_f = errors['e_freq'][-1]
        #-- possibly there are not errors defined for this fitting functions
        else:
            errors = None
        #-- optimize inside loop if necessary and if we gained prediction
        #   value:
        if optimize==2:
            params_,errors_,gain = fit.optimize(times,signal,params,model)
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
        method_kwargs['f0'] = f0
        method_kwargs['fn'] = fn
        method_kwargs['df'] = df
        #-- possibilities to escape iterative zoom in
        #print '---> {counter}/{max_loops}: freq={frequency} ({f0}-->{fn}/{df}), e_f={e_f}, freq_diff={freq_diff}'.format(**locals()),max(ampls)
        if scale_region==0 or scale_df==0:
            break
        if counter >= max_loops:
            logger.error("Frequency precision not reached in %d steps, breaking loop"%(max_loops))
            break
        if (fn-f0)/df<5:
            logger.error("Frequency precision not reached with stepsize %e , breaking loop"%(df/scale_df))
            break
        counter += 1
    #-- optimize parameters outside of loop if necessary:
    if optimize==1:
        params_,errors_,gain = fit.optimize(times,signal,params,model)
        if gain>0:
            params = params_
            logger.info('Accepted optimization (gained %g%%)'%gain)
    #-- add the errors to the parameter array if possible
    if errors is not None:
        params = numpy_ext.recarr_join(params,errors)
    # logger.info("%s model parameters via %s periodogram:\n"%(model,method)+pl.mlab.rec2txt(params,precision=8))
    # params.tofile('log_params', sep=' ', format='%s')
    # logger.info("%s model parameters via %s periodogram:\n"%(model, method) + np.fromfile('log_params'))
    logger.info("%s model parameters via %s periodogram:\n"%(model,method) + str(params))
    #-- when full output is required, return parameters, periodogram and fitting
    #   function
    if full_output:
        mymodel = getattr(evaluate,model)(times,params)
        return params,(freqs_,ampls_),mymodel
    else:
        return params



def iterative_prewhitening(times,signal,maxiter=1000,optimize=0,method='scargle',
    model='sine',full_output=False,stopcrit=None,correlation_correction=True,
    prewhiteningorder_snr=False,prewhiteningorder_snr_window=1.,**kwargs):
    """
    Fit one or more functions to a timeseries via iterative prewhitening.

    This function will use C{find_frequency} to fit the function parameters to
    the original signal (including any previously found parameters),
    optimize the parameters if needed, and remove it from the data to search
    for a new frequency in the residuals.

    It is always the original signal that is used to fit all parameters again;
    B{only the (optimized) frequency is remembered from step to step} (Vanicek's
    method).

    You best set C{maxiter} to some sensable value, to hard-limit the number of
    frequencies that will be searched for. You can additionally use a C{stopcrit}
    and stop looking for frequencies once it is reached. C{stopcrit} should be
    a tuple; the first argument is the function to call, the other arguments
    are passed to the function, after the mandatory arguments C{times,signal,
    modelfunc,allparams,pergram}. See L{stopcrit_scargle_snr} for an example of
    such a function.

    By default, the function looks for the highest (or deepest in the case of the pdm
    method) peak, but instead it is possible to go for the peak having the highest
    SNR before prewhitening by setting C{prewhiteningorder_snr} to True. In this case,
    the noise spectrum is calculated using a convolution with a
    C{prewhiteningorder_snr_window} wide box. Usage of this is strongly encouraged,
    especially combined with L{stopcrit_scargle_snr} as C{stopcrit}.

    @return: parameters, model(, model function)
    @rtype: rec array(, ndarray)
    """
    residuals = signal.copy()
    frequencies = []
    stop_criteria = []
    while maxiter:
        #-- compute the next frequency from the residuals
        params,pergram,this_fit = find_frequency(times,residuals,method=method,
                full_output=True,correlation_correction=correlation_correction,
                prewhiteningorder_snr=prewhiteningorder_snr,
                prewhiteningorder_snr_window=prewhiteningorder_snr_window,**kwargs)

        #-- do the fit including all frequencies
        frequencies.append(params['freq'][-1])
        allparams = getattr(fit,model)(times,signal,frequencies)

        #-- if there's a need to optimize, optimize the last n parameters
        if optimize>0:
            residuals_for_optimization = residuals
            if optimize<=len(params):
                model_fixed_params = getattr(evaluate,model)(times,allparams[:-optimize])
                residuals_for_optimization -= model_fixed_params
            uparams,e_uparams, gain = fit.optimize(times,residuals_for_optimization,allparams[-optimize:],model)
            #-- only accept the optimization if we gained prediction power
            if gain>0:
                allparams[-optimize:] = uparams
                logger.info('Accepted optimization (gained %g%%)'%gain)

        #-- compute the residuals to use in the next prewhitening step
        modelfunc = getattr(evaluate,model)(times,allparams)
        residuals = signal - modelfunc

        #-- exhaust the counter
        maxiter -= 1

        #-- check stop criterion
        if stopcrit is not None:
            func = stopcrit[0]
            args = stopcrit[1:]
            condition,value = func(times,signal,modelfunc,allparams,pergram,*args)
            logger.info('Stop criterion (%s): %.3g'%(func.__name__,value))
            stop_criteria.append(value)
            if condition:
                logger.info('Stop criterion reached')
                break

    #-- calculate the errors
    e_allparams = getattr(fit,'e_'+model)(times,signal,allparams,correlation_correction=correlation_correction)

    allparams = numpy_ext.recarr_join(allparams,e_allparams)
    if stopcrit is not None:
        allparams = numpy_ext.recarr_join(allparams,np.rec.fromarrays([stop_criteria],names=['stopcrit']))

    if full_output:
        return allparams,modelfunc
    else:
        return allparams


def single_prewhitening(times, signal, residuals, params, optimize=0, model='sine', full_output=False,
                        correlation_correction=True):
    """
    Fit a functions to a timeseries via a single iteration of prewhitening.
    Use this function in combination with C{find_frequency} to do step-by-step
    prewhitening.

    This function will fit the function parameters to the original signal
    (including any previously found parameters given as C{params}) and
    optimize the parameters if needed.

    @return: parameters(, model)
    @rtype: rec array(, ndarray)
    """
    # do the fit including all frequencies
    allparams = getattr(fit, model)(times, signal, params['freq'])

    # if there's a need to optimize, optimize the last n parameters
    if optimize > 0:
        residuals_for_optimization = residuals
        if optimize <= len(params):
            model_fixed_params = getattr(evaluate, model)(times, allparams[:-optimize])
            residuals_for_optimization -= model_fixed_params
        uparams, e_uparams, gain = fit.optimize(times, residuals_for_optimization, allparams[-optimize:], model)
        # only accept the optimization if we gained prediction power
        if gain > 0:
            allparams[-optimize:] = uparams
            logger.info('Accepted optimization (gained %g%%)'%gain)

    # compute the model and the errors
    modelfunc = getattr(evaluate, model)(times, allparams)
    e_allparams = getattr(fit, 'e_'+model)(times, signal, allparams, correlation_correction=correlation_correction)

    allparams = numpy_ext.recarr_join(allparams, e_allparams)

    if full_output:
        return allparams, modelfunc
    else:
        return allparams


def spectrum_2D(x,y,matrix,weights_2d=None,show_progress=False,
                subs_av=True,full_output=False,**kwargs):
    """
    Compute a 2D periodogram.

    Rows (first axis) should be spectra chronologically

    x are time points (length N)
    y are second axis units (e.g. wavelengths) (length NxM)

    If the periodogram/wavelength combination has a large range, this script
    can produce a B{ValueError} or B{MemoryError}. To solve this, you could
    iterate this function over a subset of wavelength bins yourself, and write
    the results to a file.

    This function also outputs a model, which you can use to subtract from the
    data (probably together with the average profile, which is also returned)
    and look further for any frequencies::

        output = spectrum_2D(x,y,matrix,subs_av=True)
        new_matrix = matrix - output['avprof'] - output['model'
        output2 = spectrum_2D(x,y,new_matrix,subs_av=False)

    If you want to prewhiten a model, you'd probably want to use the same
    frequency across the whole line profile. You can hack this by setting
    C{f0=frequency} and C{fn=frequency+df} with C{df} the size of the frequency
    step.

    B{Example usage}: first we generate some variable line profiles. In this case,
    this is just a radial velocity variation

    >>> times = np.linspace(0,150,100)
    >>> wavel = np.r_[4500:4600:1.0]
    >>> matrix = np.zeros((len(times),len(wavel)))
    >>> for i in xrange(len(times)):
    ...   central_wave = 5*np.sin(2*np.pi/10*times[i])
    ...   matrix[i,:] = 1 - 0.5*np.exp( -(wavel-4550-central_wave)**2/10**2)

    Once the line profiles are constructed, we can compute the Fourier transform
    of these lines:

    >>> output = spectrum_2D(times,wavel,matrix,method='scargle',model='sine',f0=0.05,fn=0.3,subs_av=True,full_output=True)

    With this output, we can find retrieve the frequency. We plot the periodograms
    for each wavelength bin on the left, and on the right the average over all
    wavelength bins:

    >>> p = pl.figure()
    >>> p = pl.subplot(121)
    >>> p = pl.imshow(output['pergram'][1][::-1],extent=[wavel[0],wavel[-1],0.05,0.3],aspect='auto')
    >>> p = pl.subplot(122)
    >>> p = pl.plot(output['pergram'][1].mean(axis=1),output['pergram'][0],'k-')
    >>> p = pl.ylim(0.05,0.3)

    We can then fix the frequency and compute all the parameters. However, we now
    choose not to subtract the average profile, but instead use the GLS periodogram
    to include the constant

    >>> output = spectrum_2D(times,wavel,matrix,method='gls',model='sine',f0=0.095,fn=0.105,full_output=True)

    >>> p = pl.figure()
    >>> p = pl.subplot(221)
    >>> p = pl.errorbar(wavel,output['pars']['const'],yerr=output['pars']['e_const'],fmt='ko-')
    >>> p = pl.subplot(222)
    >>> p = pl.errorbar(wavel,output['pars']['ampl'],yerr=output['pars']['e_ampl'],fmt='ro-')
    >>> p = pl.subplot(223)
    >>> p = pl.errorbar(wavel,output['pars']['freq'],yerr=output['pars']['e_freq'],fmt='ko-')
    >>> p = pl.subplot(224)
    >>> p = pl.errorbar(wavel,output['pars']['phase'],yerr=output['pars']['e_phase'],fmt='ro-')

    @return: dict with keys C{avprof} (2D array), C{pars} (rec array), C{model} (1D array), C{pergram} (freqs,2Darray)
    @rtype: dict
    """
    #-- compute average profile
    if subs_av:
        matrix_av = np.outer(np.ones(len(matrix)),matrix.mean(axis=0))
        matrix = matrix - matrix_av
    else:
        matrix_av = 0.

    #-- prepare output of sine-parameters
    params = []
    freq_spectrum = []
    mymodel = []
    #-- do frequency analysis
    for iwave in range(len(matrix[0])):
        signal = matrix[:,iwave]
        if weights_2d is not None:
            weights = weights_2d[:,iwave]
            kwargs['weights'] = weights

        #-- make sure output is always a tuple, in case full output was asked
        #   we don't want iterative zoom in so set scale_df=0
        out = find_frequency(x,signal,full_output=full_output,scale_df=0,**kwargs)

        #-- add the parameters of this wavelength bin to the list
        if full_output:
            params.append(out[0])
            freq_spectrum.append(out[1][1])
            mymodel.append(out[2])
        else:
            params.append(out)

    #-- prepare output
    output = {}
    output['avprof']    = matrix_av
    output['pars']      = np.hstack(params)
    if full_output:
        output['pergram']   = out[1][0],np.vstack(freq_spectrum).T
        output['model']   = np.vstack(mymodel).T

    return output

@defaults_pergram
def time_frequency(times,signal,window_width=None,n_windows=100,
         window='rectangular',detrend=None,**kwargs):
    """
    Short Time (Fourier) Transform.

    Slide a window through the timeseries, multiply the timeseries with the
    window, and perform a Fourier Transform.

    It is best to fix explicitly C{f0}, C{fn} and C{df}, to limit the computation
    time!

    Extra kwargs go to L{find_frequency}

    @param n_windows: number of slices
    @type n_windows: integer
    @param window_width: width of each slice (defaults to T/20)
    @type window_width: float
    @param detrend: detrending function, accepting times and signal as args
    @type detrend: callable
    @param window: window function to apply
    @type window: string
    @return: spectrogram, times used, parameters, errors, points used per slice
    @rtype: dict
    """
    if window_width is None:
        window_width = kwargs.get('window_width',times.ptp()/20.)
    #-- cut light curve until window fits exactly
    stft_times = np.linspace(times[0]+window_width/2.,times[-1]-window_width/2.,n_windows)

    #-- get f0,fn and df to fix in all computations. If they are not given, they
    #   are set to default values by the decorator
    f0 = kwargs.pop('f0')
    fn = kwargs.pop('fn')
    df = kwargs.pop('df')
    nyq_stat = kwargs.pop('nyq_stat',fn)

    #-- prepare arrays for parameters, points and spectrum
    pars = []
    pnts = np.zeros(n_windows)
    spec = None
    #-- compute the periodogram for each slice
    for i,t in enumerate(stft_times):
        region = (abs(times-t) <= (window_width/2.))
        times_ = times[region]
        signal_ = signal[region]
        if detrend:
            times_,signal_ = detrend(times_,signal_)
        pnts[i] = len(times_)
        if len(times_)>1:
            output = find_frequency(times_,signal_,full_output=True,f0=f0,fn=fn,df=df,nyq_stat=nyq_stat,scale_df=0,**kwargs)
            if spec is None:
                spec = np.ones((n_windows,len(output[1][1])))
            pars.append(output[0])
            spec[i,:len(output[1][1])] = output[1][1]
        else:
            nanpars = np.rec.array(np.nan*np.ones(len(pars[-1].dtype.names)),dtype=pars[-1].dtype)
            pars.append(nanpars)
            #try:
            #    pars.append(np.nan*pars[-1])
            #except:
            #    print pars
            #    print pars[-1]
            #    print pars[-1].dtype
            #    print 0*pars[-1]
            #    print np.nan*pars[-1]
            #    raise
            spec[i] = np.nan
    out = {}
    out['times']     = stft_times
    out['pars']      = np.hstack(pars)
    out['pergram']   = (output[1][0],spec)
    return out

#{ Convenience stop-criteria

def stopcrit_scargle_prob(times,signal,modelfunc,allparams,pergram,crit_value):
    """
    Stop criterium based on probability.
    """
    value = pergrams.scargle_probability(pergram[1].max(),times,pergram[0])
    print(value)
    return value>crit_value,value

def stopcrit_scargle_snr(times,signal,modelfunc,allparams,pergram,crit_value,width=6.):
    """
    Stop criterium based on signal-to-noise ratio.
    """
    width = width/2.
    argmax = np.argmax(pergram[1])
    ampls = pergram[1]
    start = max(0,pergram[0][argmax]-width)
    stop = min(pergram[0][argmax]+width,pergram[0][-1])
    if start==0:
        stop += width-pergram[0][argmax]
    if stop==pergram[0][-1]:
        start = pergram[0][-1]-pergram[0][argmax]+width
    ampls = ampls[(start<=pergram[0]) & (pergram[0]<=stop)]
    value = pergram[1][argmax]/ampls.mean()
    return value<crit_value,value
#}

def autocorrelation(frequencies,power,max_step=1.5,interval=(),threshold=None,method=1):
    """
    Compute the autocorrelation.

    The formulate to estimate the autocorrelation in the periodogram is

    R(k) = 1 / (n * s^2) * \sum_{t=1}^n (X_t - mu) * (X_{t+k} - mu)

    where n is the number of points in the interval, mu an estimator for the
    sample mean and s is an estimator for the standard deviation of the sample.

    When computed over a large interval, we have to take into account that
    we cannot sample more points than we have in the spectrum, so we have to
    average out over fewer and fewer points, the more we shift the spectrum.

    Above formula becomes

    R(k) = 1 / ((n-k) * s^2) * \sum_{t=1}^{n-k} (X_t - mu) * (X_{t+k} - mu)

    Highs and lows are cut of. The lower cut off value is the sample mean, the
    higher cut off is a user-defined multiple of the sample mean.

    This function can also be used to compute period spacings: just invert the
    frequency spacing, reverse the order and reinterpolate to be equidistant:
    Example:
        >>> periods = linspace(1./p[0][-1],1./p[0][0],2*len(p[0]))
        >>> ampl = interpol.linear_interpolation(1./p[0][::-1],p[1][::-1],periods)
        >>> ampl = where(isnan(ampl),hstack([ampl[1:],ampl[:1]]),ampl)

    @param frequencies: frequency array
    @type frequencies: numpy 1d array
    @param power: power spectrum
    @type power: numpy 1d array
    @keyword max_step: maximum time shift
    @type max_step: float
    @keyword interval: tuple of frequencies (start, end)
    @type interval: tuple of floats
    @keyword threshold: high cut off value (in units of sample mean)
    @type threshold: float
    @return: domain of autocorrelation and autocorrelation
    @rtype: (ndarray,ndarray)
    """
    #-- cut out the interesting part of the spectrum
    if interval is not ():
        cut_out = frequencies[(interval[0]<=frequencies) & (frequencies<=interval[1])]
        start_freq = cut_out[0]
        stop_freq  = cut_out[-1]
        start = np.argmin(abs(frequencies-interval[0]))
        stop = np.argmin(abs(frequencies-interval[1]))
    else:
        start =  1
        stop  = len(frequencies)-1

    #-- compute the frequency step
    Dfreq = (frequencies[start+1] - frequencies[start+0])
    max_step = int(max_step/Dfreq)
    autocorr = []
    variance = []
    mean = np.average(power)

    #-- cut of high peaks at a signal to noise level of 6
    #   cut of low values at a signal to noise level of 1
    if threshold is not None:
        power[power>=(threshold*mean)] = threshold*mean
        #power = where(less(power, mean), mean, power)

    #-- normalize power as to arrive at the AUTO-correlation function.
    mean  = np.average(power)
    power = power-mean

    #-- compute autocorrelation. If nessecary, take border effects into
    #   account.
    for i in range(2,max_step):
        end_s = min(len(power), stop+i)
        end_o = start + end_s - (start+i)
        original = power[start  :end_o]
        shifted  = power[start+i:end_s]
        if len(original) < 10:
            logger.error("AUTOCORR: too few points left in interval, breaking up.")
            break
        if method==1:
            autocorr.append(np.average(original*shifted))
        else:
            autocorr.append(np.correlate(original,shifted))
        variance.append(np.average(original*original))

    domain = np.arange(2,max_step) * Dfreq
    domain = domain[0:len(autocorr)]

    #-- normalize
    autocorr = np.array(autocorr)/np.array(variance)
    logger.info("Computed autocorrelation in interval %s with maxstep %s"%(interval,max_step))
    return domain, autocorr

if __name__=="__main__":
    import doctest
    import pylab as pl
    import sys
    from ivs.aux import argkwargparser
    from ivs.inout import ascii

    #-- if no arguments are given, we just do a test run
    if not sys.argv[1:]:
        doctest.testmod()
        pl.show()
        sys.exit()

    #-- if arguments are given, we assume the user wants to run one of the
    #   functions with arguments given in the command line
    # EXAMPLES:
    # $:> python freqanalyse.py find_frequency infile=test.dat full_output=True
    # $:> python freqanalyse.py time_frequency infile=test.dat full_output=True
    else:
        method,args,kwargs = argkwargparser.parse()
        print("Running method %s with arguments %s and keyword arguments %s"%(method,args,kwargs))
        if '--help' in args or 'help' in args or 'help' in kwargs:
            sys.exit()
        full_output = kwargs.get('full_output',False)
        times,signal = ascii.read2array(kwargs.pop('infile')).T[:2]
        out = globals()[method](times,signal, **kwargs)

        #-- when find_frequency is called
        if method=='find_frequency' and full_output:
            print(pl.mlab.rec2txt(out[0],precision=8))
            pl.figure()
            pl.subplot(211)
            pl.plot(out[1][0],out[1][1],'k-')
            pl.subplot(212)
            pl.plot(times,signal,'ko',ms=2)
            pl.plot(times,out[2],'r-',lw=2)
            pl.show()
        elif method=='find_frequency':
            print(pl.mlab.rec2txt(out))

        #-- when time_frequency is called
        elif method=='time_frequency':
            print(pl.mlab.rec2txt(out['pars'],precision=8))
            pl.figure()
            pl.imshow(out['pergram'][1].T[::-1],aspect='auto',extent=[out['times'][0],out['times'][-1],out['pergram'][0][0],out['pergram'][0][-1]])


        pl.show()

