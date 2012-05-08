"""
Fit various functions to timeseries.

Section 1. Radial velocity data
===============================

Fit orbit of massive X-ray binary LSI+65010, after Grundstrom 2007:

Necessary imports:

>>> from ivs.catalogs import vizier
>>> from ivs.timeseries import pergrams
>>> from ivs.aux import numpy_ext
>>> import pylab as pl

Read in the data, and remove the outliers:

>>> data,units,comms = vizier.search('J/ApJ/656/431/table2')
>>> times,RV = data['HJD'],data['RV']
>>> keep = RV<-30.
>>> times,RV = times[keep],RV[keep]

Find the best frequency using the Kepler periodogram, fit an orbit with that
frequency and optimize. In the latter step, we also retrieve the errors on the
parameters, and print the results to the screen:

>>> freqs,ampls = pergrams.kepler(times,RV,fn=0.2)
>>> freq = freqs[np.argmax(ampls)]
>>> pars1 = kepler(times, RV, freq)
>>> pars2,e_pars2,gain = optimize(times,RV,pars1,'kepler')
>>> print pl.mlab.rec2txt(pars1,precision=6)
           P               T0          e      omega           K        gamma
   11.581314   2451060.751789   0.190000   1.006924   11.915330   -59.178393
>>> print pl.mlab.rec2txt(numpy_ext.recarr_join(pars2,e_pars2),precision=6)
        gamma           P               T0          e      omega           K    e_gamma        e_P       e_T0        e_e    e_omega        e_K
   -59.181942   11.581472   2451060.760523   0.194192   1.015276   11.925424   0.503331   0.004101   0.573320   0.060920   0.314982   0.788034

Evaluate the orbital fits, and make phasediagrams of the fits and the data

>>> myorbit1 = evaluate.kepler(times,pars1)
>>> myorbit2 = evaluate.kepler(times,pars2)
>>> phases,phased = evaluate.phasediagram(times,RV,1/pars1['P'])
>>> phases1,phased1 = evaluate.phasediagram(times,myorbit1,1/pars1['P'])
>>> phases2,phased2 = evaluate.phasediagram(times,myorbit2,1/pars1['P'])

Now plot everything:

>>> sa1 = np.argsort(phases1)
>>> sa2 = np.argsort(phases2)
>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(freqs,ampls,'k-')
>>> p = pl.xlabel('Frequency [d$^{-1}$]')
>>> p = pl.ylabel('Statistic')
>>> p = pl.subplot(122)
>>> p = pl.plot(phases,phased,'ko')
>>> p = pl.plot(phases1[sa1],phased1[sa1],'r-',lw=2,label='Linear fit')
>>> p = pl.plot(phases2[sa2],phased2[sa2],'b--',lw=2,label='Optimization')
>>> p = pl.xlabel('Phase [$2\pi^{-1}$]')
>>> p = pl.ylabel('Amplitude [km s$^{-1}$]')
>>> p = pl.legend()

]include figure]]ivs_sigproc_fit_01.png]

Section 2. Pulsation frequency analysis
=======================================

Do a frequency analysis of the star HD129929, after Aerts 2004:

Read in the data:

>>> data,units,comms = vizier.search('J/A+A/415/241/table1')
>>> times,signal = data['HJD'],data['Umag']
>>> signal -= signal.mean()

Find the best frequency using the Scargle periodogram, fit an orbit with that
frequency and optimize. Then print the results to the screen:

>>> freqs,ampls = pergrams.scargle(times,signal,f0=6.4,fn=7)
>>> freq = freqs[np.argmax(ampls)]
>>> pars1 = sine(times, signal, freq)
>>> e_pars1 = e_sine(times,signal, pars1)
>>> pars2,e_pars2,gain = optimize(times,signal,pars1,'sine')
>>> print pl.mlab.rec2txt(numpy_ext.recarr_join(pars1,e_pars1),precision=6)
      const       ampl       freq       phase    e_const     e_ampl     e_freq    e_phase
   0.000242   0.014795   6.461705   -0.093895   0.000608   0.001319   0.000006   0.089134
>>> print pl.mlab.rec2txt(numpy_ext.recarr_join(pars2,e_pars2),precision=6)
      const       ampl       freq       phase    e_const     e_ampl     e_freq    e_phase
   0.000242   0.014795   6.461705   -0.093895   0.000609   0.001912   0.000000   0.013386

Evaluate the sines, and make phasediagrams of the fits and the data

>>> mysine1 = evaluate.sine(times,pars1)
>>> mysine2 = evaluate.sine(times,pars2)
>>> phases,phased = evaluate.phasediagram(times,signal,pars1['freq'])
>>> phases1,phased1 = evaluate.phasediagram(times,mysine1,pars1['freq'])
>>> phases2,phased2 = evaluate.phasediagram(times,mysine2,pars1['freq'])

Now plot everything:

>>> sa1 = np.argsort(phases1)
>>> sa2 = np.argsort(phases2)
>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(freqs,ampls,'k-')
>>> p = pl.xlabel('Frequency [d$^{-1}$]')
>>> p = pl.ylabel('Amplitude [mag]')
>>> p = pl.subplot(122)
>>> p = pl.plot(phases,phased,'ko')
>>> p = pl.plot(phases1[sa1],phased1[sa1],'r-',lw=2)
>>> p = pl.plot(phases2[sa2],phased2[sa2],'b--',lw=2)
>>> p = pl.xlabel('Phase [$2\pi^{-1}$]')
>>> p = pl.ylabel('Amplitude [mag]')

]include figure]]ivs_sigproc_fit_02.png]


Section 3. Exoplanet transit analysis
=====================================

Find the transits of CoRoT 8b, after Borde 2010.

>>> import urllib
>>> from ivs.io import ascii
>>> url = urllib.URLopener()
>>> filen,msg = url.retrieve('http://cdsarc.u-strasbg.fr/viz-bin/nph-Plot/Vgraph/txt?J%2fA%2bA%2f520%2fA66%2f.%2flc_white&F=white&P=0&--bitmap-size&800x400')
>>> times,signal = ascii.read2array(filen).T
>>> signal = signal / np.median(signal)
>>> url.close()

Find the best frequency using the Box Least Squares periodogram, fit a transit
model with that frequency, optimize and prewhiten.

>>> freqs,ampls = pergrams.box(times,signal,f0=0.16,fn=0.162,df=0.005/times.ptp(),qma=0.05)
>>> freq = freqs[np.argmax(ampls)]
>>> pars = box(times,signal,freq)
>>> pars = box(times,signal,freq,b0=pars['ingress'][0]-0.05,bn=pars['egress'][0]+0.05)
>>> print pl.mlab.rec2txt(pars,precision=6)
       freq      depth    ingress     egress       cont
   0.161018   0.005978   0.782028   0.799229   1.000027

Evaluate the transits, and make phasediagrams of the fits and the data

>>> transit = evaluate.box(times,pars)
>>> phases,phased = evaluate.phasediagram(times,signal,freq)
>>> phases1,phased1 = evaluate.phasediagram(times,transit,freq)

Now plot everything and print the results to the screen:

>>> sa1 = np.argsort(phases1)
>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(freqs,ampls,'k-')
>>> p = pl.xlabel('Frequency [d$^{-1}$]')
>>> p = pl.ylabel('Statistic')
>>> p = pl.subplot(122)
>>> p = pl.plot(phases,phased*100,'ko')
>>> p = pl.plot(phases1[sa1],phased1[sa1]*100,'r-',lw=2)
>>> p = pl.xlim(0.70,0.85)
>>> p = pl.xlabel('Phase [$2\pi^{-1}$]')
>>> p = pl.ylabel('Depth [%]')

]include figure]]ivs_sigproc_fit_03.png]

Section 4. Eclipsing binary fit
===============================

Splines are not a good way to fit eclipsing binaries, but just for the sake of
showing the use of the periodic spline fitting functions, we do it anyway.

We use the data on CU Cnc of Ribas, 2003:

>>> data,units,comms = vizier.search('J/A+A/398/239/table1')
>>> times,signal = data['HJD'],data['Dmag']

""" 
import time
import logging

import numpy as np
from numpy import pi,cos,sin
import numpy.linalg as la
from scipy.interpolate import splrep
import scipy.optimize

from ivs.aux import loggers
from ivs.sigproc import evaluate
import pyKEP

import re
import copy
import pylab as pl
from ivs.sigproc import lmfit

logger = logging.getLogger('TS.FIT')

#{Fitting functions


def sine(times, signal, freq, sigma=None,constant=True,error=False,t0=0):
    """
    Fit a harmonic function.
    
    This function is of the form
    
    C + \sum_j A_j sin(2pi nu_j (t-t0) + phi_j)
    
    where the presence of the constant term is an option. The error
    bars on the fitted parameters can also be requested (by Joris De Ridder).
    
    (phase in radians!)
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param freq: frequencies of the harmonics
    @type freq: numpy array or float
    @keyword sigma: standard error of observations
    @type sigma: numpy array
    @keyword constant: flag, if not None, also fit a constant
    @type constant: boolean
    @keyword error: flag, if not None, also compute the errorbars
    @type error: boolean
    @keyword t0: time zero point.
    @type t0: float
    @return: parameters
    @rtype: record array
    """
    #-- Subtract the zero point from the time points.
    times = times - t0
    
    #-- Prepare the input: if a frequency value is given, put it in a list. If
    #   an iterable is given convert it to an array
    if not hasattr(freq,'__len__'):
        freq = [freq]
    freq = np.asarray(freq)
    #   do the same for the sigmas
    if sigma is None:
        sigma = np.ones_like(signal)
    elif not hasattr(sigma,'__len__'):
        sigma = sigma * np.ones_like(signal)
    
    #-- Determine the number of fit parameters
    Ndata = len(times)
    Nfreq = len(freq)
    if not constant: Nparam = 2*Nfreq
    else:            Nparam = 2*Nfreq+1
        
    #-- The fit function used is of the form
    #      C + \sum_j a_j sin(2pi\nu_j t_i) + b_j cos(2pi\nu_j t_i)
    #   which is linear in its fit parameters. These parameters p can therefore be
    #   solved by minimizing ||b - A p||, where b are the observations, and A is the
    #   basisfunction matrix, i.e. A[i,j] is the j-th base function evaluated in 
    #   the i-th timepoint. The first Nfreq columns are the amplitudes of the sine,
    #   the second Nfreq columns are the amplitudes of the cosine, and if requested,
    #   the last column belongs to the constant
    A = np.zeros((Ndata,Nparam))
    for j in range(Nfreq):
        A[:,j]       = sin(2*pi*freq[j]*times) * sigma
        A[:,Nfreq+j] = cos(2*pi*freq[j]*times) * sigma
    if constant:
        A[:,2*Nfreq] = np.ones(Ndata) * sigma
    
    b = signal * sigma
  
    #-- Solve using SVD decomposition
    fitparam, chisq, rank, s = np.linalg.lstsq(A,b)
  
    #-- Compute the amplitudes and phases: A_j sin(2pi*\nu_j t_i + phi_j)
    amplitude = np.zeros(Nfreq)
    phase = np.zeros(Nfreq)
    for j in range(Nfreq):
        amplitude[j] = np.sqrt(fitparam[j]**2 + fitparam[Nfreq+j]**2)
        phase[j]     = np.arctan2(fitparam[Nfreq+j], fitparam[j])

    #-- If no error bars are needed, we are finished here, we collect all parameters    
    if constant:
        constn = np.zeros(len(amplitude))
        constn[0] = fitparam[2*Nfreq]
        names = ['const','ampl','freq','phase']
        fpars = [constn,amplitude,freq,phase/(2*pi)]
    else:
        names = ['ampl','freq','phase']
        fpars = [amplitude,freq,phase/(2*pi)]
        
    parameters = np.rec.fromarrays(fpars,names=names)
    
    logger.debug('SINEFIT: Calculated harmonic fit with %d frequencies through %d datapoints'%(Nfreq,Ndata))
        
    return parameters

def periodic_spline(times, signal, freq, t0=None, order=20, k=3):
    """
    Fit a periodic spline.
        
    CAUTION: this definition assumes the proposed period is either
    small compared to the total time range, or there are a lot of
    points per period.
    
    This definition basically phases all the data and constructs
    an empirical periodic function through an averaging process
    per phasebin, and then performing a splinefit through those points.
    
    In order to make the first derivative continuous, we repeat
    the first point(s) at the end and the end point(s) at the beginning.
    
    The constructed function can than be removed from the original
    data.
    
    Output is a record array containing the columns 'freq', 'knots', 'coeffs'
    and 'degree'.
    
    Example usage:    
    
    >>> myfreq = 1/20.
    >>> times = np.linspace(0,150,10000)
    >>> signal = sin(2*pi*myfreq*times+0.32*2*pi) + np.random.normal(size=len(times))
    >>> cjs = periodic_spline(times,signal,myfreq,order=30)
    >>> trend = evaluate.periodic_spline(times,cjs)
    >>> import pylab as pl
    >>> p=pl.figure()
    >>> p=pl.plot(times,signal,'ko')
    >>> p=pl.plot(times,trend,'r-')
    >>> p=pl.title("test:tfit:periodic_spline_fit")
    
    @param times: time points
    @type times: numpy 1d array
    @param signal: observation points
    @type signal: numpy 1d array
    @param freq: frequency of the periodic trend
    @type freq: float
    @keyword order: number of points to use for the spline fit.
    @type order: integer
    @keyword k: order of the spline
    @type k: integer
    @return: parameters
    @rtype: record array
    """
    #-- get keyword arguments
    if t0 is None:
        t0 = times[0]
        
    #-- Prepare the input: if a frequency value is given, put it in a list. If
    #   an iterable is given convert it to an array
    if not hasattr(freq,'__len__'):
        freq = [freq]
    freq = np.asarray(freq)
    
    N = order*3+(k-2)
    parameters = np.zeros(len(freq), dtype=[('freq','float32'),
                                    ('knots','%dfloat32'%N),
                                    ('coeffs','%dfloat32'%N),
                                    ('degree','int8')])
    for nf,ifreq in enumerate(freq):
        #-- get phased signal
        phases,phased_sig = evaluate.phasediagram(times,signal,ifreq,t0=t0)
        
        #-- construct phase domain
        x = np.linspace(0.,1.,order)[:-1]
        dx = x[1]-x[0]
        x = x+dx/2
        
        #-- make bin-averaged signal. Possible that some bins are empty, just skip
        #   them but warn the user
        found_phase = []
        found_averg = []
        for xi in x:
            this_bin = ((xi-dx/2.)<=phases) & (phases<=(xi+dx/2.))
            if np.sum(this_bin)==0:
                logger.warning("PERIODIC SPLINE: some parts of the phase are not covered")
                continue
            found_phase.append(xi)
            found_averg.append(np.average(phased_sig[this_bin]))
        
        found_phase = np.array(found_phase)
        found_averg = np.array(found_averg)
        
        #-- make circular
        found_phase = np.hstack([found_phase-1,found_phase,found_phase+1])
        found_averg = np.hstack([found_averg,  found_averg,found_averg])
        
        #-- compute spline representation
        knots,coeffs,degree = splrep(found_phase,found_averg,k=k)
        parameters['freq'][nf] = ifreq
        parameters['knots'][nf] = knots
        parameters['coeffs'][nf] = coeffs
        parameters['degree'][nf] = degree

    return parameters

def kepler(times, signal, freq, sigma=None, wexp=2., e0=0, en=0.99, de=0.01):
    """
    Fit a Kepler orbit to a time series.
    
    Example usage:
    
    >>> from ivs.timeseries import pergrams
    
    First set the parameters we want to use:
    
    >>> pars = tuple([12.456,23456.,0.37,213/180.*pi,98.76,55.])
    >>> pars = np.rec.array([pars],dtype=[('P','f8'),('T0','f8'),('e','f8'),
    ...                                ('omega','f8'),('K','f8'),('gamma','f8')])
                                       
    Then generate the signal and add some noise
    
    >>> times = np.linspace(pars[0]['T0'],pars[0]['T0']+5*pars[0]['P'],100)
    >>> signalo = evaluate.kepler(times,pars)
    >>> signal = signalo + np.random.normal(scale=20.,size=len(times))
    
    Calculate the periodogram:
    
    >>> freqs,ampls = pergrams.kepler(times,signal)
    >>> opars = kepler(times,signal,freqs[np.argmax(ampls)])
    >>> signalf = evaluate.kepler(times,opars)
    
    And make some plots
    
    >>> import pylab as pl
    >>> p = pl.figure()
    >>> p = pl.subplot(221)
    >>> p = pl.plot(times,signal,'ko')
    >>> p = pl.plot(times,signalo,'r-',lw=2)
    >>> p = pl.plot(times,signalf,'b--',lw=2)
    >>> p = pl.subplot(222)
    >>> p = pl.plot(freqs,ampls,'k-')
    
    @param times: time points
    @type times: numpy 1d array
    @param signal: observation points
    @type signal: numpy 1d array
    @param freq: frequency of kepler orbit
    @type freq: float
    @return: parameters
    @rtype: record array
    """
    if sigma is None:
        sigma = np.ones_like(times)
    T = times[-1]-times[0]
    x00 = 0.
    x0n = 359.9
    #-- initialize parameters
    f0 = freq
    fn = freq+0.05/T
    df = 0.1/T
    maxstep = int((fn-f0)/df+1)
    f1 = np.zeros(maxstep) #-- frequency
    s1 = np.zeros(maxstep) #-- power
    p1 = np.zeros(maxstep) #-- window
    l1 = np.zeros(maxstep) #-- power LS
    s2 = np.zeros(maxstep) #-- power Kepler
    k2 = np.zeros(6) #-- parameters for Kepler orbit
    
    pyKEP.kepler(times,signal,sigma,f0,fn,df,wexp,e0,en,de,x00,x0n,
         f1,s1,p1,l1,s2,k2)
    
    freq,x0,e,w,K,RV0 = k2
    pars = tuple([1/freq,x0/(2*pi*freq) + times[0], e, w, K, RV0])
    pars = np.rec.array([tuple(pars)],dtype=[('P','f8'),('T0','f8'),('e','f8'),('omega','f8'),('K','f8'),('gamma','f8')])
    return pars


  
def box(times,signal,freq,b0=0,bn=1,order=50,t0=None):
    """
    Fit box shaped transits.
    
    @param times: time points
    @type times: numpy 1d array
    @param signal: observation points
    @type signal: numpy 1d array
    @param freq: frequency of transiting signal
    @type freq: float
    @param b0: minimum start of ingress (in phase)
    @type b0: 0<float<bn
    @param bn: maximum end of egress (in phase)
    @type bn: b0<float<1
    @param order: number of phase bins
    @type order: integer
    @param t0: zeropoint of times
    @type t0: float
    @return: parameters
    @rtype: record array
    """
    if t0 is None:
        t0 = 0.
        
    #-- Prepare the input: if a frequency value is given, put it in a list. If
    #   an iterable is given convert it to an array
    if not hasattr(freq,'__len__'):
        freq = [freq]
    freq = np.asarray(freq)
    
    parameters = np.rec.fromarrays([np.zeros(len(freq)) for i in range(5)],dtype=[('freq','f8'),('depth','f8'),
                                            ('ingress','f8'),('egress','f8'),
                                            ('cont','f8')])
    bins = np.linspace(b0,bn,order)
    
    for fnr,frequency in enumerate(freq):
        parameters[fnr]['freq'] = frequency
        #-- we need to check two phase diagrams, to compensate for edge effects
        phases1,phased1 = evaluate.phasediagram(times,signal,frequency,t0=t0)
        phases2,phased2 = evaluate.phasediagram(times,signal,frequency,t0=t0+0.5/frequency)
        best_fit = np.inf
        for i,start in enumerate(bins):
            for end in bins[i+1:]:
                transit1 = (start<=phases1) & (phases1<=end)
                transit2 = (start<=phases2) & (phases2<=end)
                transit_sig1 = np.median(np.compress(  transit1,phased1))
                contini_sig1 = np.median(np.compress(1-transit1,phased1))
                depth1 = contini_sig1-transit_sig1
                transit_sig2 = np.median(np.compress(  transit2,phased2))
                contini_sig2 = np.median(np.compress(1-transit2,phased2))
                depth2 = contini_sig2-transit_sig2
                #-- check fit
                chisquare1 = sum((evaluate.box(times,[contini_sig1,frequency,depth1,start,end])-signal)**2)
                chisquare2 = sum((evaluate.box(times,[contini_sig2,frequency,depth2,start,end])-signal)**2)
                if chisquare1 < best_fit and chisquare1 < chisquare2:
                    parameters[fnr]['depth'] = depth1
                    parameters[fnr]['ingress'] = start
                    parameters[fnr]['egress'] = end
                    parameters[fnr]['cont'] = contini_sig1
                    best_fit = chisquare1
                elif chisquare2 < best_fit and chisquare2 < chisquare1:
                    start2 = start - 0.5
                    end2 = end - 0.5
                    if start2 < 0: start2 += 1
                    if end2 < 0: end2 += 1
                    parameters[fnr]['depth'] = depth2
                    parameters[fnr]['ingress'] = start2
                    parameters[fnr]['egress'] = end2
                    parameters[fnr]['cont'] = contini_sig2
                    best_fit = chisquare2
    return parameters

def gauss(x,y,threshold=0.1,constant=False,full_output=False,
          init_guess_method='analytical',window=None):
    """
    Fit a Gaussian profile to data using a polynomial fit.
    
    y = A * exp( -(x-mu)**2 / (2*sigma**2))
    
    ln(y) = ln(A) - (x-mu)**2 / (2*sigma**2)
    ln(y) = ln(A) - mu**2 / (2*sigma**2) + mu / (sigma**2) * x - x**2 / (2*sigma**2)
    
    then the parameters are given by
    
    p0 =       -  1    / (2*sigma**2)
    p1 =         mu    / (  sigma**2)
    p2 = ln(A) - mu**2 / (2*sigma**2)
    
    Note that not all datapoints are used, but only those above a certain values (namely
    10% of the maximum value), In this way, we reduce the influence of the continuum
    and concentrate on the shape of the peak itself.
    
    Afterwards, we perform a non linear least square fit with above parameters
    as starting values, but only accept it if the CHI2 has improved.
    
    If a constant has to be fitted, the nonlinear options has to be True.
    
    Example: we generate a Lorentzian curve and fit a Gaussian to it:
    
    >>> x = np.linspace(-10,10,1000)
    >>> y = evaluate.lorentz(x,[5.,0.,2.]) + np.random.normal(scale=0.1,size=len(x))
    >>> pars1,e_pars1 = gauss(x,y)
    >>> pars2,e_pars2 = gauss(x,y,constant=True)
    >>> y1 = evaluate.gauss(x,pars1)
    >>> y2 = evaluate.gauss(x,pars2)
    >>> p = pl.figure()
    >>> p = pl.plot(x,y,'k-')
    >>> p = pl.plot(x,y1,'r-',lw=2)
    >>> p = pl.plot(x,y2,'b-',lw=2)
    
    @param x: x axis data
    @type x: numpy array
    @param y: y axis data
    @type y: numpy array
    @param nl: flag for performing a non linear least squares fit
    @type nl: boolean
    @param constant: fit also a constant
    @type constant: boolean
    @rtype: tuple
    @return: A, mu, sigma (,C)
    """
    #-- if a constant needs to be determined, take the median of the 5% lowest
    #   points (but take at least three):
    if constant:
        N = len(y)
        C = np.median(np.sort(y)[:max(int(0.05*N),3)])
    else:
        C = 0.
        
    if window is not None:
        win = (window[0]<=x) & (x<=window[1])
        x,y = x[win],y[win]
    
    #-- transform to a polynomial function and perform a fit
    #   first clip where y==0
    threshold *= max(y-C)
    use_in_fit = (y-C)>=threshold
    xc = x[use_in_fit]
    yc = y[use_in_fit]
    
    lnyc = np.log(yc-C)
    p   = np.polyfit(xc, lnyc, 2)
    
    #-- determine constants
    sigma = np.sqrt(-1. / (2.*p[0]))
    mu    = sigma**2.*p[1]
    A     = np.exp(p[2] + mu**2. / (2.*sigma**2.))
    
    #-- handle NaN Exceptions
    if A!=A or mu!=mu or sigma!=sigma:
        logger.error('Initial Gaussian fit failed')
        init_success = False
        A = (yc-C).max()
        mu = xc[yc.argmax()]
        sigma = xc.ptp()/3.
    else:
        init_success = True
    
    #-- check chi2
    if constant:
        p0 = evaluate.gauss_preppars(np.asarray([A,mu,sigma,C]))
    else:
        p0 = evaluate.gauss_preppars(np.asarray([A,mu,sigma]))
    yf = evaluate.gauss(xc,p0)
    chi2 = np.sum((yc-yf)**2)
    
    #-- perform non linear least square fit
    if constant:
        pars,e_pars,gain = optimize(x,y,p0,'gauss', minimizer='leastsq')
    else:
        pars,e_pars,gain = optimize(x,y,p0,'gauss', minimizer='leastsq')
        if gain<0:
            pars = p0
    pars['sigma'] = np.abs(pars['sigma'])
    if not full_output:
        return pars,e_pars
    else:
        return pars,e_pars,p0,use_in_fit,init_success


#}
#{ Error determination

def e_sine(times,signal,parameters,correlation_correction=True,limit=10000):
    """
    Compute the errors on the parameters from a sine fit.
    
    Note: errors on the constant are only calculated when the number of datapoints
    is below 1000. Otherwise, the matrices involved become to huge.
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param parameters: record array containing the fitted parameters. This should
    have columns 'ampl','freq' and 'phase', and optionally 'const'.
    @type parameters: numpy record array
    @param correlation_correction: set to True if you want to correct for correlation effects
    @type correlation_correction: boolean
    @param limit: Calculating the error on the constant requires the calculation
    of a matrix of size Ndata x Nparam, which takes a long time for large datasets.
    The routines skips the estimation of the error on the constant if the timeseries
    is longer than C{limit} datapoints
    @type limit: integer
    @return: errors
    @rtype: Nx4 array(, Nx3 array)
    """
    #-- these are the parameters we need
    freq = parameters['freq']
    phase = parameters['phase']
    amplitude = parameters['ampl']
    chisq = np.sum((signal - evaluate.sine(times,parameters))**2)
    
    #-- for quick reference, we also need the dimensions of the data and
    #   fit parameters
    Ndata = len(times)
    Nparam = 2*len(parameters)
    Nfreq = len(freq)
    T = times.ptp()
    
    #-- do we need to include the constant?
    if 'const' in parameters.dtype.names:
        constant = True
        Nparam += 1   
    
    #-- these lists will contain the columns and their names
    errors = []
    names = []
    
    #-- If error bars on the constant are needed, we do it here. The errors
    #   on the amplitude, frequency and phase are computed below.
    if constant and Ndata<limit:
        #   First the derivative matrix: the 1st Nfreq columns the derivative w.r.t.
        #   the amplitude, the 2nd Nfreq columns w.r.t. the phase, and the if relevant
        #   the last column w.r.t. the constant. From this the covariance matrix.
        F = np.zeros((Ndata, Nparam))   
        for i in range(Ndata):
            F[i,0:Nfreq] = sin(2*pi*freq*times[i] + phase)
            F[i,Nfreq:2*Nfreq] = amplitude * cos(2*pi*freq*times[i] + phase)
            #-- and for the constant
            F[i,2*Nfreq] = 1.0 
      
        covariance = np.linalg.inv(np.dot(F.T, F))
        covariance *= chisq / (Ndata - Nparam)
  
        #error_ampl = np.sqrt(covariance.diagonal()[:Nfreq])
        #error_phase = np.sqrt(covariance.diagonal()[Nfreq:2*Nfreq])
        error_const = np.sqrt(covariance[2*Nfreq,2*Nfreq])
        error_const = np.ones(Nfreq)*error_const
        if len(error_const)>1:
            error_const[1:] = 0
        errors.append(error_const)
        names.append('e_const')
    elif constant:
        errors.append(np.zeros(Nfreq))
        names.append('e_const')
    
    #-- other variables
    errors_ = np.zeros((Nfreq,3))
    residus = signal + 0.
    for i in range(1,Nfreq+1):
        residus -= evaluate.sine(times,parameters[i-1:i])
        a       = parameters[i-1]['ampl']
        sigma_m = np.std(residus)
    
        #-- error on amplitude, frequency and phase (in that order)
        errors_[i-1,0] = np.sqrt(2./Ndata) * sigma_m
        errors_[i-1,1] = np.sqrt(6./Ndata) * 1. / (pi*T) * sigma_m / a
        errors_[i-1,2] = np.sqrt(2./Ndata) * sigma_m / a
        
        #-- correct for correlation effects
        if correlation_correction:
            rho = max(1,get_correlation_factor(residus))
            errors_[i-1,0] *= np.sqrt(rho)
            errors_[i-1,1] *= np.sqrt(rho)
            errors_[i-1,2] *= np.sqrt(rho)
    
    #-- collect results and return
    e_parameters = np.rec.fromarrays(errors+[errors_[:,0],errors_[:,1],errors_[:,2]],
                         names=names + ['e_ampl','e_freq','e_phase'])
                         
    return e_parameters


#{ Linear improvements

def diffcorr(times, signal, parameters, func_name, \
                  max_iter=100, tol=1e-6,full_output=False):
    """
    Differential corrections.
    """
    
    #-- ensure we have a flat array, and look for the functions to evaluate the
    #   fits and the coefficients
    prep_func = getattr(evaluate,func_name+'_preppars')
    diff_func = getattr(evaluate,func_name+'_diffcorr')
    eval_func = getattr(evaluate,func_name)
    if parameters.dtype.names:
        parameters = prep_func(parameters)
    
    #-- prepare arrays for coefficients and new parameters
    Deltas = np.zeros((max_iter,len(parameters)))
    params = np.zeros_like(Deltas)
    params[0] = parameters
    counter = 1
    while (counter==1) or (counter>0 and counter<max_iter and np.any(np.abs(Deltas[counter-1])>tol)):
        params[counter] = params[counter-1] + Deltas[counter-1]
        myfit = eval_func(times,params[counter])
        coeff = diff_func(times,params[counter])
        Delta,res,rank,s = la.lstsq(coeff,myfit-signal)
        Deltas[counter] = -Delta
        counter += 1
    
    #-- transform the parameters to record arrays, as well as the steps
    parameters = prep_func(params[counter-1])
    e_parameters = prep_func(Deltas[counter-1])
    e_parameters.dtype.names = ['e_'+name for name in e_parameters.dtype.names]
    
    if full_output:
        return parameters,e_parameters,0,params[:counter],Deltas[:counter]
    else:
        return parameters,e_parameters,0
    

#}


#{ Non-linear improvements

def residuals(parameters,domain,data,evalfunc,*args):
    fit = evalfunc(domain,parameters,*args)
    return data-fit

def residuals_single(parameters,domain,data,evalfunc):
    fit = evalfunc(domain,parameters)
    return sum((data-fit)**2)

def optimize(times, signal, parameters, func_name, minimizer='leastsq', args=()):
    #-- we need these function to evaluate the fit and to (un)pack the fitting
    #   parameters from and to flat arrays
    prepfunc = getattr(evaluate,func_name+'_preppars')
    evalfunc = getattr(evaluate,func_name)
    optifunc = getattr(scipy.optimize,minimizer)
    
    #-- if the initial guess of the fitting parameters aren't flat, flatten them
    #   here:
    if parameters.dtype.names:
        parameters = prepfunc(parameters)
    init_guess = parameters.copy()
    #-- keep track of the initial chi square value, to check if there is an
    #   improvement
    dof = (len(times)-len(init_guess))
    signalf_init = evalfunc(times,init_guess)
    chisq_init = np.sum((signalf_init-signal)**2)
    
    #-- optimize
    if minimizer=='leastsq':
        popt, cov, info, mesg, flag = optifunc(residuals,init_guess,
                                     args=(times,signal,evalfunc)+args,full_output=1)#,diag=[1.,10,1000,1.,100000000.])
        #-- calculate new chisquare, and check if we have improved it
        chisq = np.sum(info['fvec']*info['fvec'])
        if chisq>chisq_init or flag!=1:
            logger.error('Optimization not successful [flag=%d] (%g --> %g)'%(flag,chisq_init,chisq))
            chisq = np.inf
    
        #-- derive the errors from the nonlinear fit
        if cov is not None:
            errors = np.sqrt(cov.diagonal()) * np.sqrt(chisq/dof)
        else:
            logger.error('Error estimation via optimize not successful')
            errors = np.zeros(len(popt))
    else:
        out = optifunc(residuals_single,init_guess,
                                     args=(times,signal,evalfunc),full_output=1,disp=False)
        popt = out[0]
        #-- calculate new chisquare, and check if we have improved it
        signalf_update = evalfunc(times,popt)
        chisq = np.sum((signalf_update-signal)**2)
        errors = np.zeros(len(popt))
        if chisq>chisq_init:
            logger.error('Optimization not successful')
    
    #-- gain in chi square: if positive, we gained, if negative, we lost...
    gain = (chisq_init-chisq)/chisq_init*100.
    
    #-- transform the parameters to record arrays, as well as the errors
    parameters = prepfunc(popt)
    
    e_parameters = prepfunc(errors)
    e_parameters.dtype.names = ['e_'+name for name in e_parameters.dtype.names]
    
    return parameters,e_parameters, gain

class Model(object):
    
    def __init__(self, functions=None):
        self.functions = functions
        #-- Combine the parameters
        self.pull_parameters(functions)
    
    #{ Internal
    def pull_parameters(self, functions):
        """
        Pulls the parameter objects from the underlying functions, and combines it to 1 parameter object.
        """
        parameters = []
        for func in functions:
            parameters.append(func.parameters)
        
        #-- Create new parameter object
        new_params = lmfit.Parameters()
        pnames = []
        for i, params in enumerate(parameters):
            pname = []
            for n,par in params.items():
                pname.append(n+'_%i'%(i))
                new_params.add(n+'_%i'%(i), value=par.value, vary=par.vary, min=par.min, max=par.max, expr=par.expr)
            pnames.append(pname)
                
        self.par_names = pnames
        self.parameters = new_params
        
    def push_parameters(self, parameters=None):
        """
        Pushes the parameters in the combined parameter object to the parameter objects of the underlying 
        models or functions.
        """
        if parameters == None:
            parameters = self.parameters
        
        for pnames,function in zip(self.par_names, self.functions):
            old_parameters = function.parameters
            for name in pnames:
                old_name = re.sub('_[0123456789]*$','',name)
                old_parameters[old_name] = parameters[name]
    
    #}
    
    #{ Interaction
    
    def evaluate(self, *args):
        """
        Evaluate the model for a given values and optional a given parameter object.
        If no parameter object is given then the parameter object belonging to the model
        is used.
        
        >>> evaluate(parameters, x)
        >>> evaluate(x)
        """
        if len(args) == 1:
            #-- Use the parameters belonging to this object
            parameters = self.parameters
            x = args[0]
        elif len(args) == 2:
            #-- Use the provided parameters
            parameters = args[0]
            x = args[1]
        
        #-- Update the parameters of the individual functions before calling them
        self.push_parameters(parameters=parameters)
        
        #-- For each function, read the arguments and calculate the result
        result = np.zeros(len(x))
        for function in self.functions:
            result += function.evaluate(x)
               
        return result
        
    def setup_parameters():
        raise NotImplementedError
        
    #}    
    
    #{ Getters and setters
    
    def get_model_function(self):
        """
        Returns the model function that can be evaluated using function(parameters, x)
        """
        return self.model
        
    def get_parameters_object(self):
        """
        Return the parameter object belonging to the model
        """
        return self.parameters
    
    #}

class Function(object):
    
    def __init__(self, function=None, par_names=None, jacobian=None):
        self.function = function
        self.par_names = par_names
        self.jacobian = jacobian
        
        #create an empty parameter set based on the parameter names
        self.parameters = None
        self.setup_parameters()
    
    #{ Interaction
    
    def evaluate(self,*args):
        """
        Evaluate the function for the given values and optional the given parameter object.
        If no parameter object is given then the parameter object belonging to the function
        is used.
        
        >>> evaluate(parameters, x)
        >>> evaluate(x)
        """
        if len(args) == 1:
            #-- Use the parameters belonging to this object
            pars = []
            for name in self.par_names:
                pars.append(self.parameters[name].value)
                
            return self.function(pars,args[0])
            
        if len(args) == 2:
            #-- Use the provided parameters
            pars = []
            for name in self.par_names:
                pars.append(args[0][name].value)
                
            return self.function(pars,args[1])
            
    def setup_parameters(self,values=None, bounds=None, vary=None, exprs=None):
        """
        Create or adjust a parameter object based on the parameter names and if provided
        the values, bounds, vary and expressions.
        """
        nrpars = len(self.par_names)
        if values == None:
            values = [0 for i in range(nrpars)]
        if bounds == None:
            bounds = [(None,None) for i in range(nrpars)]
        if vary == None:
            vary = [True for i in range(nrpars)]
        if exprs == None:
            exprs = [None for i in range(nrpars)]
        
        if self.parameters == None:
            #-- Create a new parameter object
            self.parameters = lmfit.Parameters()
            for i,name in enumerate(self.par_names):
                self.parameters.add(name, value=values[i], vary=vary[i], min=bounds[i][0], max=bounds[i][1], expr=exprs[i])
        else:
            #-- Adjust an existing parameter object
            for i,name in enumerate(self.par_names):
                self.parameters[name].value = values[i]
                self.parameters[name].vary = vary[i]
                self.parameters[name].min = bounds[i][0]
                self.parameters[name].max = bounds[i][1]
                self.parameters[name].expr = exprs[i]
    
    def update_parameter(self, parameter=None, **kwargs):
        """
        Updates a specified parameter. The parameter can be given by name or by index
        """
        
        if type(parameter) == int:
            parameter = self.parameters[self.par_names[parameter]]
        elif type(parameter) == str:
            parameter = self.parameters[parameter]
        
        for key in parameter.keys():
            if key in kwargs:
                parameter[key] = kwargs[key]
    
    def get_parameters(self, full_output=False):
        """
        Returns the parameter values together with the errors if they exist. If No fitting
        has been done, or the errors could not be calculated, None is returned for the error.
        
        @param full_output: When True, also vary, the boundaries and expr are returned
        @type full_output: bool
        
        @return: the parameter values and there errors [(value, err, vary, min, max, expr),...]
        @rtype: array of tupples
        """
        out = []
        for name, par in self.parameters.items():
            if full_output:
                out.append((par.value, par.stderr, par.vary, par.min, par.max, par.expr))
            else:
                out.append((par.value, par.stderr))
        return out
        
    
    def param2str(self, full_output=False, accuracy=2):
        """
        Converts the parameter object of this function to an easy printable string
        """
        return parameters2string(self.parameters, accuracy=accuracy, full_output=full_output)
        
    #}
    
    #{ Getters and setters
    
    def get_model_function(self):
        """
        Returns the function that can be evaluated using function(parameters, x)
        """
        return self.function
        
    def get_parameters_object(self):
        """
        Return the parameter object belonging to the function
        """
        return self.parameters
        
    def get_par_names(self):
        """
        returns the names of the different parameters
        """
        return self.par_names
    
    def get_jacobian(self):
        return self.jacobian
    
    #}
    
    
class Minimizer(lmfit.Minimizer): 

    def __init__(self, x, y, model, err=None, weights=None,
             engine='leastsq', args=None, kws=None,scale_covar=True,iter_cb=None, **kwargs):
        
        self.x = x
        self.y = y
        self.model = model
        self.err = err
        self.weights = weights
        
        params = model.get_parameters_object()
        
        #-- Setup the residual function and the lmfit.minimizer object
        if weights == None:
            def residuals(params, x, y):
                return (y - model.evaluate(params,x))
            fcn_args = (x,y)
                             
        else:
            def residuals(params, x, y, weights):
                return (y - model.evaluate(params,x))*weights
            fcn_args = (x,y,weights)
        
        lmfit.Minimizer.__init__(self,residuals, params, fcn_args=fcn_args, fcn_kws=kws,
                             iter_cb=iter_cb, scale_covar=scale_covar, **kwargs)
        
        #-- Actual fitting
        if engine == 'anneal':
            self.anneal()
        elif engine == 'lbfgsb':
            self.lbfgsb()
        else:
            self.leastsq()
    
    def get_confidence_interval(self, p_names=None, sigmas=[0.65,0.95,0.99], maxiter=200, prob_func=None):
        """
        Returns the confidence intervalls of the given parameters. 
        """
        
        # if only 1 confidence intervall is asked, the output can be tupple instead of dict.
        short_output = (type(p_names)==str and type(sigmas)==float) and True or False
        if type(p_names)==str: p_names = [p_names]
        if type(sigmas)==float: sigmas = [sigmas]
        
        #Use the adjusted conf_interval() function of the lmfit package.
        out = lmfit.conf_interval(self, p_names=p_names, sigmas=sigmas, maxiter=maxiter, prob_func=prob_func, trace=False, verbose=False)
        
        if short_output:
            out = out[p_names[0]][sigmas[0]]
        return out
    
    def MC_simulation():
        raise NotImplementedError
    
    #{ Plotting Functions
    
    def plot_results(self, eval_points=1000):
        """
        Creates a basic plot with the fit results and corresponding residuals.
        
        @param eval_points: Number of points to use when plotting the best fit model.
        @type eval_points: int
        """
        
        xf = np.linspace(min(self.x),max(self.x),eval_points)
        yf = self.model.evaluate(xf)
        
        pl.subplots_adjust(left=0.10, bottom=0.1, right=0.97, top=0.95,wspace=0.0, hspace=0.0)
        
        ax = pl.subplot2grid((3,4), (0,0), rowspan=2, colspan=4)
        if self.err == None:
            pl.plot(self.x,self.y,'+b')
        else:
            pl.errorbar(self.x,self.y, yerr=self.err, ls='', marker='+', color='b')
        pl.plot(xf,yf,'-r')
        xlim = pl.xlim([min(self.x)-0.05*(max(self.x)-min(self.x)), max(self.x)+0.05*(max(self.x)-min(self.x))])
        pl.ylim([min(self.y)-0.05*(max(self.y)-min(self.y)), max(self.y)+0.05*(max(self.y)-min(self.y))])
        pl.ylabel('$y$')
        for tick in ax.axes.get_xticklabels():
            tick.set_visible(False)
            tick.set_fontsize(0.0)
            
        ax = pl.subplot2grid((3,4), (2,0), colspan=4)
        if self.err == None:
            pl.plot(self.x,self.y-self.model.evaluate(self.x), '+b')
        else:
            pl.errorbar(self.x,self.y-self.model.evaluate(self.x), yerr=self.err, ls='', marker='+', color='b')
        pl.axhline(y=0, ls=':', color='r')
        pl.xlim(xlim)
        pl.ylabel('$O-C$')
        pl.xlabel('$x$')
    
    def plot_confidence_interval(self,xname=None,yname=None, res=10, filled=True, limits=None):
        """
        Plot the confidence interval for 2 given parameters. 
        
        @param xname: The parameter on the x axis
        @param yname: The parameter on the y axis
        @param res: The resolution of the grid over which the confidence intervall is calculated
        @param filled: True for filled contour plot, False for normal contour plot
        @param limits: The upper and lower limit on the parameters for which the confidence intervall is calculated. If None, 5 times the stderr is used.
        """
        
        xn = hasattr(res,'__iter__') and res[0] or res
        yn = hasattr(res,'__iter__') and res[1] or res
        
        x, y, grid = lmfit.conf_interval2d(self,xname,yname,xn,yn, limits=limits)
        grid *= 100.
        
        pl.subplots_adjust(left=0.10, bottom=0.1, right=0.97, top=0.95,wspace=0.0, hspace=0.0)
        if filled:
            pl.contourf(x,y,grid,np.linspace(0,100,25),cmap=pl.cm.jet)
            pl.colorbar(fraction=0.08,ticks=[0,20,40,60,80,100])
        else:
            cs = pl.contour(x,y,grid,np.linspace(0,100,11),cmap=pl.cm.jet)
            cs = pl.contour(x,y,grid,[20,40,60,80,95],cmap=pl.cm.jet)
            pl.clabel(cs, inline=1, fontsize=10)
        pl.plot(self.params[xname].value, self.params[yname].value, '+r', ms=10, mew=2)
        pl.xlabel(xname)
        pl.ylabel(yname)

def minimize(x, y, model, err=None, weights=None,
             engine='leastsq', args=None, kws=None,scale_covar=True,iter_cb=None, **fit_kws):
    """
    Basic minimizer function, returns a Parameters and Minimizer object
    """
    
    fitter = Minimizer(x, y, model, err=err, weights=weights,
             engine=engine, args=args, kws=kws, scale_covar=scale_covar,iter_cb=iter_cb, **fit_kws)
    
    return fitter

#}

#{ General purpose

def get_correlation_factor(residus, full_output=False):
    """
    Calculate the correlation facor rho (Schwarzenberg-Czerny, 2003).
    
    Under white noise assumption, the residus are expected to change sign every
    2 observations (rho=1). Longer distances, 2*rho, are a sign of correlation.
    
    The errors are then underestimated by a factor 1/sqrt(rho).
    
    @param residus: residus after the fit
    @type residus: numpy array
    @param full_output: if True, the groups of data with same sign will be returned
    @type full_output: bool
    @return: rho(,same sign groups)
    @rtype: float(,list)
    """
    same_sign_groups = [1]
    
    for i in xrange(1,len(residus)):
        if np.sign(residus[i])==np.sign(residus[i-1]):
            same_sign_groups[-1] += 1
        else:
            same_sign_groups.append(0)
    
    rho = np.average(same_sign_groups)
    
    logger.debug("Correlation factor rho = %f, sqrt(rho)=%f"%(rho,np.sqrt(rho)))
    
    if full_output:
        return rho, same_sign_groups
    else:
        return rho

#}

#{ Print functions

def parameters2string(parameters, accuracy=2, full_output=False):
    #Converts a parameter object to string
    out = ""
    for name, par in parameters.items():
        if not full_output:
            out +=  '%10s = %.2f +/- %.2f \n' % (name, par.value, par.stderr)
        else:
            if par.vary:
                out +=  '%10s = %.2f +/- %.2f \t bounds = %s <-> %s \n' % (name, par.value, par.stderr, par.min, par.max)
            else:
                out +=  '%10s = %.2f +/- %.2f \t bounds = %s <-> %s (fixed) \n' % (name, par.value, par.stderr, par.min, par.max)
    return out

def confidence2string(ci, accuracy=2):
    #Converts confidence intervall dictionary to string
    out=""
    for par in ci.keys():
        out += "%s \n\t "%(par)
        sigmas = ci[par].keys()
        sigmas.sort()
        for sigma in sigmas:
            out += "%9s %% \t"%(np.round(sigma*100, decimals=1))
        out += '\n\t-'
        for sigma in sigmas:
            out += "%10s \t"%(np.round(ci[par][sigma][0], decimals=accuracy))
        out += '\n\t+'
        for sigma in sigmas:
            out += "%10s \t"%(np.round(ci[par][sigma][1], decimals=accuracy))   
        out += '\n'
    return out        

#}












if __name__=="__main__":
    import doctest
    import pylab as pl
    import sys
    doctest.testmod()
    pl.show()
    sys.exit()

    from ivs.timeseries import pergrams
    import pylab as pl

    times = np.linspace(0,150,5000)
    np.random.seed(10)
    freqs = [1.23,2.59,3.89,4.65]
    freqs += [2*freqs[0],2*freqs[0]+freqs[2],freqs[3]-freqs[0],freqs[0]+freqs[3]-freqs[2]]
    freqs = np.array(freqs)
    amplitudes = np.sort(np.random.uniform(size=len(freqs),low=1,high=10))[::-1]
    phases = np.random.uniform(size=len(freqs),low=-0.5,high=0.5)
    parins = np.rec.fromarrays([np.zeros(len(freqs)),amplitudes,freqs,phases],names=('const','ampl','freq','phase'))
    signal = evaluate.sine(times,parins)
    signal_= signal + np.random.normal(size=len(signal),scale=1)

    residus = signal_ + 0.
    frequencies = []

    for i in range(9):
        print "======== STEP %d ======"%(i)
        
        
        pergram = pergrams.scargle(times,residus,threads=2)
        frequency = pergram[0][np.argmax(pergram[1])]
        frequencies.append(frequency)
        parameters = sine(times,signal_,frequencies)
        e_parameters = e_sine(times,signal_,parameters)
        parameters,e_parameters,gain = optimize(times,signal_,parameters,'sine')
        frequencies = list(parameters['freq'])
        
        signalf = evaluate.sine(times,parameters)
        
        
        if i<len(parins):
            print ''
            for par in ['const','ampl','freq','phase']:
                print par,parins[i][par],parameters[par],e_parameters['e_%s'%(par)]
            print 'S/N',parameters['ampl']/e_parameters['e_ampl']            
            
        else:
            print ''
            for par in ['const','ampl','freq','phase']:
                print par,parameters[par],e_parameters['e_%s'%(par)]
            print 'S/N',parameters['ampl']/e_parameters['e_ampl']
        
        
        
        print ''
        
        
        
        residus = signal_-signalf


    parameters   = sine(times,signal_,frequencies)
    signalf = evaluate.sine(times,parameters)

    pl.subplot(121)
    pl.plot(times,signal_,'ko')
    pl.plot(times,signal,'r-')
    pl.plot(times,signalf,'b-')

    pl.subplot(122)
    pl.plot(*pergram)
    pl.show()
