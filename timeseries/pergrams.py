# -*- coding: utf-8 -*-
"""
Contains many different periodogram calculations

Section 1. Basic usage
======================

Given a time series of the form

>>> times = np.linspace(0,1,1000)
>>> signal = np.sin(times)

The basic interface is

>>> freq,ampl = scargle(times,signal)

If you give no extra information, default values for the start, end and step
frequency will be chosen. See the 'defaults_pergram' decorator for a list of
keywords common to all periodogram calculations.

Many periodograms can be computed in parallel, by supplying an extra keyword
'threads':

>>> freq,ampl = scargle(times,signal,threads=2)

B{Warning}: the timeseries must be B{sorted in time} and B{cannot contain the
same timepoint twice}. Otherwise, a 'ValueError, concatenation problem' can
occur.

If something goes wrong in the periodogram computation, be sure to run
L{check_input} on your input data. This will print out some basic diagnostics
to see if your data are valid.

Section 2. Nyquist frequency
============================

The periodogram functions are written such that they never exceed the value
of the Nyquist frequency. This behaviour can be changed.

By default, the Nyquist frequency is defined as half of the inverse of the smallest
time step in the data. That is, the C{nyq_stat} variable is set to the function
C{np.min}. If you prefer the Nyquist to be defined as the median, set the
C{nyq_stat} variable for any periodogram to C{np.median}. If you want a more
complex or self-defined function, that is also acceptable. If you give a number
to C{nyq_stat}, nothing will be computed but that value will be considered the
nyquist frequency.

Section 3. Speed comparison
===========================

>>> import time

The timeseries is generated for N=500,1000,2000,5000,10000,20000,50000

>>> techniques = ['scargle','deeming','gls','fasper',
...               'schwarzenberg_czerny','pdm','box']
>>> Ns = [1000,2000,5000,10000,20000,30000]
>>> for tech in techniques[:0]:
...     freqs = []
...     clock = []
...     for N in Ns:
...         times = np.linspace(0,1,N)
...         signal = np.sin(times) + np.random.normal(size=N)
...         c0 = time.time()
...         freq,ampl = locals()[tech](times,signal)
...         clock.append(time.time()-c0)
...         freqs.append(len(freq))
...     p = pl.plot(freqs,clock,'o-',label=tech)
...     print tech,np.polyfit(np.log10(np.array(freqs)),np.log10(np.array(clock)),1)
>>> #p = pl.legend(loc='best',fancybox=True)
>>> #q = p.get_frame().set_alpha(0.5)
>>> #p,q = pl.xlabel('Number of frequencies'),pl.ylabel('Seconds')

]]include figure]]ivs_timeseries_pergram_speeds.png]

Section 2. Periodogram comparison
=================================

We generate a sinusoidal signal, add some noise and compute the periodogram with
the different methods.

>>> times = np.sort(np.random.uniform(size=500,low=0,high=100))
>>> signal = np.sin(2*np.pi/10.*times)
>>> signal += np.random.normal(size=500)

Fourier based periodgrams suited for unequidistant data: L{scargle}, L{deeming},
L{fasper} and L{schwarzenberg_czerny}.

>>> f1,a1 = scargle(times,signal,fn=0.35)
>>> f2,a2 = fasper(times,signal,fn=0.35)
>>> f3,a3 = deeming(times,signal,fn=0.35)
>>> f4,a4 = scargle(times,signal,norm='power',fn=0.35)
>>> f5,a5 = schwarzenberg_czerny(times,signal,fn=0.35,nh=1,mode=2)

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(f1,a1,'k-',lw=4,label='scargle')
>>> p = pl.plot(f2,a2,'r--',lw=4,label='fasper')
>>> p = pl.plot(f3,a3,'b-',lw=2,label='deeming')
>>> p = pl.xlim(f1[0],f1[-1])
>>> leg = pl.legend(fancybox=True,loc='best')
>>> leg.get_frame().set_alpha(0.5)
>>> p = pl.subplot(122)
>>> p = pl.plot(f4,a4,'k-',lw=2,label='scargle')
>>> p = pl.plot(f5,a5/2.,'r--',lw=2,label='schwarzenberg-czerny')
>>> p = pl.xlim(f1[0],f1[-1])
>>> leg = pl.legend(fancybox=True,loc='best')
>>> leg.get_frame().set_alpha(0.5)

]]include figure]]ivs_timeseries_pergrams_fourier.png]

Least-square fitting: L{gls} and L{kepler}.

>>> f1,a1 = gls(times,signal,fn=0.35)
>>> f2,a2 = kepler(times,signal,fn=0.35)

>>> p = pl.figure()
>>> p = pl.plot(f1,a1,'k-',lw=2,label='gls')
>>> p = pl.plot(f2,a2,'r--',lw=2,label='kepler')
>>> p = pl.xlim(f1[0],f1[-1])
>>> leg = pl.legend(fancybox=True,loc='best')
>>> leg.get_frame().set_alpha(0.5)

]]include figure]]ivs_timeseries_pergrams_lsq.png]

Phase folding techniques: L{box} and L{pdm}

>>> f1,a1 = box(times,signal,fn=0.35)
>>> f2,a2 = pdm(times,signal,fn=0.35)

>>> p = pl.figure()
>>> p = pl.plot(f1,a1,'k-',lw=2,label='box')
>>> p = pl.plot(f2,a2,'r--',lw=2,label='pdm')
>>> p = pl.xlim(f1[0],f1[-1])
>>> leg = pl.legend(fancybox=True,loc='best')
>>> leg.get_frame().set_alpha(0.5)

]]include figure]]ivs_timeseries_pergrams_phase.png]

"""
import logging
import numpy as np
from numpy import cos,sin,pi
from scipy.special import jn
from ivs.aux.decorators import make_parallel
from ivs.aux import loggers
from ivs.aux import termtools
from ivs.timeseries.decorators import parallel_pergram,defaults_pergram,getNyquist

import pyscargle
import pyscargle_single
import pyfasper
import pyfasper_single
import pyclean
import pyGLS
import pyKEP
import pydft
import multih
import deeming as fdeeming
import eebls

logger = logging.getLogger("TS.PERGRAMS")


#{ Periodograms


@defaults_pergram
@parallel_pergram
@make_parallel
def scargle(times, signal, f0=None, fn=None, df=None, norm='amplitude',
            weights=None, single=False):
    """
    Scargle periodogram of Scargle (1982).
    
    Several options are available (possibly combined):
        1. weighted Scargle
        2. Amplitude spectrum
        3. Distribution power spectrum
        4. Traditional power spectrum
        5. Power density spectrum (see Kjeldsen, 2005 or Carrier, 2010)
    
    This definition makes use of a Fortran-routine written by Jan Cuypers, Conny
    Aerts and Peter De Cat. A slightly adapted version is used for the weighted
    version (adapted by Pieter Degroote).
    
    Through the option "norm", it's possible to norm the periodogram as to get a
    periodogram that has a known statistical distribution. Usually, this norm is
    the variance of the data (NOT of the noise or residuals, see Schwarzenberg-
    Czerny 1998!).
    
    Also, it is possible to retrieve the power density spectrum in units of
    [ampl**2/frequency]. In this routine, the normalisation constant is taken
    to be the total time span T. Kjeldsen (2005) chooses to multiply the power
    by the 'effective length of the observing run', which is calculated as the
    reciprocal of the area under spectral window (in power, and take 2*Nyquist
    as upper frequency value).
    
    REMARK: this routine does B{not} automatically remove the average. It is the
    user's responsibility to do this adequately: e.g. subtract a B{weighted}
    average if one computes the weighted periodogram!!
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param weights: weights of the datapoints
    @type weights: numpy array
    @param norm: type of normalisation
    @type norm: str
    @param f0: start frequency
    @type f0: float
    @param fn: stop frequency
    @type fn: float
    @param df: step frequency
    @type df: float
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """ 
    if single: pyscargle_ = pyscargle_single
    else:
        pyscargle_ = pyscargle
    #-- initialize variables for use in Fortran routine
    sigma=0.;xgem=0.;xvar=0.;n=len(times)
    T = times.ptp()
    nf=int((fn-f0)/df+0.001)+1
    f1=np.zeros(nf,'d');s1=np.zeros(nf,'d')
    ss=np.zeros(nf,'d');sc=np.zeros(nf,'d');ss2=np.zeros(nf,'d');sc2=np.zeros(nf,'d')
    
    #-- run the Fortran routine
    if weights is None:
        f1,s1=pyscargle_.scar2(signal,times,f0,df,f1,s1,ss,sc,ss2,sc2)
    else:
        w=np.array(weights,'float')
        logger.debug('Weighed scargle')
        f1,s1=pyscargle_.scar3(signal,times,f0,df,f1,s1,ss,sc,ss2,sc2,w)
    
    #-- search for peaks/frequencies/amplitudes    
    if not s1[0]: s1[0]=0. # it is possible that the first amplitude is a none-variable
    fact  = np.sqrt(4./n)
    if norm =='distribution': # statistical distribution
        s1 /= np.var(signal)
    elif norm == "amplitude": # amplitude spectrum
        s1 = fact * np.sqrt(s1)
    elif norm == "density": # power density
        s1 = fact**2 * s1 * T    
    return f1, s1



@defaults_pergram
@parallel_pergram
@make_parallel
def fasper(times,signal, f0=None, fn=None, df=None, single=True, norm='amplitude'):
    """
    Fasper periodogram from Numerical Recipes.
    
    Normalisation here is not correct!!
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param f0: start frequency
    @type f0: float
    @param fn: stop frequency
    @type fn: float
    @param df: step frequency
    @type df: float
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """
    #-- average nyquist frequency and oversampling rate
    nyq = 1./(2*np.diff(times).mean())
    mynyq = 1./(2*np.diff(times).min())
    T = times.ptp()
    ofac = 1./(df*T)
    hifac = fn/mynyq*mynyq/nyq
    #-- prepare input for fasper
    n = len(times)
    nout = int(4*ofac*hifac*n*4.)
    wk1 = np.zeros(nout)
    wk2 = np.zeros(nout)
    jmax,prob = 0,0.
    #import pyfasper2
    if not single:
        wk1,wk2,nwk,nout,jmax,prob = pyfasper.fasper(times,signal,ofac,hifac,wk1,wk2,nout,jmax,prob)
    else:
        wk1,wk2,nwk,nout,jmax,prob = pyfasper_single.fasper(times,signal,ofac,hifac,wk1,wk2,nout,jmax,prob)
    #wk1,wk2,nout,jmax,prob = fasper_py(times,signal,ofac,hifac)
    wk1,wk2 = wk1[:nout],wk2[:nout]*1.5
    fact  = np.sqrt(4./n)
    if norm =='distribution': # statistical distribution
        wk2 /= np.var(signal)
    elif norm == "amplitude": # amplitude spectrum
        wk2 = fact * np.sqrt(wk2)
    elif norm == "density": # power density
        wk2 = fact**2 * wk2 * T     
    if f0 is not None:
        keep = f0<wk1
        wk1,wk2 = wk1[keep],wk2[keep]
    return wk1,wk2









@defaults_pergram
@parallel_pergram
@make_parallel
def deeming(times,signal, f0=None, fn=None, df=None, norm='amplitude'):
    """
    Deeming periodogram of Deeming et al (1975).
    
    Thanks to Jan Cuypers
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param norm: type of normalisation
    @type norm: str
    @param f0: start frequency
    @type f0: float
    @param fn: stop frequency
    @type fn: float
    @param df: step frequency
    @type df: float
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """        
    #-- initialize variables for use in Fortran routine
    nf=int((fn-f0)/df+0.001)+1
    n = len(times)
    T = times.ptp()
    f1,s1 = fdeeming.deeming1(times,signal,f0,df,nf)
    s1 /= n
    fact  = np.sqrt(4./n)
    fact  = np.sqrt(4./n)
    if norm =='distribution': # statistical distribution
        s1 /= np.var(signal)
    elif norm == "amplitude": # amplitude spectrum
        s1 = fact * np.sqrt(s1)
    elif norm == "density": # power density
        s1 = fact**2 * s1 * T
    
    return f1,s1
    

@defaults_pergram
@parallel_pergram
@make_parallel
def gls(times,signal, f0=None, fn=None, df=None, errors=None, wexp=2):
    """
    Generalised Least Squares periodogram of Zucher et al (2010).
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param f0: start frequency
    @type f0: float
    @param fn: stop frequency
    @type fn: float
    @param df: step frequency
    @type df: float
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """
    T = times.ptp()
    n = len(times)
    if errors is None:
        errors = np.ones(n)
    maxstep = int((fn-f0)/df+1)
    
    #-- initialize parameters
    f1 = np.zeros(maxstep) #-- frequency
    s1 = np.zeros(maxstep) #-- power
    p1 = np.zeros(maxstep) #-- window
    l1 = np.zeros(maxstep) #-- power LS
    
    #-- calculate generalized least squares
    pyGLS.gls(times+0.,signal+0.,errors,f0,fn,df,wexp,f1,s1,p1,l1)
    return f1,s1





@defaults_pergram
@parallel_pergram
@make_parallel
def clean(times,signal, f0=None, fn=None, df=None, freqbins=None, niter=10.,
          gain=1.0):
    """
    Cleaned Fourier periodogram of Roberts et al (1987)
        
    Parallization probably isn't such a good idea here because of the frequency
    bins.
    
    Fortran module probably from John Telting.
    
    Should always start from zero, so f0 is not an option
    
    Generate some signal with heavy window side lobes:
    
    >>> times_ = np.linspace(0,150,1000)
    >>> times = np.array([times_[i] for i in xrange(len(times_)) if (i%10)>7])
    >>> signal = np.sin(2*pi/10*times) + np.random.normal(size=len(times))
    
    Compute the scargle periodogram as a reference, and compare with the
    the CLEAN periodogram with different gains.
    
    >>> niter,freqbins = 10,[0,1.2]
    >>> p1 = scargle(times,signal,fn=1.2,norm='amplitude',threads=2)
    >>> p2 = clean(times,signal,fn=1.2,gain=1.0,niter=niter,freqbins=freqbins)
    >>> p3 = clean(times,signal,fn=1.2,gain=0.1,niter=niter,freqbins=freqbins)
    
    Make a figure of the result:
    
    >>> p=pl.figure()
    >>> p=pl.plot(p1[0],p1[1],'k-',label="Scargle")
    >>> p=pl.plot(p2[0],p2[1],'r-',label="Clean (g=1.0)")
    >>> p=pl.plot(p3[0],p3[1],'b-',label="Clean (g=0.1)")
    >>> p=pl.legend()
    
    ]]include figure]]ivs_timeseries_pergrams_clean.png]
    
    @keyword freqbins: frequency bins for clean computation
    @type freqbins: list or array
    @keyword niter: number of iterations
    @type niter: integer
    @keyword gain: gain for clean computation
    @type gain: float between 0 (no cleaning) and 1 (full cleaning)
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """
    T = times.ptp()
    n = len(times)
    if freqbins is None:
        freqbins = [f0,fn]
    
    startfreqs = np.array(freqbins[0::2])
    endfreqs = np.array(freqbins[1::2])
    nbins = len(freqbins)-1
    
    nf = int(fn/df)
    
    #-- do clean computation, seems not so straightforward to thread cleaning
    f,wpow,wpha = pyclean.main_clean(times,signal,fn,nf,gain,niter,nbins,\
                    startfreqs,endfreqs)
    
    return f,wpow









@defaults_pergram
@parallel_pergram
@make_parallel
def schwarzenberg_czerny(times, signal, f0=None, fn=None, df=None, nh=2, mode=1):
    """
    Multi harmonic periodogram of Schwarzenberg-Czerny (1996).
    
    This periodogram follows an F-distribution, so it is possible to perform
    hypothesis testing.
    
    If the number of the number of harmonics is 1, then this peridogram reduces
    to the Lomb-Scargle periodogram except for its better statistic behaviour.
    This script uses a Fortran procedure written by Schwarzenberg-Czerny.
    
    Modes:
        - mode=1: AoV Fisher-Snedecor F(df1,df2) statistic
        - mode=2: total power fitted in all harmonics for mode=2
    
    @param times: list of observations times
    @type times: numpy 1d array
    @param signal: list of observations
    @type signal: numpy 1d array
    @keyword f0: start frequency (cycles per day) (default: 0.)
    @type f0: float
    @keyword fn: stop frequency (cycles per day) (default: 10.)
    @type fn: float
    @keyword df: step frequency (cycles per day) (default: 0.001)
    @type df: float
    @keyword nh: number of harmonics to take into account
    @type nh: integer
    @return: frequencies, f-statistic
    @rtype: array,array
    """
    T = times.ptp()
    n = len(times)
    frequencies = np.arange(f0, fn+df, df)
    ll   = len(frequencies)
    th   = np.zeros(len(frequencies))
    #-- use Fortran subroutine
    th  = multih.sfou(n,times,signal,ll,f0,df,nh,mode,th)
    
    # th *= 0.5 seemed necessary to fit the F-distribution
        
    return frequencies,th
    
def DFTpower(time, signal, f0=None, fn=None, df=None,full_output=False):

    """
    Computes the modulus square of the fourier transform. 
    
    Unit: square of the unit of signal. Time points need not be equidistant.
    The normalisation is such that a signal A*sin(2*pi*nu_0*t)
    gives power A^2 at nu=nu_0
    
    @param time: time points [0..Ntime-1] 
    @type time: ndarray
    @param signal: signal [0..Ntime-1]
    @type signal: ndarray
    @param f0: the power is computed for the frequencies
                      freq = arange(startfreq,stopfreq,stepfreq)
    @type f0: float
    @param fn: see startfreq
    @type fn: float
    @param df: see startfreq
    @type df: float
    @return: power spectrum of the signal
    @rtype: array 
    """
    freqs = np.arange(f0,fn,df)
    Ntime = len(time)
    Nfreq = int(np.ceil((fn-f0)/df))
  
    A = np.exp(1j*2.*pi*f0*time) * signal
    B = np.exp(1j*2.*pi*df*time)
    ft = np.zeros(Nfreq, complex) 
    ft[0] = A.sum()
    for k in range(1,Nfreq):
        A *= B
        ft[k] = np.sum(A)
    
    if full_output:
        return freqs,ft**2*4.0/Ntime**2
    else:
        return freqs,(ft.real**2 + ft.imag**2) * 4.0 / Ntime**2    


def DFTpower2(time, signal, freqs):

    """
    Computes the power spectrum of a signal using a discrete Fourier transform.

    @param time: time points, not necessarily equidistant
    @type time: ndarray
    @param signal: signal corresponding to the given time points
    @type signal: ndarray
    @param freqs: frequencies for which the power spectrum will be computed. Unit: inverse of 'time'.
    @type freqs: ndarray
    @return: power spectrum. Unit: square of unit of 'signal'
    @rtype: ndarray
    """
    
    powerSpectrum = np.zeros(len(freqs))

    for i, freq in enumerate(freqs):
        arg = 2.0 * np.pi * freq * time
        powerSpectrum[i] = np.sum(signal * np.cos(arg))**2 + np.sum(signal * np.sin(arg))**2

    powerSpectrum = powerSpectrum * 4.0 / len(time)**2
    return(powerSpectrum)


    
def DFTscargle(times, signal,f0,fn,df):
    
    """
    Compute Discrete Fourier Transform for unevenly spaced data ( Scargle, 1989).
    
    Doesn't work yet!
    
    It is recommended to start f0 at 0.
    It is recommended to stop fn  at the nyquist frequency
    
    This makes use of a FORTRAN algorithm written by Scargle (1989).
    
    @param times: observations times
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param f0: start frequency
    @type f0: float
    @param fn: end frequency
    @type fn: float
    @param df: frequency step
    @type df: float
    @return: frequencies, dft, Re(dft), Im(dft)
    """
    f0 = 0.
    #fn *= 2*np.pi
    #df *= 2*np.pi
    #-- initialize
    nfreq = int((fn-f0)/(df))
    print "Nfreq=",nfreq
    tzero = times[0]
    si = 1.
    lfreq = 2*nfreq+1
    mm = 2*nfreq
    ftrx = np.zeros(mm)
    ftix = np.zeros(mm)
    om = np.zeros(mm)
    w = np.zeros(mm)
    wz = df #pi/(nn+dt)
    nn = len(times)
    
    #-- calculate DFT
    ftrx,ftix,om,w = pydft.ft(signal,times,wz,nfreq,si,lfreq,tzero,df,ftrx,ftix,om,w,nn,mm)
    
    if f0==0:
        ftrx[1:] *= np.sqrt(2)
        ftix[1:] *= np.sqrt(2)
        w[1:] *= np.sqrt(2)
        
    cut_off = len(w)
    for i in range(0,len(w))[::-1]:
        if w[i] != 0: cut_off = i+1;break
    
    om = om[:cut_off]
    ftrx = ftrx[:cut_off]
    ftix = ftix[:cut_off]
    w = w[:cut_off]
    
    # norm amplitudes for easy inversion
    T = times[-1]-times[0]
    N = len(times)
    
    w *= T/(2.*N)
    ftrx *= T/(2.*N)
    ftix *= T/(2.*N)
    
    return om*2*np.pi,w,ftrx,ftix
    
    
def FFTpower(signal, timestep):

    """
    Computes power spectrum of an equidistant time series 'signal'
    using the FFT algorithm. The length of the time series need not
    be a power of 2 (zero padding is done automatically). 
    Normalisation is such that a signal A*sin(2*pi*nu_0*t)
    gives power A^2 at nu=nu_0  (IF nu_0 is in the 'freq' array)
    
    @param signal: the time series [0..Ntime-1]
    @type signal: ndarray
    @param timestep: time step fo the equidistant time series
    @type timestep: float
    @return: frequencies and the power spectrum
    @rtype: array,array
    
    """
  
    # Compute the FFT of a real-valued signal. If N is the number 
    # of points of the original signal, 'Nfreq' is (N/2+1).
  
    fourier = np.fft.rfft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)
  
    # Compute the power
  
    power = np.abs(fourier)**2 * 4.0 / Ntime**2
  
    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep
  
    freq = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep
  
    # That's it!
  
    return (freq, power)





  
  
  
def FFTpowerdensity(signal, timestep):
  
    """
    Computes the power density of an equidistant time series 'signal',
    using the FFT algorithm. The length of the time series need not
    be a power of 2 (zero padding is done automatically). 

    @param signal: the time series [0..Ntime-1]
    @type signal: ndarray
    @param timestep: time step fo the equidistant time series
    @type timestep: float
    @return: frequencies and the power density spectrum
    @rtype: array,array

    """
  
    # Compute the FFT of a real-valued signal. If N is the number 
    # of points of the original signal, 'Nfreq' is (N/2+1).
  
    fourier = np.fft.rfft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)
  
    # Compute the power density
  
    powerdensity = np.abs(fourier)**2 / Ntime * timestep
  
    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep
  
    freq = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep
  
    # That's it!
  
    return (freq, powerdensity)

  


  
  



def weightedpower(time, signal, weight, freq):

    """
    Compute the weighted power spectrum of a time signal.
    For each given frequency a weighted sine fit is done using
    chi-square minimization.
    
    @param time: time points [0..Ntime-1] 
    @type time: ndarray
    @param signal: observations [0..Ntime-1]
    @type signal: ndarray
    @param weight: 1/sigma_i^2 of observation
    @type weight: ndarray
    @param freq: frequencies [0..Nfreq-1] for which the power 
                 needs to be computed
    @type freq: ndarray
    @return: weighted power [0..Nfreq-1]
    @rtype: array
    
    """

    result = np.zeros(len(freq))

    for i in range(len(freq)):
        if (freq[i] != 0.0):
            sine   = np.sin(2.0*pi*freq[i]*time)
            cosine = np.cos(2.0*pi*freq[i]*time)
            a11= np.sum(weight*sine*sine)
            a12 = np.sum(weight*sine*cosine)
            a21 = a12
            a22 = np.sum(weight*cosine*cosine)
            b1 = np.sum(weight*signal*sine)
            b2 = np.sum(weight*signal*cosine)
            denominator = a11*a22-a12*a21
            A = (b1*a22-b2*a12)/denominator
            B = (b2*a11-b1*a21)/denominator
            result[i] = A*A+B*B
        else:
            result[i] = np.sum(signal)/len(signal)

    return(result)    
    
    
    
    






@defaults_pergram
@parallel_pergram
@make_parallel
def pdm(times, signal,f0=None,fn=None,df=None,Nbin=5,Ncover=2,
         D=0,forbit=None,asini=None,e=None,omega=None,nmax=10):
    """
    Phase Dispersion Minimization of Jurkevich-Stellingwerf (1978)
    
    This definition makes use of a Fortran routine written by Jan Cuypers and
    Conny Aerts.
    
    Inclusion of linear frequency shift by Pieter Degroote (see Cuypers 1986)
    
    Inclusion of binary orbital motion by Pieter Degroote (see Shibahashi &
    Kurtz 2012). When orbits are added, times must be in days, then asini is
    in AU.
    
    For circular orbits, give only forbit and asini.
    
    @param times: time points
    @type times: numpy array
    @param signal: observations
    @type signal: numpy array
    @param f0: start frequency
    @type f0: float
    @param fn: stop frequency
    @type fn: float
    @param df: step frequency
    @type df: float
    @param Nbin: number of bins (default: 5)
    @type Nbin: int
    @param Ncover: number of covers (default: 1)
    @type Ncover: int
    @param D: linear frequency shift parameter
    @type D: float
    @return: frequencies, theta statistic
    @rtype: array,array
    """
    T = times.ptp()
    n  = len(times)
    
    #-- initialize variables
    xvar     = signal.std()**2.
    xx = (n-1) * xvar
    nf = int((fn-f0) / df + 0.001) + 1
    f1 = np.zeros(nf,'d')
    s1 = np.zeros(nf,'d')
    
    #-- use Fortran subroutine
    #-- Normal PDM
    if D is None and asini is None:
        f1, s1 = pyscargle.justel(signal,times,f0,df,Nbin,Ncover,xvar,xx,f1,s1,n,nf)
    #-- PDM with linear frequency shift
    elif asini is None:
        f1, s1 = pyscargle.justel2(signal,times,f0,df,Nbin,Ncover,xvar,xx,D,f1,s1,n,nf)
    #-- PDM with circular binary orbit
    elif asini is not None and (e is None or e==0):
        f1, s1 = pyscargle.justel3(signal,times,f0,df,Nbin,Ncover,xvar,xx,asini,
                  forbit,f1,s1,n,nf)
    #-- PDM with eccentric binary orbit
    elif e>0:
        forbit = 2*pi*forbit
        ans,bns = np.array([[__ane__(n,e),__bne__(n,e)] for n in range(1,nmax+1)]).T
        ksins = np.sqrt(ans**2*np.cos(omega)**2+bns**2*np.sin(omega)**2)
        thns = np.arctan(bns/ans*np.tan(omega))
        tau = -np.sum(bns*np.sin(omega))
        f1, s1 = pyscargle.justel4(signal,times,f0,df,Nbin,Ncover,xvar,xx,asini,
        forbit,e,omega,ksins,thns,tau,f1,s1,n,nf,nmax)
        
    
    #-- it is possible that the first computed value is a none-variable
    if not s1[0]: s1[0] = 1. 
    
    return f1, s1














@defaults_pergram
@parallel_pergram
@make_parallel
def box(times, signal, f0=None, fn=None, df=None, Nbin=10, qmi=0.005, qma=0.75 ):
    """
    Box-Least-Squares spectrum of Kovacs et al (2002).

    [ see Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369 ]

    This is the slightly modified version of the original BLS routine 
    by considering Edge Effect (EE) as suggested by 
    Peter R. McCullough [ pmcc@stsci.edu ].

    This modification was motivated by considering the cases when 
    the low state (the transit event) happened to be devided between 
    the first and last bins. In these rare cases the original BLS 
    yields lower detection efficiency because of the lower number of 
    data points in the bin(s) covering the low state.

    For further comments/tests see  www.konkoly.hu/staff/kovacs.html
    
    Transit fraction and precision are given by nb,qmi and qma
    
    Remark: output parameter parameter contains:
    [frequency,depth,transit fraction width,fractional start, fraction end]
    
    @param times: observation times
    @type times: numpy 1D array
    @param signal: observations
    @type signal: numpy 1D array
    @param f0: start frequency
    @type f0: float
    @param fn: end frequency
    @type fn: float
    @param df: frequency step
    @type df: float
    @param Nbin: number of bins in the folded time series at any test period
    @type Nbin: integer
    @param qmi: minimum fractional transit length to be tested
    @type qmi: 0<float<qma<1
    @param qma: maximum fractional transit length to be tested
    @type qma: 0<qmi<float<1
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """
    #-- initialize some variables needed in the FORTRAN module
    n = len(times)
    T = times.ptp()
    u = np.zeros(n)
    v = np.zeros(n)
    
    #-- frequency vector and variables
    nf = (fn-f0)/df
    if f0<2./T: f0=2./T
    
    #-- calculate EEBLS spectrum and model parameters
    power,depth,qtran,in1,in2 = eebls.eebls(times,signal,u,v,nf,f0,df,Nbin,qmi,qma,n)
    frequencies = np.linspace(f0,fn,nf)
    
    #-- to return parameters of fit, do this:
    # pars = [max_freq,depth,qtran+(1./float(nb)),(in1-1)/float(nb),in2/float(nb)]
    return frequencies,power





@defaults_pergram
@parallel_pergram
@make_parallel
def kepler(times,signal, f0=None, fn=None, df=None, e0=0., en=0.91, de=0.1,
           errors=None, wexp=2, x00=0.,x0n=359.9):
    """
    Keplerian periodogram of Zucker et al (2010).
    
    @param times: observation times
    @type times: numpy 1D array
    @param signal: observations
    @type signal: numpy 1D array
    @param f0: start frequency
    @type f0: float
    @param fn: end frequency
    @type fn: float
    @param df: frequency step
    @type df: float
    @param e0: start eccentricity
    @type e0: float
    @param en: end eccentricity
    @type en: float
    @param de: eccentricity step
    @type de: float
    @param x00: start x0
    @type x00: float
    @param x0n: end x0
    @type x0n: float
    @return: frequencies, amplitude spectrum
    @rtype: array,array
    """
    T = times.ptp()
    n = len(times)
    if errors is None:
        errors = np.ones(n)
    maxstep = int((fn-f0)/df+1)
    
    #-- initialize parameters
    f1 = np.zeros(maxstep) #-- frequency
    s1 = np.zeros(maxstep) #-- power
    p1 = np.zeros(maxstep) #-- window
    l1 = np.zeros(maxstep) #-- power LS
    s2 = np.zeros(maxstep) #-- power Kepler
    k2 = np.zeros(6) #-- parameters for Kepler orbit
    
    #-- calculate Kepler periodogram
    pyKEP.kepler(times+0,signal+0,errors,f0,fn,df,wexp,e0,en,de,\
          x00,x0n,f1,s1,p1,l1,s2,k2)
    return f1,s2





def Zwavelet(time, signal, freq, position, sigma=10.0):

    """
    Weighted Wavelet Z-transform of Foster (1996)
    
    Computes "Weighted Wavelet Z-Transform"
    which is a type of time-frequency diagram more suitable
    for non-equidistant time series.
    See: G. Foster, 1996, Astronomical Journal, 112, 1709 
    
    The result can be plotted as a colorimage (imshow in pylab).
    
    Mind the max, min order for position. It's often useful to try log(Z)
    and/or different sigmas.
    
    @param time: time points [0..Ntime-1]
    @type time: ndarray
    @param signal: observed data points [0..Ntime-1]
    @type signal: ndarray
    @param freq: frequencies (omega/2pi in Foster's paper) [0..Nfreq-1]
                 the array should not contain 0.0
    @type freq: ndarray
    @param position: time parameter: tau in Foster's paper [0..Npos-1]
    @type position: ndarray
    @param sigma: smoothing parameter in time domain: sigma in Foster's paper
    @type sigma: float
    @return: Z[0..Npos-1, 0..Nfreq-1]: the Z-transform: time-freq diagram
    @rtype: array
    
    """
  
    y = np.zeros(3)
    S = np.zeros([3,3])
    Z = np.zeros([len(position),len(freq)])
  
    for i in range(len(position)):
        
        tau = position[i]
        arg = 2.0*pi*(time-tau)
 
        for j in range(len(freq)):
        
            nu = freq[j]
      
            # Compute statistical weights akin the Morlet wavelet
      
            weight = np.exp(-(time-tau)**2 * (nu / 2./sigma)**2)
            W = np.sum(weight)
      
            # Compute the base functions. A 3rd base function is the constant 1.
      
            cosine = np.cos(arg*nu)
            sine = np.sin(arg*nu)
          
            print arg, nu
            print sine, cosine
          
            # Compute the innerproduct of the base functions
            # phi_0 = 1 (constant), phi_1 = cosine, phi_2 = sine
      
            S[0,0] = 1.0
            S[0,1] = S[1,0] = np.sum(weight * 1.0 * cosine) / W
            S[0,2] = S[2,0] = np.sum(weight * 1.0 * sine) / W
            S[1,1] = np.sum(weight * cosine * cosine) / W
            S[1,2] = S[2,1] = np.sum(weight * cosine * sine) / W
            S[2,2] = np.sum(weight * sine * sine) / W
            print S
            invS = np.linalg.inv(S)
      
            # Determine the best-fit coefficients y_k of the base functions
      
            for k in range(3):
                y[k] =   np.sum(weight * 1.0 * signal) / W * invS[k,0]     \
                       + np.sum(weight * cosine * signal) / W * invS[k,1]  \
                       + np.sum(weight * sine * signal) / W * invS[k,2]      
      
            # Compute the best-fit model
      
            model = y[0] + y[1] * cosine + y[2] * sine
      
            # Compute the weighted variation of the signal and the model functions
      
            Vsignal = np.sum(weight * signal**2) / W -  (np.sum(weight * signal) / W)**2
            Vmodel = np.sum(weight * model**2) / W -  (np.sum(weight * model) / W)**2
      
            # Calculate the weighted Wavelet Z-Transform
      
            Neff = W**2 / np.sum(weight**2)
            Z[i,j] = (Neff - 3) * Vmodel / 2. / (Vsignal - Vmodel)
      
    # That's it!
  
    return Z
    
#}

#{ Pure Python versions

@defaults_pergram
@parallel_pergram
@make_parallel
def pdm_py(time, signal, f0=None, fn=None, df=None, Nbin=10, Ncover=5, D=0.):

    """
    Computes the theta-statistics to do a Phase Dispersion Minimisation.
    See Stellingwerf R.F., 1978, ApJ, 224, 953)
    
    Joris De Ridder
    
    Inclusion of linear frequency shift by Pieter Degroote (see Cuypers 1986)
    
    @param time: time points  [0..Ntime-1]
    @type time: ndarray       
    @param signal: observed data points [0..Ntime-1]
    @type signal: ndarray
    @param f0: start frequency
    @type f0: float
    @param fn: stop frequency
    @type fn: float
    @param df: step frequency
    @type df: float
    @param Nbin: the number of phase bins (with length 1/Nbin)
    @type Nbin: integer
    @param Ncover: the number of covers (i.e. bin shifts)
    @type Ncover: integer
    @param D: linear frequency shift parameter
    @type D: float
    @return: theta-statistic for each given frequency [0..Nfreq-1]
    @rtype: array
    """
    freq = np.arange(f0,fn+df,df)
    
    Ntime = len(time)
    Nfreq = len(freq)
  
    binsize = 1.0 / Nbin
    covershift = 1.0 / (Nbin * Ncover)
  
    theta = np.zeros(Nfreq)
  
    for i in range(Nfreq):
  
        # Compute the phases in [0,1[ for all time points
        phase = np.fmod((time - time[0]) * freq[i] + D/2.*time**2, 1.0)
    
        # Reset the number of (shifted) bins without datapoints
    
        Nempty = 0
    
        # Loop over all Nbin * Ncover (shifted) bins
    
        for k in range(Nbin):
            for n in range(Ncover):
        
                # Determine the left and right boundary of one such bin
                # Note that due to the modulo, right may be < left. So instead
                # of 0-----+++++------1, the bin might be 0++-----------+++1 .
        
                left = np.fmod(k * binsize + n * covershift, 1.0) 
                right = np.fmod((k+1) * binsize + n * covershift, 1.0) 

                # Select all data points in that bin
        
                if (left < right):
                    bindata = np.compress((left <= phase) & (phase < right), signal)
                else:
                    bindata = np.compress(~((right <= phase) & (phase < left)), signal)

                # Compute the contribution of that bin to the theta-statistics  
          
                if (len(bindata) != 0):
                    theta[i] += (len(bindata) - 1) * bindata.var()
                else:
                    Nempty += 1
  
        # Normalize the theta-statistics        

        theta[i] /= Ncover * Ntime - (Ncover * Nbin - Nempty)     
    
    # Normalize the theta-statistics again
  
    theta /= signal.var()  
    
    # That's it!
 
    return freq,theta




def fasper_py(x,y,ofac,hifac, MACC=4):
    """ function fasper
        Given abscissas x (which need not be equally spaced) and ordinates
        y, and given a desired oversampling factor ofac (a typical value
        being 4 or larger). this routine creates an array wk1 with a
        sequence of nout increasing frequencies (not angular frequencies)
        up to hifac times the "average" Nyquist frequency, and creates
        an array wk2 with the values of the Lomb normalized periodogram at
        those frequencies. The arrays x and y are not altered. This
        routine also returns jmax such that wk2(jmax) is the maximum
        element in wk2, and prob, an estimate of the significance of that
        maximum against the hypothesis of random noise. A small value of prob
        indicates that a significant periodic signal is present.
    
    Reference: 
        Press, W. H. & Rybicki, G. B. 1989
        ApJ vol. 338, p. 277-280.
        Fast algorithm for spectral analysis of unevenly sampled data
        (1989ApJ...338..277P)

    Arguments:
        X   : Abscissas array, (e.g. an array of times).
        Y   : Ordinates array, (e.g. corresponding counts).
        Ofac : Oversampling factor.
        Hifac : Hifac * "average" Nyquist frequency = highest frequency
            for which values of the Lomb normalized periodogram will
            be calculated.
        
    Returns:
        Wk1 : An array of Lomb periodogram frequencies.
        Wk2 : An array of corresponding values of the Lomb periodogram.
        Nout : Wk1 & Wk2 dimensions (number of calculated frequencies)
        Jmax : The array index corresponding to the MAX( Wk2 ).
        Prob : False Alarm Probability of the largest Periodogram value
        MACC : Number of interpolation points per 1/4 cycle
                of highest frequency

    History:
        02/23/2009, v1.0, MF
        Translation of IDL code (orig. Numerical recipies)
    """
    #Check dimensions of input arrays
    n = long(len(x))
    if n != len(y):
        print 'Incompatible arrays.'
        return

    nout  = 0.5*ofac*hifac*n
    nfreqt = long(ofac*hifac*n*MACC)   #Size the FFT as next power
    nfreq = 64L             # of 2 above nfreqt.

    while nfreq < nfreqt: 
        nfreq = 2*nfreq

    ndim = long(2*nfreq)
    
    #Compute the mean, variance
    ave = y.mean()
    ##sample variance because the divisor is N-1
    var = ((y-y.mean())**2).sum()/(len(y)-1) 
    # and range of the data.
    xmin = x.min()
    xmax = x.max()
    xdif = xmax-xmin

    #extirpolate the data into the workspaces
    wk1 = np.zeros(ndim, dtype='complex')
    wk2 = np.zeros(ndim, dtype='complex')

    fac  = ndim/(xdif*ofac)
    fndim = ndim
    ck  = ((x-xmin)*fac) % fndim
    ckk  = (2.0*ck) % fndim

    for j in range(0L, n):
        __spread__(y[j]-ave,wk1,ndim,ck[j],MACC)
        __spread__(1.0,wk2,ndim,ckk[j],MACC)

    #Take the Fast Fourier Transforms
    wk1 = np.fft.ifft( wk1 )*len(wk1)
    wk2 = np.fft.ifft( wk2 )*len(wk1)

    wk1 = wk1[1:nout+1]
    wk2 = wk2[1:nout+1]
    rwk1 = wk1.real
    iwk1 = wk1.imag
    rwk2 = wk2.real
    iwk2 = wk2.imag
    
    df  = 1.0/(xdif*ofac)
    
    #Compute the Lomb value for each frequency
    hypo2 = 2.0 * abs( wk2 )
    hc2wt = rwk2/hypo2
    hs2wt = iwk2/hypo2

    cwt  = np.sqrt(0.5+hc2wt)
    swt  = np.sign(hs2wt)*(np.sqrt(0.5-hc2wt))
    den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2
    cterm = (cwt*rwk1+swt*iwk1)**2./den
    sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)

    wk1 = df*(np.arange(nout, dtype='float')+1.)
    wk2 = (cterm+sterm)/(2.0*var)
    pmax = wk2.max()
    jmax = wk2.argmax()


    #Significance estimation
    #expy = exp(-wk2)          
    #effm = 2.0*(nout)/ofac       
    #sig = effm*expy
    #ind = (sig > 0.01).nonzero()
    #sig[ind] = 1.0-(1.0-expy[ind])**effm

    #Estimate significance of largest peak value
    expy = np.exp(-pmax)          
    effm = 2.0*(nout)/ofac       
    prob = effm*expy

    if prob > 0.01: 
        prob = 1.0-(1.0-expy)**effm

    return wk1,wk2,nout,jmax,prob
#}

#{ Helper functions

def windowfunction(time, freq):

    """
    Computes the modulus square of the window function of a set of 
    time points at the given frequencies. The time point need not be 
    equidistant. The normalisation is such that 1.0 is returned at 
    frequency 0.
    
    @param time: time points  [0..Ntime-1]
    @type time: ndarray       
    @param freq: frequency points. Units: inverse unit of 'time' [0..Nfreq-1]
    @type freq: ndarray       
    @return: |W(freq)|^2      [0..Nfreq-1]
    @rtype: array
    
    """
  
    Ntime = len(time)
    Nfreq = len(freq)
    winkernel = np.empty_like(freq)

    for i in range(Nfreq):
        winkernel[i] = np.sum(np.cos(2.0*pi*freq[i]*time))**2     \
                     + np.sum(np.sin(2.0*pi*freq[i]*time))**2

    # Normalise such that winkernel(nu = 0.0) = 1.0 

    return winkernel/Ntime**2


def check_input(times,signal,**kwargs):
    """
    Check the input arguments for periodogram calculations for mistakes.
    
    If you get an error when trying to compute a periodogram, and you don't
    understand it, just feed the input you gave to this function, and it will
    perform some basic checks.
    """
    #-- check if the input are arrays and have the same 1D shape
    is_array0 = isinstance(times,np.ndarray)
    is_array1 = isinstance(signal,np.ndarray)
    if not is_array0: print(termtools.red('ERROR: time input is not an array'))
    if not is_array1: print(termtools.red('ERROR: signal input is not an array'))
    if not is_array0 or not is_array1:
        times = np.asarray(times)
        signal = np.asarray(signal)
        print(termtools.green("---> FIXED: inputs are arrays"))
    print(termtools.green("OK: inputs are arrays"))
    onedim = (len(times.shape)==1) & (len(signal.shape)==1)
    same_shape = times.shape==signal.shape
    if not onedim or not same_shape:
        print(termtools.red('ERROR: input is not 1D or not of same length'))
        return False
    print(termtools.green("OK: inputs are 1D and have same length"))
    #-- check if the signal constains nans or infs:
    isnan0 = np.sum(np.isnan(times))
    isnan1 = np.sum(np.isnan(signal))
    isinf0 = np.sum(np.isinf(times))
    isinf1 = np.sum(np.isinf(signal))
    if isnan0: print(termtools.red('ERROR: time array contains nans'))
    if isnan1: print(termtools.red('ERROR: signal array contains nans'))
    if isinf0: print(termtools.red('ERROR: time array contains infs'))
    if isinf1: print(termtools.red('ERROR: signal array contains infs'))
    if not isnan0 and not isnan1 and not isinf0 and not isinf1:
        print(termtools.green('OK: no infs or nans'))
    else:
        keep = -np.isnan(times) & -np.isnan(signal) & -np.isinf(times) & -np.isinf(signal)
        times,signal = times[keep],signal[keep]
        print(termtools.green('---> FIXED: infs and nans removed'))
    #-- check if the timeseries is sorted
    is_sorted = np.all(np.diff(times)>0)
    if not is_sorted:
        print(termtools.red('ERROR: time array is not sorted'))
        sa = np.argsort(times)
        times,signal = times[sa],signal[sa]
        print(termtools.green('---> FIXED: time array is sorted'))
    else:
        print(termtools.green("OK: time array is sorted"))
    print(termtools.green("No inconsistencies found or inconsistencies are fixed"))
    
    #-- check keyword arguments:
    fnyq = getNyquist(times,nyq_stat=np.min)
    print("Default Nyquist frequency: {}".format(fnyq))
    if 'nyq_stat' in kwargs:
        fnyq = getNyquist(times,nyq_stat=kwargs['nyq_stat'])
        print("Nyquist value manually set to {}".format(fnyq))
    if 'fn' in kwargs and kwargs['fn']>fnyq:
        print(termtools.red("Final frequency 'fn' is larger than the Nyquist frequency"))
    return times,signal


def __spread__(y, yy, n, x, m):
    """
    Given an array yy(0:n-1), extirpolate (spread) a value y into
    m actual array elements that best approximate the "fictional"
    (i.e., possible noninteger) array element number x. The weights
    used are coefficients of the Lagrange interpolating polynomial
    Arguments:
        y : 
        yy : 
        n : 
        x : 
        m : 
    Returns:
        
    """
    nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]
    if m > 10. :
        print 'factorial table too small in spread'
        return

    ix=long(x)
    if x == float(ix): 
        yy[ix]=yy[ix]+y
    else:
        ilo = long(x-0.5*float(m)+1.0)
        ilo = min( max( ilo , 1 ), n-m+1 ) 
        ihi = ilo+m-1
        nden = nfac[m]
        fac=x-ilo
        for j in range(ilo+1,ihi+1): fac = fac*(x-j)
        yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))
        for j in range(ihi-1,ilo-1,-1):
            nden=(nden/(j+1-ilo))*(j-ihi)
            yy[j] = yy[j] + y*fac/(nden*(x-j))    

def __ane__(n,e):
    return 2.*np.sqrt(1-e**2)/e/n*jn(n,n*e)
    
def __bne__(n,e):
    return 1./n*(jn(n-1,n*e)-jn(n+1,n*e))

    
def getSignificance(wk1, wk2, nout, ofac):
    """
    Returns the peak false alarm probabilities
    
    Hence the lower is the probability and the more significant is the peak
    """
    expy = np.exp(-wk2)          
    effm = 2.0*(nout)/ofac       
    sig = effm*expy
    ind = (np.sig > 0.01).nonzero()
    sig[ind] = 1.0-(1.0-expy[ind])**effm
    return sig
  

def scargle_probability(peak_value,times,freqs,correct_for_frange=False,**kwargs):
    """
    Compute the probability to observe a peak in the Scargle periodogram.
    
    If C{correct_for_frange=True}, the Bonferroni correction will be applied
    to a smaller number of frequencies (i.e. the independent number of
    frequencies in C{freqs}). To be conservative, set C{correct_for_frange=False}.
    
    Example simulation:
    
    >>> times = np.linspace(0,1,5000)
    >>> N = 500
    >>> probs = np.zeros(N)
    >>> peaks = np.zeros(N)
    >>> for i in range(N):
    ...     signal = np.random.normal(size=len(times))
    ...     f,s = scargle(times,signal,threads='max',norm='distribution')
    ...     peaks[i] = s.max()
    ...     probs[i] = scargle_probability(s.max(),times,f)
    
    Now make a plot:
    
    >>> p = pl.figure()
    >>> p = pl.subplot(131)
    >>> p = pl.plot(probs,'ko')
    >>> p = pl.plot([0,N],[0.01,0.01],'r-',lw=2)
    >>> p = pl.subplot(132)
    >>> p = pl.plot(peaks[np.argsort(peaks)],probs[np.argsort(peaks)],'ro')
    >>> p = pl.plot(peaks[np.argsort(peaks)],1-(1-np.exp(-np.sort(peaks)))**(10000.),'g-')
    >>> #p = pl.plot(peaks[np.argsort(peaks)],1-(1-(1-np.sort(peaks)/2500.)**2500.)**(10000.),'b--')
    >>> p = pl.subplot(133)
    >>> for i in np.logspace(-3,0,100):
    ...     p = pl.plot([i*100],[np.sum(probs<i)/float(N)*100],'ko')
    >>> p = pl.plot([1e-6,100],[1e-6,100],'r-',lw=2)
    >>> p = pl.xlabel('Should observe this many points below threshold')
    >>> p = pl.ylabel('Observed this many points below threshold')
    
    ]]include figure]]ivs_timeseries_pergrams_prob.png]
    
    """
    #-- independent frequencies
    nr_obs = len(times)
    ni = 2*nr_obs
    #-- correct the nr of independent frequencies for the frequency range
    #   that is tested, but only if it is requested
    if correct_for_frange:
        nyqstat = kwargs.pop('nyqstat',np.min)
        nyquist = getNyquist(times,nyqstat=nyqstat)
        ni = int(freqs.ptp()/nyquist*ni)
    #p_value = 1. - (1.- (1-2*peak_value/nr_obs)**(nr_obs/2))**ni
    p_value = 1. - (1.- np.exp(-peak_value))**ni
    return p_value


if __name__=="__main__":
    import sys
    import pylab as pl
    from ivs.aux import loggers
    logger = loggers.get_basic_logger()
    
    #-- run tests
    if '--test' in sys.argv[1] or '-t' in sys.argv[1]:
        import doctest
        doctest.testmod()
        pl.show()
    #-- command line interface    
    else:
    
        method,args,kwargs = argkwargparser.parse()
        print "Running method %s with arguments %s and keyword arguments %s"%(method,args,kwargs)
        if '--help' in args or 'help' in args or 'help' in kwargs:
            sys.exit()
        times,signal = ascii.read2array(kwargs.pop('infile')).T[:2]
        freqs,ampls = globals()[method](times,signal, **kwargs)
        pl.plot(freqs,ampls,'k-')
        pl.show()