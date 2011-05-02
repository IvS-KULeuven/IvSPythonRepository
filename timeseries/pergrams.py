# -*- coding: utf-8 -*-
"""
Contains many different periodogram calculations

The basic interface is

>>> freq,ampl = scargle(times,signal)

If you give no extra information, default values for the start, end and step
frequency will be chosen. See the 'defaults_pergram' decorator for a list of
keywords common to all periodogram calculations.

All periodograms can be computed in parallel, by supplying an extra keyword
'threads'.
"""
import numpy as np
from ivs.misc.decorators import make_parallel
from ivs.timeseries.decorators import parallel_pergram,defaults_pergram

import pyscargle
import pyclean
import pyGLS
import pyKEP
import multih
import deeming as fdeeming
import eebls

@defaults_pergram
@parallel_pergram
@make_parallel
def scargle(times, signal, f0=None, fn=None, df=None, norm='amplitude', weights=None):
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
    the variance of the data (NOT of the noise or residuals!).
    
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
    #-- initialize variables for use in Fortran routine
    sigma=0.;xgem=0.;xvar=0.;n=len(times)
    T = times.ptp()
    nf=int((fn-f0)/df+0.001)+1
    f1=np.zeros(nf,'d');s1=np.zeros(nf,'d')
    ss=np.zeros(nf,'d');sc=np.zeros(nf,'d');ss2=np.zeros(nf,'d');sc2=np.zeros(nf,'d')
    
    #-- run the Fortran routine
    if weights is None:
        f1,s1=pyscargle.scar2(signal,times,f0,df,f1,s1,ss,sc,ss2,sc2)
    else:
        w=np.array(weights,'float')
        logger.debug('Weighed scargle')
        f1,s1=pyscargle.scar3(signal,times,f0,df,f1,s1,ss,sc,ss2,sc2,w)
    
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
def deeming(times,signal, f0=None, fn=None, df=None, norm='amplitude'):
    """
    Deeming periodogram of Deeming et al. (1975).
    
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
    Generalised Least Squares periodogram of Zucher et al. (2010).
    
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
    Cleaned Fourier periodogram of Roberts et al. (1987)
        
    Parallization probably isn't such a good idea here because of the frequency
    bins.
    
    Fortran module probably from John Telting.
    
    Should always start from zero, so f0 is not an option
        >>> times_ = np.linspace(0,150,1000)
        >>> times = np.array([times_[i] for i in xrange(len(times_)) if (i%10)>7])
        >>> signal = np.sin(2*pi/10*times) + np.random.normal(size=len(times))
        >>> niter,freqbins = 10,[0,1.2]
        >>> p1 = scargle(times,signal,fn=1.2,norm='amplitude')
        >>> p2 = clean(times,signal,fn=1.2,gain=1.0,niter=niter,freqbins=freqbins)
        >>> p3 = clean(times,signal,fn=1.2,gain=0.1,niter=niter,freqbins=freqbins)
        >>> from pylab import figure,plot,legend
        >>> p=figure()
        >>> p=plot(p1[0],p1[1],'k-',label="Scargle")
        >>> p=plot(p2[0],p2[1],'r-',label="Clean (g=1.0)")
        >>> p=plot(p3[0],p3[1],'b-',label="Clean (g=0.1)")
        >>> p=legend()
    
    @keyword freqbins: frequency bins for clean computation
    @type freqbins: list or array
    @keyword niter: number of iterations
    @type niter: integer
    @keyword gain: gain for clean computation
    @type gain: float between 0 (no cleaning) and 1 (full cleaning)
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
    @return: frequencies, f-statistics.
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









@defaults_pergram
@parallel_pergram
@make_parallel
def pdm(times, signal, f0=None, fn=None, df=None, nb=5, nc=2, D=0):
    """
    Phase Dispersion Minimization of Jurkevich-Stellingwerf (1978)
    
    This definition makes use of a Fortran routine written by Jan Cuypers and
    Conny Aerts.
    
    Inclusion of linear frequency shift by Pieter Degroote
    
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
    @param nb: number of bins (default: 5)
    @type nb: int
    @param nc: number of covers (default: 1)
    @type nc: int
    @param D: linear frequency shift parameter
    @type D: float
    @return: frequencies, theta statistic
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
    if D is None:
        f1, s1 = pyscargle.justel(signal,times,f0,df,nb,nc,xvar,xx,f1,s1,n,nf)
    else:
        f1, s1 = pyscargle.justel2(signal,times,f0,df,nb,nc,xvar,xx,D,f1,s1,n,nf)
    
    #-- it is possible that the first computed value is a none-variable
    if not s1[0]: s1[0] = 1. 
    
    return f1, s1







@defaults_pergram
@parallel_pergram
@make_parallel
def bls(times, signal, f0=None, fn=None, df=None, nb=50, qmi=0.005, qma=0.75 ):
    """
    Box-Least-Squares spectrum of Kovacs et al. (2002).

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
    @param nb: number of bins in the folded time series at any test period
    @type nb: integer
    @param qmi: minimum fractional transit length to be tested
    @type qmi: 0<float<qma<1
    @param qma: maximum fractional transit length to be tested
    @type qma: 0<qmi<float<1
    @return: [frequency,power], [frequencies, power spectrum], parameters
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
    power,depth,qtran,in1,in2 = eebls.eebls(times,signal,u,v,nf,f0,df,nb,qmi,qma,n)
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
    Keplerian periodogram of Zucker et al. (2010).
    
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




if __name__=="__main__":
    import pylab as pl
    from ivs.misc import loggers
    logger = loggers.get_basic_logger()
    
    x = np.linspace(0,100,100)
    y = np.sin(2*np.pi/10.*x) + np.random.normal(size=len(x),scale=0.2)
    for i,norm in enumerate(['power','amplitude','distribution','density']):
        f1,s1 = scargle(x,y,norm=norm)
        f2,s2 = deeming(x,y,norm=norm)
        pl.subplot(2,2,i+1)
        pl.plot(f1,s1,lw=2)
        pl.plot(f2,s2)
    
    pl.figure()
    f1,s1 = gls(x,y)
    f2,s2 = clean(x,y)
    f3,s3 = schwarzenberg_czerny(x,y,nh=2)
    f4,s4 = pdm(x,y)
    f5,s5 = bls(x,y)
    f6,s6 = kepler(x,y)
    
    pl.subplot(2,3,1)
    pl.plot(f1,s1)
    pl.subplot(2,3,2)
    pl.plot(f2,s2)
    pl.subplot(2,3,3)
    pl.plot(f3,s3)
    pl.subplot(2,3,4)
    pl.plot(f4,s4)
    pl.subplot(2,3,5)
    pl.plot(f5,s5)
    pl.subplot(2,3,6)
    pl.plot(f6,s6)
    
    
    pl.show()