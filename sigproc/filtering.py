"""
Filter non-equidistant signals (spectra, timeseries...).

Currently implemented:
    - Gaussian
    - Sinc
    - Box
    - Iterative Nonlinear filter
    
Example usage:

Generate some data:

>>> x = np.linspace(0,150,10000)
>>> x_= np.linspace(0,150,2000)
>>> y = np.sin(2*pi/20.*x)+np.random.normal(size=len(x))

Apply some filter:

>>> x1,y1,pnts = filter_signal(x,y,"gauss",sigma=1)
>>> x2,y2,pnts = filter_signal(x,y,"pijpers",delta=1)
>>> x3,y3,pnts = filter_signal(x,y,"box",window_width=10)
>>> x3_,y3_,pnts = filter_signal(x,y,"box",window_width=10,x_template=x_,threads=2)

>>> import pylab as pl
>>> p = pl.plot(x,y,'ko')
>>> p = pl.plot(x1,y1,'ro',mec='r')
>>> p = pl.plot(x2,y2,'bo',mec='b')
>>> p = pl.plot(x3,y3,'go',mec='g')
>>> p = pl.plot(x3_,y3_,'gs')

"""
import copy
import logging
import numpy as np
from numpy import sqrt,exp,pi,cos,sinc,trapz,average

from scipy.integrate import trapz
from multiprocessing import Process, Manager,cpu_count

from ivs.aux.decorators import make_parallel
from ivs.timeseries.decorators import parallel_pergram,defaults_filtering

logger = logging.getLogger("TS.FILTERING")

#{ Main filter program for convolution-based filters

@defaults_filtering
@parallel_pergram
@make_parallel
def filter_signal(x,y,ftype,f0=None,fn=None,step=1,x_template=None,**kwargs):
    """
    Filter a signal.
    
    Example usage:
    
    #Generate some data:
        #>>> import time
        #>>> x = linspace(0,150,10000)
        #>>> y = sin(2*pi/20.*x)+random.normal(size=len(x))
    
    #Apply some filter:
        #>>> c0 = time.time()
        #>>> y1,pnts = filter_signal(x,y,"gauss",sigma=1)
        #>>> print time.time()-c0
        #>>> c0 = time.time()
        #>>> y1,pnts = filter_signal(x,y,"gauss",sigma=1,threads=2)
        #>>> print time.time()-c0

    @param ftype: one of 'gauss','pijpers','box','inl'
    @type ftype: string
    @rtype: tuple
    @return: output from the used filter
    """
    if x_template is None:
        x_template = x + 0.
    if f0 is None: f0 = 0
    if fn is None: fn = len(x_template)
    
    #-- set window and kernel function
    ftype = ftype.lower()    
    window = globals()[ftype+'_window']
    kernel = globals()[ftype+'_kernel']
    
    logger.info("Applying filter: %s"%(ftype))
    
    #-- set window width
    lower_window,higher_window = window(**kwargs)
    
    #-- start running through timeseries (make searchsorted local for speedup)
    logger.debug("FILTER between index %d-%d with step %d"%(f0,fn,step))
    x_template = x_template[f0:fn:step]
    
    searchsorted = x.searchsorted
    out = [kernel(x,y,t,index=searchsorted([t-1e-30-lower_window,t+1e-30+higher_window]),
                        **kwargs) for t in x_template]
    #-- that's it!
    return tuple([x_template] + list(np.array(out).T))    
    
#}

#{ Basic functions for convolution-based filters
    
def _fixed_window(window_width=0):
    """
    Define general window
    @rtype: (float,float)
    @return: window limits
    """
    lower_window = window_width/2.
    higher_window = window_width/2.
    logger.debug("Window: width=%f"%(window_width))
    return lower_window,higher_window
#}
#{ Convolution-based filter kernels and Windows
    
def gauss_window(sigma=1.,limit=4.):
    """
    Define Gaussian_window
    @rtype: (float,float)
    @return: window limits
    """
    return _fixed_window(window_width=2*sigma*limit)
    

def gauss_kernel(times, signal, t,sigma=1.,index=(0,-1),norm_weights=True):
    """
    Define Gaussian kernel.
    
    #Import necessary modules:
        #>>> from pylab import plot,figure,title,subplot
        #>>> from sigproc.timeseries import tfreq
    
    #Generate some data:
        #>>> x = linspace(0,150,10000)
        #>>> y = random.normal(size=len(x))
    
    #Apply Gaussian filter twice, subtract and calculate periodogram:
        #>>> y1,pnts = filter_signal(x,y,"gauss",sigma=1)
        #>>> y2,pnts = filter_signal(x,y,"gauss",sigma=3)
    
    #Plot the results:
        #>>> p=figure(figsize=(8,18))
        #>>> p=subplot(311);p=title("GAUSSIAN FILTER")
        #>>> p=plot(x,y,'ko')
        #>>> p=plot(x,y1,'r-',linewidth=2)
        #>>> p=plot(x,y2,'b-',linewidth=2)
        #>>> p=subplot(312)
        #>>> p=plot(x,y1-y2,'ko')
        #>>> p=subplot(313)
        #>>> p=plot(per1[0],per1[1],'k-')
        #>>> p=plot(per2[0],per2[1],'r-')
    
    @rtype: (ndarray,)
    """
    #-- compute convolution kernel -- should we normalize? so that sum(weights)=1
    index0,indexn = index
    times_ = times[index0:indexn]
    signal_ = signal[index0:indexn]
    weights = 1./(sqrt(2.*pi)*sigma) * exp( -(times_-t)**2./(2.*sigma**2.))
    if norm_weights and len(times)>1:
        convolved = trapz(weights*signal_,x=times_-t)
        norm_fact = trapz(weights,x=times_-t)
    else:
        convolved = average(signal_,weights=weights)
        norm_fact = 1
    convolved_signal = convolved/norm_fact
    return convolved_signal,len(times_)

def mad(a, c=0.6745):
    """
    Median Absolute Deviation of a 1D array:
    
    median(abs(a - median(a))) / c
    
    @rtype: float
    """
    return np.median(np.abs(a-np.median(a)))/c

def inl_window(**kwargs):
    """
    Defines the INL window
    
    @rtype: (float,float)
    @return: window limits
    """
    return _fixed_window(**kwargs)

def inl_kernel(times_in,signal_in,t,c=0.6745,sig_level=3.,tolerance=0.01,index=(0,-1)):
    """
    Iterative Nonlinear Filter (Aigrain).
    
    Example usage:
    
    #Import necessary modules:
        #>>> from pylab import plot,figure,title,subplot
        #>>> from sigproc.timeseries import tfreq
    
    #Generate some data:
        #>>> x = linspace(0,150,10000)
        #>>> y = sin(2*pi*x/20.*x)+random.normal(size=len(x))
    
    #Apply INL filter twice, subtract and calculate periodogram:
        #>>> y1,pnts1,mads = filter_signal(x,y,"inl",window_width=1)
        #>>> y2,pnts2,mads = filter_signal(x,y,"inl",window_width=10)
        #>>> freq1,per1,stat1 = tfreq.scargle(x,y,norm='amplitude',fn=4)
        #>>> freq2,per2,stat2 = tfreq.scargle(x,y1-y2,norm='amplitude',fn=4)
    
    #Plot the results:
        #>>> p=figure(figsize=(8,18))
        #>>> p=subplot(311);p=title("INL FILTER")
        #>>> p=plot(x,y,'ko')
        #>>> p=plot(x,y1,'r-',linewidth=2)
        #>>> p=plot(x,y2,'b-',linewidth=2)
        #>>> p=subplot(312)
        #>>> p=plot(x,y1-y2,'ko')
        #>>> p=subplot(313)
        #>>> p=plot(per1[0],per1[1],'k-')
        #>>> p=plot(per2[0],per2[1],'r-')
    
    @rtype: (ndarray,ndarray)
    @return: continuum,sigma
    """
    index0,indexn = index
    times = times_in[index0:indexn]
    signal = signal_in[index0:indexn]
    sigma = mad(signal,c=c)
    continuum = np.median(signal)
    outliers = (np.abs(continuum-signal)>sig_level*sigma)
    while 1:
        signal_ = np.compress(outliers==0,signal)
        sigma = mad(signal_,c=c)
        continuum_ = np.median(signal_)
        if abs((continuum_-continuum)/continuum_)<tolerance:
            break
        else:
            continuum = continuum_
            outliers = (np.abs(continuum - signal_)>sig_level*sigma)
        
    return continuum_, sigma,len(times)

def pijpers_window(delta=1.):
    """
    Defines the window for the Pijpers filter.
    
    @rtype: (float,float)
    @return: window limits
    """
    return _fixed_window(window_width=100*delta)
    return _fixed_window(**kwargs)
    

def pijpers_kernel(times_in, signal_in,t,limit=0.000001,delta=1.,sigma=1.,gamma=0,
                         r=0.0001,index0=0,indexn=-1,index=(0,-1),norm_weights=True):
    """
    Defines the Pijpers (2006) filter kernel.
    
    Equals box-multiplying in the frequency domain.
    
    Example usage:
    
    #Import necessary modules:
        #>>> from pylab import plot,figure,title,subplot
        #>>> from sigproc.timeseries import tfreq
    
    #Generate some data:
        #>>> x = linspace(0,150,10000)
        #>>> y = random.normal(size=len(x))
    
    #Apply Gaussian filter twice, subtract and calculate periodogram:
        #>>> y1,pnts1 = filter_signal(x,y,"pijpers",delta=1)
        #>>> y2,pnts2 = filter_signal(x,y,"pijpers",delta=3)
        #>>> freq1,per1,stat1 = tfreq.scargle(x,y,norm='amplitude',fn=4)
        #>>> freq2,per2,stat2 = tfreq.scargle(x,y1-y2,norm='amplitude',fn=4)
    
    #Plot the results:
        #>>> p=figure(figsize=(8,18))
        #>>> p=subplot(311);p=title("PIJPERS FILTER")
        #>>> p=plot(x,y,'ko')
        #>>> p=plot(x,y1,'r-',linewidth=2)
        #>>> p=plot(x,y2,'b-',linewidth=2)
        #>>> p=subplot(312)
        #>>> p=plot(x,y1-y2,'ko')
        #>>> p=subplot(313)
        #>>> p=plot(per1[0],per1[1],'k-')
        #>>> p=plot(per2[0],per2[1],'r-')
    
    @rtype: (ndarray,)
    @return: convolved signal
    """
    index0,indexn = index
    times = times_in[index0:indexn]
    signal = signal_in[index0:indexn]
    
    delta = delta*2*pi
    
    x = delta * (times-t)
    #-- compute convolution kernel -- should we normalize? so that sum(weights)=1
    weights = delta/pi * sinc(x/pi) * exp(-1/4. * r**2 * x**2) * cos(gamma*x)
    #-- compensate for severe unequidistant sampling
    convolved = np.trapz(weights*signal,x=times-t)
    norm_fact = np.trapz(weights,x=times-t)
    convolved_signal = convolved/norm_fact
    return convolved_signal,len(times)

def box_window(**kwargs):
    """
    Defines the window for the box-filter.
    
    @rtype: (float,float)
    @return: window limits
    """
    return _fixed_window(**kwargs)

def box_kernel(times_in,signal_in,t,window_width=None,index=(0,-1),norm_weights=True):
    """
    Defines the Box filter kernel (moving average).
    
    Equals sinc-multiplication of frequency domain.
    
    Example usage:

    #Import necessary modules:
        #>>> from pylab import plot,figure,title,subplot
        #>>> from sigproc.timeseries import tfreq
    
    #Generate some data:
        #>>> x = linspace(0,150,10000)
        #>>> y = random.normal(size=len(x))
    
    #Apply Gaussian filter twice, subtract and calculate periodogram:
        #>>> y1,pnts = filter_signal(x,y,"box",window_width=1)
        #>>> y2,pnts = filter_signal(x,y,"box",window_width=3)
    
    @rtype: (ndarray,)
    @return: convolved signal
    """
    index0,indexn = index
    times = times_in[index0:indexn]
    signal = signal_in[index0:indexn]
    weights = np.ones(len(times))
    if norm_weights: weights /= np.sum(weights)
    
    convolved_signal = np.sum(weights*signal)
    return convolved_signal,len(times)

#}

def test():
    """
    
    >>> from pylab import show
    >>> p=show()
    
    """
    import doctest
    doctest.testmod()


if __name__ == "__main__":
    test()