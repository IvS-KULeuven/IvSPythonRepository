# -*- coding: utf-8 -*-
"""
Various decorator functions for time series analysis
    - Parallel periodogram
    - Autocompletion of default arguments
"""
import functools
import logging
from multiprocessing import Manager,Process,cpu_count
import numpy as np
from ivs.aux import loggers

logger = logging.getLogger("TS.DEC")
logger.addHandler(loggers.NullHandler)

def parallel_pergram(fctn):
    """
    Run periodogram calculations in parallel.

    This splits up the frequency range between f0 and fn in 'threads' parts.

    This must decorate a 'make_parallel' decorator.
    """
    @functools.wraps(fctn)
    def globpar(*args,**kwargs):
        #-- construct a manager to collect all calculations
        manager = Manager()
        arr = manager.list([])
        all_processes = []
        #-- get information on frequency range
        f0 = kwargs['f0']
        fn = kwargs['fn']
        threads = kwargs.pop('threads',1)
        if threads=='max':
            threads = cpu_count()
        elif threads=='safe':
            threads = cpu_count()-1
        else:
            threads = float(threads)

        #-- extend the arguments to include the parallel array
        myargs = tuple(list(args) + [arr] )
        #-- however, some functions cannot be parallelized
        if fctn.__name__ in ['fasper']:
            threads = 1

        #-- distribute the periodogram calcs over different threads, and wait
        for i in range(int(threads)):
            #-- define new start and end frequencies
            kwargs['f0'] = f0 + i*(fn-f0) / threads
            kwargs['fn'] = f0 +(i+1)*(fn-f0) / threads
            logger.debug("parallel: starting process %s: f=%.4f-%.4f"%(i,kwargs['f0'],kwargs['fn']))
            p = Process(target=fctn, args=myargs, kwargs=kwargs)
            p.start()
            all_processes.append(p)

        for p in all_processes: p.join()

        logger.debug("parallel: all processes ended")

        #-- join all periodogram pieces
        freq = np.hstack([output[0] for output in arr])
        ampl = np.hstack([output[1] for output in arr])
        sort_arr = np.argsort(freq)
        ampl = ampl[sort_arr]
        freq = freq[sort_arr]
        ampl[np.isnan(ampl)] = 0.

        if len(arr[0])>2:
            rest = []
            for i in range(2,len(arr[0])):
                rest.append(np.hstack([output[i] for output in arr]))
            rest = np.array(rest).T
            rest = rest[sort_arr].T
            return tuple([freq,ampl]+list(rest))
        else:
            return freq,ampl

    return globpar



def defaults_pergram(fctn):
    """
    Set default parameters common to all periodograms.
    """
    @functools.wraps(fctn)
    def globpar(*args,**kwargs):
        #-- this is the information we need the compute everything
        times = args[0]
        signal = args[1]
        T = times.ptp()

        #-- get information on frequency range. If it is not given, compute the
        #   start (0.1/T) and stop (Nyquist) frequency.
        #   Also compute the frequency step as 0.1/T
        nyq_stat = kwargs.pop('nyq_stat',np.min)
        nyquist = getNyquist(times,nyq_stat=nyq_stat)
        f0 = kwargs.get('f0',0.01/T)
        fn = kwargs.get('fn',nyquist)
        df = kwargs.get('df',0.1/T)
        if df==0: df = 0.1/T
        if f0==0: f0 = 0.01/T
        if fn>nyquist:
            fn = nyquist
        kwargs['f0'] = f0
        kwargs['df'] = df
        kwargs['fn'] = fn

        #-- maybe the data needs to be windowed
        window = kwargs.pop('window',None)
        if window is not None:
            signal = signal*windowfunctions.getWindowFunction(window)(times)
            signal -= signal.mean()
            logger.debug('Signal is windowed with %s'%(window))

        #-- normalise weights if they are given
        weights = kwargs.get('weights',None)
        if weights is not None:
            if weights.sum() != len(weights):
                weights = weights / float(weights.sum()) * len(weights)
                logger.debug("Weights were initially not normalized: normalization performed.")
                kwargs['weights'] = weights
        return fctn(times,signal,*args[2:],**kwargs)

    return globpar


def defaults_filtering(fctn):
    """
    Set default parameters common to all filtering functions.
    """
    @functools.wraps(fctn)
    def globpar(*args,**kwargs):
        #-- this is the information we need the compute everything
        #   it is possible a template x-axis has been given. If so, the parallelization
        #   depends on the length of the template, not of the original thing.
        fn = ('x_template' in kwargs) and len(kwargs['x_template']) or len(args[0])
        f0 = kwargs.get('f0',0)
        fn = kwargs.get('fn',fn)
        kwargs['f0'] = f0
        kwargs['fn'] = fn
        return fctn(*args,**kwargs)
    return globpar





def getNyquist(times,nyq_stat=np.inf):
    """
    Calculate Nyquist frequency.

    Typical use is minimum or median of time points differences.

    If C{nyq_stat} is not callable, it is assumed to be a number and that number
    will just be returned: this you can do to search for frequencies above the
    nyquist frequency

    @param times: sorted array containing time points
    @type times: numpy array
    @param nyq_stat: statistic to use or absolute value of the Nyquist frequency
    @type nyq_stat: callable or float
    @return: Nyquist frequency
    @rtype: float
    """
    if not hasattr(nyq_stat,'__call__'):
        nyquist = nyq_stat
    else:
        nyquist = 1/(2.*nyq_stat(np.diff(times)))
    return nyquist
