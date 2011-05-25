# -*- coding: utf-8 -*-
"""
Decorators specifically for SED package
"""
import functools
import logging
import numpy as np
from multiprocessing import Manager,Process,cpu_count

logger = logging.getLogger('SED.DEC')

def parallel_gridsearch(fctn):
    """
    Decorator to run SED grid fitting in parallel.
    
    This splits up the effective temperature range between teffrange[0] and
    teffrange[1] in 'threads' parts.
    
    This must decorate a 'make_parallel' decorator.
    """
    @functools.wraps(fctn)
    def globpar(*args,**kwargs):
        #-- construct a manager to collect all calculations 
        manager = Manager() 
        arr = manager.list([]) 
        all_processes = []
        #-- get information on threading
        threads = kwargs.pop('threads',1)
        if threads=='max':
            threads = cpu_count()
        elif threads=='half':
            threads = cpu_count()/2
        elif threads=='safe':
            threads = cpu_count()-1
        threads = int(threads)
        index = np.arange(len(args[0]))
        
        #-- distribute the periodogram calcs over different threads, and wait
        for i in range(threads):
            #-- extend the arguments to include the parallel array, and split
            #   up the first three input arrays
            myargs = tuple([args[0][i::threads],args[1][i::threads],args[2][i::threads],args[3][i::threads]] + list(args[4:]) + [arr] )
            kwargs['index'] = index[i::threads]
            logger.debug("parallel: starting process %s"%(i))
            p = Process(target=fctn, args=myargs, kwargs=kwargs) 
            p.start()
            all_processes.append(p)
        
        for p in all_processes: p.join() 
        
        logger.debug("parallel: all processes ended") 
        
        #-- join all periodogram pieces
        chisqs = np.hstack([output[0] for output in arr]) 
        scales = np.hstack([output[1] for output in arr]) 
        e_scales = np.hstack([output[2] for output in arr])
        lumis = np.hstack([output[3] for output in arr])
        index = np.hstack([output[4] for output in arr])
        sa = np.argsort(index)
        return chisqs[sa],scales[sa],e_scales[sa],lumis[sa],index[sa]
        
    return globpar

def iterate_gridsearch(fctn):
    """
    Decorator to run SED iteratively and zooming in on the minimum.
    
    iterations: number of iterative zoom-ins
    increase: increase in number of grid points in each search (1 means no increase)
    size: speed of zoomin: the higher, the slower
    """
    @functools.wraps(fctn)
    def globpar(*args,**kwargs):
        iterations = kwargs.pop('iterations',1)
        increase = kwargs.pop('increase',1)
        speed = kwargs.pop('speed',2)
        
        N = 0
        for nr_iter in range(iterations):
            data_ = fctn(*args,**kwargs)
            #-- append results to record array
            if nr_iter == 0:
                data = data_
                startN = len(data)
            else:
                data = np.core.records.fromrecords(data.tolist()+data_.tolist(),dtype=data.dtype)
            
            #-- select next stage
            best = np.argmin(data['chisq'])
            limit = data['chisq'][best]+speed*0.5**nr_iter*data['chisq'][best]
            
            kwargs['teffrange'] = (data['teff'][data['chisq']<=limit]).min(),(data['teff'][data['chisq']<=limit]).max()
            kwargs['loggrange'] = (data['logg'][data['chisq']<=limit]).min(),(data['logg'][data['chisq']<=limit]).max()
            kwargs['ebvrange'] = (data['ebv'][data['chisq']<=limit]).min(),(data['ebv'][data['chisq']<=limit]).max()
            kwargs['zrange'] = (data['z'][data['chisq']<=limit]).min(),(data['z'][data['chisq']<=limit]).max()
            kwargs['points'] = increase**(nr_iter+1)*startN
            
            logger.info('Best parameters (stage %d/%d): teff=%.0f logg=%.3f E(B-V)=%.3f Z=%.2f (CHI2=%g, cutoff=%g)'\
                     %(nr_iter+1,iterations,data['teff'][best],data['logg'][best],\
                       data['ebv'][best],data['z'][best],data['chisq'][best],limit))
        
        return data
        
    return globpar
