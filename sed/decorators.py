# -*- coding: utf-8 -*-
"""
Decorators specifically for SED package
"""
import functools
import logging
import numpy as np
import pylab as pl
from multiprocessing import Manager,Process,cpu_count
from . import model
from ivs.units import conversions
from ivs.units import constants

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
        index = np.arange(len(args[-1]))

        #-- distribute the periodogram calcs over different threads, and wait
        for i in range(threads):
            #-- extend the arguments to include the parallel array, and split
            #   up the first four input arrays
            #myargs = tuple([args[0][i::threads],args[1][i::threads],args[2][i::threads],args[3][i::threads]] + list(args[4:]) + [arr] )
            myargs = tuple(list(args[:3]) + [args[j][i::threads] for j in range(3,len(args))] +  [arr] )
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
        return chisqs[sa],scales[sa],e_scales[sa],lumis[sa]#,index[sa]

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



def standalone_figure(fctn):
    """
    Accept 'savefig' as an extra keyword. If it is given, start a new figure and
    save it to the filename given, and close it!
    """
    @functools.wraps(fctn)
    def dofig(*args,**kwargs):
        savefig = kwargs.pop('savefig',None)
        colors = kwargs.get('colors',False)
        #-- if savefig is not a string but a boolean, make name: it consists
        #   of the ID, the model used and the function name
        if savefig is True:
            savefig = args[0].ID + '__%s_'%(args[0].label) + model.defaults2str()
            for char in ['/','*']:
                savefig = savefig.replace(char,'_')
            savefig = savefig.replace('.','p')
            savefig = savefig + '_' + fctn.__name__
            if colors:
                savefig = savefig + '_colors'
        #-- start figure
        if savefig:
            pl.figure()
        out = fctn(*args,**kwargs)
        #-- end figure
        if savefig:
            pl.savefig(savefig)
            pl.close()
        return out

    return dofig



def blackbody_input(fctn):
    """
    Prepare input and output for blackbody-like functions.

    If the user gives wavelength units and Flambda units, we only need to convert
    everything to SI (and back to the desired units in the end).

    If the user gives frequency units and Fnu units, we only need to convert
    everything to SI ( and back to the desired units in the end).

    If the user gives wavelength units and Fnu units, we need to convert
    the wavelengths first to frequency.
    """
    @functools.wraps(fctn)
    def dobb(x,T,**kwargs):
        wave_units = kwargs.get('wave_units','AA')
        flux_units = kwargs.get('flux_units','erg/s/cm2/AA')
        #-- prepare input
        #-- what kind of units did we receive?
        curr_conv = constants._current_convention
        # X: wavelength/frequency
        x_unit_type = conversions.get_type(wave_units)
        x = conversions.convert(wave_units,curr_conv,x)
        # T: temperature
        if isinstance(T,tuple):
            T = conversions.convert(T[1],'K',T[0])
        # Y: flux
        y_unit_type = conversions.change_convention('SI',flux_units)
        #-- if you give Jy vs micron, we need to first convert wavelength to frequency
        if y_unit_type=='kg1 rad-1 s-2' and x_unit_type=='length':
            x = conversions.convert(conversions._conventions[curr_conv]['length'],'rad/s',x)
            x_unit_type = 'frequency'
        elif y_unit_type=='kg1 m-1 s-3' and x_unit_type=='frequency':
            x = conversions.convert('rad/s',conversions._conventions[curr_conv]['length'],x)
            x_unit_type = 'length'
        #-- correct for rad
        if x_unit_type=='frequency':
            x /= (2*np.pi)
        print(y_unit_type)
        #-- run function
        I = fctn((x,x_unit_type),T)

        #-- prepare output
        disc_integrated = kwargs.get('disc_integrated',True)
        ang_diam = kwargs.get('ang_diam',None)
        if disc_integrated:
            I *= np.sqrt(2*np.pi)
            if ang_diam is not None:
                scale = conversions.convert(ang_diam[1],'sr',ang_diam[0]/2.)
                I *= scale
        I = conversions.convert(curr_conv,flux_units,I)
        return I

    return dobb
