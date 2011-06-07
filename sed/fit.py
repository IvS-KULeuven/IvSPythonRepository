# -*- coding: utf-8 -*-
"""
Fit SED models to observed data using various approaches.
"""
import logging
import sys
import pyfits
import itertools

import numpy as np
from scipy.interpolate import Rbf

from ivs.statistics import pca
from ivs.sed import model
from ivs.sed import filters
from ivs.sed.decorators import iterate_gridsearch,parallel_gridsearch
from ivs.misc import numpy_ext
from ivs.misc import progressMeter
from ivs.misc.decorators import make_parallel


logger = logging.getLogger("SED.FIT")

#{ PCA functions

def get_PCA_grid(colors,res=3,teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf),
                 ebvrange=(-np.inf,np.inf),**kwargs):
    """
    Transform a flux grid to a colour grid.
    
    The resulting grid is actually log10(flux ratio) instead of flux ratio! This
    works better in the PCA.
    
    Extra keyword arguments can be used to set the atmosphere model.
    
    @parameter colors: list of desired colours (['GENEVA.U-B','STROMGREN.M1'...])
    @type colors: list of strings
    @keyword res: resolution in E(B-V)
    @type res: integer
    @keyword teffrange: range of Teffs to use
    @type teffrange: tuple
    @keyword loggrange: range of loggs to use
    @type loggrange: tuple
    @keyword ebvrange: range of E(B-V)s to use
    @type ebvrange: tuple
    @return: log10(colorgrid), (teffs,loggs,ebvs)
    @rtype: N-D numpy array,(1Darray,1Darray,1Darray)
    """
    #-- read in the color parameters from the FITS file
    gridfile = model.get_file(integrated=True,**kwargs)
    ff = pyfits.open(gridfile)
    ext = ff[1]
    
    teff = ext.data.field('teff')
    logg = ext.data.field('logg')
    ebv = ext.data.field('ebv')
    keep = (ebvrange[0]<=ebv) & (ebv<=ebvrange[1])
    keep = keep & (teffrange[0]<=teff) & (teff<=teffrange[1])
    keep = keep & (loggrange[0]<=logg) & (logg<=loggrange[1])
    teff,logg,ebv = teff[keep],logg[keep],ebv[keep]
    
    A = model._get_flux_from_table(ext,colors)
    A = A[keep]
    
    #-- order and set the resolution of the grid
    B = np.vstack([teff.T,logg.T,ebv.T,A.T]).T
    B = numpy_ext.sort_order(B,[0,1,2])
    B = B[::res]
    logger.info('Calculated color grid for PCA (DIM=%s, using %s)'%(B[:,3:].shape,' '.join(colors)))
    return np.log10(B[:,3:]),(B[:,0],B[:,1],B[:,2])








def get_PCA(A):
    """
    Find the principal components of a color grid.
    
    @param A: color grid obtained via C{get_PCA_grid}.
    @type A. numpy N-D array
    @return: PCA loadings, PCA scores, (column means, column standard devs)
    @rtype: N-D array, N-D array, (1D array, 1D array)
    """
    #-- relocate/standardize the color grid
    means = A.mean(axis=0)
    stds  = A.std(axis=0)
    A_st = (A-means)/stds
    #-- find the principal components
    # T = PCA scores
    # P = PCA loadings
    T,P,explained_var = pca.PCA_svd(A_st,standardize=False)
    #-- report the explained variance per axis (and a maximum of four to show)
    ev = explained_var*100
    logger.info("PCA: Explained variance: %s"%(' '.join(['Axis%d=%.2f%%'%(i,j) for i,j in enumerate(ev[:4])])))
    return P,T,(means,stds)









def calibrate_PCA(T,pars,function='linear'):
    """
    Define the interpretations of the principal components.
    
    T are the PCA scores, pars a list of axes (e.g. teff, logg and ebv).
    
    @param T: PCA scores
    @type T: N-D array
    @param pars: real-life axes
    @type pars: list of 1D arrays
    @return: Rbf interpolating function
    @rtype: n-tuple
    """
    D = len(pars) # Dimension of parameter space
    calib = []
    for i in range(D):
        args = [T[:,j] for j in range(D)] + [pars[i],function]
        calib.append(Rbf(function='linear',*args[:-1]))
    logger.info('Calibrated first %d axes of PCA (of total %d)'%(len(pars),T.shape[1]))
    return tuple(calib)







def get_PCA_parameters(obsT,calib,P,means,stds,e_obsT=None,mc=None):
    """
    Derive fundamental parameters of a sample of observations given a PCA.
    
    Monte Carlo simulations are only available when you give 1 target.
    
    @param obsT: observed colours of the sample
    @type obsT: list of lists
    @param calib: PCA calibration obtained via C{calibrate_PCA}
    @type calib: list of Rbf functions
    @param P: PCA loadings
    @type P: N-D array
    @param means: means of PCA grid
    @type means: numpy array
    @param stds: standard deviations of PCA grid
    @type stds: numpy array
    @param mc: number of MonteCarlo simulations
    @type mc: integer
    @return: list of fundamental parameters
    @rtype: list of lists
    """
    #-- put obsT in same format as is used for PCA
    obsT = np.asarray(obsT)
    obsT = np.log10(obsT)
    obsT = np.dot((obsT-means)/stds,P.T)
    
    #-- this function is made for a sample of observations: if only one
    #   observations is given, make it appear to be part of a sample
    if not len(obsT.shape)==2:
        obsT = np.array([obsT])
    
    #-- if we want MC simulations, expand the array
    if mc is not None:
        if e_obsT is None:
            e_obsT = 0.01*obsT[0]
        obsT_ = np.array([obsT[0]+np.random.normal(size=len(obsT[0]),scale=e_obsT) for i in xrange(mc)])
        obsT_[0] = obsT[0]
        obsT = obsT_
    
    #-- prepare output of parameters
    pars = np.zeros((len(obsT),len(calib)))
    for i in range(len(calib)):
        args = [obsT[:,j] for j in range(len(calib))]
        pars[:,i] = calib[i](*args)
    return pars


#}

#{ Grid search

def stat_chi2(meas,e_meas,colors,syn,full_output=False):
    """
    Calculate Chi2 and compute angular diameter.
    
    Colors and absolute fluxes are used to compute the Chi2, only absolute
    fluxes are used to compute angular diameter. If no absolute fluxes are
    given, the angular diameter is set to 0.
    
    @param meas: array of measurements
    @type meas: 1D array
    @param e_meas: array containing measurements errors
    @type e_meas: 1D array
    @param colors: boolean array separating colors (True) from absolute fluxes (False)
    @type colors: 1D boolean array
    @param syn: synthetic fluxes and colors
    @type syn: 1D array
    @param full_output: set to True if you want individual chisq
    @type full_output: boolean
    @return: chi-square, scale, e_scale
    @rtype: float,float,float
    """
    if sum(-colors) > 0:
        ratio = (meas/syn)[-colors]
        weights = (meas/e_meas)[-colors]
        #-- weighted average and standard deviation
        scale = np.average(ratio,weights=weights)
        e_scale = np.sqrt(np.dot(weights, (ratio-scale)**2)/weights.sum())
    else:
        scale,e_scale = 0,0
    #-- we don't need to scale the colors, only the absolute fluxes
    chisq = np.where(colors, (syn-meas)**2/e_meas**2, (syn*scale-meas)**2/e_meas**2)
    if full_output:
        return chisq,meas/syn,meas/e_meas
    else:
        return chisq.sum(),scale,e_scale

@iterate_gridsearch
def igrid_search(meas,e_meas,photbands,teffrange=(-np.inf,np.inf),
                 loggrange=(-np.inf,np.inf),ebvrange=(-np.inf,np.inf),
                 zrange=(-np.inf,np.inf),
                 points=None,res=None,**kwargs):
    """
    Fit measured photometry to an interpolated grid of SEDs.
    
    Usage of this function is extend via the decorator.
    
    If you give an extra keyword C{iterations}, the function will iteratively
    zoom in on the minimum found in the previous stage.
    
    If you give an extra keyword C{threads}, the gridsearch will be threaded.
    
    If C{points=None}, the grid search is done on the predefined grid points.
    Otherwise, C{points} grid points will be generated, uniformly distributed
    between the ranges defined by C{teffrange}, C{loggrange} and C{ebvrange},
    and afterwards all the points outside the grid are removed (i.e., you could
    end up with less points that given). It is probably safest to go for
    C{points=None} and adapt the resolution of the grid (C{res}=2). If you set
    the resolution to '2', one out of every two points will be selected.
    
    Extra keyword arguments can be used to give more details on the atmosphere
    models to use.
    
    Colors are automatically detected.
    
    You can fix one parameter e.g. via setting teffrange=(10000,10000), but make
    sure to either fix it on a grid point, or encompassing grid points!
    
    @param meas: array of measurements
    @type meas: 1D array
    @param e_meas: array containing measurements errors
    @type e_meas: 1D array
    @param photbands: a list of photometric passbands, corresponding each
    measurement
    @type photbands: list of strings
    @param teffrange: range of temperatures to use
    @type teffrange: 2-tuple
    @param loggrange: range of surface gravities to use
    @type loggrange: 2-tuple
    @param ebvrange: range of reddenings to use
    @type ebvrange: 2-tuple
    @param points: points to sample (when None, predefined grid points are used)
    @type points: int
    @param res: resolution of the original grid (the higher, the coarser)
    @type res: int
    @keyword threads: number of threads to use (defaults to 1)
    @type threads: int
    @keyword iterations: number of iterations for zoom-in
    @type iterations: int
    @return: record array containing the searched grid, chi-squares and scale
    factors
    @rtype: record array
    """
    colors = np.array([filters.is_color(photband) for photband in photbands],bool)
    logger.info('Grid search with parameters teffrange=%s, loggrange=%s, ebvrange=%s, zrange=%s, points=%s'%(teffrange,loggrange,ebvrange,zrange,points))
    
    #-- we first get/set the grid. Calling this function means it will be
    #   memoized, so that we can safely thread (and don't have to memoize for
    #   each thread)
    markers,(unique_teffs,unique_loggs,unique_ebvs,unique_zs),gridpnts,flux = \
                 model._get_itable_markers(photbands,ebvrange=(-np.inf,np.inf),
                        zrange=(-np.inf,np.inf),include_Labs=True)
    teffs,loggs,ebvs,zs = gridpnts.T
    
    unique_teffs = unique_teffs[(teffrange[0]<=unique_teffs) & (unique_teffs<=teffrange[1])]
    unique_loggs = unique_loggs[(loggrange[0]<=unique_loggs) & (unique_loggs<=loggrange[1])]
    
    #-- if we gave a certain number of points, we need to choose our points
    #   randomly in the grid
    if points:
        #-- generate appropriate evaluation points uniformly within the predefined
        #   edges
        #-- first list all effective temperatures and their minimum logg in the grid
        teff_min_logg = np.zeros((len(unique_teffs),2))
        for i,iteff in enumerate(unique_teffs):
            teff_min_logg[i] = iteff,(loggs[teffs==iteff]).min()
        #-- we have one square per logg: calculate their sizes
        unique_min_loggs = sorted(list(set(teff_min_logg[:,1])))
        limits_and_sizes = []
        for index,unique_min_logg in enumerate(unique_min_loggs):
            min_teff = teff_min_logg[:,0][teff_min_logg[:,1]==unique_min_logg].min()
            #-- we need to avoid having gaps in the grid:
            if index>0:
                min_teff = max_teff
            max_teff = teff_min_logg[:,0][teff_min_logg[:,1]==unique_min_logg].max()
            min_logg = unique_min_logg
            max_logg = loggs.max()
            #-- we're at too low temperatures
            if max_teff<teffrange[0]: continue
            else:
                min_teff = max(teffrange[0],min_teff)
                max_teff = min(teffrange[1],max_teff)
            #-- we're at too low surface gravities:
            min_logg = max(loggrange[0],min_logg)
            max_logg = min(loggrange[1],max_logg)
            #-- make sure there are points defined even if some range in parameters
            #   equals zero
            if (max_teff-min_teff)>1 and (max_logg-min_logg)>1:
                size = (max_teff-min_teff)*(max_logg-min_logg)
            elif (max_teff-min_teff)>1:
                size = (max_teff-min_teff)
            elif (max_logg-min_logg)>1:
                size = (max_logg-min_logg)
            else:
                size = int(float(points)/(len(unique_min_loggs)))
            #print 'size',size
            limits_and_sizes.append([(min_teff,max_teff),(min_logg,max_logg),size])
        total_size = sum([row[-1] for row in limits_and_sizes])
        #print total_size
        #print limits_and_sizes
        teffs,loggs,ebvs,zs = np.hstack([np.random.uniform(low=[lims[0][0],lims[1][0],ebvrange[0],zrange[0]],
                                                       high=[lims[0][1],lims[1][1],ebvrange[1],zrange[1]],
                                                       size=(int(lims[-1]/total_size*points),4)).T for lims in limits_and_sizes])
        #print teffs.shape
    keep = (teffrange[0]<=teffs) & (teffs<=teffrange[1]) &\
            (loggrange[0]<=loggs) & (loggs<=loggrange[1]) &\
            (ebvrange[0]<=ebvs) & (ebvs<=ebvrange[1]) &\
            (zrange[0]<=zs) & (zs<=zrange[1])
    #print sum(keep)
    teffs,loggs,ebvs,zs = teffs[keep],loggs[keep],ebvs[keep],zs[keep]
    
    if res:
        teffs,loggs,ebvs,zs = teffs[::res],loggs[::res],ebvs[::res],zs[::res]
    logger.info('Evaluating %d points in parameter space'%(len(teffs)))
    
    #-- run over the grid and calculate chisq and scale factors for each point
    @parallel_gridsearch
    @make_parallel
    def do_grid_search(teffs,loggs,ebvs,zs,meas,e_meas,photbands,colors,**kwargs):
        """
        index is to trace the results after parallelization
        """
        index = kwargs.pop('index',None)
        chisqs = np.zeros_like(teffs)
        scales = np.zeros_like(teffs)
        e_scales = np.zeros_like(teffs)
        lumis = np.zeros_like(teffs)
        if index is None:
            p = progressMeter.ProgressMeter(total=len(teffs))
        for n,(teff,logg,ebv,z) in enumerate(itertools.izip(teffs,loggs,ebvs,zs)):
            if index is None: p.update(1)
            syn_flux,Labs = model.get_itable(teff=teff,logg=logg,ebv=ebv,z=z,photbands=photbands)
            chisqs[n],scales[n],e_scales[n] = stat_chi2(meas,e_meas,colors,syn_flux)
            lumis[n] = Labs
        return chisqs,scales,e_scales,lumis,index
    
    #-- then we do the grid search
    chisqs,scales,e_scales,lumis,index = do_grid_search(teffs,loggs,ebvs,zs,meas,e_meas,photbands,colors,**kwargs)
    #-- transform output to record array
    data_rec = np.rec.fromarrays([teffs,loggs,ebvs,zs,chisqs,scales,e_scales,lumis],
                   dtype=[('teff','f8'),('logg','f8'),('ebv','f8'),('z','f8'),
                          ('chisq','f8'),('scale','f8'),('e_scale','f8'),('Labs','f8')])
    return data_rec

#}
def do_grid_search2(teffs,loggs,ebvs,zs,meas,e_meas,photbands,colors,**kwargs):
        """
        index is to trace the results after parallelization
        """
        index = kwargs.pop('index',None)
        chisqs = np.zeros_like(teffs)
        scales = np.zeros_like(teffs)
        e_scales = np.zeros_like(teffs)
        lumis = np.zeros_like(teffs)
        if index is None:
            p = progressMeter.ProgressMeter(total=len(teffs))
        for n,(teff,logg,ebv,z) in enumerate(itertools.izip(teffs,loggs,ebvs,zs)):
            if index is None: p.update(1)
            syn_flux,Labs = model.get_itable(teff=teff,logg=logg,ebv=ebv,z=z,photbands=photbands)
            chisqs[n],scales[n],e_scales[n] = stat_chi2(meas,e_meas,colors,syn_flux)
            lumis[n] = Labs
        return chisqs,scales,e_scales,lumis,index



if __name__=="__main__":
    from ivs.misc import loggers
    from pylab import *
    from numpy import *
    logger = loggers.get_basic_logger()
    random.seed(1111)
    A,grid = get_PCA_grid(['GENEVA.U-B','GENEVA.B1-B','GENEVA.B2-B','GENEVA.V-B','GENEVA.V1-B','GENEVA.G-B','2MASS.J-H','2MASS.KS-H'],ebvrange=(0,0.5),res=10)
    P,T,(means,stds) = get_PCA(A)
    calib = calibrate_PCA(T,grid,function='linear')
    sample_index = int(np.random.uniform(high=len(A)))
    sample = 10**A[sample_index]
    print [bla[sample_index] for bla in grid]
    print get_PCA_parameters(sample,calib,P,means,stds)