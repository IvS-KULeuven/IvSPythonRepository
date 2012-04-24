# -*- coding: utf-8 -*-
"""
Fit SED models to observed data using various approaches.
"""
import logging
import sys
import pyfits
import itertools

import numpy as np
from numpy import inf
from scipy.interpolate import Rbf

from ivs.statistics import pca
from ivs.sed import model
from ivs.sed import filters
from ivs.sed.decorators import iterate_gridsearch,parallel_gridsearch
from ivs.aux import numpy_ext
from ivs.aux import progressMeter
from ivs.aux.decorators import make_parallel
from ivs.units import constants

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

def generate_grid_single(photbands,teffrange=(-inf,inf),loggrange=(-inf,inf),
                  ebvrange=(-inf,inf),zrange=(-inf,inf),
                  points=None,res=None,clear_memory=True,**kwargs):
    """
    Generate grid points at which to fit an interpolated grid of SEDs.
    
    If C{points=None}, the points are chosen on the predefined grid points.
    Otherwise, C{points} grid points will be generated, uniformly distributed
    between the ranges defined by C{teffrange}, C{loggrange} and C{ebvrange}. If
    you set the resolution to C{2}, one out of every two points will be selected.
    
    Extra keyword arguments can be used to give more details on the atmosphere
    models to use.
    
    Colors are automatically detected.
    
    You can fix one parameter e.g. via setting teffrange=(10000,10000).
    
    >>> photbands = ['GENEVA.G','GENEVA.B-V']
    
    Start a figure:
    
    >>> p = pl.figure()
    >>> rows,cols = 2,4
    
    On the grid points, but only one in every 100 points (otherwise we have over
    a million points):
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,res=100)
    
    >>> p = pl.subplot(rows,cols,1)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,1+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    Randomly distributed over the grid's ranges:
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,points=10000)
    
    >>> p = pl.subplot(rows,cols,2)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,2+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    Confined to a small area in the grid's range:
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,teffrange=(8000,10000),loggrange=(4.1,4.2),zrange=(0,inf),ebvrange=(1.,2),points=10000)
    
    >>> p = pl.subplot(rows,cols,3)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,3+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    Confined to a small area in the grid's range with some parameters fixed:
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,teffrange=(8765,8765),loggrange=(4.1,4.2),zrange=(0,0),ebvrange=(1,2),points=10000)
    
    >>> p = pl.subplot(rows,cols,4)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,4+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    ]include figure]]ivs_sed_fit_grids.png]
    
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
    @keyword clear_memory: flag to clear memory from previously loaded SED tables.
    If you set it to False, you can easily get an overloaded memory!
    @type clear_memory: boolean
    @return: record array containing the searched grid, chi-squares and scale
    factors
    @rtype: record array
    """
    #test
    logger.info('Grid search with parameters teffrange=%s, loggrange=%s, ebvrange=%s, zrange=%s, points=%s'%(teffrange,loggrange,ebvrange,zrange,points))
    
    #-- we first get/set the grid. Calling this function means it will be
    #   memoized, so that we can safely thread (and don't have to memoize for
    #   each thread). We also have an exact view of the size of the grid here...
    markers,(unique_teffs,unique_loggs,unique_ebvs,unique_zs),gridpnts,flux = \
                 model._get_itable_markers(photbands,ebvrange=(-np.inf,np.inf),
                        zrange=(-np.inf,np.inf),include_Labs=True,
                        clear_memory=clear_memory,**kwargs)
    if not kwargs:
        logger.info('Received grid (%s)'%model.defaults2str())
    else:
        logger.info('Received custom grid (%s)'%kwargs)
    teffs,loggs,ebvs,zs = gridpnts.T
    
    #-- We need to avoid having only one grid point! If nessessary the grid needs to be 
    #   broader to get points in the entire intervall. 
    index1 = teffrange[0] in unique_teffs and unique_teffs.searchsorted(teffrange[0]) or \
            max(0,unique_teffs.searchsorted(teffrange[0])-1)
    index2 = teffrange[1] in unique_teffs and unique_teffs.searchsorted(teffrange[1]) or \
            min(len(unique_teffs),unique_teffs.searchsorted(teffrange[1]))
    unique_teffs = unique_teffs[index1:index2+1]
    index1 = max(0,unique_teffs.searchsorted(loggrange[0])-1)
    index2 = unique_teffs.searchsorted(loggrange[1])+1
    unique_loggs = unique_loggs[index1:index2+1]
    
    #-- if we gave a certain number of points, we need to choose our points
    #   randomly in the grid: the grid is usually not a square, so we have to
    #   subdivide it
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
            if (max_teff-min_teff)>1 and (max_logg-min_logg)>0.01:
                size = (max_teff-min_teff)*(max_logg-min_logg)
            elif (max_teff-min_teff)>1:
                size = (max_teff-min_teff)
            elif (max_logg-min_logg)>1:
                size = (max_logg-min_logg)
            else:
                size = int(float(points)/(len(unique_min_loggs)))
            if size==0: size=2
            #-- sizes of ebv and z:
            zrange_   = max(  zrange[0],min(unique_zs)),   min(  zrange[1],max(unique_zs))
            ebvrange_ = max(ebvrange[0],min(unique_ebvs)), min(ebvrange[1],max(unique_ebvs))
            limits_and_sizes.append([(min_teff,max_teff),(min_logg,max_logg),ebvrange_,zrange_,size])
        total_size = sum([row[-1] for row in limits_and_sizes])
        #-- in the following case, we fall in between the grid points. We correct
        #   for this
        if len(limits_and_sizes)==0:
            total_size = points
            limits_and_sizes = [[teffrange,loggrange,ebvrange,zrange,points]]
        logger.debug('Limits and sizes of boxes:'+str(limits_and_sizes))
        teffs,loggs,ebvs,zs = np.hstack([np.random.uniform(low=[lims[0][0],lims[1][0],lims[2][0],lims[3][0]],
                                                       high=[lims[0][1],lims[1][1],lims[2][1],lims[3][1]],
                                                       size=(int(lims[-1]/total_size*points),4)).T for lims in limits_and_sizes])
    #-- override grid generation... seems to still be troubled....                                                   
    #if points:
    #    teffs,loggs,ebvs,zs = np.random.uniform(low=[teffrange[0],loggrange[0],ebvrange_[0],zrange_[0]],
    #                                                   high=[teffrange[1],loggrange[1],ebvrange_[1],zrange_[1]],
    #                                                   size=(points,4)).T
    keep = (teffrange[0]<=teffs) & (teffs<=teffrange[1]) &\
            (loggrange[0]<=loggs) & (loggs<=loggrange[1]) &\
            (ebvrange[0]<=ebvs) & (ebvs<=ebvrange[1]) &\
            (zrange[0]<=zs) & (zs<=zrange[1])
    teffs,loggs,ebvs,zs = teffs[keep],loggs[keep],ebvs[keep],zs[keep]
    
    if res:
        teffs,loggs,ebvs,zs = teffs[::res],loggs[::res],ebvs[::res],zs[::res]
    logger.info('Evaluating %d points in parameter space'%(len(teffs)))
    
    ##-- run over the grid and calculate chisq and scale factors for each point
    ##-- then we do the grid search
    #chisqs,scales,e_scales,lumis,index = do_grid_search(teffs,loggs,ebvs,zs,meas,e_meas,photbands,colors,**kwargs)
    ##-- transform output to record array
    #data_rec = np.rec.fromarrays([teffs,loggs,ebvs,zs,chisqs,scales,e_scales,lumis],
                   #dtype=[('teff','f8'),('logg','f8'),('ebv','f8'),('z','f8'),
                          #('chisq','f8'),('scale','f8'),('e_scale','f8'),('Labs','f8')])
    return teffs,loggs,ebvs,zs

def generate_grid(photbands,teffrange=((-inf,inf),(-inf,inf)),
                  loggrange=((-inf,inf),(-inf,inf)),ebvrange=(-inf,inf),
                  zrange=((-inf,inf),(-inf,inf)),
                  radiusrange=((1,1),(0.1,10.)),grids=None,
                  points=None,res=None,clear_memory=False,
                  type='single', **kwargs):                    
    """
    Generate grid points at which to fit an interpolated grid of SEDs.
    
    If C{points=None}, the points are chosen on the predefined grid points.
    Otherwise, C{points} grid points will be generated, uniformly distributed
    between the ranges defined by C{teffrange}, C{loggrange} and C{ebvrange}. If
    you set the resolution to C{2}, one out of every two points will be selected.
    
    Extra keyword arguments can be used to give more details on the atmosphere
    models to use.
    
    Colors are automatically detected.
    
    You can fix one parameter e.g. via setting teffrange=(10000,10000).
    
    >>> photbands = ['GENEVA.G','GENEVA.B-V']
    
    Start a figure:
    
    >>> p = pl.figure()
    >>> rows,cols = 2,4
    
    On the grid points, but only one in every 100 points (otherwise we have over
    a million points):
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,res=100)
    
    >>> p = pl.subplot(rows,cols,1)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,1+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    Randomly distributed over the grid's ranges:
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,points=10000)
    
    >>> p = pl.subplot(rows,cols,2)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,2+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    Confined to a small area in the grid's range:
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,teffrange=(8000,10000),loggrange=(4.1,4.2),zrange=(0,inf),ebvrange=(1.,2),points=10000)
    
    >>> p = pl.subplot(rows,cols,3)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,3+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    Confined to a small area in the grid's range with some parameters fixed:
    
    >>> teffs,loggs,ebvs,zs = generate_grid(photbands,teffrange=(8765,8765),loggrange=(4.1,4.2),zrange=(0,0),ebvrange=(1,2),points=10000)
    
    >>> p = pl.subplot(rows,cols,4)
    >>> p = pl.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlim(pl.xlim()[::-1]);p = pl.ylim(pl.ylim()[::-1])
    >>> p = pl.xlabel('Teff');p = pl.ylabel('Logg')
    
    >>> p = pl.subplot(rows,cols,4+cols)
    >>> p = pl.scatter(ebvs,zs,c=teffs,s=(loggs+2)*10,edgecolors='none',cmap=pl.cm.spectral)
    >>> p = pl.xlabel('E(B-V)');p = pl.ylabel('Z')
    
    ]include figure]]ivs_sed_fit_grids.png]
    
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
    @keyword clear_memory: flag to clear memory from previously loaded SED tables.
    If you set it to False, you can easily get an overloaded memory!
    @type clear_memory: boolean
    @return: record array containing the searched grid, chi-squares and scale
    factors
    @rtype: record array
    """                  
    
    #-- Select the grid
    #   but remove metallicity, as it will be fitted!
    if type=='single':
        #--Single grid, uses the basic function
        teffs,loggs,ebvs,zs = generate_grid_single(photbands,teffrange=teffrange,
                      loggrange=loggrange,ebvrange=ebvrange,
                      zrange=zrange,points=points)
        radii = [1 for i in teffs]
        return teffs,loggs,ebvs,zs,radii
    
    #-- first collect the effetive temperatures, loggs, ebvs, zs for the
    #   different stars in the multiple system
    pars = []
    if grids is None:
        grids = model.defaults_multiple
        #-- but remove metallicity, as it will be fitted!
        for grid in grids:
            if 'z' in grid:
                thrash = grid.pop('z')
    for i,grid in enumerate(grids):
        #-- it is possible that we want certain parameters to be the same for
        #   all components
        teffrange_ = hasattr(teffrange[0],'__iter__') and teffrange[i] or teffrange
        loggrange_ = hasattr(loggrange[0],'__iter__') and loggrange[i] or loggrange
        ebvrange_ = hasattr(ebvrange[0],'__iter__') and ebvrange[i] or ebvrange
        zrange_ = hasattr(zrange[0],'__iter__') and zrange[i] or zrange
        pars += list(generate_grid_single(photbands,teffrange=teffrange_,
                      loggrange=loggrange_,ebvrange=ebvrange_,
                      zrange=zrange_,points=points,**grid))
    #-- the L{generate_grid} method does not guarantee the number of points.
    #   We have to strip some points if the arrays don't have the same shape
    nmin = np.min([len(i) for i in pars])
    pars = [i[:nmin] for i in pars]
    pars = np.array(pars)
    #-- permute parameters so that the different blocks from the generate_grid
    #   are not clustered together
    for i in range(0,len(pars),4):
        permutation = np.random.permutation(len(pars[0]))
        pars[i:i+4] = pars[i:i+4,permutation]
    #-- make arrays of the output parameters
    teffs,loggs,ebvs,zs = pars[0::4].T,pars[1::4].T,pars[2::4].T,pars[3::4].T
    
    #-- keep in mind that we probably want all the members in the system to have
    #   the same value for the interstellar reddening and metallicity, though
    #   this is not prerequisitatory
    if not hasattr(teffrange[0],'__iter__'): teffs = np.column_stack([teffs[:,0]]*len(grids))
    if not hasattr(loggrange[0],'__iter__'): loggs = np.column_stack([loggs[:,0]]*len(grids))
    #if not hasattr(ebvrange[0],'__iter__'): ebvs = np.column_stack([ebvs[:,0]]*len(grids))
    #if not hasattr(zrange[0],'__iter__'): zs = np.column_stack([zs[:,0]]*len(grids))
    ebvs = np.column_stack([ebvs[:,0]]*len(grids))
    zs = np.column_stack([zs[:,0]]*len(grids))
    
    if type=='binary':
        #-- The radius of the stars is calculated bassed on logg and the provided masses
        masses = 'masses' in kwargs and  kwargs['masses'] or (1,1)
        G = constants.GG_cgs
        Msol = constants.Msol_cgs
        Rsol = constants.Rsol_cgs
        radius1 = np.sqrt(G*masses[0]*Msol/10**loggs[:,0])/Rsol
        radius2 = np.sqrt(G*masses[1]*Msol/10**loggs[:,1])/Rsol
        #radii = radius2/radius1
        
        #radii = np.array([np.ones(len(radii)),radii]).T
        radii = np.array([radius1,radius2]).T
    elif type=='multiple':
        #-- We have random different radii for the stars
        radii = np.array([10**np.random.uniform(low=radiusrange[0][0], high=radiusrange[0][1], size=(len(teffs))), 
                    10**np.random.uniform(low=radiusrange[1][0], high=radiusrange[1][1], size=(len(teffs)))]).T
        #radii = 10**np.random.uniform(low=[np.log10(i[0]) for i in radiusrange],
                              #high=[np.log10(i[1]) for i in radiusrange],size=(len(teffs),2))                       
    
    return teffs,loggs,ebvs,zs,radii                     
    
@parallel_gridsearch
@make_parallel
def igrid_search(meas,e_meas,photbands,*args,**kwargs):
        """
        Run over gridpoints and evaluate model C{model_func} via C{stat_func}.
        
        The measurements are defined via C{meas, e_meas, photbands, colors} and
        should be 1d arrays of equal length. C{colors} should be a boolean
        array, C{photbands} should be a string array.
        
        The grid points are defined via C{args}. C{args} should be a tuple of 
        1x? dimensional arrays of equal length. For single stars, this is
        typically effective temperatures, loggs, reddenings and metallicities.
        For multiple systems, (at least some of the) previously mentioned
        parameters are typically doubled, and radius ratios are added. Remember
        to specify the C{model_func} to match single or multiple systems.
        
        At each grid point, the pre-calculated photometry will be retrieved via
        the keyword C{model_func} and compared to the measurements via the function
        definded via C{stat_func}. This function should be of the same form as
        L{stat_chi2}.
        
        Extra arguments are passed to L{parallel_gridsearch} for parallelization
        and to {model_func} for further specification of grids etc.
        
        The index array is returned to trace the results after parallelization.
        
        @param meas: the measurements that have to be compared with the models
        @type meas: 1D numpy array of floats
        @param e_meas: errors on the measurements
        @type e_meas: 1D numpy array of floats
        @param photbands: names of the photometric passbands
        @type photbands: 1D numpy array of strings
        @keyword model_func: function to translate parameters to synthetic (model) data
        @type model_func: function
        @keyword stat_func: function to evaluate the fit
        @type stat_func: function
        @return: (chi squares, scale factors, error on scale factors, absolute
        luminosities (R=1Rsol), index
        @rtype: 4/5X1d array
        """
        model_func = kwargs.pop('model_func',model.get_itable)
        stat_func = kwargs.pop('stat_func',stat_chi2)
        index = kwargs.pop('index',None)
        N = len(args[0])
        #-- prepare output arrays
        chisqs = np.zeros(N)
        scales = np.zeros(N)
        e_scales = np.zeros(N)
        lumis = np.zeros(N)
        colors = np.array([filters.is_color(photband) for photband in photbands],bool)
        #-- show a progressMeter when not parallelized
        if index is None:
            p = progressMeter.ProgressMeter(total=N)
        #-- run over the grid, retrieve synthetic fluces and compare with
        #   observations.
        for n,pars in enumerate(itertools.izip(*args)):
            if index is None: p.update(1)
            syn_flux,Labs = model_func(*pars,photbands=photbands,**kwargs)
            chisqs[n],scales[n],e_scales[n] = stat_func(meas,e_meas,colors,syn_flux)
            lumis[n] = Labs
        #-- return results
        if index is not None:
            return chisqs,scales,e_scales,lumis,index
        else:
            return chisqs,scales,e_scales,lumis

def residual_single(params, meas, photbands, kwargs):
    teff = params['teff'].value
    logg = params['logg'].value
    ebv = params['ebv'].value
    z = params['z'].value

    iflux, Labs = model.get_itable(teff=teff, logg=logg, ebv=ebv, z=z, photbands=photbands, **kwargs)
    scale = np.average( meas / iflux )
    
    print np.array(meas-iflux*scale).shape
    return (meas-iflux*scale)#**2 / e_meas**2
        
def residual_multiple(params, meas, e_meas, photbands, kwargs):
    teff = params['teff'].value
    logg = params['logg'].value
    rad = params['rad'].value
    teff2 = params['teff2'].value
    logg2 = params['logg2'].value
    rad2 = params['rad2'].value
    ebv = params['ebv'].value
    z = params['z'].value
    
    print teff,logg,rad,teff2,logg2,rad2
    
    if logg > 5.00: logg=5.00
    if logg2 > 6.00: logg2=6.00
    if logg2 < 5.00: logg2=5.00

    iflux, Labs = model.get_itable_multiple(teff=(teff,teff2), logg=(logg,logg2), ebv=(ebv,ebv),
                                      z=(z,z), radius=(rad,rad2), photbands=photbands, **kwargs)
    scale = np.average( meas / iflux )
    
    return (meas-iflux*scale)**2 / e_meas**2

def residual_binary(params, meas, e_meas, photbands, kwargs):
    teff = params['teff'].value
    logg = params['logg'].value
    rad = params['rad'].value
    teff2 = params['teff2'].value
    logg2 = params['logg2'].value
    rad2 = params['rad2'].value
    ebv = params['ebv'].value
    z = params['z'].value

    iflux, Labs = model.get_itable_multiple(teff=(teff,teff2), logg=(logg,logg2), ebv=(ebv,ebv),
                                      z=(z,z), radius=(rad,rad2), photbands=photbands, **kwargs)
    scale = np.average( meas / iflux )
    
    return (meas-iflux*scale)**2 / e_meas**2

def dfun(params, meas, photbands, kwargs):
    dt = (model.get_itable(teff=5200, logg=4.25, ebv=0.0, z=0.0, photbands=photbands, **kwargs)[0] - 
            model.get_itable(teff=5000, logg=4.25, ebv=0.0, z=0.0, photbands=photbands, **kwargs)[0] ) / (5200-5000)
            
    dg = (model.get_itable(teff=5000, logg=4.50, ebv=0.0, z=0.0, photbands=photbands, **kwargs)[0] - 
            model.get_itable(teff=5000, logg=4.25, ebv=0.0, z=0.0, photbands=photbands, **kwargs)[0] ) / (4.50-4.25)
            
    de = (model.get_itable(teff=5000, logg=4.25, ebv=0.1, z=0.0, photbands=photbands, **kwargs)[0] - 
            model.get_itable(teff=5000, logg=4.25, ebv=0.0, z=0.0, photbands=photbands, **kwargs)[0] ) / (0.1-0.0)
            
    dz = np.ones(len(meas))
    
    print np.array([dt,dg,de,dz]).shape
    
    return [dt,dg,de,dz]

def iminimize(meas,e_meas,photbands,**kwargs):
    """
    Using http://newville.github.com/lmfit-py/
    """
    from ivs.lmfit import minimize, Parameters
    model_func = kwargs.pop('model_func',model.get_itable)
    engine = kwargs.pop('engine','leastsq')
    type = kwargs.pop('type','single')
    
    if type == 'single':
        residual = residual_single
        
        params = Parameters()
        params.add('teff', value=kwargs.pop('teff'), min=15000, max=25000)
        params.add('logg', value=kwargs.pop('logg'), min=4.0, max=4.99)
        params.add('ebv', value=kwargs.pop('ebv'), min=0.0, max=0.2, vary=False)
        params.add('z', value=kwargs.pop('z'), vary=False)
    
    if type == 'multiple':
        residual = residual_multiple
        params = Parameters()
        params.add('teff', value=kwargs['teff'][0], min=5000, max=7000)
        params.add('logg', value=kwargs['logg'][0], min=4.0, max=4.99)
        params.add('rad', value=kwargs['rad'][0], min=0.5, max=2.0, vary=False)
        params.add('teff2', value=kwargs.pop('teff')[1], min=20000, max=40000)
        params.add('logg2', value=kwargs.pop('logg')[1], min=5.0, max=6.0)
        params.add('rad2', value=kwargs.pop('rad')[1], min=0.05, max=0.5, vary=False)
        params.add('ebv', value=kwargs.pop('ebv'), min=0.0, max=0.01, vary=False)
        params.add('z', value=kwargs.pop('z'), vary=False)
        
    if type == 'binary':
        residual = residual_binary
        params = Parameters()
        params.add('mass', value=kwargs['masses'][0], vary=False)
        params.add('mass2', value=kwargs.pop('masses')[1], vary=False)
        params.add('teff', value=kwargs['teff'][0], min=5000, max=7000)
        params.add('logg', value=kwargs['logg'][0], min=4.0, max=4.99)
        params.add('rad', expr='sqrt(mass / 10**logg)')
        params.add('teff2', value=kwargs.pop('teff')[1], min=20000, max=40000)
        params.add('logg2', value=kwargs.pop('logg')[1], min=5.0, max=6.0)
        params.add('rad2', expr='sqrt(mass2 / 10**logg2)')
        params.add('ebv', value=kwargs.pop('ebv'), min=0.0, max=0.01, vary=True)
        params.add('z', value=kwargs.pop('z'), vary=False)
    
    out = minimize(residual, params, args=(meas, photbands, kwargs), engine=engine, Dfun=dfun, col_deriv=1)#, ftol=0.0005)#, diag=(500,0.5,0.01,0.01), factor=1/5.0)
    
    return out, params

#}


if __name__=="__main__":
    from ivs.aux import loggers
    import time
    import pylab as plt
    logger = loggers.get_basic_logger(clevel='DEBUG')
    
    photbands = ['GENEVA.G','GENEVA.B-V']
    
    c0 = time.time()
    teffs,loggs,ebvs,zs,radii = generate_grid(photbands,teffrange=(5000,5800),loggrange=(4.20,4.70),zrange=(0,0),ebvrange=(0.05,0.08), grid='kurucz',points=10000)
    print 'Time: %i'%(time.time()-c0)
    
    plt.figure(2)
    plt.scatter(teffs,loggs,c=ebvs,s=(zs+5)*10,edgecolors='none',cmap=plt.cm.spectral)
    plt.xlim(plt.xlim()[::-1])
    plt.ylim(plt.ylim()[::-1])
    plt.xlabel('Teff')
    plt.ylabel('Logg')
    plt.show()
    
    sys.exit()
    
    
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()
    
    
    sys.exit()
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