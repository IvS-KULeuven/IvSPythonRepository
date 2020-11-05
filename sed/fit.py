# -*- coding: utf-8 -*-
"""
Fit SED models to observed data using various approaches.
"""
import logging
import sys
import astropy.io.fits as pf
import itertools
import re
import copy

import numpy as np
from numpy import inf
import scipy
from scipy.interpolate import Rbf
from scipy.optimize import fmin,fmin_powell
from ivs.statistics import pca
from ivs.sed import model
from ivs.sed import filters
from ivs.sed.decorators import parallel_gridsearch
from ivs.sigproc import fit as sfit
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
    ff = pf.open(gridfile)
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
        obsT_ = np.array([obsT[0]+np.random.normal(size=len(obsT[0]),scale=e_obsT) for i in range(mc)])
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

def stat_chi2(meas,e_meas,colors,syn,full_output=False, **kwargs):
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
    #-- if syn represents only one measurement
    #syn*=2
    if len(syn.shape)==1:
        if sum(~colors) > 0:
            if 'distance' in kwargs:
                scale = 1/kwargs['distance']**2
                e_scale = scale / 100
            else:
                ratio = (meas/syn)[~colors]
                weights = (meas/e_meas)[~colors]
                #-- weighted average and standard deviation
                scale = np.average(ratio,weights=weights)
                #print 'bla',weights.shape,ratio.shape,scale
                e_scale = np.sqrt(np.dot(weights, (ratio-scale)**2)/weights.sum())
        else:
            scale,e_scale = 0,0
        #-- we don't need to scale the colors, only the absolute fluxes
        chisq = np.where(colors, (syn-meas)**2/e_meas**2, (syn*scale-meas)**2/e_meas**2)
        if full_output:
            return chisq,meas/syn,meas/e_meas
        else:
            return chisq.sum(),scale,e_scale
    #-- if syn is many measurements, we need to vectorize this:
    else:
        if sum(~colors) > 0:
            if 'distance' in kwargs:
                scale = 1/kwargs['distance']**2
                scale = np.ones_like(syn[~colors][0,:]) * scale
                e_scale = scale / 100
            else:
                ratio = (meas/syn)[~colors]
                weights = (meas/e_meas)[~colors]
                #-- weighted average and standard deviation
                scale = np.average(ratio,weights=weights.reshape(-1),axis=0)
                e_scale = np.sqrt(np.dot(weights.T, (ratio-scale)**2)/weights.sum(axis=0))[0]
                #scale = np.average(ratio,axis=0)
                #e_scale = np.zeros_like(scale)
        else:
            scale,e_scale = np.zeros(syn.shape[1]),np.zeros(syn.shape[1])
        #-- we don't need to scale the colors, only the absolute fluxes
        chisq = np.where(colors.reshape(-1,1), (syn-meas)**2/e_meas**2, (syn*scale-meas)**2/e_meas**2)
        if full_output:
            return chisq,meas/syn,meas/e_meas
        else:
            return chisq.sum(axis=0),scale,e_scale


def generate_grid_single_pix(photbands, points=None, clear_memory=True, **kwargs):
    """
    Generate a grid of parameters.
    """
    #-- Find the parameters provided and store them separately.
    ranges, parameters = {}, []
    for key in list(kwargs.keys()):
        if re.search('range$', key):
            ranges[key] = kwargs.pop(key)
            parameters.append(re.sub('range$', '', key))
    #-- report on the received grid
    if not kwargs:
        logger.info('Received grid (%s)'%model.defaults2str())
    else:
        logger.info('Received custom grid (%s)'%kwargs)

    #-- get the pixelgrid
    axis_values,gridpnts,flux,colnames = \
                 model._get_pix_grid(photbands,teffrange=(-inf,inf),
                 loggrange=(-inf,inf),ebvrange=(-inf,inf),
                 zrange=(-inf,inf),rvrange=(-inf,inf),vradrange=(0,0),
                 include_Labs=True,clear_memory=clear_memory,
                 variables=parameters, **kwargs)

    #-- we first generate random teff-logg coordinates, since the grid is
    #   not exactly convex in these parameters. We assume it is for all the
    #   other parameters. We need to extract the teff/logg points, but make them
    #   unique
    colnames = list(colnames)
    teff_index = colnames.index('teff')
    logg_index = colnames.index('logg')
    teffs,loggs = gridpnts[:,teff_index],gridpnts[:,logg_index]

    #-- get ranges for teff and logg
    teffrange = ranges.pop('teffrange', (-inf,inf))
    loggrange = ranges.pop('loggrange', (-inf,inf))
    correctTeff, correctLogg = False, False
    if teffrange[0] == teffrange[1]:
        teffrange = [teffrange[0], teffrange[0]+1]
        correctTeff = True
    if loggrange[0] == loggrange[1]:
        loggrange = [loggrange[0], loggrange[0]+0.01]
        correctLogg = True

    #-- we need to cut the grid to fit the teff and logg range: we replace the
    #   values for the upper and lower limit in the grid with those from the
    #   given ranges. This is a bit elaborate, but I don't see a better way
    #   of doin' it.
    teffl_index = max(np.searchsorted(axis_values[teff_index],teffrange[0])-1,0)
    teffu_index = min(np.searchsorted(axis_values[teff_index],teffrange[1]),len(axis_values[teff_index])-1)
    teff_lower = axis_values[teff_index][teffl_index]
    teff_upper = axis_values[teff_index][teffu_index]
    cut = (teffs<teff_lower) | (teff_upper<teffs)
    if teff_lower<teffrange[0]: teffs[teffs==teff_lower] = teffrange[0]
    if teff_upper>teffrange[1]: teffs[teffs==teff_upper] = teffrange[1]

    loggl_index = max(np.searchsorted(axis_values[logg_index],loggrange[0])-1,0)
    loggu_index = min(np.searchsorted(axis_values[logg_index],loggrange[1]),len(axis_values[logg_index])-1)
    logg_lower = axis_values[logg_index][loggl_index]
    logg_upper = axis_values[logg_index][loggu_index]
    cut = cut | (loggs<logg_lower) | (logg_upper<loggs)
    if logg_lower<loggrange[0]: loggs[loggs==logg_lower] = loggrange[0]
    if logg_upper>loggrange[1]: loggs[loggs==logg_upper] = loggrange[1]
    teffs = teffs[~cut]
    loggs = loggs[~cut]

    #-- Generate a grid in logg/teff keeping in mind that this is not a rectangular space
    gridpnts_ = numpy_ext.unique_arr(np.column_stack([teffs,loggs]))

    #-- now we can generate random points:
    sample1 = numpy_ext.random_rectangular_grid(gridpnts_,points)
    if correctTeff: sample1[:,0] = np.array([teffrange[0] for i in sample1[:,0]])
    if correctLogg: sample1[:,1] = np.array([loggrange[0] for i in sample1[:,1]])

    for name in colnames:
        if not name+'range' in ranges: ranges[name+'range'] = (-inf, inf)
    sample2 = np.random.uniform(low =[max(ax.min(),ranges[name+'range'][0]) for ax,name in zip(axis_values,colnames) if not name in ['teff','logg']],\
                                high=[min(ax.max(),ranges[name+'range'][1]) for ax,name in zip(axis_values,colnames) if not name in ['teff','logg']],\
                                size=((len(sample1),len(colnames)-2)))

    colnames.remove('teff')
    colnames.remove('logg')
    #-- return grid and column names
    out_dict_ = {}
    for col,name in zip(np.column_stack([sample1,sample2]).T,['teff','logg']+colnames):
        out_dict_[name] = col

    #-- Check if all collumns that were provided are also returned
    out_dict = {}
    for name in parameters:
        if name in out_dict_:
            out_dict[name] = out_dict_[name]
        else:
            out_dict[name] = np.array([ranges[name+'range'][0]
                                      for i in out_dict_['teff']])
    return out_dict


def generate_grid_pix(photbands, points=None, clear_memory=False,**kwargs):
    """
    Generate a grid of parameters for 1 or more stars. Based on the generate_grid_single_pix
    method. The radius of the components is based on the masses if given, otherwise on the
    radrange arguments. If masses are given the radrange arguments are ignored, meaning that
    the returned radius can be outside those limits. If only 1 component is provided no radius
    will be returned.

    This function can handle multiple components in the same way as a single component.
    parameter ranges are provided as <parname><component>range=(...,...). fx:
    teffrange = (5000, 10000), loggrange = (3.5, 4.5), teff2range = (10000, 20000),
    logg2range = (5.5, 6.0)

    For the first (or only) component both no component number (teffrange) or 1 as component
    number (teff1range) can be used.

    Returns a dictionary with for each parameter that has a range provided, an array of values.
    In the case of multiple components, the radius will be returned even is no radius ranges
    are provided.
    """

    #-- Find all ranges and the number of components
    radiusrange = []
    ranges, parameters, components = {}, set(), set()
    for key in list(kwargs.keys()):
        if re.search('range$', key) and not re.search('^rad\d?',key):
            ranges[key] = kwargs.pop(key)
            name, comp = re.findall('(.*?)(\d?)range$', key)[0]
            parameters.add(name)
            components.add(comp)
        elif re.search('range$', key) and re.search('^rad\d?',key):
            radiusrange.append(kwargs.pop(key))

    #-- If only one component we can directly return the grid
    if len(components) == 1:
        kwargs.update(ranges)
        return generate_grid_single_pix(photbands, points=points, clear_memory=clear_memory, **kwargs)

    #-- For each component get the grid from grid_single_pix
    pars, npoints = {}, +inf
    for i, (comp, grid) in enumerate(zip(components, model.defaults_multiple)):
        ranges_ = {}
        for par in parameters:
            ranges_[par+'range'] = ranges[par+comp+'range'] if par+comp+'range' in ranges else ranges[par+'range']

        kwargs.update(ranges_)
        kwargs.update(grid)
        grid_ = generate_grid_single_pix(photbands, points=points, clear_memory=clear_memory, **kwargs)

        #-- prepare a permutation so different blocks are not clustered together
        permutation = np.random.permutation(len(grid_['teff']))

        for key in list(grid_.keys()):
            npoints = min(npoints,len(grid_[key]))
            pars[key+comp] = grid_[key][permutation]

    #-- The generate_grid_single_pix method does not guarantee the number of points.
    #   So force that all arrays have the same length.
    for key in list(pars.keys()):
        pars[key] = pars[key][:npoints]

    #-- Check that ebv, z and rv is the same for each component
    for comp in components:
        if 'ebv' in parameters: pars['ebv'+comp] = pars['ebv']
        if 'z' in parameters: pars['z'+comp] = pars['z']
        if 'rv' in parameters: pars['rv'+comp] = pars['rv']

    #-- Check if we are dealing with a binary or not and set the radii accordingly
    if 'masses' in kwargs and kwargs['masses'] != None:
        #-- The radius of the stars is calculated based on logg and the provided masses
        masses = kwargs['masses']
        G, Msol, Rsol = constants.GG_cgs, constants.Msol_cgs, constants.Rsol_cgs
        for i, comp in enumerate(components):
            pars['rad'+comp] = np.sqrt(G*masses[i]*Msol/10**pars['logg'+comp])/Rsol
    else:
        #-- We have random different radii for the stars
        if radiusrange == []: radiusrange = [(0.1,10) for i in components]
        for i, comp in enumerate(components):
            pars['rad'+comp] = np.random.uniform(low=radiusrange[i][0], high=radiusrange[i][1], size=npoints)
    return pars


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
    #   broader to get points in the entire interval.
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
    keep = (teffrange[0]<=teffs) & (teffs<=teffrange[1]) &\
            (loggrange[0]<=loggs) & (loggs<=loggrange[1]) &\
            (ebvrange[0]<=ebvs) & (ebvs<=ebvrange[1]) &\
            (zrange[0]<=zs) & (zs<=zrange[1])
    teffs,loggs,ebvs,zs = teffs[keep],loggs[keep],ebvs[keep],zs[keep]

    if res:
        teffs,loggs,ebvs,zs = teffs[::res],loggs[::res],ebvs[::res],zs[::res]
    logger.info('Evaluating %d points in parameter space'%(len(teffs)))

    return teffs,loggs,ebvs,zs

def generate_grid(photbands,teffrange=((-inf,inf),(-inf,inf)),
                  loggrange=((-inf,inf),(-inf,inf)),ebvrange=(-inf,inf),
                  zrange=((-inf,inf),(-inf,inf)),
                  radiusrange=((1,1),(0.1,10.)),grids=None,
                  points=None,res=None,clear_memory=False,
                  type='single',primary_hottest=False, **kwargs):
    """
    Generate grid points at which to fit an interpolated grid of SEDs.

    If C{points=None}, the points are chosen on the predefined grid points.
    Otherwise, C{points} grid points will be generated, uniformly distributed
    between the ranges defined by C{teffrange}, C{loggrange} and C{ebvrange}. If
    you set the resolution to C{2}, one out of every two points will be selected.

    Setting C{primary_hottest} makes sure that, when doing a binary fit, that
    the the effective temperature of the first star is hotter than the second.

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
                      zrange=zrange,points=points,clear_memory=clear_memory)
        radii = np.ones(len(teffs))#[1 for i in teffs]
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
    #   this is not obligatory
    if not hasattr(teffrange[0],'__iter__'): teffs = np.column_stack([teffs[:,0]]*len(grids))
    if not hasattr(loggrange[0],'__iter__'): loggs = np.column_stack([loggs[:,0]]*len(grids))
    #if not hasattr(ebvrange[0],'__iter__'): ebvs = np.column_stack([ebvs[:,0]]*len(grids))
    #if not hasattr(zrange[0],'__iter__'): zs = np.column_stack([zs[:,0]]*len(grids))
    ebvs = np.column_stack([ebvs[:,0]]*len(grids))
    zs = np.column_stack([zs[:,0]]*len(grids))

    if type=='binary':
        #-- The radius of the stars is calculated based on logg and the provided masses
        masses = 'masses' in kwargs and  kwargs['masses'] or (1,1)
        G = constants.GG_cgs
        Msol = constants.Msol_cgs
        Rsol = constants.Rsol_cgs
        radius1 = np.sqrt(G*masses[0]*Msol/10**loggs[:,0])/Rsol
        radius2 = np.sqrt(G*masses[1]*Msol/10**loggs[:,1])/Rsol
        #radii = radius2/radius1

        #radii = np.array([np.ones(len(radii)),radii]).T
        radii = np.array([radius1,radius2]).T

        #-- maybe we need to lower the temperatures of the secondary, so that
        #   it is not hotter than the primary effective temperature
        if primary_hottest:
            wrong = teffs[:,1]>teffs[:,0]
            teffs[wrong,1] = np.random.uniform(low=teffs[:,1].min()*np.ones(sum(wrong)),\
                                               high=teffs[wrong,0])
            logger.info('Ensured primary is the hottest (%d/%d corrections)'%(sum(wrong),len(wrong)))


    elif type=='multiple':
        #-- We have random different radii for the stars
        radii = np.array([10**np.random.uniform(low=radiusrange[0][0], high=radiusrange[0][1], size=(len(teffs))),
                    10**np.random.uniform(low=radiusrange[1][0], high=radiusrange[1][1], size=(len(teffs)))]).T
        #radii = 10**np.random.uniform(low=[np.log10(i[0]) for i in radiusrange],
                              #high=[np.log10(i[1]) for i in radiusrange],size=(len(teffs),2))

    return teffs,loggs,ebvs,zs,radii

#{ Fitting: grid search

def igrid_search_pix(meas,e_meas,photbands, constraints={},**kwargs):
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
        luminosities (R=1Rsol)
        @rtype: array
        """
        model_func = kwargs.pop('model_func',model.get_itable_pix)
        stat_func = kwargs.pop('stat_func',stat_chi2)
        colors = np.array([filters.is_color(photband) for photband in photbands],bool)
        #-- run over the grid, retrieve synthetic fluces and compare with
        #   observations.
        syn_flux,lumis = model_func(photbands=photbands,**kwargs)
        chisqs,scales,e_scales = stat_func(meas.reshape(-1,1),\
                                           e_meas.reshape(-1,1),\
                                           colors,syn_flux, **constraints)
        #-- return results
        return chisqs,scales,e_scales,lumis

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
    fitkws = {}
    if 'distance' in kwargs and kwargs['distance'] != None:
        fitkws = {'distance':kwargs['distance']}
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
    for n,pars in enumerate(zip(*args)):
        if index is None: p.update(1)
        syn_flux,Labs = model_func(*pars,photbands=photbands,**kwargs)
        chisqs[n],scales[n],e_scales[n] = stat_func(meas,e_meas,colors,syn_flux, **fitkws)
        lumis[n] = Labs
    #-- return results
    if index is not None:
        return chisqs,scales,e_scales,lumis,index
    else:
        return chisqs,scales,e_scales,lumis

#}

#{ Fitting: minimizer

def create_parameter_dict(**pars):

    #-- Find all the parameters first
    parnames = set()
    atributes = set()
    for key in list(pars.keys()):
        name, att = re.findall("(.*)_([a-zA-Z]+)$", key)[0]
        parnames.add(name)
        atributes.add(att)

    #-- create dictionary with the attributes
    parnames = np.array(list(parnames))
    result = dict(names=parnames)
    for att in atributes:
        result[att] = np.array([None for i in parnames])

    #-- read the attributes
    for key in list(pars.keys()):
        name, att = re.findall("(.*)_([a-zA-Z]+)$", key)[0]
        result[att][parnames == name] = pars[key]

    return result

def _get_info_from_minimizer(minimizers, startpars, photbands, meas, e_meas, **fitkws):
    scales, lumis, chisqrs, nfevs, allpars = [], [], [], [], {}
    for n in fitkws['pnames']:
        allpars[n] = np.array([])
        allpars[n+'start'] = np.array([])
    for mini, start in zip(minimizers, startpars):
        chisqrs.append(mini.chisqr)
        nfevs.append(mini.nfev)

        val, err = mini.model.get_parameters(full_output=False)
        for n, v in zip(fitkws['pnames'], val):
            allpars[n] = np.append(allpars[n], [v])
            allpars[n+'start'] = np.append(allpars[n+'start'], [start[n].value])

        synth, lum = mini.model.evaluate(photbands, **fitkws)
        distance = fitkws['distance'] if 'distance' in fitkws else None
        if distance != None:
            scale = 1/distance**2
        else:
            ratio = (meas/synth[:,0])
            weights = (meas/e_meas)
            scale = np.average(ratio,weights=weights)
        lumis.append(lum[0])
        scales.append(scale)

    return np.array(chisqrs), np.array(nfevs), np.array(scales), np.array(lumis), allpars

def _iminimize_model(varlist, x, *args, **kws):
    pnames = kws.pop('pnames')
    pars = {}
    for n, v in zip(pnames, varlist):
        pars[n] = np.array([v])
    pars.update(kws)
    #print pars
    return model.get_itable_pix(wave_units=None, photbands=x, **pars)

def _iminimize_residuals(synth, meas, weights=None, **kwargs):
    synth = synth[0][:,0] #select the flux.
    e_meas = 1 / weights
    if 'distance' in kwargs:
        scale = 1/kwargs['distance']**2
    else:
        ratio = (meas/synth)
        weights = (meas/e_meas)
        scale = np.average(ratio,weights=weights)
    return (meas - synth*scale)/e_meas

def iminimize(meas,e_meas,photbands, points=None, return_minimizer=False,**kwargs):
    """
    minimizer based on the sigproc.fit lmfit minimizer.
    provide the observed data, the fitting model, residual function and parameter
    information about the variables, and this function will return the best fit
    parameters together with extra information about the fit.

    if the fitkws keyword is supplied, this dict will be made available to the
    model_func (fit model) during the fitting process. The order of the parameters
    will also be made available as the 'pnames' keyword.
    """

    kick_list = kwargs.pop('kick_list', None)
    constraints = kwargs.pop('constraints', dict())
    fitmodel = kwargs.pop('model_func',_iminimize_model)
    residuals = kwargs.pop('res_func',_iminimize_residuals)
    epsfcn = kwargs.pop('epsfcn', 0.0005)# using ~3% step to derive jacobian.
    fitkws_old = kwargs.pop('fitkws', None)
    #-- get the parameters
    parameters = create_parameter_dict(**kwargs)

    #-- setup the fitting model
    pnames = parameters.pop('names')
    fmodel = sfit.Function(function=fitmodel, par_names=pnames)
    fmodel.setup_parameters(**parameters)

    #-- fit the model to the data
    fitkws = dict(pnames=pnames)
    fitkws.update(constraints)
    if points == None:
        startpars = [copy.deepcopy(fmodel.parameters)]
        minimizer = sfit.minimize(photbands,meas, fmodel, weights=1/e_meas, kws=fitkws, \
                                      resfunc=residuals, engine='leastsq', epsfcn=epsfcn,\
                                      ftol=0.001, xtol=0.001)
        minimizer = [minimizer]
    else:
        minimizer, startpars, newmodels, chisqr = sfit.grid_minimize(photbands, meas, fmodel, \
                           weights=1/e_meas, kws=fitkws, resfunc=residuals, engine='leastsq', \
                           epsfcn=epsfcn, points=points, parameters=kick_list, return_all=True)

    if return_minimizer:
        #-- return the actual minimizer used by the calculate ci methods
        return minimizer[0]
    else:
        #-- for all other users return the actual results
        chisqr, nfev, scale, lumis, grid = _get_info_from_minimizer(minimizer, startpars,\
                                                        photbands, meas, e_meas, **fitkws)
        ##-- collect all parameter info
        #val, err, vary, min, max, expr = fmodel.get_parameters(full_output=True)
        #pars = dict(name=pnames, value=val, cilow=min, cihigh=max)

        return grid, chisqr, nfev, scale, lumis

def calculate_iminimize_CI(meas,e_meas,photbands, CI_limit=0.66, **kwargs):
    """
    Calculate the confidence intervalls for every parameter.
    When ci fails, the boundaries are returned as ci.
    """

    maxiter = kwargs.pop('maxiter', 100)

    mini = iminimize(meas,e_meas,photbands, return_minimizer=True, **kwargs)
    val, err, vary, low, high, expr = mini.model.get_parameters(full_output=True)
    pnames = mini.model.par_names

    #-- setup stat function
    def prob_func(Ndata, Nparas, new_chi, best_chi, Nfix=1.):
        k = Ndata-Nparas
        factor = max(best_chi/k,1)
        return scipy.stats.distributions.chi2.cdf(new_chi/factor,k)

    cilow, cihigh = low.copy(), high.copy()
    for i,p in enumerate(pnames):
        try:
            ci = mini.calculate_CI(parameters=[p], short_output=True,\
                    sigma=CI_limit, maxiter=maxiter, prob_func=prob_func)
            logger.debug('Calculated ci for parameter %s: %s'%(p,ci) )
            if ci[0] != None:  cilow[i] = ci[0]
            if ci[1] != None:  cihigh[i] = ci[1]
        except Exception:
            logger.warning('Failed to calculate CI for parameter: %s'%(p))

    return dict(name=pnames, value=val, cilow=cilow, cihigh=cihigh)

def calculate_iminimize_CI2D(meas,e_meas,photbands, xpar, ypar, df=None, **kwargs):
    """
    Calculate a 2D confidence grid for the given parameters.
    You can provide limits on the parameters yourself, if none are
    provided, the function will take care of it, but it might crash
    if the limits are outside the model range.
    returns: x vals, y vals, ci values
    """
    maxiter = kwargs.pop('maxiter', 25)
    res = kwargs.pop('res', (10,10))
    if type(res) == int: res = (res,res)
    limits = kwargs.pop('limits', None)
    citype = kwargs.pop('type', 'prob')

    mini = iminimize(meas,e_meas,photbands, return_minimizer=True, **kwargs)
    val, err, vary, low, high, expr = mini.model.get_parameters(full_output=True)
    val, low, high = np.array(val, dtype=float), np.array(low, dtype=float), np.array(high, dtype=float)

    #-- set some keywords so the process speeds up a bit:
    mini.userkws.update({'epsfcn':0.001, 'xtol':0.005, 'ftol':0.005})


    name = mini.model.par_names
    xval, yval = val[name==xpar], val[name==ypar]

    logger.debug('Starting calculations of CI2D with parameters:\n%s'%\
                               (mini.model.param2str(full_output=True)))

    #-- if not provided choose own limits
    if limits != None:
        xmin, xmax = limits[0]
        ymin, ymax = limits[1]
    else:
        xmin = low[name==xpar] if low[name==xpar] != None else xval - 5*err[name==xpar]
        xmax = high[name==xpar] if high[name==xpar] != None else xval + 5*err[name==xpar]
        ymin = low[name==ypar] if low[name==ypar] != None else yval - 5*err[name==ypar]
        ymax = high[name==ypar] if high[name==ypar] != None else yval + 5*err[name==ypar]

    #-- fix limits to include best fit in CI
    xstep = (xmax - xmin) / float(res[0])
    ystep = (ymax - ymin) / float(res[1])
    xmin = xval - np.floor((xval - xmin) / xstep) * xstep
    xmax = xval + (res[0] - 1 - np.floor((xval - xmin) / xstep)) * xstep
    ymin = yval - np.floor((yval - ymin) / ystep) * ystep
    ymax = yval + (res[1] - 1 - np.floor((yval - ymin) / ystep)) * ystep
    limits = [(xmin, xmax), (ymin, ymax)]

    logger.debug('Recalculated grid limits to: \n\t%s: %s <-> %s \n\t%s %s <-> %s'%(xpar,
                 limits[0][0], limits[0][1], ypar, limits[1][0], limits[1][1]))

    x,y,ci_chi2 = mini.calculate_CI_2D(xpar=xpar, ypar=ypar, res=res, limits=limits,
                                       ctype='chi2')

    #-- do the statistics
    N = mini.ndata
    if df == None: df = mini.nvarys
    k = N-df
    factor = max(mini.chisqr/k,1)
    ci_raw = scipy.stats.distributions.chi2.cdf(ci_chi2,k)
    ci_red = scipy.stats.distributions.chi2.cdf(ci_chi2/factor,k)

    return x, y, ci_chi2, ci_raw, ci_red

def iminimize2(meas,e_meas,photbands,*args,**kwargs):
    model_func = kwargs.pop('model_func',model.get_itable)
    #res_func = kwargs.pop('res_func',residual_single) # not yet defined!
    method = kwargs.pop('fitmethod','fmin')
    stat_func = kwargs.pop('stat_func',stat_chi2)

    colors = np.array([filters.is_color(photband) for photband in photbands],bool)

    # the help function which returns the chisquare  NOTE: metallicity is not yet included in the fitting!
    def residual_single(parameters):
        syn_flux,Labs = model_func(*parameters,photbands=photbands,**kwargs)
        # in case any of the parameters goes out of its bounds
        #if isnan(syn_flux).any():
        chisq,scale,e_scale = stat_func(meas,e_meas,colors,syn_flux,full_output=False)
        return chisq
    res_func = residual_single
    # calling the fitting function NOTE: the initial metallicity is returned in the output!

    if method=='fmin': # fmin
        optpars,fopt,niter,funcalls,warnflag = fmin(res_func,np.array(args),xtol=0.0001,disp=0,full_output=True)
    elif method=='fmin_powell': #fmin_powell
        optpars = fmin_powell(res_func,np.array(args))
    else:
        raise NotImplementedError
    logger.debug("Optimization finished")
    syn_flux,Labs = model_func(*optpars,photbands=photbands,**kwargs)

    # when any of the parameters goes out of the bounds of the grid, syn_flux contains NaN
    if np.isnan(syn_flux).any():
        warnflag = 3

    stats = stat_func(meas,e_meas,colors,syn_flux,full_output=False)
    optpars = np.hstack([optpars,stats])
    # stats: chisq, scale, e_scale
    return optpars,warnflag
#}


if __name__=="__main__":
    from ivs.aux import loggers
    import time
    import pylab as plt
    logger = loggers.get_basic_logger(clevel='DEBUG')

    photbands = ['GENEVA.G','GENEVA.B-V']

    c0 = time.time()
    teffs,loggs,ebvs,zs,radii = generate_grid(photbands,teffrange=(5000,5800),loggrange=(4.20,4.70),zrange=(0,0),ebvrange=(0.05,0.08), grid='kurucz',points=10000)
    print('Time: %i'%(time.time()-c0))

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
    from ivs.aux import loggers
    from pylab import *
    from numpy import *
    logger = loggers.get_basic_logger()
    random.seed(1111)
    A,grid = get_PCA_grid(['GENEVA.U-B','GENEVA.B1-B','GENEVA.B2-B','GENEVA.V-B','GENEVA.V1-B','GENEVA.G-B','2MASS.J-H','2MASS.KS-H'],ebvrange=(0,0.5),res=10)
    P,T,(means,stds) = get_PCA(A)
    calib = calibrate_PCA(T,grid,function='linear')
    sample_index = int(np.random.uniform(high=len(A)))
    sample = 10**A[sample_index]
    print([bla[sample_index] for bla in grid])
    print(get_PCA_parameters(sample,calib,P,means,stds))
