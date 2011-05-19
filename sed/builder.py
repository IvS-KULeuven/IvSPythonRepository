# -*- coding: utf-8 -*-
"""
SED builder program.

Keep this preferentially out of the IVS repository: this is more of a program
that a module!
"""
import sys
import time
import os
import logging

import pylab as pl
from matplotlib import mlab
import Image
import numpy as np
import scipy.stats
from scipy.interpolate import Rbf

from ivs import config
from ivs.misc import numpy_ext
from ivs.misc.decorators import memoized
from ivs.io import ascii
from ivs.io import fits
from ivs.sed import model
from ivs.sed import filters
from ivs.sed import fit
from ivs.sed import distance
from ivs.sed import extinctionmodels
from ivs.catalogs import crossmatch
from ivs.catalogs import vizier
from ivs.catalogs import mast
from ivs.catalogs import sesame
from ivs.units import conversions
from ivs.units import constants


logger = logging.getLogger("SED.BUILD")
#logger.setLevel(10)

def fix_master(master,e_default=None):
    """
    Clean/extend/fix record array received from C{get_photometry}.
    
    This function does a couple of things:
    
        1. Adds common but uncatalogized colors like 2MASS.J-H if not already
           present. WARNING: these colors can mix values from the same system
           but from different catalogs!
        2. Removes photometry for which no calibration is available
        3. Adds a column 'color' with flag False denoting absolute flux measurement
        and True denoting color.
        4. Adds a column 'include' with flag True meaning the value will be
        included in the fit
        5. Sets some default errors to photometric values, for which we know
        that the catalog values are not trustworthy.
        6. Sets a lower limit to the allowed errors of 1%. Errors below this
        value are untrostworthy because the calibration error is larger than that.
        7. USNOB1 and ANS photometry are set the have a minimum error of 30%
        due to uncertainties in response curves.
        
    To do BUG FIX: only make colors from measurements from the same source
    catalog!
    
    @param master: record array containing photometry. This should have the fields
    'meas','e_meas','unit','photband','source','_r','_RAJ2000','DEJ2000',
    'cmeas','e_cmeas','cwave','cunit'
    @type master: record array
    @param e_default: default error for measurements without errors
    @type e_default: float
    @return: record array extended with fields 'include' and 'color', and with
    rows added (see above description)
    @rtype: numpy recard array
    """
    #-- we recognize uncalibrated stuff as those for which no absolute flux was
    #   obtained, and remove them:
    master = master[-np.isnan(master['cmeas'])]
    
    #-- add common uncatalogized colors:
    columns = list(master.dtype.names)
    for color in ['2MASS.J-H','2MASS.KS-H','TD1.1565-1965','TD1.2365-1965','TD1.2365-2740',
                  'JOHNSON.J-H','JOHNSON.K-H','JOHNSON.B-V','JOHNSON.U-B','GALEX.FUV-NUV',
                  'TYCHO2.BT-VT','ANS.15N-15W','ANS.15W-18','ANS.18-22','ANS.22-25','ANS.25-33']:
        #-- get the filter system and the separate bands for these colors
        system,band = color.split('.')
        band0,band1 = band.split('-')
        band0,band1 = '%s.%s'%(system,band0),'%s.%s'%(system,band1)
        
        #-- if the separate bands are available (and not the color itself),
        #   calculate the colors here
        if band0 in master['photband'] and band1 in master['photband'] and not color in master['photband']:
            #-- where are the bands located?
            index0 = list(master['photband']).index(band0)
            index1 = list(master['photband']).index(band1)
            
            #-- start a new row to add the color
            row = list(master[index1])
            row[columns.index('photband')] = color
            row[columns.index('cwave')] = np.nan
            
            #-- it could be a magnitude difference
            if master['unit'][index0]=='mag':
                row[columns.index('meas')] = master['meas'][index0]-master['meas'][index1]
                row[columns.index('e_meas')] = np.sqrt(master['e_meas'][index0]**2+master['e_meas'][index1]**2)
                #-- error will not always be available...
                try:
                    row[columns.index('cmeas')],row[columns.index('e_cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],row[columns.index('e_meas')],photband=color)
                except AssertionError:
                    row[columns.index('cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],photband=color)
                    row[columns.index('e_cmeas')] = np.nan
            #-- or it could be a flux ratio
            else:
                row[columns.index('meas')] = master['meas'][index0]/master['meas'][index1]
                row[columns.index('e_meas')] = np.sqrt(((master['e_meas'][index0]/master['meas'][index0])**2+(master['e_meas'][index1]/master['meas'][index1])**2)*row[columns.index('meas')]**2)
                row[columns.index('cmeas')],row[columns.index('e_cmeas')] = row[columns.index('meas')],row[columns.index('e_meas')]
            row[columns.index('cunit')] = 'flux_ratio'
            master = numpy_ext.recarr_addrows(master,[tuple(row)])
            
    #-- add an extra column with a flag to distinguish colors from absolute
    #   fluxes, and a column with flags to include/exclude photometry
    #   By default, exclude photometry if the effective wavelength is above
    #   150 mum or if the measurement is nan
    extra_cols = [[filters.is_color(photband) for photband in master['photband']],
                  [(cwave<150e4) or np.isnan(meas) for cwave,meas in zip(master['cwave'],master['cmeas'])]]
    dtypes = [('color',np.bool),('include',np.bool)]
    master = numpy_ext.recarr_addcols(master,extra_cols,dtypes)
    
    #-- set default errors if no errors are available and set really really
    #   small errors to some larger default value
    if e_default is not None:
        no_error = np.isnan(master['e_cmeas'])
        master['e_cmeas'][no_error] = np.abs(e_default*master['cmeas'][no_error])
        small_error = master['e_cmeas']<(0.01*np.abs(master['cmeas']))
        master['e_cmeas'][small_error] = 0.01*np.abs(master['cmeas'][small_error])
    
    #-- the measurements from USNOB1 are not that good because the response
    #   curve is approximated by the JOHNSON filter. Set the errors to min 30%
    #   The same holds for ANS: here, the transmission curves are very uncertain
    for i,photband in enumerate(master['photband']):
        if 'USNOB1' in photband or 'ANS' in photband:
            master['e_cmeas'][i] = max(0.30*master['cmeas'][i],master['e_cmeas'][i])
    
    return master


def decide_phot(master,names=None,wrange=None,ptype='all',include=False):
    """
    Exclude/include photometric passbands containing one of the strings listed in
    photbands.
    
    This function will set the flags in the 'include' field of the extended
    master record array.
    
    Colours also have a wavelength in this function, and it is defined as the
    average wavelength of the included bands (doesn't work for stromgren C1 and M1).
    
    include=False excludes photometry based on the input parameters.
    include=True includes photometry based on the input parameters.
    
    Some examples:
    
    1. Exclude all measurements:
    >>> decide_phot(master,wrange=(-np.inf,+np.inf),ptype='all',include=False)
    
    2. Include all TD1 fluxes and colors:
    >>> decide_phot(master,names=['TD1'],ptype='all',include=True)
    
    3. Include all V band measurements from all systems (but not the colors):
    >>> decide_phot(master,names=['.V'],ptype='abs',include=True)
    
    4. Include all Geneva colors and exclude Geneva magnitudes:
    >>> decide_phot(master,names=['GENEVA'],ptype='col',include=True)
    >>> decide_phot(master,names=['GENEVA'],ptype='abs',include=False)
    
    5. Exclude all infrared measurements beyond 1 micron:
    >>> decide_phot(master,wrange=(1e4,np.inf),ptype='all',include=False)
    
    6. Include all AKARI measurements below 10 micron:
    >>> decide_phot(master,names=['AKARI'],wrange=(-np.inf,1e5),ptype='all',include=True)
    
    @param master: record array containing all photometry
    @type master: numpy record array
    @param names: strings excerpts to match filters
    @type names: list of strings
    @param wrange: wavelength range (most likely angstrom) to include/exclude
    @type wrange: 2-tuple (start wavelength,end wavelength)
    @param ptype: type of photometry to include/exclude: absolute values, colors
    or both
    @type ptype: string, one of 'abs','col','all'
    @param include: flag setting exclusion or inclusion
    @type include: boolean
    @return: master record array with adapted 'include' flags
    @rtype: numpy record array
    """ 
    #-- exclude/include passbands based on their names
    if names is not None:
        for index,photband in enumerate(master['photband']):
            for name in names:
                if name in photband:
                    if ptype=='all' or (ptype=='abs' and -master['color'][index]) or (ptype=='col' and master['color'][index]):
                        master['include'][index] = include
                        break
    #-- exclude/include colors based on their wavelength
    if wrange is not None:
        for index,photband in enumerate(master['photband']):
            system,color = photband.split('.')
            if not '-' in color or ptype=='abs':
                continue
            band0,band1 = color.split('-')
            cwave0 = filters.eff_wave('%s.%s'%(system,band0))
            cwave1 = filters.eff_wave('%s.%s'%(system,band1))
            if (wrange[0]<cwave0<wrange[1]) or (wrange[0]<cwave0<wrange[1]):
                master['include'][index] = include
        #-- exclude/include passbands based on their wavelength
        if not ptype=='col':
            master['include'][(wrange[0]<master['cwave']) & (master['cwave']<wrange[1])] = include    


def photometry2str(master):
    """
    String representation of master record array
    
    @param master: master record array containing photometry
    @type master: numpy record array
    """
    master = master[np.argsort(master['photband'])]
    print '%20s %12s %12s %12s %10s %12s %12s %11s %s'%('PHOTBAND','MEAS','E_MEAS','UNIT','CWAVE','CMEAS','E_CMEAS','UNIT','SOURCE')
    print '=========================================================================================================================='
    for i,j,k,l,m,n,o,p in zip(master['photband'],master['meas'],master['e_meas'],master['unit'],master['cwave'],master['cmeas'],master['e_cmeas'],master['source']):
        print '%20s %12g %12g %12s %10.0f %12g %12g erg/s/cm2/A %s'%(i,j,k,l,m,n,o,p)
    

@memoized
def get_parallax(ID):
    """
    Retrieve a star's parallax and galactic coordinates in degrees.
    
    @param ID: target's name
    @type ID: string
    @return: parallax and galactic coordinates
    @rtype: (plx,e_plx),(long,lat)
    """
    data,units,comms = vizier.search('I/311/hip2',ID=ID)
    if data is None or not len(data):
        logger.warning('No parallax found')
        plx = None
    else:
        plx = (data['Plx'][0],data['e_Plx'][0])
    info = sesame.search(ID)
    ra,dec = info['jpos'].split()
    gal = conversions.convert('equ','gal',(str(ra),str(dec)),epoch='2000')
    gal = float(gal[0])/np.pi*180,float(gal[1])/np.pi*180
    return plx,gal

@memoized
def get_schaller_grid():
    """
    Download Schaller 1992 evolutionary tracks and return an Rbf interpolation
    function.
    
    @return: Rbf interpolation function
    @rtype: Rbf interpolation function
    """
    #-- translation between table names and masses
    masses = [1,1.25,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,40,60][:-1]
    tables = ['table20','table18','table17','table16','table15','table14',
              'table13','table12','table11','table10','table9','table8',
              'table7','table6','table5','table4','table3'][:-1]
    #-- read in all the tables and compute radii and luminosities.
    all_teffs,all_loggs,all_radii = [],[],[]
    for mass,table in zip(masses,tables):
        data,comms,units = vizier.search('J/A+AS/96/269/%s'%(table))
        all_teffs.append(10**data['logTe'])
        all_radii.append(np.sqrt((10**data['logL']*constants.Lsol_cgs)/(10**data['logTe'])**4/(4*np.pi*constants.sigma_cgs)))
        all_loggs.append(np.log10(constants.GG_cgs*mass*constants.Msol_cgs/(all_radii[-1]**2)))
    all_teffs = np.hstack(all_teffs)
    all_radii = np.hstack(all_radii)/constants.Rsol_cgs
    all_loggs = np.hstack(all_loggs)
    #-- remove low temperature models, the evolutionary tracks are hard to
    #   interpolate there.
    keep = all_teffs>5000
    all_teffs = all_teffs[keep]
    all_radii = all_radii[keep]
    all_loggs = all_loggs[keep]
    #-- make linear interpolation model between all modelpoints
    mygrid = Rbf(np.log10(all_teffs),all_loggs,all_radii,function='linear')
    logger.info('Interpolation of Schaller 1992 evolutionary tracks to compute radii')
    return mygrid


def get_radii(teffs,loggs):
    """
    Retrieve radii from stellar evolutionary tracks from Schaller 1992.
    
    @param teffs: model effective temperatures
    @type teffs: numpy array
    @param loggs: model surface gravities
    @type loggs: numpy array
    @return: model radii (solar units)
    @rtype: numpy array
    """
    mygrid = get_schaller_grid()
    radii = mygrid(np.log10(teffs),loggs)
    return radii


def calculate_distance(plx,gal,teffs,loggs,scales,n=75000):
    """
    Calculate distances and radii of a target given its parallax and location
    in the galaxy.
    
    
    """
    #-- compute distance up to 25 kpc, which is about the maximum distance from
    #   earth to the farthest side of the Milky Way galaxy
    #   rescale to set the maximum to 1
    d = np.logspace(np.log10(0.1),np.log10(25000),100000)
    if plx is not None:
        dprob = distance.distprob(d,gal[1],plx)
        dprob = dprob / dprob.max()
    else:
        dprob = np.ones_like(d)
    #-- compute the radii for the computed models, and convert to parsec
    radii = get_radii(teffs[-n:],loggs[-n:])
    radii = conversions.convert('Rsol','pc',radii)
    d_models = radii/np.sqrt(scales[-n:])
    #-- we set out of boundary values to zero
    if plx is not None:
        dprob_models = np.interp(d_models,d[-n:],dprob[-n:],left=0,right=0)
    else:
        dprob_models = np.ones_like(d_models)
    #-- reset the radii to solar units for return value
    radii = conversions.convert('pc','Rsol',radii)
    return (d_models,dprob_models,radii),(d,dprob)
    


class SED(object):

    def __init__(self,ID,photfile=None,plx=None):
        """
        Initialize SED class.
        
        @param plx: parallax (and error) of the object
        @type plx: tuple (plx,e_plx)
        """
        self.ID = ID
        #-- the file containing photometry should have the following name. We
        #   save photometry to a file to avoid having to download photometry
        #   each time from the internet
        if photfile is None:
            self.photfile = '%s.phot'%(ID).replace(' ','_')
        else:
            self.photfile = photfile
        
        #-- keep information on the star from SESAME, but override parallax
        #   and set new galactic coordinates
        self.info = sesame.search(ID)
        plx_,gal_ = get_parallax(ID)
        if plx is not None:
            self.info['plx'] = plx
        elif plx_ is not None:
            self.info['plx'] = plx_
        self.info['gal'] = gal_
        
        #-- prepare for information on fitting processes
        self.results = {}
    
    #{ Handling photometric data
    def get_photometry(self,radius=None):
        """
        Search photometry on the net or from the phot file if it exists.
        
        For bright stars, you can set radius a bit higher...
        
        @param radius: search radius (arcseconds)
        @type radius: float.
        """
        if radius is None:
            if 'mag.V.v' in self.info and self.info['mag.V.v']<6.:
                radius = 60.
            else:
                radius = 10.
        if not os.path.isfile(self.photfile):
            #-- get and fix photometry. Set default errors to 1%, and set
            #   USNOB1 errors to 3%
            master = crossmatch.get_photometry(ID=self.ID,radius=radius,extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            master['_RAJ2000'] -= self.info['jradeg']
            master['_DEJ2000'] -= self.info['jdedeg']
            
            #-- fix the photometry: set default errors to 2% and print it to the
            #   screen
            self.master = fix_master(master,e_default=0.1)
            photometry2str(master)
            
            #-- write to file
            self.save_photometry()
        else:
            self.load_photometry()
    
    def exclude(self,names=None,wrange=None):
        """
        Exclude photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=False,ptype='all')
    
    def exclude_colors(self,names=None,wrange=None):
        """
        Exclude photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=False,ptype='col')
    
    def exclude_abs(self,names=None,wrange=None):
        """
        Exclude photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=False,ptype='abs')
    
    
    def include(self,names=None,wrange=None):
        """
        Include photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=True,ptype='all')
    
    def include_colors(self,names=None,wrange=None):
        """
        Include photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=True,ptype='col')
    
    def include_abs(self,names=None,wrange=None):
        """
        Include photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=True,ptype='abs')
    
    def set_photometry_scheme(self,scheme):
        """
        Scheme = 'abs': means excluding all colors, including all absolute values
        scheme = 'color': means including all colors, excluding all absolute values
        scheme = 'combo': means inculding all colors, and one absolute value per
          system (the one with the smallest relative error)
        """
        if 'abs' in scheme.lower():
            self.master['include'][self.master['color']] = False
            self.master['include'][-self.master['color']] = True
            logger.info('Fitting procedure will use only absolute fluxes (%d)'%(sum(self.master['include'])))
        elif 'col' in scheme.lower():
            self.master['include'][self.master['color']] = True
            self.master['include'][-self.master['color']] = False
            logger.info('Fitting procedure will use only colors (%d)'%(sum(self.master['include'])))
        elif 'com' in scheme.lower():
            self.master['include'][self.master['color']] = True
            self.master['include'][-self.master['color']] = False
            systems = np.array([photband.split('.')[0] for photband in self.master['photband']])
            set_systems = sorted(list(set(systems)))
            for isys in set_systems:
                keep = -self.master['color'] & (isys==systems)
                if not sum(keep): continue
                index = np.argmax(self.master['e_cmeas'][keep]/self.master['cmeas'][keep])
                index = np.arange(len(self.master))[keep][index]
                self.master['include'][index] = True
            logger.info('Fitting procedure will use colors + one absolute flux for each system (%d)'%(sum(self.master['include'])))
            
        
    
    def save_photometry(self,photfile=None):
        #-- write to file
        if photfile is not None:
            self.photfile = photfile
        logger.info('Save photometry to file %s'%(self.photfile))
        ascii.write_array(self.master,self.photfile,header=True,auto_width=True,use_float='%g')
    
    def load_photometry(self):
        logger.info('Load photometry from file %s'%(self.photfile))
        self.master = ascii.read2recarray(self.photfile)
    
    def add_photometry_fromarrays(self,meas,e_meas,units,photbands,source):
        """
        Add photometry from a file or arrays
        
        By default unflagged, and at the ra and dec of the star. No color and
        included.
        """
        extra_rows = []
        for m,e_m,u,p,s in zip(meas,e_meas,units,photbands,source):
            cm,e_cm = conversions.convert(u,self.master['cunit'][0],m,e_m,photband=p)
            extra_rows.append((m,e_m,np.nan,u,p,s,0.,self.info['jradeg'],info['jdedeg'],eff_wave,cm,e_cm,self.master['cunit'][0],0,1))
        self.master = numpy_ext.recarr_addrows(self.master,extra_rows)

    #}
    
    #{ Fitting routines
    def igrid_search(self,teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf),
                          ebvrange=(-np.inf,np.inf),zrange=(-np.inf,np.inf),
                          threads='max',iterations=3,increase=1,speed=1,res=1):
        #-- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are include in the fitting process:')
        photometry2str(self.master[include_grid])
        grid_results = fit.igrid_search(self.master['cmeas'][include_grid],
                                        self.master['e_cmeas'][include_grid],
                                        self.master['photband'][include_grid],
                                    teffrange=teffrange,loggrange=loggrange,
                                    ebvrange=ebvrange,zrange=zrange,
                                    threads=threads,iterations=iterations,
                                    increase=increase,speed=speed,res=res)
        
        #-- exclude failures
        grid_results = grid_results[-np.isnan(grid_results['chisq'])]
        
        #-- inverse sort according to chisq: this means the best models are at
        #   the end (mainly for plotting reasons, so that the best models
        #   are on top).
        sa = np.argsort(grid_results['chisq'])[::-1]
        grid_results = grid_results[sa]
        
        #-- do the statistics
        #   degrees of freedom: teff,logg,E(B-V),theta,Z
        N,df = sum(include_grid),5
        k = N-df
        if k<=0:
            logger.warning('Not enough data to compute CHI2: it will not make sense')
            k = 1
        #   rescale if needed and compute confidence intervals
        factor = max(grid_results['chisq'][-1]/k,1)
        logger.warning('CHI2 rescaling factor equals %g'%(factor))
        CI_raw = scipy.stats.distributions.chi2.cdf(grid_results['chisq'],k)
        CI_red = scipy.stats.distributions.chi2.cdf(grid_results['chisq']/factor,k)
        
        #-- add the results to the record array and to the results dictionary
        grid_results = mlab.rec_append_fields(grid_results, 'CI_raw', CI_raw)
        grid_results = mlab.rec_append_fields(grid_results, 'CI_red', CI_red)
        if not 'igrid_search' in self.results:
            self.results['igrid_search'] = {}
        
        self.results['igrid_search']['grid'] = grid_results
        self.results['igrid_search']['factor'] = factor
        self.results['igrid_search']['CI'] = {}
        
        start_CI = np.argmin(np.abs(CI_red-0.95))
        for name in grid_results.dtype.names:
            self.results['igrid_search']['CI'][name+'L'] = grid_results[name][start_CI:].min()
            self.results['igrid_search']['CI'][name] = grid_results[name][-1]
            self.results['igrid_search']['CI'][name+'U'] = grid_results[name][start_CI:].max()
            logger.info('95%% CI %s: %g <= %g <= %g'%(name,self.results['igrid_search']['CI'][name+'L'],
                                                           self.results['igrid_search']['CI'][name],
                                                           self.results['igrid_search']['CI'][name+'U']))
        
        self.set_best_model()
        self.set_distance()
    #}
    
    #{ Interfaces
    
    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004'):
        """
        Get reddenend and unreddened model
        """
        #-- get reddened and unreddened model
        scale = self.results[mtype]['CI']['scale']
        wave,flux = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                    logg=self.results[mtype]['CI']['logg'],
                                    ebv=self.results[mtype]['CI']['ebv'],
                                    law=law)
        wave_ur,flux_ur = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                          logg=self.results[mtype]['CI']['logg'],
                                          ebv=0,
                                          law=law)
        flux,flux_ur = flux*scale,flux_ur*scale
        
        #-- synthetic flux
        include = self.master['include']
        synflux = np.zeros(len(self.master['photband']))
        keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
        synflux_,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
                                  logg=self.results[mtype]['CI']['logg'],
                                  ebv=self.results[mtype]['CI']['ebv'],
                                  photbands=self.master['photband'][keep])
        synflux[keep] = synflux_
        
        #synflux,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
        #                          logg=self.results[mtype]['CI']['logg'],
        #                          ebv=self.results[mtype]['CI']['ebv'],
        #                          photbands=self.master['photband'])
        synflux[-self.master['color']] *= scale
        chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
        #-- calculate effective wavelengths of the photometric bands via the model
        #   values
        eff_waves = filters.eff_wave(self.master['photband'],model=(wave,flux))
        self.results['model'] = wave,flux,flux_ur
        self.results['synflux'] = eff_waves,synflux,self.master['photband']
        self.results['chi2'] = chi2
    
    def set_distance(self):
        #-- calculate the distance
        cutlogg = (self.results['igrid_search']['grid']['logg']<=4.4) & (self.results['igrid_search']['grid']['CI_red']<=0.95)
        gal = self.info['gal']
        plx = self.info.get('plx',None)
        (d_models,dprob_models,radii),(d,dprob)\
                   = calculate_distance(plx,self.info['gal'],self.results['igrid_search']['grid']['teff'][cutlogg],\
                                                   self.results['igrid_search']['grid']['logg'][cutlogg],\
                                                   self.results['igrid_search']['grid']['scale'][cutlogg])
        #-- calculate model extinction
        res = 100
        self.results['igrid_search']['drimmel'] = np.ravel(np.array([extinctionmodels.findext(gal[0], gal[1], model='drimmel', distance=myd) for myd in d[::res]]))
        self.results['igrid_search']['marshall'] = np.ravel(np.array([extinctionmodels.findext(gal[0], gal[1], model='marshall', distance=myd) for myd in d[::res]]))
        self.results['igrid_search']['d_mod'] = (d_models,dprob_models,radii)
        self.results['igrid_search']['d'] = (d,dprob)
    
    #}
    
    #{ Input and output
    def plot_grid(self,ptype='CI_red',mtype='igrid_search'):
        """
        Plot grid as scatter plot
        
        Possible plot types: 'CI_red','z','ebv'
        """
        region = self.results[mtype]['grid']['CI_red']<0.95
        #-- get the colors and the color scale
        colors = self.results[mtype]['grid'][ptype][region]
        vmin,vmax = colors.min(),colors.max()
        if 'CI' in ptype:
            colors *= 100.
            vmax = 95.
            vmin = colors.min()
        #-- grid scatter plot
        pl.scatter(self.results[mtype]['grid']['teff'][region],
                   self.results[mtype]['grid']['logg'][region],
             c=colors,edgecolors='none',cmap=pl.cm.spectral,vmin=vmin,vmax=vmax)
        #-- mark best value
        pl.plot(self.results[mtype]['grid']['teff'][-1],
                self.results[mtype]['grid']['logg'][-1],'r+',ms=40,mew=3)
        #-- set the limits to only include the 95 interval
        pl.xlim(self.results[mtype]['grid']['teff'][region].max(),
                self.results[mtype]['grid']['teff'][region].min())
        pl.ylim(self.results[mtype]['grid']['logg'][region].max(),
                self.results[mtype]['grid']['logg'][region].min())
        cbar = pl.colorbar()
        pl.xlabel('log (effective temperature [K]) [dex]')
        pl.ylabel(r'log (surface gravity [cm s$^{-2}$]) [dex]')
        #-- set the colorbar label
        if 'CI'in ptype:
            cbar.set_label('Probability [%]')
        elif ptype=='z':
            cbar.set_label('Metallicity Z ($Z_\odot$)')
        elif ptype=='ebv':
            cbar.set_label('E(B-V) [mag]')
        logger.info('Plotted teff-logg diagram of %s'%(ptype))
            
    def plot_sed(self,colors=False,mtype='igrid_search'):
        
        
        def plot_sed_getcolors(master,color_dict=None):
            myphotbands = [iphotb.split('.')[1] for iphotb in master['photband'][master['color']]]
            if color_dict is None:
                color_dict = {myphotbands[0]:0}
                for mycol in myphotbands[1:]:
                    if not mycol in color_dict:
                        max_int = max([color_dict[entry] for entry in color_dict])
                        color_dict[mycol] = max_int+1
            x = [color_dict[mycol] for mycol in myphotbands]
            y = master['cmeas']
            e_y = master['e_cmeas']
            return x,y,e_y,color_dict
            
        
        x__,y__,e_y__,color_dict = plot_sed_getcolors(self.master)
        
        systems = np.array([system.split('.')[0] for system in self.master['photband']],str)
        set_systems = sorted(list(set(systems)))
        color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))]
        pl.gca().set_color_cycle(color_cycle)
        
        #-- for plotting reasons, we translate every color to an integer
        
        for system in set_systems:
            keep = (systems==system) & (self.master['color']==colors)
            #-- synthetic:
            pl.gca()._get_lines.count += 1
            if sum(keep):
                pl.gca()._get_lines.count -= 1
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                    y = self.results['synflux'][1][keep]
                else:
                    x = self.results['synflux'][0][keep]
                    y = x*self.results['synflux'][1][keep]
                pl.plot(x,y,'x',ms=10,mew=2,alpha=0.75)
            #-- include:
            keep = (systems==system) & (self.master['color']==colors) & self.master['include']
            if sum(keep):
                pl.gca()._get_lines.count -= 1
                if colors:
                    #-- translate every color to an integer
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results['synflux'][0][keep]
                    y = x*self.master['cmeas'][keep]
                    e_y = x*self.master['e_cmeas'][keep]
                p = pl.errorbar(x,y,yerr=e_y,fmt='o',label=system,
                                capsize=10,ms=7)
            
            #-- exclude:
            label = np.any(keep) and '_nolegend_' or system
            keep = (systems==system) & (self.master['color']==colors) & -self.master['include']
            if sum(keep):
                pl.gca()._get_lines.count -= 1
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results['synflux'][0][keep]
                    y = x*self.master['cmeas'][keep]
                    e_y = x*self.master['e_cmeas'][keep]
                    
                pl.errorbar(x,y,yerr=e_y,fmt='o',label=label,
                            capsize=10,ms=7,mew=2,mfc='1',mec=color_cycle[pl.gca()._get_lines.count%len(color_cycle)])
        
        #-- only set logarithmic scale if absolute fluxes are plotted
        #   and only plot the real model then
        if not colors:
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            pl.gca().set_autoscale_on(False)
        
            #-- the model
            wave,flux,flux_ur = self.results['model']
            pl.plot(wave,wave*flux,'r-')
            pl.plot(wave,wave*flux_ur,'k-')
            pl.ylabel(r'$\lambda F_\lambda$ [erg/s/cm$^2$]')
            pl.xlabel('wavelength [$\AA$]')
        else:
            xlabels = color_dict.keys()
            xticks = [color_dict[key] for key in xlabels]
            pl.xticks(xticks,xlabels,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Color name')
            pl.xlim(min(xticks)-0.5,max(xticks)+0.5)
        
        pl.grid()
        if colors:
            pl.legend(loc='best',prop=dict(size='x-small'))
        else:
            pl.legend(loc='upper right',prop=dict(size='x-small'))
        angdiam = 0.8#2*np.arctan(np.sqrt(scale))/np.pi*180*3600*1000
        loc = (0.05,0.05)
        teff = self.results[mtype]['grid']['teff'][-1]
        logg = self.results[mtype]['grid']['logg'][-1]
        ebv = self.results[mtype]['grid']['ebv'][-1]
        pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\n$\Theta$=%.3f'%(teff,logg,ebv,0),
                    loc,xycoords='axes fraction')
        logger.info('Plotted SED as %s'%(colors and 'colors' or 'absolute fluxes'))
        
    def plot_chi2(self,colors=False,mtype='igrid_search'):
        
        include_grid = self.master['include']
        eff_waves,synflux,photbands = self.results['synflux']
        chi2 = self.results['chi2']
        systems = np.array([system.split('.')[0] for system in self.master['photband'][include_grid]],str)
        set_systems = sorted(list(set(systems)))
        pl.gca().set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0, 1.0, len(set_systems))])
        for system in set_systems:
            pl.gca()._get_lines.count += 1
            keep = systems==system
            if sum(keep) and not colors:
                pl.gca()._get_lines.count -= 1
                pl.loglog(eff_waves[include_grid][keep],chi2[include_grid][keep],'o',label=system)
            elif sum(keep) and colors:
                pl.gca()._get_lines.count -= 1
                pl.semilogy(range(len(eff_waves[include_grid][keep])),chi2[include_grid][keep],'o',label=system)
        pl.legend(loc='upper right',prop=dict(size='x-small'))
        pl.grid()
        k = 5.
        pl.annotate('Exp. val. = %d'%(k),(0.69,0.12),xycoords='axes fraction',color='b')
        pl.annotate('Total $\chi^2$ = %.1f'%(self.results[mtype]['grid']['chisq'][-1]),(0.69,0.065),xycoords='axes fraction',color='r')
        pl.annotate('Error scale = %.2f'%(np.sqrt(self.results[mtype]['factor'])),(0.69,0.03),xycoords='axes fraction',color='k')
        xlims = pl.xlim()
        pl.plot(xlims,[k,k],'b-',lw=2)
        pl.plot(xlims,[self.results[mtype]['grid']['chisq'][-1],self.results[mtype]['grid']['chisq'][-1]],'r-',lw=2)
        pl.xlim(xlims)
        pl.xlabel('wavelength [$\AA$]')
        pl.ylabel(r'Statistic ($\chi^2$)')
        logger.info('Plotted CHI2 of %s'%(colors and 'colors' or 'absolute fluxes'))
        
    def plot_distance(self,mtype='igrid_search'):
        #-- necessary information
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        
        ax_d = pl.gca()
        
        
        gal = self.info['gal']
        #-- the plot
        dzoom = dprob>1e-10
        pl.plot(d,dprob,'k-')
        pl.grid()
        pl.xlabel('Distance [pc]')
        pl.ylabel('Probability [unnormalized]')
        pl.xlim(d[dzoom].min(),d[dzoom].max())
        xlims = pl.xlim()
        pl.twiny(ax_d)
        pl.xlim(xlims)
        xticks = pl.xticks()
        pl.xticks(xticks[0],['%.2f'%(conversions.convert('pc','Rsol',np.sqrt(self.results['igrid_search']['grid']['scale'][-1])*di)) for di in xticks[0]])
        pl.xlabel('Radius [$R_\odot$]')
        pl.twinx(ax_d)
        res = 100
        d_keep = (xlims[0]<=d[::res]) & (d[::res]<=xlims[1])
        if len(self.results['igrid_search']['drimmel']):
            pl.plot(d[::res][d_keep],self.results['igrid_search']['drimmel'].ravel()[d_keep],'b-',label='Drimmel')
        if len(self.results['igrid_search']['marshall']):
            pl.plot(d[::res][d_keep],self.results['igrid_search']['marshall'].ravel()[d_keep],'b--',label='Marshall')
        ebv = self.results[mtype]['grid']['ebv'][-1]
        pl.plot(xlims,[ebv*3.1,ebv*3.1],'r--',lw=2,label='measured')
        pl.ylabel('Visual extinction $A_v$ [mag]')
        pl.legend(loc='lower right',prop=dict(size='x-small'))
        pl.xlim(xlims)
        logger.info('Plotted distance/reddening')
    
    def plot_grid_model(self,ptype='prob'):
        """
        Grid of models
        """
        pl.title(self.info['spType'])
        cutlogg = (self.results['igrid_search']['grid']['logg']<=4.4) & (self.results['igrid_search']['grid']['CI_red']<=0.95)
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        gal = self.info['gal']
        
        n = 75000
        region = self.results['igrid_search']['grid']['CI_red']<0.95
        total_prob = 100-(1-self.results['igrid_search']['grid']['CI_red'][cutlogg][-n:])*d_prob_models*100
        tp_sa = np.argsort(total_prob)[::-1]
        if ptype=='prob':
            pl.scatter(self.results['igrid_search']['grid']['teff'][cutlogg][-n:][tp_sa],self.results['igrid_search']['grid']['logg'][cutlogg][-n:][tp_sa],
                c=total_prob[tp_sa],edgecolors='none',cmap=pl.cm.spectral,
                vmin=total_prob.min(),vmax=total_prob.max())
        elif ptype=='radii':
            pl.scatter(self.results['igrid_search']['grid']['teff'][cutlogg][-n:][tp_sa],self.results['igrid_search']['grid']['logg'][cutlogg][-n:][tp_sa],
                c=radii,edgecolors='none',cmap=pl.cm.spectral,
                vmin=radii.min(),vmax=radii.max())
        pl.xlim(self.results['igrid_search']['grid']['teff'][region].max(),self.results['igrid_search']['grid']['teff'][region].min())
        pl.ylim(self.results['igrid_search']['grid']['logg'][region].max(),self.results['igrid_search']['grid']['logg'][region].min())
        cbar = pl.colorbar()
        pl.xlabel('log (effective temperature [K]) [dex]')
        pl.ylabel(r'log (surface gravity [cm s$^{-2}$]) [dex]')
        
        if ptype=='prob':
            cbar.set_label('Probability (incl. plx) [%]')
        elif ptype=='radii':
            cbar.set_label('Model radii [$R_\odot$]')
        logger.info('Plotted teff-logg diagram of models (%s)'%(ptype))
        
        
    
    def plot_MW_side(self):
        im = Image.open(config.get_datafile('images','ESOmilkyway.tif'))
        left,bottom,width,height = 0.0,0.0,1.0,1.0
        startm,endm = 183,-177
        startv,endv = -89,91

        xwidth = startm-endm
        ywidth = 90.
        ratio = ywidth/xwidth

        gal = list(self.info['gal'])
        if gal[0]>180:
            gal[0] = gal[0] - 360.
        #-- boundaries of ESO image
        pl.imshow(im,extent=[startm,endm,startv,endv],origin='lower')
        pl.plot(gal[0],gal[1],'rx',ms=15,mew=2)
        pl.xlim(startm,endm)
        pl.ylim(startv,endv)
        pl.xticks([])
        pl.yticks([])
    
    def plot_MW_top(self):
        im = Image.open(config.get_datafile('images','topmilkyway.jpg'))
        pl.imshow(im,origin='lower')
        pl.box(on=False)
        pl.xticks([])
        pl.yticks([])
        xlims = pl.xlim()
        ylims = pl.ylim()
        gal = self.info['gal']
        pl.plot(2800,1720,'ro',ms=10)
        pl.plot([2800,-5000*np.sin(gal[0]/180.*np.pi)+2800],[1720,5000*np.cos(gal[0]/180.*np.pi)+1720],'r-',lw=2)
        
        #-- necessary information
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        
        x = d/10.*1.3
        y = dprob*1000.
        theta = gal[0]/180.*np.pi + np.pi/2.
        x_ = np.cos(theta)*x - np.sin(theta)*y + 2800
        y_ = np.sin(theta)*x + np.cos(theta)*y + 1720
    
        pl.plot(x_,y_,'r-',lw=2)
        index = np.argmax(y)
        pl.plot(np.cos(theta)*x[index] + 2800,np.sin(theta)*x[index] + 1720,'rx',ms=15,mew=2)
    
        pl.xlim(xlims)
        pl.ylim(ylims)
    
    def plot_finderchart(self):
        
        try:
            data,coords,size = mast.get_dss_image(self.ID)
            pl.imshow(data[::-1],extent=[-size[0]/2*60,size[0]/2*60,
                                        -size[1]/2*60,size[1]/2*60],cmap=pl.cm.RdGy_r)#Greys
            xlims,ylims = pl.xlim(),pl.ylim()
            keep_this = -self.master['color'] & (self.master['cmeas']>0)
            toplot = self.master[keep_this]
            systems = np.array([system.split('.')[0] for system in toplot['photband']],str)
            set_systems = sorted(list(set(systems)))
            pl.gca().set_color_cycle([pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))])    
            for system in set_systems:
                keep = systems==system
                if sum(keep):
                    pl.plot(toplot['_RAJ2000'][keep][0]*60,
                            toplot['_DEJ2000'][keep][0]*60,'x',label=system,
                            mew=2.5,ms=15,alpha=0.5)
                else:
                    pl.gca()._get_lines.count += 1
            pl.legend(numpoints=1,prop=dict(size='x-small'))
            pl.xlim(xlims)
            pl.ylim(ylims)
            pl.xlabel(r'Right ascension $\alpha$ [arcmin]')
            pl.ylabel(r'Declination $\delta$ [arcmin]')    
        except:
            pass
    
    
        
        
    def make_plots(self):
        """
        Make all available plots
        """
        pl.figure(figsize=(22,12))
        rows,cols = 3,4
        pl.subplots_adjust(left=0.04, bottom=0.07, right=0.97, top=0.96,
                wspace=0.17, hspace=0.24)
        pl.subplot(rows,cols,1);self.plot_grid(ptype='CI_red')
        pl.subplot(rows,cols,2);self.plot_grid(ptype='ebv')
        pl.subplot(rows,cols,3);self.plot_grid(ptype='z')
        pl.subplot(rows,cols,4);self.plot_distance()
        
        pl.subplot(3,2,3);self.plot_sed(colors=False)
        pl.subplot(3,2,5);self.plot_sed(colors=True)
        
        pl.subplot(rows,cols,7);self.plot_chi2(colors=False)
        pl.subplot(rows,cols,11);self.plot_chi2(colors=True)
        
        pl.subplot(rows,cols,8);self.plot_grid_model(ptype='prob')
        pl.subplot(rows,cols,12);self.plot_grid_model(ptype='radii')
        
        pl.figure(figsize=(12,12))
        pl.axes([0,0.0,1.0,0.5]);self.plot_MW_side()
        pl.axes([0,0.5,0.5,0.5]);self.plot_MW_top()
        pl.axes([0.5,0.5,0.5,0.5]);self.plot_finderchart()
    
    def save2fits(self,filename=None,overwrite=True):
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.fits'
        if overwrite:
            if os.path.isfile(filename): os.remove(filename)
            
        eff_waves,synflux,photbands = self.results['synflux']
        chi2 = self.results['chi2']
        
        master = mlab.rec_append_fields(self.master, 'synflux',synflux)
        master = mlab.rec_append_fields(master, 'mod_eff_wave',eff_waves)
        master = mlab.rec_append_fields(master, 'chi2',chi2)
        
        results_dict = {}
        keys = sorted(self.results['igrid_search'])
        for key in keys:
            if 'CI' in key:
                for ikey in self.results['igrid_search'][key]:
                    results_dict[ikey] = self.results['igrid_search'][key][ikey]
                
            
        fits.write_recarray(master,filename)
        fits.write_array(list(self.results['model']),filename,
                                names=('wave','flux','dered_flux'),
                                units=('A','erg/s/cm2/A','erg/s/cm2/A'),
                                header_dict=results_dict)
        
        fits.write_recarray(self.results['igrid_search']['grid'],filename,header_dict=results_dict)
    
    #}


if __name__ == "__main__":
    #-- PCA analysis
    master['include'] = True
    exclude(master,names=['STROMGREN.HBN-HBW','USNOB1'],wrange=(1.5e4,np.inf))
    do_pca = False
    include_pca = master['include']

    if sum(keep)>2:
        do_pca = True
        logger.info("Start of PCA analysis to find fundamental parameters")
        colors,index = np.unique1d(master['photband'][include_pca],return_index=True)
        A,grid = fit.get_PCA_grid(colors,ebvrange=(0,0.5),res=10)
        P,T,(means,stds) = fit.get_PCA(A)
        calib = fit.calibrate_PCA(T,grid,function='linear')
        obsT,e_obsT = master['cmeas'][include_pca][index], master['e_cmeas'][include_pca][index]
        pars = fit.get_PCA_parameters(obsT,calib,P,means,stds,e_obsT=e_obsT,mc=mc)
        teff,logg,ebv = pars[0]
        logger.info("PCA result: teff=%.0f, logg=%.2f, E(B-V)=%.3f"%(teff,logg,ebv))

        #-- find angular diameter
        logger.info('Estimation of angular diameter')
        iflux = model.get_itable(teff=teff,logg=logg,ebv=ebv,photbands=master['photband'][include_grid])
        scale_pca = np.median(master['cmeas'][include_grid]/iflux)
        angdiam = 2*np.arctan(np.sqrt(scale_pca))/np.pi*180*3600*1000
        logger.info('Angular diameter = %.3f mas'%(angdiam))

        #-- get model
        wave_pca,flux_pca = model.get_table(teff=teff,logg=logg,ebv=ebv,law='fitzpatrick1999')
        wave_ur_pca,flux_ur_pca = model.get_table(teff=teff,logg=logg,ebv=0,law='fitzpatrick1999')

    logger.info('...brought to you in %.3fs'%(time.time()-c0))

    pl.figure()
    pl.title(ID)
    toplot = master[-master['color']]
    systems = np.array([system.split('.')[0] for system in toplot['photband']],str)
    set_systems = sorted(list(set(systems)))
    pl.gca().set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0, 1.0, len(set_systems))])
    for system in set_systems:
        keep = systems==system
        pl.errorbar(master['cwave'][keep],master['cmeas'][keep],
                yerr=master['e_cmeas'][keep],fmt='o',label=system)
    pl.gca().set_xscale('log',nonposx='clip')
    pl.gca().set_yscale('log',nonposy='clip')
    pl.gca().set_autoscale_on(False)
    pl.plot(wave_pca,flux_pca*scale_pca,'r-')
    pl.plot(wave_ur_pca,flux_ur_pca*scale_pca,'k-')
    pl.grid()
    pl.legend(loc='lower left',prop=dict(size='x-small'))
    pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\n'%(teff,logg,ebv)+r'$\theta$=%.3f mas'%(angdiam),(0.75,0.80),xycoords='axes fraction')
    pl.xlabel('wavelength [$\mu$m]')
    pl.ylabel(r'$F_\nu$ [Jy]')
    #### -- END SOME FIGURES -- ####


    if mc and do_pca:
        pl.figure()
        pl.subplot(131)
        pl.hexbin(pars[:,0],pars[:,1])
        pl.xlim(pars[:,0].max(),pars[:,0].min())
        pl.ylim(pars[:,1].max(),pars[:,1].min())
        pl.colorbar()
        pl.subplot(132)
        pl.hexbin(pars[:,0],pars[:,2])
        pl.xlim(pars[:,0].max(),pars[:,0].min())
        pl.ylim(pars[:,2].min(),pars[:,2].max())
        pl.colorbar()
        pl.subplot(133)
        pl.hexbin(pars[:,1],pars[:,2])
        pl.xlim(pars[:,1].max(),pars[:,1].min())
        pl.ylim(pars[:,2].min(),pars[:,2].max())
        pl.colorbar()
    pl.show()









    ################################
