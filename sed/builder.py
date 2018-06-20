# -*- coding: utf-8 -*-
"""
SED builder program.
"""

import re
import sys
import time
import os
import logging
import itertools
import glob
import json

import pylab as pl
from matplotlib import mlab
try:
    from PIL import Image
except ImportError:
    print("The PIL package is discontinued")
import numpy as np
import scipy.stats
import astropy.io.fits as pf

from ivs import config
from ivs.aux import numpy_ext
from ivs.aux import termtools
from ivs.inout import ascii
from ivs.inout import fits
from ivs.inout import hdf5
from ivs.sed import model
from ivs.sed import filters
from ivs.sed import fit
from ivs.sed import distance
from ivs.sed import extinctionmodels
from ivs.sed.decorators import standalone_figure
from ivs.spectra import tools
from ivs.catalogs import crossmatch
from ivs.catalogs import vizier
from ivs.catalogs import mast
from ivs.catalogs import sesame
from ivs.catalogs import corot
from ivs.units import conversions
from ivs.units import constants
from uncertainties import unumpy, ufloat
from uncertainties.unumpy import sqrt as usqrt
from uncertainties.unumpy import tan as utan
from ivs.sigproc import evaluate
try:
    # This module has now been removed,
    # perhaps future re-implementation if demanded
    from ivs.stellar_evolution import evolutionmodels
except ImportError:
    print("Warning: The ivs.stellar_evolution module has been removed from the"
          "repository as of 03.11.2017.")
    print("Stellar evolution models are no longer available to use.")

logger = logging.getLogger("SED.BUILD")
# logger.setLevel(10)


def fix_master(master, e_default=None):
    """
    Clean/extend/fix record array received from C{get_photometry}.
    This function does a couple of things:

    1.  Adds common but uncatalogized colors like 2MASS.J-H if not
        already present. WARNING: these colors can mix values from the same
        system but from different catalogs!
    2.  Removes photometry for which no calibration is available
    3.  Adds a column 'color' with flag False denoting absolute flux
        measurement and True denoting color.
    4.  Adds a column 'include' with flag True meaning the value will be
        included in the fit
    5.  Sets some default errors to photometric values, for which we know
        that the catalog values are not trustworthy.
    6.  Sets a lower limit to the allowed errors of 1%. Errors below this
        value are untrustworthy because the calibration error is larger than
        that.
    7.  USNOB1 and ANS photometry are set the have a minimum error of 30%
        due to uncertainties in response curves.

    :param master: record array containing photometry. This should have the
    fields ``meas``, ``e_meas``, ``unit``, ``photband``, ``source``, ``_r``,
    ``_RAJ2000``, ``DEJ2000``, ``cmeas``, ``e_cmeas``, ``cwave``, ``cunit``
    :type master: record array
    :param e_default: default error for measurements without errors
    :type e_default: float
    :return: record array extended with fields 'include' and 'color', and with
    rows added (see above description)
    :rtype: numpy recard array
    """
    # -- we recognize uncalibrated stuff as those for which no absolute flux
    # was obtained, and remove them:
    master = master[~np.isnan(master['cmeas'])]
    cats = np.array([ph.split('.')[0] for ph in master['source']])
    set_cats = sorted(set(cats))
    # -- add common uncatalogized colors:
    columns = list(master.dtype.names)
    # -- if the separate bands are available (and not the color itself),
    #   calculate the colors here. We need to do this for every separate
    #   system!
    add_rows = []
    for cat in set_cats:
        master_ = master[cats == cat]

        for color in ['2MASS.J-H', '2MASS.KS-H', 'TD1.1565-1965',
                      'TD1.2365-1965', 'TD1.2365-2740', 'JOHNSON.J-H',
                      'JOHNSON.K-H', 'JOHNSON.B-V', 'JOHNSON.U-B',
                      'JOHNSON.R-V', 'JOHNSON.I-V', 'GALEX.FUV-NUV',
                      'TYCHO2.BT-VT', 'WISE.W1-W2', 'WISE.W3-W2', 'WISE.W4-W3',
                      'SDSS.U-G', 'SDSS.G-R', 'SDSS.R-I', 'SDSS.R-Z',
                      'UVEX.U-G', 'UVEX.G-R', 'UVEX.R-I', 'UVEX.HA-I',
                      'DENIS.I-J', 'DENIS.KS-J', 'ANS.15N-15W', 'ANS.15W-18',
                      'ANS.18-22', 'ANS.22-25', 'ANS.25-33']:

            # -- get the filter system and the separate bands for these colors
            system, band = color.split('.')
            band0, band1 = band.split('-')
            band0, band1 = '%s.%s' % (system, band0), '%s.%s' % (system, band1)

            if band0 in master_['photband'] and band1 in master_['photband'] and not color in master_['photband']:
                # -- where are the bands located?
                index0 = list(master_['photband']).index(band0)
                index1 = list(master_['photband']).index(band1)

                # -- start a new row to add the color
                row = list(master_[index1])
                row[columns.index('photband')] = color
                row[columns.index('cwave')] = np.nan

                # -- it could be a magnitude difference
                if master_['unit'][index0]=='mag':
                    row[columns.index('meas')] = master_['meas'][index0]-master_['meas'][index1]
                    row[columns.index('e_meas')] = np.sqrt(master_['e_meas'][index0]**2+master_['e_meas'][index1]**2)
                    # -- error will not always be available...
                    try:
                        row[columns.index('cmeas')],row[columns.index('e_cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],row[columns.index('e_meas')],photband=color)
                    except AssertionError:
                        row[columns.index('cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],photband=color)
                        row[columns.index('e_cmeas')] = np.nan
                # -- or it could be a flux ratio
                else:
                    row[columns.index('meas')] = master_['meas'][index0]/master_['meas'][index1]
                    row[columns.index('e_meas')] = np.sqrt(((master_['e_meas'][index0]/master_['meas'][index0])**2+(master_['e_meas'][index1]/master_['meas'][index1])**2)*row[columns.index('meas')]**2)
                    row[columns.index('cmeas')],row[columns.index('e_cmeas')] = row[columns.index('meas')],row[columns.index('e_meas')]
                row[columns.index('cunit')] = 'flux_ratio'
                add_rows.append(tuple(row))
    master = numpy_ext.recarr_addrows(master,add_rows)

    # -- add an extra column with a flag to distinguish colors from absolute
    #   fluxes, and a column with flags to include/exclude photometry
    #   By default, exclude photometry if the effective wavelength is above
    #   150 mum or if the measurement is nan
    if not 'color' in master.dtype.names:
        extra_cols = [[filters.is_color(photband) for photband in master['photband']],
                    [(cwave<150e4) or np.isnan(meas) for cwave,meas in zip(master['cwave'],master['cmeas'])]]
        dtypes = [('color',np.bool),('include',np.bool)]
        master = numpy_ext.recarr_addcols(master,extra_cols,dtypes)
    else:
        iscolor = [filters.is_color(photband) for photband in master['photband']]
        master['color'] = iscolor
    # -- add an extra column with indices to distinguish different data sets, e.g. based
    #   on photometric phase.
    if not 'phase' in master.dtype.names:
        extra_cols = [[0]*len(master['meas'])]
        dtypes = [('phase',np.int)]
        master = numpy_ext.recarr_addcols(master,extra_cols,dtypes)
    # -- set default errors if no errors are available and set really really
    #   small errors to some larger default value
    if e_default is not None:
        no_error = np.isnan(master['e_cmeas'])
        master['e_cmeas'][no_error] = np.abs(e_default*master['cmeas'][no_error])
        small_error = master['e_cmeas']<(0.01*np.abs(master['cmeas']))
        master['e_cmeas'][small_error] = 0.01*np.abs(master['cmeas'][small_error])

    # -- the measurements from USNOB1 are not that good because the response
    #   curve is approximated by the JOHNSON filter. Set the errors to min 30%
    #   The same holds for ANS: here, the transmission curves are very uncertain
    for i,photband in enumerate(master['photband']):
        if 'USNOB1' in photband or 'ANS' in photband:
            master['e_cmeas'][i] = max(0.30*master['cmeas'][i],master['e_cmeas'][i])

    # -- remove negative fluxes
    master = master[master['cmeas']>0]

    return master


def decide_phot(master,names=None,wrange=None,sources=None,indices=None,ptype='all',include=False):
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

    1 Exclude all measurements::

      >> decide_phot(master,wrange=(-np.inf,+np.inf),ptype='all',include=False)

    2 Include all TD1 fluxes and colors::

      >> decide_phot(master,names=['TD1'],ptype='all',include=True)

    3 Include all V band measurements from all systems (but not the colors)::

      >> decide_phot(master,names=['.V'],ptype='abs',include=True)

    4 Include all Geneva colors and exclude Geneva magnitudes::

      >> decide_phot(master,names=['GENEVA'],ptype='col',include=True)
      >> decide_phot(master,names=['GENEVA'],ptype='abs',include=False)

    5 Exclude all infrared measurements beyond 1 micron::

      >> decide_phot(master,wrange=(1e4,np.inf),ptype='all',include=False)

    6 Include all AKARI measurements below 10 micron::

      >> decide_phot(master,names=['AKARI'],wrange=(-np.inf,1e5),ptype='all',include=True)

    :param master: record array containing all photometry
    :type master: numpy record array
    :param names: strings excerpts to match filters
    :type names: list of strings
    :param wrange: wavelength range (most likely angstrom) to include/exclude
    :type wrange: 2-tuple (start wavelength,end wavelength)
    :param sources: list of sources
    :type sources: list of strings
    :param indices: list of indices (integers)
    :type indices: list of integers
    :param ptype: type of photometry to include/exclude: absolute values, colors
    or both
    :type ptype: string, one of 'abs','col','all'
    :param include: flag setting exclusion or inclusion
    :type include: boolean
    :return: master record array with adapted 'include' flags
    :rtype: numpy record array
    """
    # -- exclude/include passbands based on their names
    if names is not None:
        logger.info('%s photometry based on photband containining one of %s'%((include and 'Include' or "Exclude"),names))
        for index,photband in enumerate(master['photband']):
            for name in names:
                if name in photband:
                    if ptype=='all' or (ptype=='abs' and -master['color'][index]) or (ptype=='col' and master['color'][index]):
                        master['include'][index] = include
                        break
    # -- exclude/include colors based on their wavelength
    if wrange is not None:
        logger.info('%s photometry based on photband wavelength between %s'%((include and 'Include' or "Exclude"),wrange))
        for index,photband in enumerate(master['photband']):
            system,color = photband.split('.')
            if not '-' in color or ptype=='abs':
                continue
            band0,band1 = color.split('-')
            cwave0 = filters.eff_wave('%s.%s'%(system,band0))
            cwave1 = filters.eff_wave('%s.%s'%(system,band1))
            if (wrange[0]<cwave0<wrange[1]) or (wrange[0]<cwave0<wrange[1]):
                master['include'][index] = include
        # -- exclude/include passbands based on their wavelength
        if not ptype=='col':
            master['include'][(wrange[0]<master['cwave']) & (master['cwave']<wrange[1])] = include
    # -- exclude/include passbands based on their source
    if sources is not None:
        logger.info('%s photometry based on source catalog in %s'%((include and 'Include' or "Exclude"),sources))
        for index,msource in enumerate(master['source']):
            for source in sources:
                if source==msource:
                    if ptype=='all' or (ptype=='abs' and -master['color'][index]) or (ptype=='col' and master['color'][index]):
                        master['include'][index] = include
                        break
    # -- exclude/include passbands based on their index
    if indices is not None:
        logger.info('%s photometry based on index'%((include and 'Include' or "Exclude")))
        indices = np.array(indices,int)
        if not indices.shape: indices = indices.reshape(1)
        master['include'][indices] = include


def photometry2str(master, comment='', sort='photband', color=False,
                   index=False):
    """
    String representation of master record array.

    Sorting is disabled when C{index=True}.

    :param master: master record array containing photometry
    :type master: numpy record array
    """
    if sort and not index:
        master = master[np.argsort(master[sort])]

    templateh = '{:20s} {:>12s} {:>12s} {:12s} {:>10s} {:>12s} {:>12s} {:12s} {:30s}'
    templated = '{:20s} {:12g} {:12g} {:12s} {:10.0f} {:12g} {:12g} {:12s} {:30s}'
    header = ['PHOTBAND', 'MEAS', 'E_MEAS', 'UNIT', 'CWAVE', 'CMEAS',
              'E_CMEAS', 'CUNIT', 'SOURCE']
    if 'comments' in master.dtype.names:
        templateh += ' {:s}'
        templated += ' {:s}'
        header += ['COMMENTS']
    if index:
        templateh = '{:3s} '+templateh
        templated = '{:3d} '+templated
        header = ['NR']+header

    txt = [comment+templateh.format(*header)]
    txt.append(comment+'='*170)
    columns = [master[col.lower()] for col in header if not col == 'NR']
    for nr, contents in enumerate(zip(*columns)):
        contents = list(contents)
        if 'comments' in master.dtype.names:
            contents[-1] = contents[-1].replace('_', ' ')
        if index:
            contents = [nr] + contents
        line = templated.format(*contents)
        if color:
            mycolor = termtools.green if master['include'][nr] else termtools.red
            line = mycolor(line)
        txt.append(comment + line)
    return "\n".join(txt)

#@memoized
#def get_schaller_grid():
    #"""
    #Download Schaller 1992 evolutionary tracks and return an Rbf interpolation
    #function.

    #:return: Rbf interpolation function
    #:rtype: Rbf interpolation function
    #"""
    ## -- translation between table names and masses
    ##masses = [1,1.25,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,40,60][:-1]
    ##tables = ['table20','table18','table17','table16','table15','table14',
    ##          'table13','table12','table11','table10','table9','table8',
    ##          'table7','table6','table5','table4','table3'][:-1]
    ## -- read in all the tables and compute radii and luminosities.
    #data,comms,units = vizier.search('J/A+AS/96/269/models')
    #all_teffs = 10**data['logTe']
    #all_radii = np.sqrt((10**data['logL']*constants.Lsol_cgs)/(10**data['logTe'])**4/(4*np.pi*constants.sigma_cgs))
    #all_loggs = np.log10(constants.GG_cgs*data['Mass']*constants.Msol_cgs/(all_radii**2))
    #all_radii /= constants.Rsol_cgs
    ## -- remove low temperature models, the evolutionary tracks are hard to
    ##   interpolate there.
    #keep = all_teffs>5000
    #all_teffs = all_teffs[keep]
    #all_radii = all_radii[keep]
    #all_loggs = all_loggs[keep]
    ## -- make linear interpolation model between all modelpoints
    #mygrid = Rbf(np.log10(all_teffs),all_loggs,all_radii,function='linear')
    #logger.info('Interpolation of Schaller 1992 evolutionary tracks to compute radii')
    #return mygrid







#def get_radii(teffs,loggs):
    #"""
    #Retrieve radii from stellar evolutionary tracks from Schaller 1992.

    #:param teffs: model effective temperatures
    #:type teffs: numpy array
    #:param loggs: model surface gravities
    #:type loggs: numpy array
    #:return: model radii (solar units)
    #:rtype: numpy array
    #"""
    #mygrid = get_schaller_grid()
    #radii = mygrid(np.log10(teffs),loggs)
    #return radii


#def calculate_distance(plx,gal,teffs,loggs,scales,n=75000):
    #"""
    #Calculate distances and radii of a target given its parallax and location
    #in the galaxy.


    #"""
    ## -- compute distance up to 25 kpc, which is about the maximum distance from
    ##   earth to the farthest side of the Milky Way galaxy
    ##   rescale to set the maximum to 1
    #d = np.logspace(np.log10(0.1),np.log10(25000),100000)
    #if plx is not None:
        #dprob = distance.distprob(d,gal[1],plx)
        #dprob = dprob / dprob.max()
    #else:
        #dprob = np.ones_like(d)
    ## -- compute the radii for the computed models, and convert to parsec
    ##radii = np.ones(len(teffs[-n:]))
    #radii = get_radii(teffs[-n:],loggs[-n:])
    #radii = conversions.convert('Rsol','pc',radii)
    #d_models = radii/np.sqrt(scales[-n:])
    ## -- we set out of boundary values to zero
    #if plx is not None:
        #dprob_models = np.interp(d_models,d[-n:],dprob[-n:],left=0,right=0)
    #else:
        #dprob_models = np.ones_like(d_models)
    ## -- reset the radii to solar units for return value
    #radii = conversions.convert('pc','Rsol',radii)
    #return (d_models,dprob_models,radii),(d,dprob)



class SED(object):
    """
    Class that facilitates the use of the ivs.sed module.

    This class is meant to be an easy interface to many of the ivs.sed module's
    functionality.

    The most important attributes of SED are:

    1 C{sed.ID}: star's identification (str)
    2 C{sed.photfile}: name of the file containing all photometry (str)
    3 C{sed.info}: star's information from Simbad (dict)
    4 C{sed.master}: photometry data (record array)
    5 C{sed.results}: results and summary of the fitting process (dict)

    After fitting, e.g. via calling :class:`igrid_search`, you can call
    :class:`get_model` to retrieve the full SED matching the best fitting
    parameters (or, rather, closely matching them, see the documentation).

    """

    def __init__(self, ID=None, photfile=None, plx=None, load_fits=True,
                 load_hdf5=True, label=''):
        """
        Initialize SED class.

        Different ways to initialize:

        B{1. If no previous data are saved:}

        >>> mysed = SED('HD129929')

        B{2. If previous data exists:}

        >>> #mysed = SED(photfile='HD129929.phot') # will set ID with 'oname' field from SIMBAD
        >>> #mysed = SED(ID='bla', photfile='HD129929.phot') # Sets custom ID

        The C{ID} variable is used internally to look up data, so it should be
        something SIMBAD understands and that designates the target.

        :param plx: parallax (and error) of the object
        :type plx: tuple (plx,e_plx)
        """
        self.ID = ID
        self.label = label
        self.info = {}
        # -- the file containing photometry should have the following name. We
        #   save photometry to a file to avoid having to download photometry
        #   each time from the internet
        if photfile is None:
            self.photfile = '%s.phot' % (ID).replace(' ', '_')
            # -- keep information on the star from SESAME, but override parallax with
            #   the value from Van Leeuwen's new reduction. Set the galactic
            #   coordinates, which are not available from SESAME.
        else:
            self.photfile = photfile

        # -- load information from the photometry file if it exists
        if not os.path.isfile(self.photfile):
            try:
                self.info = sesame.search(os.path.basename(ID),fix=True)
            except KeyError:
                logger.warning('Star %s not recognised by SIMBAD'%(os.path.basename(ID)))
                try:
                    self.info = sesame.search(os.path.basename(ID),db='N',fix=True)
                    logger.info('Star %s recognised by NED'%(os.path.basename(ID)))
                except KeyError:
                    logger.warning('Star %s not recognised by NED'%(os.path.basename(ID)))
            # -- final attempt: if it's a CoRoT star, try the catalog
            if 'corot' in ID.lower() and not 'galpos' in self.info:
                corot_id = int("".join([char for char in ID if char.isdigit()]))
                try:
                    self.info['jradeg'],self.info['jdedeg'] = ra,dec = corot.resolve(corot_id)
                    gal = conversions.convert('equatorial','galactic',(ra,dec),epoch='2000')
                    self.info['galpos'] = float(gal[0])/np.pi*180,float(gal[1])/np.pi*180
                    logger.info('Resolved position via CoRoT EXO catalog')
                except:
                    logger.warning("Star %s not recognised by CoRoT"%(os.path.basename(ID)))

            if plx is not None:
                if not 'plx' in self.info:
                    self.info['plx'] = {}
                self.info['plx']['v'] = plx[0]
                self.info['plx']['e'] = plx[1]
        else:
            self.load_photometry()
            logger.info('Photometry loaded from file')
            # -- if no ID was given, set the official name as the ID.
            if self.ID is None:
                self.ID = os.path.splitext(os.path.basename(self.photfile))[0]#self.info['oname']
                logger.info('Name from file used to set ID of object')
        # --load information from the FITS file if it exists
        self.results = {}
        self.constraints = {}
        if load_fits:
            self.load_fits()
        if load_hdf5:
            self.load_hdf5()

        # -- prepare for information on fitting processes
        self.CI_limit = 0.95

    def __repr__(self):
        """
        Machine readable string representation of an SED object.
        """
        return "SED('{}')".format(self.ID)

    def __str__(self):
        """
        Human readable string representation of an SED object.
        """
        txt = []
        # -- object designation
        txt.append("Object identification: {:s}".format(self.ID))
        if hasattr(self,'info') and self.info and 'oname' in self.info:
            txt.append("Official designation: {:s}".format(self.info['oname']))
        # -- additional info
        for key in sorted(self.info.keys()):
            if isinstance(self.info[key],dict):
                txt.append(" {:10s} = ".format(key)+", ".join(["{}: {}".format(i,j) for i,j in self.info[key].items()]))
            else:
                txt.append(" {:10s} = {}".format(key,self.info[key]))

        #txt.append("Included photometry:")
        #if hasattr(self,'master') and self.master is not None:
            #include_grid = self.master['include']
            #txt.append(photometry2str(self.master[include_grid]))
        #txt.append("Excluded photometry:")
        #if hasattr(self,'master') and self.master is not None:
            #include_grid = self.master['include']
            #txt.append(photometry2str(self.master[-include_grid]))
        if hasattr(self,'master'):
            txt.append(photometry2str(self.master,color=True))
        return "\n".join(txt)

    #{ Handling photometric data
    def get_photometry(self,radius=None,ra=None,dec=None,
                       include=None,exclude=None,force=False,
                       units='erg/s/cm2/AA'):
        """
        Search photometry on the net or from the phot file if it exists.

        For bright stars, you can set radius a bit higher...

        :param radius: search radius (arcseconds)
        :type radius: float.
        """
        if radius is None:
            if 'mag.V.v' in self.info and self.info['mag.V.v']<6.:
                radius = 60.
            else:
                radius = 10.
        if not os.path.isfile(self.photfile) or force:
            # -- get and fix photometry. Set default errors to 1%, and set
            #   USNOB1 errors to 3%
            if ra is None and dec is None:
                master = crossmatch.get_photometry(ID=os.path.basename(self.ID),radius=radius,
                                       include=include,exclude=exclude,to_units=units,
                                       extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            else:
                master = crossmatch.get_photometry(ID=os.path.basename(self.ID),ra=ra,dec=dec,radius=radius,
                                       include=include,exclude=exclude,to_units=units,
                                       extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            if 'jradeg' in self.info:
                master['_RAJ2000'] -= self.info['jradeg']
                master['_DEJ2000'] -= self.info['jdedeg']

            # -- fix the photometry: set default errors to 2% and print it to the
            #   screen
            self.master = fix_master(master,e_default=0.1)
            logger.info('\n'+photometry2str(master))

            # -- write to file
            self.save_photometry()

    def get_spectrophotometry(self,directory=None,force_download=False):
        """
        Retrieve and combine spectra.

        B{WARNING:} this function creates FUSE and DIR directories!
        """
        if directory is None and os.path.dirname(self.ID) == '':
            directory = os.getcwd()
        elif os.path.dirname(self.ID) != '':
            directory = os.path.dirname(self.ID)

        if not os.path.isdir(directory):
            os.mkdir(directory)

        # -- add spectrophotometric filters to the set
        photbands = filters.add_spectrophotometric_filters(R=200,lambda0=950,lambdan=3350)
        if hasattr(self,'master') and self.master is not None and not force_download:
            if any(['BOXCAR' in photband for photband in self.master['photband']]):
                return None
        # -- FUSE spectra
        fuse_direc = os.path.join(directory,'FUSE')
        iue_direc = os.path.join(directory,'IUE')
        if not os.path.isdir(fuse_direc) or force_download:
            out1 = mast.get_FUSE_spectra(ID=os.path.basename(self.ID),directory=fuse_direc,select=['ano'])
            if out1 is None:
                out1 = []
        else:
            out1 = glob.glob(fuse_direc+'/*')
        # -- IUE spectra
        if not os.path.isdir(iue_direc) or force_download:
            out2 = vizier.get_IUE_spectra(ID=os.path.basename(self.ID),directory=iue_direc,select='lo')
            if out2 is None:
                out2 = []
        else:
            out2 = glob.glob(iue_direc+'/*')
        # -- read them in to combine
        list_of_spectra = [fits.read_fuse(ff)[:3] for ff in out1]
        list_of_spectra+= [fits.read_iue(ff)[:3] for ff in out2]
        # -- combine
        wave,flux,err,nspec = tools.combine(list_of_spectra)
        # -- add to the master
        N = len(flux)
        units = ['erg/s/cm2/AA']*N
        source = ['specphot']*N
        self.add_photometry_fromarrays(flux,err,units,photbands,source,flags=nspec)

        # -- write to file
        self.master = fix_master(self.master,e_default=0.1)
        logger.info('\n'+photometry2str(self.master))
        self.save_photometry()




    def exclude(self,names=None,wrange=None,sources=None,indices=None):
        """
        Exclude (any) photometry from fitting process.

        If called without arguments, all photometry will be excluded.
        """
        if names is None and wrange is None and sources is None and indices is None:
            wrange = (-np.inf,np.inf)
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=False,ptype='all')

    def exclude_colors(self,names=None,wrange=None,sources=None,indices=None):
        """
        Exclude (color) photometry from fitting process.
        """
        if names is None and wrange is None and sources is None and indices is None:
            wrange = (-np.inf,0)
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=False,ptype='col')

    def exclude_abs(self,names=None,wrange=None,sources=None,indices=None):
        """
        Exclude (absolute) photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=False,ptype='abs')


    def include(self,names=None,wrange=None,sources=None,indices=None):
        """
        Include (any) photometry in fitting process.
        """
        if names is None and wrange is None and sources is None and indices is None:
            wrange = (-np.inf,np.inf)
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=True,ptype='all')

    def include_colors(self,names=None,wrange=None,sources=None,indices=None):
        """
        Include (color) photometry in fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=True,ptype='col')

    def include_abs(self,names=None,wrange=None,sources=None,indices=None):
        """
        Include (absolute) photometry in fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=True,ptype='abs')

    def set_photometry_scheme(self,scheme,infrared=(1,'micron')):
        """
        Set a default scheme of colors/absolute values to fit the SED.
        Possible values:

        1 ``abs``: means excluding all colors, including all absolute values
        2 ``color``: means including all colors, excluding all absolute values
        3 ``combo``: means inculding all colors, and one absolute value per
        system (the one with the smallest relative error)
        4 ``irfm``: means mimic infrared flux method: choose absolute values
        in the infrared (define with C{infrared}), and colours in the optical

        :param infrared: definition of start of infrared for infrared flux method
        :type infrared: tuple (value <float>, unit <str>)
        """
        # -- only absolute values: real SED fitting
        if 'abs' in scheme.lower():
            self.master['include'][self.master['color']] = False
            self.master['include'][~self.master['color']] = True
            logger.info('Fitting procedure will use only absolute fluxes (%d)'%(sum(self.master['include'])))
        # -- only colors: color fitting
        elif 'col' in scheme.lower():
            self.master['include'][self.master['color']] = True
            self.master['include'][~self.master['color']] = False
            logger.info('Fitting procedure will use only colors (%d)'%(sum(self.master['include'])))
        # -- combination: all colors and one absolute value per system
        elif 'com' in scheme.lower():
            self.master['include'][self.master['color'] & self.master['include']] = True
            systems = np.array([photband.split('.')[0] for photband in self.master['photband']])
            set_systems = sorted(list(set(systems)))
            for isys in set_systems:
                keep = -self.master['color'] & (isys==systems) & self.master['include']
                if not sum(keep): continue
                index = np.argmax(self.master['e_cmeas'][keep]/self.master['cmeas'][keep])
                index = np.arange(len(self.master))[keep][index]
                self.master['include'][keep] = False
                self.master['include'][index] = True
            logger.info('Fitting procedure will use colors + one absolute flux for each system (%d)'%(sum(self.master['include'])))
        # -- infrared flux method scheme
        elif 'irfm' in scheme.lower():
            # -- check which measurements are in the infrared
            for i,meas in enumerate(self.master):
                # -- if measurement is in the infrared and it is a color, remove it (only use absolute values in IR)
                # -- for colors, make sure *both* wavelengths are checked
                if meas['color']:
                    bands = filters.get_color_photband(meas['photband'])
                    eff_waves = filters.eff_wave(bands)
                    is_infrared = any([(conversions.Unit(*infrared)<conversions.Unit(eff_wave,'AA')) for eff_wave in eff_waves])
                else:
                    is_infrared = conversions.Unit(*infrared)<conversions.Unit(meas['cwave'],'AA')
                if is_infrared and meas['color']:
                    self.master['include'][i] = False
                # -- if measurement is in the infrared and it is not color, keep it (only use absolute values in IR)
                elif is_infrared:
                    self.master['include'][i] = True
                # -- if measurement is not in the infrared and it is a color, keep it (only use colors in visible)
                elif not is_infrared and meas['color']:
                    self.master['include'][i] = True
                # -- if measurement is not in the infrared and it is not a color, remove it (only use colors in visible)
                elif not is_infrared:
                    self.master['include'][i] = False
            use_colors = sum(self.master['include'] & self.master['color'])
            use_abs = sum(self.master['include'] & ~self.master['color'])
            logger.info('Fitting procedure will use colors (%d) < %.3g%s < absolute values (%d)'%(use_colors,infrared[0],infrared[1],use_abs))


    def add_photometry_fromarrays(self,meas,e_meas,units,photbands,source,flags=None,phases=None):
        """
        Add photometry from numpy arrays

        By default unflagged, and at the ra and dec of the star. No color and
        included.

        :param meas: original measurements (fluxes, magnitudes...)
        :type meas: array
        :param e_meas: error on original measurements in same units as C{meas}
        :type e_meas: array
        :param units: units of original measurements
        :type units: array of strings
        :param photbands: photometric passbands of original measurements
        :type photbands: array of strings
        :param source: source of original measurements
        :type source: array of strings
        """
        if flags is None:
            flags = np.nan*np.ones(len(meas))
        if phases is None:
            phases = np.zeros(len(meas))
        # -- if master record array does not exist, make a new one
        if not hasattr(self, 'master') or self.master is None:
            dtypes = [('meas', 'f8'), ('e_meas', 'f8'), ('flag', 'U20'),
                      ('unit', 'U30'), ('photband', 'U30'), ('source', 'U50'),
                      ('_r', 'f8'), ('_RAJ2000', 'f8'), ('_DEJ2000', 'f8'),
                      ('cwave', 'f8'), ('cmeas', 'f8'), ('e_cmeas', 'f8'),
                      ('cunit', 'U50'), ('color', bool), ('include', bool),
                      ('phase', int), ('bibcode', 'U20'), ('comments', 'U200')]
            logger.info('No previous measurements available, initialising master record')
            self.master = np.rec.fromarrays(np.array([ [] for i in dtypes]), dtype=dtypes)
            _to_unit = 'erg/s/cm2/AA'
        else:
            _to_unit = self.master['cunit'][0]

        extra_master = np.zeros(len(meas),dtype=self.master.dtype)

        for i,(m,e_m,u,p,s,f) in enumerate(zip(meas,e_meas,units,photbands,source,flags)):
            photsys,photband = p.split('.')
            if filters.is_color(p):
                to_unit = 'flux_ratio'
                color = True
            else:
                to_unit = _to_unit
                color = False
            if e_m>0:
                cm,e_cm = conversions.convert(u,to_unit,m,e_m,photband=p)
            else:
                cm,e_cm = conversions.convert(u,to_unit,m,photband=p),np.nan
            eff_wave = filters.eff_wave(p)
            extra_master['cmeas'][i] = cm
            extra_master['e_cmeas'][i] = e_cm
            extra_master['cwave'][i] = eff_wave
            extra_master['cunit'][i] = to_unit
            extra_master['color'][i] = color
            extra_master['include'][i] = True
            extra_master['meas'][i] = meas[i]
            extra_master['e_meas'][i] = e_meas[i]
            extra_master['unit'][i] = units[i]
            extra_master['photband'][i] = photbands[i]
            extra_master['source'][i] = source[i]
            extra_master['flag'][i] = flags[i]
            extra_master['phase'][i] = phases[i]
            if 'bibcode' in extra_master.dtype.names:
                extra_master['bibcode'][i] = '-'
            if 'comments' in extra_master.dtype.names:
                extra_master['comments'][i] = '-'

        logger.info('Original measurements:\n%s'%(photometry2str(self.master)))
        logger.info('Appending:\n%s'%(photometry2str(extra_master)))
        self.master = fix_master(np.hstack([self.master,extra_master]),e_default=0.1)
        #self.master = np.hstack([self.master,extra_array])
        logger.info('Final measurements:\n%s'%(photometry2str(self.master)))

    #}
    #{ Additional information

    def is_target(self,name):
        """
        Check if this SED represents the object ``name``.

        Purpose: solve alias problems. Maybe the ID is 'HD129929', and you are
        checking for "V* V836 Cen", which is the same target.

        :param name: object name
        :type name: str
        :return: True if this object instance represent the target "name".
        :rtype: bool
        """
        try:
            info = sesame.search(name)
            oname = info['oname']
        except:
            logger.warning('Did not find {:s} on Simbad}'.format(name))
            return False
        if oname==self.info['oname']:
            return True

    def has_photfile(self):
        """
        Check if this SED has a phot file.

        :return: True if this object instance has a photfile
        :rtype: bool
        """
        return os.path.isfile(self.photfile)

    def get_distance_from_plx(self,plx=None,lutz_kelker=True,unit='pc'):
        """
        Get the probability density function for the distance, given a parallax.

        If no parallax is given, the catalogue value will be used. If that one
        is also not available, a ValueError will be raised.

        Parallax must be given in mas.

        If the parallax is a float without uncertainty, an uncertainty of 0 will
        be assumed.

        If C{lutz_kelker=True}, a probability density function will be given for
        the distance, otherwise the distance is computed from simply inverting
        the parallax.

        Distance is returned in parsec (pc).

        :return: distance
        :rtype: (float,float)
        """
        # -- we need parallax and galactic position
        if plx is None and 'plx' in self.info:
            plx = self.info['plx']['v'],self.info['plx']['e']
        elif plx is None:
            raise ValueError('distance cannot be computed from parallax (no parallax)')
        if not 'galpos' in self.info:
            raise ValueError('distance cannot be computed from parallax (no position)')
        gal = self.info['galpos']
        # -- if the parallax has an uncertainty and the Lutz-Kelker bias needs
        #   to be taken into account, compute probability density function
        if lutz_kelker and hasattr(plx,'__iter__'):
            d = np.logspace(np.log10(0.1),np.log10(25000),100000)
            dprob = distance.distprob(d,gal[1],plx)
            dprob = dprob / dprob.max()
            logger.info('Computed distance to target with Lutz-Kelker bias')
            return d,dprob
        # -- just invert (with or without error)
        else:
            dist = (conversions.Unit(plx,'kpc-1')**-1).convert(unit)
            logger.info('Computed distance to target via parallax inversion: {:s}'.format(dist))
            return dist.get_value()

    def get_interstellar_reddening(self,distance=None, Rv=3.1):
        """
        Under construction.
        """
        gal = self.info['galpos']
        if distance is None:
            distance = self.get_distance_from_plx(lutz_kelker=False,unit='pc')[0]
        output = {}
        for model in ['arenou','schlegel','drimmel','marshall']:
            ext = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance)
            if ext is not None:
                output[model] = ext/Rv
        return output

    def get_angular_diameter(self):
        """
        Under construction.
        """
        raise NotImplementedError

    def compute_distance(self,mtype='igrid_search'):
        """
        Compute distance from radius and angular diameter (scale factor).

        The outcome is added to the C{results} attribute::

            distance = r/sqrt(scale)

        This particularly useful when you added constraints from solar-like
        oscillations (:class:`add_constraint_slo`).

        :return: distance,uncertainty (pc)
        :rtype: (float,float)
        """
        grid = self.results[mtype]['grid']
        radius = grid['radius']
        e_radius = grid['e_radius']
        scale = grid['scale']
        e_scale = grid['escale']

        radius = conversions.unumpy.uarray([radius,e_radius])
        scale = conversions.unumpy.uarray([scale,e_scale])

        distance = radius/scale**0.5
        distance = conversions.convert('Rsol','pc',distance)
        distance,e_distance = conversions.unumpy.nominal_values(distance),\
                      conversions.unumpy.std_devs(distance)
        self.results[mtype]['grid'] = pl.mlab.rec_append_fields(grid,\
                        ['distance','e_distance'],[distance,e_distance])
        d,e_d = distance.mean(),distance.std()
        logger.info("Computed distance from R and scale: {0}+/-{1} pc".format(d,e_d))
        return d,e_d




    #}

    def generate_ranges(self,start_from='igrid_search',distribution='uniform',**kwargs):   #type='single',
        """
        Generate sensible search range for each parameter.
        """
        limits = {}
        exist_previous = (start_from in self.results and 'CI' in self.results[start_from])

        # -- run over all parnames, and generate the parameters ranges for each
        #   component:
        #   (1) if a previous parameter search is done and the user has not set
        #       the parameters range, derive a sensible parameter range from the
        #       previous results.
        #   (2) if no previous search is done, use infinite upper and lower bounds
        #   (3) if ranges are given, use those. Cycle over them, and use the
        #       strategy from (1) if no parameter ranges are given for a specific
        #       component.
        for par_range_name in kwargs:
            # -- when "teffrange" is given, par_range will be the value of the
            #   keyword "teffrange", and parname will be "teff".
            parname = par_range_name.rsplit('range')[0]
            parrange = kwargs.get(par_range_name,None)
            # -- the three cases described above:
            if exist_previous and parrange is None:
                lkey = parname+'_l'
                ukey = parname+'_u'
                # -- if the parameters was not used in the fit, stick to the
                #   default given value
                if not lkey in self.results[start_from]['CI']:
                    parrange = kwargs[par_range_name]
                # -- else we can derive a better parameter range
                else:
                    parrange = (self.results[start_from]['CI'][lkey],
                                        self.results[start_from]['CI'][ukey])
            elif parrange is None:
                parrange = (-np.inf,np.inf)

            # -- now, the ranges are (lower,upper) for the uniform distribution,
            #   and (mu,scale) for the normal distribution
            if distribution=='normal':
                parrange = ((i[1]-i[0])/2.,(i[1]-i[0])/6.)
            elif distribution!='uniform':
                raise NotImplementedError('Any distribution other than "uniform" and "normal" has not been implemented yet!')
            limits[par_range_name] = parrange
        # -- this returns the kwargs but with filled in limits, and confirms
        #   the type if it was given, or gives the type when it needed to be derived
        logger.info('Parameter ranges calculated starting from {0:s} and using distribution {1:s}.'.format(start_from,distribution))
        return limits

    #{ Fitting routines

    def clip_grid(self,mtype='igrid_search',CI_limit=None):
        """
        Clip grid on CI limit, to save memory.

        :param mtype: type or results to clip
        :type mtype: str
        :param CI_limit: confidence limit to clip on
        :type CI_limit: float (between 0 (clips everything) and 1 (clips nothing))
        """
        if CI_limit is None:
            CI_limit = self.CI_limit
        new_grid = self.results[mtype]['grid']
        new_grid = new_grid[new_grid['ci_red']<=CI_limit]
        self.results[mtype]['grid'] = new_grid
        logger.info("Clipped grid at {:.6f}%".format(CI_limit*100))

    def collect_results(self, grid=None, fitresults=None, mtype='igrid_search', selfact='chisq', **kwargs):
        """creates record array(s) of all fit results and removes the failures"""

        gridnames = sorted(grid.keys())
        fitnames = sorted(fitresults.keys())
        pardtypes = [(name,'f8') for name in gridnames]+[(name,'f8') for name in fitnames]
        pararrays = [grid[name] for name in gridnames]+[fitresults[name] for name in fitnames]

        grid_results = np.rec.fromarrays(pararrays,dtype=pardtypes)

        # -- exclude failures
        failures = np.isnan(grid_results[selfact])
        if sum(failures):
            logger.info('Excluded {0} failed results (nan)'.format(sum(failures)))
            grid_results = grid_results[~failures]

        # -- make room for chi2 statistics
        grid_results = mlab.rec_append_fields(grid_results, 'ci_raw', np.zeros(len(grid_results)))
        grid_results = mlab.rec_append_fields(grid_results, 'ci_red', np.zeros(len(grid_results)))

        # -- take the previous results into account if they exist:
        if not mtype in self.results:
            self.results[mtype] = {}
        elif 'grid' in self.results[mtype]:
            logger.info('Appending previous results ({:d}+{:d})'.format(len(self.results[mtype]['grid']),len(grid_results)))
            ex_names = grid_results.dtype.names
            ex_grid = np.rec.fromarrays([self.results[mtype]['grid'][exname] for exname in ex_names],
                                        names=ex_names)
            grid_results = np.hstack([ex_grid,grid_results])

        # -- inverse sort according to chisq: this means the best models are at the end
        #   (mainly for plotting reasons, so that the best models are on top).
        sa = np.argsort(grid_results[selfact])[::-1]
        grid_results = grid_results[sa]

        self.results[mtype]['grid'] = grid_results

        logger.info('Total of %d grid points, best chisq=%s'%(len(grid_results), grid_results[-1]))

    def calculateDF(self, **ranges):
        """ Calculates the degrees of freedom from the given ranges"""

        df, df_info = 1, ['theta']
        for range_name in ranges:
            if re.search('ebv\d?range$', range_name):
                if not 'ebv' in df_info:
                    df += 1
                    df_info.append('ebv')
                else:
                    continue
            elif not np.allclose(ranges[range_name][0],ranges[range_name][1]):
                df += 1
                df_info.append(range_name[0:-5])
        df_info.sort()
        logger.info('Degrees of freedom = {} ({})'.format(df,', '.join(df_info)))

        return df, df_info

    def calculate_statistics(self, df=None, ranges=None, mtype='igrid_search', selfact='chisq'):
        """
        Calculates the Chi2 and reduced Chi2 based on nr of observations and degrees of
        freedom. If df is not provided, tries to calculate df from the fitting ranges.
        """

        # -- If nessessary calculate degrees of freedom from the ranges
        if df == None and ranges != None:
            df,df_info = self.calculateDF(**ranges)
        elif df == None:
            logger.warning('Cannot compute degrees of freedom!!! CHI2 might not make sense. (using df=5)')
            df = 5

        # -- Do the statistics
        N = sum(self.master['include'])
        k = N-df
        if k<=0:
            logger.warning('Not enough data to compute CHI2: it will not make sense')
            k = 1
        logger.info('Calculated statistics based on df={0} and Nobs={1}'.format(df,N))

        # -- Rescale if needed and compute confidence intervals
        results = self.results[mtype]['grid']
        factor = max(results[selfact][-1]/k,1)
        logger.warning('CHI2 rescaling factor equals %g'%(factor))
        results['ci_raw'] = scipy.stats.distributions.chi2.cdf(results[selfact],k)
        results['ci_red'] = scipy.stats.distributions.chi2.cdf(results[selfact]/factor,k)

        # -- Store the results
        self.results[mtype]['grid'] = results
        self.results[mtype]['factor'] = factor

    def calculate_confidence_intervals(self,mtype='igrid_search',chi2_type='red',CI_limit=None):
        """
        Compute confidence interval of all columns in the results grid.

        :param mtype: type of results to compute confidence intervals of
        :type mtype: str
        :param chi2_type: type of chi2 (raw or reduced)
        :type chi2_type: str ('raw' or 'red')
        :param CI_limit: confidence limit to clip on
        :type CI_limit: float (between 0 (clips everything) and 1 (clips nothing))
        """
        # -- get some info
        grid_results = self.results[mtype]['grid']
        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit
        # -- the chi2 grid is ordered: where is the value closest to the CI limit?
        region = self.results[mtype]['grid']['ci_'+chi2_type]<=CI_limit
        if sum(region)==0:
            raise ValueError("No models in the sample have a chi2_{} below the limit {}. Try increasing the number of models, increasing the CI_limit or choose different photometry.".format(chi2_type,CI_limit))
        # -- now compute the confidence intervals
        cilow, cihigh, value = [],[],[]
        for name in grid_results.dtype.names:
            cilow.append(grid_results[name][region].min())
            value.append(grid_results[name][-1])
            cihigh.append(grid_results[name][region].max())

        return dict(name=grid_results.dtype.names, value=value, cilow=cilow, cihigh=cihigh)

    def store_confidence_intervals(self, mtype='igrid_search', name=None, value=None, cilow=None, \
                                   cihigh=None, **kwargs):
        """
        Saves the provided confidence intervals in the result dictionary of self.
        The provided confidence intervals will be saved in a dictionary:
        self.results[mtype]['CI'] with as key for the value the name, for cilow:
        name_l and for cihigh: name_u.

        :param mtype: the search type
        :type mtype: str
        :param name: names of the parameters
        :type name: array
        :param value: best fit values
        :type value: array
        :param cilow: lower CI limit
        :type cilow: array
        :param cihigh: upper CI limit
        :type cihigh: array
        """
        if not 'CI' in self.results[mtype]:
            self.results[mtype]['CI'] = {}
        for name, val, cil, cih in zip(name, value, cilow, cihigh):
            self.results[mtype]['CI'][name+'_l'] = cil
            self.results[mtype]['CI'][name] = val
            self.results[mtype]['CI'][name+'_u'] = cih
            try:
                logger.info('CI %s: %g <= %g <= %g'%(name,cil,val,cih))
            except Exception:
                logger.info('CI %s: nan <= %g <= nan'%(name,val))


    def igrid_search(self,points=100000,teffrange=None,loggrange=None,ebvrange=None,
                          zrange=(0,0),rvrange=(3.1,3.1),vradrange=(0,0),
                          df=None,CI_limit=None,set_model=True,exc_interpolpar=[],**kwargs):
        """
        Fit fundamental parameters using a (pre-integrated) grid search.

        If called consecutively, the ranges will be set to the CI_limit of previous
        estimations, unless set explicitly.

        If called for the first time, the ranges will be +/- np.inf by defaults,
        unless set explicitly.

        There is an option to exclude a certain parameter from the interpolation. This is done by
        including it in a list attached to the keyword 'exc_interpolpar', e.g. exc_interpolpar = ['z'].
        """
        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit
        # -- set defaults limits
        ranges = self.generate_ranges(teffrange=teffrange,\
                             loggrange=loggrange,ebvrange=ebvrange,\
                             zrange=zrange,rvrange=rvrange,vradrange=vradrange)
        # -- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'%(photometry2str(self.master[include_grid])))
        # -- an initial check on the conditions for exclusion of a parameter from interpolation
        if not exc_interpolpar == None:
            for var in exc_interpolpar:
                if ranges[var+'range'][0] != ranges[var+'range'][1]:
                    raise IOError('Exclusion of parameters from interpolation is only possible if the lower and upper ranges of those ranges are equal to an actual grid point.')
        # -- build the grid, run over the grid and calculate the CHI2
        pars = fit.generate_grid_pix(self.master['photband'][include_grid],points=points,**ranges)
        pars['exc_interpolpar'] = exc_interpolpar
        #for var in ['teff','logg','ebv','z','rv','vrad']:
            #if ranges[var+'range'][0] == ranges[var+'range'][1]:
                #insertpars[var] = ranges[var+'range'][0]

        chisqs,scales,e_scales,lumis = fit.igrid_search_pix(self.master['cmeas'][include_grid],
                             self.master['e_cmeas'][include_grid],
                             self.master['photband'][include_grid],**pars)
        fitres = dict(chisq=chisqs, scale=scales, escale=e_scales, labs=lumis)
        pars.pop('exc_interpolpar')
        # -- collect all the results in a record array
        self.collect_results(grid=pars, fitresults=fitres, mtype='igrid_search')
        # -- do the statistics
        self.calculate_statistics(df=df, ranges=ranges, mtype='igrid_search')
        # -- compute the confidence intervals
        ci = self.calculate_confidence_intervals(mtype='igrid_search',chi2_type='red',CI_limit=CI_limit)
        self.store_confidence_intervals(mtype='igrid_search', **ci)
        # -- remember the best model
        if set_model: self.set_best_model()
    def generate_fit_param(self, start_from='igrid_search', **pars):
        """
        generates a dictionary with parameter information that can be handled by fit.iminimize
        """
        result = dict()
        for key in list(pars.keys()):
            if re.search("range$", key):
                result[key[0:-5]+"_min"] =  pars[key][0]
                result[key[0:-5]+"_max"] =  pars[key][1]
                # -- if min == max: fix the parameter and force value = min
                if np.allclose([pars[key][0]], [pars[key][1]]):
                    result[key[0:-5]+"_vary"] = False
                    result[key[0:-5]+"_value"] = pars[key][0]
                else:
                    result[key[0:-5]+"_vary"] = True
            else:
                # -- Store the startvalue. If None, look in start_from for a new startvalue.
                if pars[key] != None and not key+"_value" in result:
                    result[key+"_value"] = pars[key]
                elif pars[key] == None and not key+"_value" in result:
                    result[key+"_value"] = self.results[start_from]['CI'][key]

        return result

    def calculate_iminimize_CI(self, mtype='iminimize', CI_limit=0.66, **kwargs):

        # -- Get the best fit parameters and ranges
        pars = {}
        skip = ['scale', 'chisq', 'nfev', 'labs', 'ci_raw', 'ci_red', 'scale', 'escale']
        for name in list(self.results[mtype]['CI'].keys()):
            name = re.sub('_[u|l]$', '', name)
            if not name in pars and not name in skip:
                pars[name] = self.results[mtype]['CI'][name]
                pars[name+"range"] = [self.results[mtype]['CI'][name+"_l"],\
                                       self.results[mtype]['CI'][name+"_u"]]
        pars.update(kwargs)
        pars = self.generate_fit_param(**pars)

        # -- calculate the confidence intervalls
        include_grid = self.master['include']
        ci = fit.calculate_iminimize_CI(self.master['cmeas'][include_grid],
                             self.master['e_cmeas'][include_grid],
                             self.master['photband'][include_grid],
                             CI_limit=CI_limit,constraints=self.constraints, **pars)

        # -- Add the scale factor
        ci['name'] = np.append(ci['name'], 'scale')
        ci['value'] = np.append(ci['value'], self.results[mtype]['grid']['scale'][-1])
        ci['cilow'] = np.append(ci['cilow'], min(self.results[mtype]['grid']['scale']))
        ci['cihigh'] = np.append(ci['cihigh'], max(self.results[mtype]['grid']['scale']))

        logger.info('Calculated %s%% confidence intervalls for all parameters'%(CI_limit))

        self.store_confidence_intervals(mtype='iminimize', **ci)

    def calculate_iminimize_CI2D(self,xpar, ypar, mtype='iminimize', limits=None, res=10, **kwargs):

        # -- get the best fitting parameters
        pars = {}
        skip = ['scale', 'chisq', 'nfev', 'labs', 'ci_raw', 'ci_red', 'scale', 'escale']
        for name in list(self.results[mtype]['CI'].keys()):
            name = re.sub('_[u|l]$', '', name)
            if not name in pars and not name in skip:
                pars[name] = self.results[mtype]['CI'][name]
                pars[name+"range"] = [self.results[mtype]['CI'][name+"_l"],\
                                       self.results[mtype]['CI'][name+"_u"]]
        pars.update(kwargs)
        pars = self.generate_fit_param(**pars)

        logger.info('Calculating 2D confidence intervalls for %s--%s'%(xpar,ypar))

        # -- calculate the confidence intervalls
        include_grid = self.master['include']
        x,y,chi2, raw, red = fit.calculate_iminimize_CI2D(self.master['cmeas'][include_grid],
                             self.master['e_cmeas'][include_grid],
                             self.master['photband'][include_grid],
                             xpar, ypar, limits=limits, res=res,
                             constraints=self.constraints, **pars)

        # -- store the CI
        if not 'CI2D' in self.results[mtype]: self.results[mtype]['CI2D'] = {}
        self.results[mtype]['CI2D'][xpar+"-"+ypar] = dict(x=x, y=y, ci_chi2=chi2,\
                                                           ci_raw=raw, ci_red=red)


    def _get_imin_ci(self, mtype='iminimize',**ranges):
        """ returns ci information for store_confidence_intervals """
        names, values, cil, cih = [],[],[],[]
        for key in list(ranges.keys()):
            name = key[0:-5]
            names.append(name)
            values.append(self.results[mtype]['grid'][name][-1])
            cil.append(ranges[key][0])
            cih.append(ranges[key][1])
        names.append('scale')
        values.append(self.results[mtype]['grid']['scale'][-1])
        cil.append(min(self.results[mtype]['grid']['scale']))
        cih.append(max(self.results[mtype]['grid']['scale']))
        return dict(name=names, value=values, cilow=cil, cihigh=cih)

    def iminimize(self, teff=None, logg=None, ebv=None, z=0, rv=3.1, vrad=0, teffrange=None,
                     loggrange=None, ebvrange=None, zrange=None, rvrange=None, vradrange=None,
                     points=None, distance=None, start_from='igrid_search',df=None, CI_limit=None,
                     calc_ci=False, set_model=True, **kwargs):
        """
        Basic minimizer method for SED fitting implemented using the lmfit library from sigproc.fit
        """

        # -- set defaults limits and starting values
        ranges = self.generate_ranges(teffrange=teffrange,loggrange=loggrange,\
                             ebvrange=ebvrange,zrange=zrange,rvrange=rvrange,\
                             vradrange=vradrange, start_from=start_from)
        pars = self.generate_fit_param(teff=teff, logg=logg, ebv=ebv, z=z, \
                            rv=rv, vrad=vrad, start_from=start_from, **ranges)
        fitkws = dict(distance=distance) if distance != None else dict()

        # -- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'% \
        (photometry2str(self.master[include_grid])))

        # -- pass all ranges and starting values to the fitter
        grid, chisq, nfev, scale, lumis = fit.iminimize(self.master['cmeas'][include_grid],
                                            self.master['e_cmeas'][include_grid],
                                            self.master['photband'][include_grid],
                                            fitkws=fitkws, points=points,**pars)

        logger.info('Minimizer Succes with startpoints=%s, chi2=%s, nfev=%s'%(len(chisq), chisq[0], nfev[0]))
        # -- handle the results
        fitres = dict(chisq=chisq, nfev=nfev, scale=scale, labs=lumis)
        self.collect_results(grid=grid, fitresults=fitres, mtype='iminimize', selfact='chisq')

        # -- Do the statistics
        self.calculate_statistics(df=5, ranges=ranges, mtype='iminimize', selfact='chisq')

        # -- Store the confidence intervals
        ci = self._get_imin_ci(mtype='iminimize',**ranges)
        self.store_confidence_intervals(mtype='iminimize', **ci)
        if calc_ci:
            self.calculate_iminimize_CI(mtype='iminimize', CI_limit=CI_limit)

        # -- remember the best model
        if set_model: self.set_best_model(mtype='iminimize')


    def imc(self,teffrange=None,loggrange=None,ebvrange=None,zrange=None,start_from='igrid_search',\
                 distribution='uniform',points=None,fitmethod='fmin',disturb=True):

        limits = self.generate_ranges(teffrange=teffrange,loggrange=loggrange,\
                                      ebvrange=ebvrange,zrange=zrange,distribution=distribution,\
                                      start_from=start_from)

        # -- grid search on all include data: extract the best CHI2
        include = self.master['include']
        meas = self.master['cmeas'][include]
        if disturb:
            emeas = self.master['e_cmeas'][include]
        else:
            emeas = self.master['e_cmeas'][include]*1e-6
        photbands = self.master['photband'][include]
        logger.info('The following measurements are included in the fitting process:\n%s'%(photometry2str(self.master[include])))

        # -- generate initial guesses
        teffs,loggs,ebvs,zs,radii = fit.generate_grid(self.master['photband'][include],type='single',points=points+25,**limits)
        NrPoints = len(teffs)>points and points or len(teffs)
        firstoutput = np.zeros((len(teffs)-NrPoints,9))
        output = np.zeros((NrPoints,9))

        # -- fit the original data a number of times
        for i,(teff,logg,ebv,z) in enumerate(zip(teffs[NrPoints:],loggs[NrPoints:],ebvs[NrPoints:],zs[NrPoints:])):
            try:
                fittedpars,warnflag = fit.iminimize2(meas,emeas,photbands,teff=teff,logg=logg,ebv=ebv,z=z,fitmethod=fitmethod)
                firstoutput[i,1:] = fittedpars
                firstoutput[i,0] = warnflag
            except IOError:
                firstoutput[i,0] = 3

        logger.info("{0}/{1} fits on original data failed (max func call)".format(sum(firstoutput[:,0]==1),firstoutput.shape[0]))
        logger.info("{0}/{1} fits on original failed (max iter)".format(sum(firstoutput[:,0]==2),firstoutput.shape[0]))
        logger.info("{0}/{1} fits on original data failed (outside of grid)".format(sum(firstoutput[:,0]==3),firstoutput.shape[0]))

        # -- retrieve the best fitting result and make it the first entry of output
        keep = (firstoutput[:,0]==0) & (firstoutput[:,1]>0)
        best = firstoutput[keep,-2].argmin()
        output[-1,:] = firstoutput[keep][best,:]

        # calculate the factor with which to multiply the scale
        #factor = np.sqrt(output[-1,5]/len(meas))
        #print factor

        # -- now do the actual Monte Carlo simulation
        for i,(teff,logg,ebv,z) in enumerate(zip(teffs[:NrPoints-1],loggs[:NrPoints-1],ebvs[:NrPoints-1],zs[:NrPoints-1])):
            newmeas = meas + np.random.normal(scale=emeas) #*factor)
            try:
                fittedpars,warnflag = fit.iminimize2(newmeas,emeas,photbands,teff,logg,ebv,z,fitmethod=fitmethod)
                output[i,1:] = fittedpars
                output[i,0] = warnflag
            except IOError:
                output[i,0] = 3

        logger.info("{0}/{1} MC simulations failed (max func call)".format(sum(output[:,0]==1),NrPoints))
        logger.info("{0}/{1} MC simulations failed (max iter)".format(sum(output[:,0]==2),NrPoints))
        logger.info("{0}/{1} MC simulations failed (outside of grid)".format(sum(output[:,0]==3),NrPoints))

        # -- remove nonsense results
        keep = (output[:,0]==0) & (output[:,1]>0)
        output = output[keep]
        output = np.rec.fromarrays(output[:,1:-1].T,names=['teff','logg','ebv','z','labs','chisq','scale'])
        # -- derive confidence intervals and median values
        #print np.median(output[:,1:],axis=0)

        self.results['imc'] = {}
        self.results['imc']['CI'] = {}
        self.results['imc']['grid'] = output
        CI_limit = 0.95
        for name in output.dtype.names:
            sortarr = np.sort(output[name])
            trimarr = scipy.stats.trimboth(sortarr,(1-CI_limit)/2.) # trim 2.5% of top and bottom, to arrive at 95% CI
            self.results['imc']['CI'][name+'_l'] = trimarr.min()
            self.results['imc']['CI'][name] = output[name][-1]#np.median(output[name])
            self.results['imc']['CI'][name+'_u'] = trimarr.max()
            logger.info('%i%% CI %s: %g <= %g <= %g'%(CI_limit*100,name,self.results['imc']['CI'][name+'_l'],
                                                           self.results['imc']['CI'][name],
                                                           self.results['imc']['CI'][name+'_u']))
        self.set_best_model(mtype='imc')


    #}

    #{ Interfaces

    def clear(self):
        """
        Clear the results.
        """
        self.results = {}

    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004',**kwargs):
        """
        Get reddenend and unreddened model
        """
        logger.info('Interpolating approximate full SED of best model')

        # -- synthetic flux
        include = self.master['include']
        synflux = np.zeros(len(self.master['photband']))
        keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
        keep = keep & include

        if ('igrid_search' in mtype) | ('iminimize' in mtype): #mtype in ['igrid_search', 'iminimize']:
            # -- get the metallicity right
            files = model.get_file(z='*')
            if type(files) == str: files = [files] #files needs to be a list!
            metals = np.array([pf.getheader(ff)['Z'] for ff in files])
            metals = metals[np.argmin(np.abs(metals-self.results[mtype]['CI']['z']))]
            scale = self.results[mtype]['CI']['scale']
            # -- get (approximated) reddened and unreddened model
            wave,flux = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                    logg=self.results[mtype]['CI']['logg'],
                                    ebv=self.results[mtype]['CI']['ebv'],
                                    z=metals,
                                    law=law)
            wave_ur,flux_ur = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                        logg=self.results[mtype]['CI']['logg'],
                                        ebv=0,
                                        z=metals,
                                        law=law)
            # -- get synthetic photometry
            synflux_,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
                               logg=self.results[mtype]['CI']['logg'],
                               ebv=self.results[mtype]['CI']['ebv'],
                               z=self.results[mtype]['CI']['z'],
                               photbands=self.master['photband'][keep])

            flux,flux_ur = flux*scale,flux_ur*scale

            synflux[keep] = synflux_

            #synflux,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
            #                          logg=self.results[mtype]['CI']['logg'],
            #                          ebv=self.results[mtype]['CI']['ebv'],
            #                          photbands=self.master['photband'])
            synflux[~self.master['color']] *= scale
            chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
            # -- calculate effective wavelengths of the photometric bands via the model
            #   values
            eff_waves = filters.eff_wave(self.master['photband'],model=(wave,flux))
            self.results[mtype]['model'] = wave,flux,flux_ur
            self.results[mtype]['synflux'] = eff_waves,synflux,self.master['photband']
            self.results[mtype]['chi2'] = chi2

    def set_model(self,wave,flux,label='manual',unique_phase=None):
        """
        Manually set the best SED model.

        This is particularly useful when working with calibrators for which there
        is an observed SED.

        The label will be used to store the model in the C{results} attribute.

        :param wave: wavelength array (angstrom)
        :type wave: ndarray
        :param flux: flux array (erg/s/cm2/AA)
        :type flux: ndarray
        :param label: key used to store the model and synthetic photometry in C{results}
        :type label:s str
        """
        # -- necessary information
        photbands = self.master['photband']
        is_color = self.master['color']
        if unique_phase != None:
            include = self.master['include'] & (self.master['phase'] == unique_phase)
        else:
            include = self.master['include']
        synflux = np.zeros(len(photbands))

        # -- compute synthetic photometry
        synflux[(~is_color) & include] = model.synthetic_flux(wave,flux,photbands[(~is_color) & include])
        synflux[is_color & include] = model.synthetic_color(wave,flux,photbands[is_color & include])
        chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
        eff_waves = filters.eff_wave(photbands)

        if not label in self.results:
            self.results[label] = {}
        self.results[label]['model'] = wave,flux,flux
        self.results[label]['synflux'] = eff_waves,synflux,photbands
        self.results[label]['chi2'] = chi2
        logger.debug('Stored model SED in {}'.format(label))

    def get_model(self,label='igrid_search'):
        """
        Retrieve the best SED model.

        B{Warning}: the grid search interpolation is also done with interpolation
        in metallicity, while this is not the case for the best full SED model.
        Interpolation of the full SED is only done in teff/logg, and the
        metallicity is chosen to be equal to the closest grid point. On the
        other hand, while the reddening value is interpolated in the grid search,
        this is not the case for the best model (the model is simply reddened
        according the found best value). So don't be surprised if you find
        differences between the values of C{self.results[label]['synflux']},
        which are the value sfor the interpolated photometric points of the
        best model, and the synthetic photometry obtained by manual integration
        of the returned full SED model. Those differences should be incredible
        small and mainly due to the metallicity.
        """
        wave,flux,urflux = self.results[label]['model']
        return wave,flux,urflux

    def sample_gridsearch(self,mtype='igrid_search',NrSamples=1,df=None,selfact='chisq'):
        """
        Retrieve an element from the results of a grid search according to the derived probability.

        An array of length "NrSamples" containing 'grid-indices' is returned, so the actual parameter values
        of the corresponding model can be retrieved from the results dictionary.

        :param NrSamples: the number of samples you wish to draw
        :type NrSamples: int
        """
        # -- this function is only checked to work with the results of an igrid_search
        if not 'igrid_search' in mtype:
            return None

        ranges = self.generate_ranges(teffrange=None,loggrange=None,ebvrange=None,\
                            zrange=None,rvrange=None,vradrange=None,start_from=mtype)

        # -- If nessessary calculate degrees of freedom from the ranges
        if df == None and ranges != None:
            df,df_info = self.calculateDF(**ranges)
        elif df == None:
            logger.warning('Cannot compute degrees of freedom!!! CHI2 might not make sense. (using df=5)')
            df = 5

        # -- Do the statistics
        N = sum(self.master['include'])
        k = N-df
        if k<=0:
            logger.warning('Not enough data to compute CHI2: it will not make sense')
            k = 1
        logger.info('Statistics based on df={0} and Nobs={1}'.format(df,N))
        factor = max(self.results[mtype]['grid'][selfact][-1]/k,1)

        # -- Compute the pdf and cdf
        probdensfunc = scipy.stats.distributions.chi2.pdf(self.results[mtype]['grid'][selfact]/factor,k)
        cumuldensfunc = probdensfunc.cumsum()

        # -- Uniformly sample the cdf, to get a sampling according to the pdf
        sample = pl.uniform(cumuldensfunc[0],cumuldensfunc[-1],NrSamples)
        indices = np.zeros(NrSamples,int)
        for i in range(NrSamples):
            indices[i] = (abs(cumuldensfunc-sample[i])).argmin()
        return indices

    def chi2(self,select=None,reduced=False,label='igrid_search'):
        """
        Calculate chi2 of best model.

        TEMPORARY!!! API WILL CHANGE!!!

        This function is not called anywhere in the builder.
        """
        # -- select photometry
        master = self.master.copy()
        keep = np.ones(len(master),bool)
        if isinstance(select,dict):
            for key in select:
                keep = keep & (master[key]==select[key])
        master = master[keep]
        # -- compute synthetic phtoometry
        photbands = master['photband']
        is_color = master['color']
        synflux = np.zeros(len(photbands))

        wave,flux,urflux = self.get_model(label=label)
        synflux[~is_color] = model.synthetic_flux(wave,flux,photbands[~is_color])
        synflux[is_color] = model.synthetic_color(wave,flux,photbands[is_color])
        x2 = (master['cmeas']-synflux)**2/master['e_cmeas']**2
        x2 = x2.mean()
        return x2

    #}
    #{ Add constraints

    def add_constraint_distance(self,distance=None,mtype='igrid_search',**kwargs):
        """
        Use the distance to compute additional information.

        Compute radii, absolute luminosities and masses, and add them to the
        results.

        Extra kwargs go to :class:`get_distance_from_plx`.

        B{Warning:} after calling this function, the C{labs} column in the grid
        is actually absolutely calibrated and reflects the true absolute
        luminosity instead of the absolute luminosity assuming 1 solar radius.

        :param distance: distance in solar units and error
        :type distance: tuple (float,float)
        :param mtype: type of results to add the information to
        :type mtype: str
        """
        if distance is None:
            kwargs['lutz_kelker'] = False # we can't handle asymmetric error bars
            kwargs['unit'] = 'Rsol'
            distance = self.get_distance_from_plx(**kwargs)

        # -- compute the radius, absolute luminosity and mass: note that there is
        #   an uncertainty on the scaling factor *and* on the distance!
        scale = conversions.unumpy.uarray([self.results[mtype]['grid']['scale'],\
                                           self.results[mtype]['grid']['escale']])
        distance = conversions.ufloat(distance)
        radius = distance*np.sqrt(self.results[mtype]['grid']['scale']) # in Rsol
        labs = self.results[mtype]['grid']['labs']*radius**2
        mass = conversions.derive_mass((self.results[mtype]['grid']['logg'],'[cm/s2]'),\
                                      (radius,'Rsol'),unit='Msol')

        # -- store the results in the grid
        mass,e_mass = conversions.unumpy.nominal_values(mass),\
                      conversions.unumpy.std_devs(mass)
        radius,e_radius = conversions.unumpy.nominal_values(radius),\
                      conversions.unumpy.std_devs(radius)
        labs,e_labs = conversions.unumpy.nominal_values(labs),\
                      conversions.unumpy.std_devs(labs)
        self.results[mtype]['grid']['labs'] = labs # Labs is already there, overwrite
        # -- check if the others are already in there or not:
        labels,data = ['e_labs','radius','e_radius','mass','e_mass'],\
                      [e_labs,radius,e_radius,mass,e_mass]
        if 'e_labs' in self.results[mtype]['grid'].dtype.names:
            for idata,ilabel in zip(data,labels):
                self.results[mtype]['grid'][ilabel] = idata
        else:
            self.results[mtype]['grid'] = pl.mlab.rec_append_fields(self.results[mtype]['grid'],\
                     labels,data)

        # -- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: distance (improved luminosity, added info on radius and mass)')

    def add_constraint_slo(self,numax,Deltanu0,mtype='igrid_search',chi2_type='red'):
        """
        Use diagnostics from solar-like oscillations to put additional constraints on the parameters.

        If the results are constrained with the distance before :class:`add_constraint_distance`,
        then these results are combined with the SLO constraints.

        :param numax: frequency of maximum amplitude
        :type numax: 3-tuple (value,error,unit)
        :param Deltanu0: large separation (l=0)
        :type Deltanu0: 3-tuple (value,error,unit)
        :param chi2_type: type of chi2 (raw or reduced)
        :type chi2_type: str ('raw' or 'red')
        """
        grid = self.results[mtype]['grid']
        # -- we need the teffs, so that we can compute the logg and radius using
        #   the solar-like oscillations information.
        teff = grid['teff']
        logg = grid['logg']
        pvalues = [1-grid['ci_'+chi2_type]]
        names = ['ci_'+chi2_type+'_phot']

        # -- we *always* have a logg, that's the way SEDs work:
        logg_slo = conversions.derive_logg_slo((teff,'K'),numax)
        logg_slo,e_logg_slo = conversions.unumpy.nominal_values(logg_slo),\
                      conversions.unumpy.std_devs(logg_slo)
        #   compute the probabilities: scale to standard normal
        logg_prob = scipy.stats.distributions.norm.sf(abs(logg_slo-logg)/e_logg_slo)

        pvalues.append(logg_prob)
        names.append('ci_logg_slo')

        # -- but we only sometimes have a radius (viz. when the distance is
        #   known). There is an uncertainty on that radius!
        radius_slo = conversions.derive_radius_slo(numax,Deltanu0,(teff,'K'),unit='Rsol')
        radius_slo,e_radius_slo = conversions.unumpy.nominal_values(radius_slo),\
                                  conversions.unumpy.std_devs(radius_slo)
        if 'radius' in grid.dtype.names:
            radius = grid['radius']
            e_radius = grid['e_radius']
            #total_error = np.sqrt( (e_radius_slo**2+e_radius**2))
            total_error = np.sqrt( (e_radius_slo**2+e_radius**2)/2. + (radius-radius_slo)**2/4.)
            radius_prob = scipy.stats.distributions.norm.sf(abs(radius_slo-radius)/total_error)
            pvalues.append(radius_prob)
            names.append('ci_radius_slo')
            # -- combined standard deviation and mean for two populations with
            #   possibly zero intersection (wiki: standard deviation)
            #total_error = np.sqrt( (e_radius_slo**2+e_radius**2)/2. + (radius-radius_slo)**2/4.)
            total_mean = (radius+radius_slo)/2.
            grid['radius'] = total_mean
            grid['e_radius'] = total_error
            logger.info('Added contraint: combined existing radius estimate with slo estimate')
        #   otherwise, we just add info on the radius
        else:
            labs = grid['labs']*conversions.unumpy.uarray([radius_slo,e_radius_slo])**2
            labs,e_labs = conversions.unumpy.nominal_values(labs),\
                          conversions.unumpy.std_devs(labs)
            grid['labs'] = labs
            mass = conversions.derive_mass((self.results[mtype]['grid']['logg'],'[cm/s2]'),\
                                           (radius_slo,e_radius_slo,'Rsol'),unit='Msol')
            mass,e_mass = conversions.unumpy.nominal_values(mass),\
                          conversions.unumpy.std_devs(mass)
            # -- now we can also derive the distance:
            scale = conversions.unumpy.uarray([self.results[mtype]['grid']['scale'],\
                                           self.results[mtype]['grid']['escale']])
            distance = radius_slo/conversions.sqrt(scale)
            distance,e_distance = conversions.unumpy.nominal_values(distance),\
                                  conversions.unumpy.std_devs(distance)
            grid = pl.mlab.rec_append_fields(grid,['radius','e_radius','e_labs','mass','e_mass','distance','e_distance'],\
                                        [radius_slo,e_radius_slo,e_labs,mass,e_mass,distance,e_distance])
            logger.info('Added constraint: {0:s} via slo, improved luminosity'.format(', '.join(['radius','e_radius','e_labs','mass','e_mass'])))

        # -- combine p values using Fisher's method
        combined_pvals = evaluate.fishers_method(pvalues)

        # -- add the information to the grid
        self.results[mtype]['grid'] = pl.mlab.rec_append_fields(grid,\
                     names,pvalues)

        # -- and replace the original confidence intervals, and re-order
        self.results[mtype]['grid']['ci_'+chi2_type] = 1-combined_pvals
        sa = np.argsort(self.results[mtype]['grid']['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = self.results[mtype]['grid'][sa]

        # -- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: {0:s} via slo and replaced ci_{1:s} with combined CI'.format(', '.join(names),chi2_type))

    def add_constraint_reddening(self,distance=None,ebv=None,e_ebv=0.1,Rv=3.1,\
                model=None,mtype='igrid_search',chi2_type='red',upper_limit=False):
        """
        Use reddening maps to put additional constraints on the parameters.

        This constraint assumes that, if C{upper_limit=False}, given the
        distance to target, the reddening from the reddening maps is similar to
        the one derived from the SED. All models where this is not the case will
        be deemed improbable.

        If C{upper_limit=True}, all models with an E(B-V) value above C{ebv}
        will be considered improbable.

        When you don't set C{ebv}, the value will, by default, be derived from
        Drimmel maps if C{upper_limit=False} and from Schlegel if
        C{upper_limit=True}. You can change this behaviour by setting C{model}
        manually.

        :param distance: distance and uncertainty in parsec
        :type distance: tuple (float,float)
        :param ebv: E(B-V) reddening in magnitude
        :type ebv: float
        :param e_ebv: error on the reddening in percentage
        :type e_ebv: float
        :param model: model reddening maps
        :type model: str
        :param mtype: type of results to add the information to
        :type mtype: str
        :param upper_limit: consider the E(B-V) value as an upper limit
        :type upper_limit: bool
        """
        # -- for upper limits on E(B-V), we best use Schlegel maps by default,
        #   otherwise we best use Drimmel maps.
        if model is None and upper_limit:
            model = 'schlegel'
        elif model is None:
            model = 'drimmel'
        # -- if I need to figure out the reddening myself, I need the distance!
        if ebv is None and distance is None:
            distance = self.get_distance_from_plx(lutz_kelker=False,unit='pc')
        # -- let me figure out the reddening if you didn't give it:
        if ebv is None:
            gal = self.info['galpos']
            ebv = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance[0])/Rv
            ebv_u = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance[0]-distance[1])/Rv
            ebv_l = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance[0]+distance[1])/Rv
            e_ebv = max(ebv-ebv_l,ebv_u-ebv)
        else:
            e_ebv = 0.1*ebv
        grid = self.results[mtype]['grid']
        ebvs = grid['ebv']
        # -- probabilities are differently calculated depending on upper limit.
        if not upper_limit:
            ebv_prob = scipy.stats.distributions.norm.cdf(abs(ebv-ebvs)/e_ebv)
            # -- combine p values using Fisher's method
            combined_pvals = evaluate.fishers_method([1-grid['ci_'+chi2_type],ebv_prob])
            grid['ci_'+chi2_type] = 1-combined_pvals
        else:
            grid['ci_'+chi2_type] = np.where(ebvs<=ebv,grid['ci_'+chi2_type],1.)

        # -- and replace the original confidence intervals, and re-order
        sa = np.argsort(grid['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = grid[sa]

        # -- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: E(B-V)={0}+/-{1}'.format(ebv,e_ebv))

    def add_constraint_angular_diameter(self,angdiam):
        raise NotImplementedError

    def add_constraint_mass(self,mass,mtype='igrid_search',chi2_type='red'):
        """
        Add constraints on the mass.

        C{mass} must be a tuple, if the second element is smaller than the first,
        it is assumed to be (mu,sigma) from a normal distribution. If the second
        element is larger than the first, it is assumed to be (lower, upper)
        from a uniform distribution.
        """
        normal = mass[0]>mass[1]
        grid = self.results[mtype]['grid']
        masses = grid['mass']
        if normal:
            mass_prob = scipy.stats.distributions.norm.cdf(abs(mass[0]-masses)/mass[1])
            # -- combine p values using Fisher's method
            combined_pvals = evaluate.fishers_method([1-grid['ci_'+chi2_type],mass_prob])
            grid['ci_'+chi2_type] = 1-combined_pvals
        else:
            grid['ci_'+chi2_type] = np.where((masses<=mass[0]) | (mass[1]<=masses),1.,grid['ci_'+chi2_type])

        # -- and replace the original confidence intervals, and re-order
        sa = np.argsort(grid['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = grid[sa]

        # -- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: mass={0}+/-{1}'.format(*mass))

    def add_constraint_evolution_models(self,models='siess2000',\
               ylabels=['age','labs','radius'],e_y=None,
               function='linear',mtype='igrid_search',chi2_type='red'):
        """
        Use stellar evolutionary models to put additional constraints on the parameters.
        """
        grid = self.results[mtype]['grid']

        # -- make sure y's and ylabels are iterable
        if isinstance(ylabels,str):
            ylabels = [ylabels]
        # -- cycle over all yvalues and compute the pvalues
        pvalues = [1-grid['ci_'+chi2_type]]
        add_info = []
        # -- create the evolutionary grid and interpolate the stellar evolutinary
        #   grid on the SED-integrated grid points
        output = np.array([evolutionmodels.get_itable(iteff,ilogg,iz) for iteff,ilogg,iz in zip(grid['teff'],grid['logg'],grid['z'])])
        output = np.rec.fromarrays(output.T,names=['age','labs','radius'])
        for label in ylabels:
            y_interpolated = output[label]
            add_info.append(y_interpolated)
            # -- only add the contraint when it is possible to do so.
            if not label in grid.dtype.names:
                logger.info("Cannot put constraint on {} (not a parameter)".format(label))
                continue
            y_computed = grid[label]
            # -- error on the y-value: either it is computed, it is given or it is
            #   assumed it is 10% of the value
            if e_y is None and ('e_'+label) in grid.dtype.names:
                e_y = np.sqrt(grid['e_'+label]**2+(0.1*y_interpolated)**2)
            elif e_y is None:
                e_y = np.sqrt((0.1*y_computed)**2+(0.1*y_interpolated)**2)
            y_prob = scipy.stats.distributions.norm.sf(abs(y_computed-y_interpolated)/e_y)
            pl.figure()
            pl.subplot(221)
            pl.title(label)
            pl.scatter(y_computed,y_interpolated,c=y_prob,edgecolors='none',cmap=pl.cm.spectral)
            pl.plot([pl.xlim()[0],pl.xlim()[1]],[pl.xlim()[0],pl.xlim()[1]],'r-',lw=2)
            pl.xlim(pl.xlim())
            pl.ylim(pl.xlim())
            pl.xlabel('Computed')
            pl.ylabel('Interpolated')
            pl.colorbar()
            pl.subplot(223)
            pl.title('y_interpolated')
            pl.scatter(grid['teff'],grid['logg'],c=y_interpolated,edgecolors='none',cmap=pl.cm.spectral)
            pl.colorbar()
            pl.subplot(224)
            pl.title('y_computed')
            pl.scatter(grid['teff'],grid['logg'],c=y_computed,edgecolors='none',cmap=pl.cm.spectral)
            pl.colorbar()
            pvalues.append(y_prob)

        pl.figure()
        pl.subplot(221)
        pl.title('p1')
        sa = np.argsort(pvalues[0])
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=pvalues[0][sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        pl.subplot(222)
        pl.title('p2')
        sa = np.argsort(pvalues[1])
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=pvalues[1][sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        pl.subplot(223)
        pl.title('p3')
        sa = np.argsort(pvalues[2])
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=pvalues[2][sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')

        # -- combine p values using Fisher's method
        combined_pvals = evaluate.fishers_method(pvalues)
        pl.subplot(224)
        pl.title('pcombined')
        sa = np.argsort(combined_pvals)
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=combined_pvals[sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        pl.show()

        # -- add the information to the grid (assume that if there one label
        #   not in there already, none of them are
        if not 'c_'+ylabels[0] in grid.dtype.names:
            self.results[mtype]['grid'] = pl.mlab.rec_append_fields(grid,\
                        ['c_'+ylabel for ylabel in ylabels],add_info)
        # -- if they are already in there, overwrite!
        else:
            for info,ylabel in zip(add_info,ylabels):
                self.results[mtype]['grid']['c_'+ylabel] = add_info

        # -- and replace the original confidence intervals, and re-order
        self.results[mtype]['grid']['ci_'+chi2_type] = 1-combined_pvals

        sa = np.argsort(self.results[mtype]['grid']['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = self.results[mtype]['grid'][sa]

        # -- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: {0:s} via stellar models and replaced ci_{1:s} with combined CI'.format(', '.join(ylabels),chi2_type))


    #}

    #{ Plotting routines

    def _label_dict(self, param):
        """
        Returns the label belonging to a certain parameter

        If the label is not present, the function will just return param.
        """
        #split parameter in param name and componentent number
        param, component = re.findall('(.*?)(\d?$)', param)[0]

        ldict = dict(teff='Effective temperature [K]',\
                    z='log (Metallicity Z [$Z_\odot$]) [dex]',\
                    logg=r'log (surface gravity [cm s$^{-2}$]) [dex]',\
                    ebv='E(B-V) [mag]',\
                    ci_raw='Raw probability [%]',\
                    ci_red='Reduced probability [%]',\
                    #labs=r'log (Absolute Luminosity [$L_\odot$]) [dex]',\
                    labs=r'Absolute Luminosity [$L_\odot$]',\
                    radius=r'Radius [$R_\odot$]',\
                    mass=r'Mass [$M_\odot$]',
                    mc=r'MC [Nr. points in hexagonal bin]',
                    rv=r'Extinction parameter $R_v$')
        if param in ldict:
            param = ldict[param]
        if component != '':
            return param #+ " - " + component
        else:
            return param

    @standalone_figure
    def plot_grid(self,x='teff',y='logg',ptype='ci_red',mtype='igrid_search',limit=0.95,d=None,**kwargs):
        """
        Plot grid as scatter plot

        PrameterC{ptype} sets the colors of the scattered points (e.g., 'ci_red','z','ebv').

        Example usage:

        First set the SED:

        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()

        Then make the plots:

        >>> p = pl.figure()
        >>> p = pl.subplot(221);mysed.plot_grid(limit=None)
        >>> p = pl.subplot(222);mysed.plot_grid(x='ebv',y='z',limit=None)
        >>> p = pl.subplot(223);mysed.plot_grid(x='teff',y='ebv',limit=None)
        >>> p = pl.subplot(224);mysed.plot_grid(x='logg',y='z',limit=None)

        ]]include figure]]ivs_sed_builder_plot_grid_01.png]

        """
        # -- if no distance is set, derive the most likely distance from the plx:
        #if d is None:
            #try:
                #d = self.get_distance_from_plx()
            #except ValueError:
                #d = None
                #logger.info('Distance to {0} unknown'.format(self.ID))
            #if isinstance(d,tuple):
                #d = d[0][np.argmax(d[1])]
        ## -- it's possible that we still don't have a distance
        #if d is not None:
            #logger.info('Assumed distance to {0} = {1:.3e} pc'.format(self.ID,d))
            #radius  = d*np.sqrt(self.results[mtype]['grid']['scale'])
            #radius  = conversions.convert('pc','Rsol',radius) # in Rsol
            #labs = np.log10(self.results[mtype]['grid']['labs']*radius**2) # in [Lsol]
            #mass = conversions.derive_mass((self.results[mtype]['grid']['logg'].copy(),'[cm/s2]'),\
                                           #(radius,'Rsol'),unit='Msol')
        # -- compute angular diameter
        theta = 2*conversions.convert('sr','mas',self.results[mtype]['grid']['scale'])

        if limit is not None:
            region = self.results[mtype]['grid']['ci_red']<=limit
        else:
            region = self.results[mtype]['grid']['ci_red']<np.inf
        # -- get the colors and the color scale
        if d is not None and ptype=='labs':
            colors = locals()[ptype][region]
        elif ptype in self.results[mtype]['grid'].dtype.names:
            colors = self.results[mtype]['grid'][ptype][region]
        elif isinstance(ptype,str):
            colors = locals()[ptype][region]
        else:
            colors = ptype[region]
            ptype = 'custom variable'

        if 'ci_' in ptype.lower():
            colors *= 100.
            vmin = colors.min()
            vmax = 95.
        else:
            vmin = kwargs.pop('vmin',colors.min())
            vmax = kwargs.pop('vmax',colors.max())

        # -- define abbrevation for plotting
        if x in self.results[mtype]['grid'].dtype.names:
            X = self.results[mtype]['grid'][x]
        else:
            X = locals()[x]

        if y in self.results[mtype]['grid'].dtype.names:
            Y = self.results[mtype]['grid'][y]
        else:
            Y = locals()[y]

        # -- make the plot
        if mtype == 'imc':
            pl.hexbin(X,Y,mincnt=1,cmap=pl.cm.spectral)  #bins='log'
            ptype = 'mc'

            # -- set the limits
            pl.xlim(X.max(),X.min())
            pl.ylim(Y.max(),Y.min())
            cbar = pl.colorbar()
        else:
            #if limit is not None:
                #region = self.results[mtype]['grid']['ci_red']<limit
            #else:
                #region = self.results[mtype]['grid']['ci_red']<np.inf
            ## -- get the colors and the color scale
            #if d is not None and ptype=='labs':
                #colors = locals()[ptype][region]
            #elif ptype in self.results[mtype]['grid'].dtype.names:
                #colors = self.results[mtype]['grid'][ptype][region]
            #else:
                #colors = locals()[ptype][region]

            #if 'ci' in ptype:
                #colors *= 100.
                #vmin = colors.min()
                #vmax = 95.
            #else:
                #vmin = kwargs.pop('vmin',colors.min())
                #vmax = kwargs.pop('vmax',colors.max())

            # -- grid scatter plot
            pl.scatter(X[region],Y[region],
                 c=colors,edgecolors='none',cmap=pl.cm.spectral,vmin=vmin,vmax=vmax)
            # -- set the limits to only include the 95 interval
            pl.xlim(X[region].max(),X[region].min())
            pl.ylim(Y[region].max(),Y[region].min())
            cbar = pl.colorbar()

        # -- mark best value
        pl.plot(X[-1],Y[-1],'r+',ms=40,mew=3)

        pl.xlabel(self._label_dict(x))
        pl.ylabel(self._label_dict(y))
        cbar.set_label(self._label_dict(ptype))

        logger.info('Plotted %s-%s diagram of %s'%(x,y,ptype))

    @standalone_figure
    def plot_CI2D(self,xpar='teff',ypar='logg',mtype='iminimize', ptype='ci_red',**kwargs):
        """
        Plot a 2D confidence intervall calculated using the CI2D computation
        from calculate_iminimize_CI2D. Make sure you first calculate the grid
        you want using the calculate_iminimize_CI2D function.
        """

        grid = self.results[mtype]['CI2D'][xpar+"-"+ypar][ptype]
        x = self.results[mtype]['CI2D'][xpar+"-"+ypar]['x']
        y = self.results[mtype]['CI2D'][xpar+"-"+ypar]['y']

        if ptype == 'ci_red' or ptype == 'ci_raw':
            grid = grid*100.0
            levels = np.linspace(0,100,25)
            ticks = [0,20,40,60,80,100]
        elif ptype == 'ci_chi2':
            grid = np.log10(grid)
            levels = np.linspace(np.min(grid), np.max(grid), 25)
            ticks = np.round(np.linspace(np.min(grid), np.max(grid), 11), 2)

        pl.contourf(x,y,grid,levels,**kwargs)
        pl.xlabel(self._label_dict(xpar))
        pl.ylabel(self._label_dict(ypar))
        cbar = pl.colorbar(fraction=0.08,ticks=ticks)
        cbar.set_label(ptype!='ci_chi2' and r'Probability' or r'$^{10}$log($\chi^2$)')





    @standalone_figure
    def plot_data(self,colors=False, plot_unselected=True,
                  unit_wavelength='angstrom',unit_flux=None,**kwargs):
        """
        Plot only the SED data.

        Extra kwargs are passed to plotting functions.

        The decorator provides an keyword C{savefig}. When set to C{True}, an
        image name is generated automatically (very long!), and the figure is
        closed. When set to a string, the image is saved with that name, and
        the figure is closed. When C{savefig} is not given, the image will
        stay open so that the user can still access all plot elements and do
        further enhancements.

        :param colors: if False, plot absolute values, otherwise plot colors
        (flux ratios)
        :type colors: boolean
        :param plot_unselected: if True, all photometry is plotted, otherwise
        only those that are selected
        """
        if not plot_unselected:
            master = self.master[self.master['include']]
        else:
            master = self.master
        if unit_flux is None:
            unit_flux = master['cunit'][0]
        wave,flux,e_flux = master['cwave'],master['cmeas'],master['e_cmeas']
        sources = master['source']
        iscolor = np.array(master['color'],bool)
        photbands = master['photband']
        indices = np.arange(len(master))

        allsystems = np.array([i.split('.')[0] for i in photbands])
        systems = sorted(set(allsystems))
        color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(systems))]

        if not colors:
            color_cycle = itertools.cycle(color_cycle)
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            wave = conversions.convert('angstrom',unit_wavelength,wave)
            flux,e_flux = conversions.convert(master['cunit'][0],unit_flux,flux,e_flux,wave=(wave,unit_wavelength))
            mf = []
            # plot photometric systems in different colors
            for system in systems:
                keep = (allsystems==system) & ~iscolor
                if keep.sum():
                    # plot each photometric points separately, so that we could
                    # use it interactively. Label them all with a unique ID
                    # and make them pickable.
                    color = next(color_cycle)
                    #for i in range(sum(keep)):
                        #label = system if i==0 else '_nolegend_'
                        #pltlin,caplins,barlincs = pl.errorbar(wave[keep][i],flux[keep][i],yerr=e_flux[keep][i],fmt='o',label=label,ms=7,picker=5,color=color,**kwargs)
                        #pltlin.sed_index = indices[keep][i]
                        #caplins[0].sed_index = indices[keep][i]
                        #caplins[1].sed_index = indices[keep][i]
                        #barlincs[0].sed_index = indices[keep][i]
                    pl.errorbar(wave[keep],flux[keep],yerr=e_flux[keep],fmt='o',label=system,ms=7,color=color,**kwargs)
                    mf.append(flux[keep])
            if keep.sum():
                pl.ylabel(conversions.unit2texlabel(unit_flux,full=True))
            pl.xlabel('Wavelength [{0}]'.format(conversions.unit2texlabel(unit_wavelength)))
            # -- scale y-axis (sometimes necessary for data with huge errorbars)
            mf = np.log10(np.hstack(mf))
            lmin,lmax = np.nanmin(mf),np.nanmax(mf)
            lrange = np.abs(lmin-lmax)
            pl.ylim(10**(lmin-0.1*lrange),10**(lmax+0.1*lrange))
        else:
            pl.gca().set_color_cycle(color_cycle)
            names = []
            start_index = 1
            for system in systems:
                keep = (allsystems==system) & iscolor
                if keep.sum():
                    pl.errorbar(list(range(start_index,start_index+keep.sum())),flux[keep],yerr=e_flux[keep],fmt='o',label=system,ms=7,**kwargs)
                    names += [ph.split('.')[1] for ph in photbands[keep]]
                start_index += keep.sum()
            pl.xticks(list(range(1,len(names)+1)),names,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Index')
        leg = pl.legend(prop=dict(size='small'),loc='best',fancybox=True) #,numpoints=1) #prop=dict(size='small'),loc='best',fancybox=True)
        leg.get_frame().set_alpha(0.5)
        pl.grid()
        pl.title(self.ID)


    @standalone_figure
    def plot_sed(self,colors=False,mtype='igrid_search',plot_redded=True,plot_deredded=False,
            plot_unselected=True,wave_units='AA',flux_units='erg/s/cm2',**kwargs):
        """
        Plot a fitted SED together with the data.

        Example usage:

        First set the SED:

        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()

        Then make the plots:

        >>> p = pl.figure()
        >>> p = pl.subplot(121)
        >>> mysed.plot_sed(colors=False)
        >>> p = pl.subplot(122)
        >>> mysed.plot_sed(colors=True)

        ]]include figure]]ivs_sed_builder_plot_sed_01.png]
        """
        annotation = kwargs.pop('annotation',True)

        def plot_sed_getcolors(master,color_dict=None):
            myphotbands = [iphotb.split('.')[1] for iphotb in master['photband'][master['color']]]
            if not myphotbands:  # -- If there are no colours none can be returned (added by joris 30-01-2012)
                return [],[],[],None
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

        # -- get the color cycle
        systems = np.array([system.split('.')[0] for system in self.master['photband']],str)
        set_systems = sorted(list(set(systems)))
        if ('absolutesymbolcolor' in kwargs) and (kwargs.pop('absolutesymbolcolor') == True):
            sortedphotsystems,plotcolorvalues = filters.get_plotsymbolcolorinfo()
            selectedcolorinds = np.array([np.where(sortedphotsystems == system)[0][0] for system in set_systems])
            color_cycle = itertools.cycle([pl.cm.spectral(j) for j in plotcolorvalues[selectedcolorinds]])
        else:
            color_cycle = itertools.cycle([pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))])

        # -- for plotting reasons, we translate every color to an integer
        for system in set_systems:
            color = next(color_cycle)
            keep = (systems==system) & (self.master['color']==colors)
            if not plot_unselected:
                keep = keep & self.master['include']
            # -- synthetic:
            if sum(keep) and mtype in self.results and 'synflux' in self.results[mtype]:
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                    y = self.results[mtype]['synflux'][1][keep]
                else:
                    x = self.results[mtype]['synflux'][0][keep]
                    y = self.results[mtype]['synflux'][1][keep]
                    # -- convert to correct units
                    y = conversions.convert('erg/s/cm2/AA',flux_units,y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)
                pl.plot(x,y,'x',ms=10,mew=2,alpha=0.75,color=color,**kwargs)
            # -- include:
            keep = (systems==system) & (self.master['color']==colors) & self.master['include']
            if sum(keep):
                if colors:
                    # -- translate every color to an integer
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    if mtype in self.results and 'synflux' in self.results[mtype]:
                        x = self.results[mtype]['synflux'][0][keep]
                        for ind in range(len(x)):
                           if np.isnan(x[ind]):
                               try:
                                   x[ind] = self.master['cwave'][keep][ind]
                               except:
                                   pass
                    else:
                        x = self.master['cwave'][keep]
                    y = self.master['cmeas'][keep]
                    e_y = self.master['e_cmeas'][keep]
                    y,e_y = conversions.convert('erg/s/cm2/AA',flux_units,y,e_y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)

                p = pl.errorbar(x,y,yerr=e_y,fmt='o',label=system,
                                capsize=10,ms=7,color=color,**kwargs)

            # -- exclude:
            label = np.any(keep) and '_nolegend_' or system
            keep = (systems==system) & (self.master['color']==colors) & ~self.master['include']
            if sum(keep) and plot_unselected:
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results[mtype]['synflux'][0][keep]
                    if np.any(np.isnan(x)):
                        x = self.master['cwave'][keep]
                    y = self.master['cmeas'][keep]
                    e_y = self.master['e_cmeas'][keep]
                    y,e_y = conversions.convert('erg/s/cm2/AA',flux_units,y,e_y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)
                pl.errorbar(x,y,yerr=e_y,fmt='o',label=label,
                            capsize=10,ms=7,mew=2,color=color,mfc='1',mec=color,**kwargs)

        # -- only set logarithmic scale if absolute fluxes are plotted
        #   and only plot the real model then
        if not colors:
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            pl.gca().set_autoscale_on(False)

            # -- the model
            if mtype in self.results and 'model' in self.results[mtype]:
                wave,flux,flux_ur = self.results[mtype]['model']

                flux = conversions.convert('erg/s/cm2/AA',flux_units,flux,wave=(wave,'AA'))
                flux_ur = conversions.convert('erg/s/cm2/AA',flux_units,flux_ur,wave=(wave,'AA'))
                wave = conversions.convert('AA',wave_units,wave)

                if plot_redded:
                    pl.plot(wave,flux,'r-',**kwargs)
                if plot_deredded:
                    pl.plot(wave,flux_ur,'k-',**kwargs)
            pl.ylabel(conversions.unit2texlabel(flux_units,full=True))
            pl.xlabel('wavelength [{0}]'.format(conversions.unit2texlabel(wave_units)))
        else:
            xlabels = list(color_dict.keys())
            xticks = [color_dict[key] for key in xlabels]
            pl.xticks(xticks,xlabels,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Color name')
            pl.xlim(min(xticks)-0.5,max(xticks)+0.5)

        pl.grid()
        if colors:
            leg = pl.legend(loc='best',prop=dict(size='x-small'),numpoints=1)
        else:
            leg = pl.legend(loc='upper right',prop=dict(size='x-small'),numpoints=1)
            leg.get_frame().set_alpha(0.5)
        loc = (0.35,0.05)
        if mtype in self.results and 'grid' in self.results[mtype] and annotation:
            teff = self.results[mtype]['grid']['teff'][-1]
            logg = self.results[mtype]['grid']['logg'][-1]
            ebv = self.results[mtype]['grid']['ebv'][-1]
            met = self.results[mtype]['grid']['z'][-1]
            scale = self.results[mtype]['grid']['scale'][-1]
            angdiam = 2*conversions.convert('sr','mas',scale)
            try:
                teff2 = self.results[mtype]['CI']['teff2']
                logg2 = self.results[mtype]['CI']['logg2']
                radii = self.results[mtype]['CI']['rad2']/self.results[mtype]['CI']['rad']
                pl.annotate('Teff=%i   %i K\nlogg=%.2f   %.2f cgs\nE(B-V)=%.3f mag\nr2/r1=%.2f\n$\Theta$=%.3g mas'%(teff,teff2,logg,logg2,ebv,radii,angdiam),
                        loc,xycoords='axes fraction')
            except:
                pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\nlogZ=%.3f dex\n$\Theta$=%.3g mas'%(teff,logg,ebv,met,angdiam),
                        loc,xycoords='axes fraction')
                pass
        logger.info('Plotted SED as %s'%(colors and 'colors' or 'absolute fluxes'))
        teff = "%d" % teff
        logg = "%.2f" % logg
        ebv = "%.3f" % ebv
        metallicity = "%.3f" % met
        #wave,flux = model.get_table(teff=teff,logg=logg,ebv=ebv,z=z,grid='kurucz2')
        '''f = open('/home/anae/python/SEDfitting/Bastars_info.txt', 'a')
        f.writelines(str(teff)+'\t'+str(logg)+'\t'+str(ebv)+'\t'+str(metallicity)+'\n')
        f.closed'''



    @standalone_figure
    def plot_chi2(self,colors=False,mtype='igrid_search',**kwargs):
        """
        Plot chi2 statistic for every datapoint included in the fit.

        To plot the statistic from the included absolute values, set
        C{colors=False}.

        To plot the statistic from the included colors, set C{colors=True}.

        Example usage:

        First set the SED:

        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()

        Then make the plots:

        >>> p = pl.figure()
        >>> p = pl.subplot(121)
        >>> mysed.plot_chi2(colors=False)
        >>> p = pl.subplot(122)
        >>> mysed.plot_chi2(colors=True)

        ]]include figure]]ivs_sed_builder_plot_chi2_01.png]

        :param colors: flag to distinguish between colors and absolute values
        :type colors: boolean
        """
        if 'phase' in kwargs:
            uniquephase = kwargs.pop('phase')
            phase = True
        else:
            uniquephase = None
            phase = False

        include_grid = self.master['include']
        systems = np.array([system.split('.')[0] for system in self.master['photband'][include_grid]],str)
        set_systems = sorted(list(set(systems)))
        color_cycle = itertools.cycle([pl.cm.spectral(i) for i in np.linspace(0, 1.0, len(set_systems))])
        if mtype in self.results:
            eff_waves,synflux,photbands = self.results[mtype]['synflux']
            chi2 = self.results[mtype]['chi2']
            for system in set_systems:
                color = next(color_cycle)
                keep = (systems==system)
                if phase:
                    keep = keep & (self.master['phase'][include_grid] == uniquephase)
                if sum(keep) and not colors:
                    try:
                        pl.loglog(eff_waves[include_grid][keep],chi2[include_grid][keep],'o',label=system,color=color)
                    except:
                        logger.critical('Plotting of CHI2 of absolute values failed')
                elif sum(keep) and colors:
                    pl.semilogy(list(range(len(eff_waves[include_grid][keep]))),chi2[include_grid][keep],'o',label=system,color=color)
            pl.legend(loc='upper right',prop=dict(size='x-small'))
            pl.grid()
            pl.annotate('Total $\chi^2$ = %.1f'%(self.results[mtype]['grid']['chisq'][-1]),(0.59,0.120),xycoords='axes fraction',color='r')
            pl.annotate('Total Reduced $\chi^2$ = %0.2f'%(sum(chi2[include_grid][keep])),(0.59,0.075),xycoords='axes fraction',color='r')
            if 'factor' in self.results[mtype]:
                pl.annotate('Error scale = %.2f'%(np.sqrt(self.results[mtype]['factor'])),(0.59,0.030),xycoords='axes fraction',color='k')
            xlims = pl.xlim()
            pl.plot(xlims,[self.results[mtype]['grid']['chisq'][-1],self.results[mtype]['grid']['chisq'][-1]],'r-',lw=2)
            pl.xlim(xlims)
            pl.xlabel('wavelength [$\AA$]')
            pl.ylabel(r'Reduced ($\chi^2$)')
            logger.info('Plotted CHI2 of %s'%(colors and 'colors' or 'absolute fluxes'))
        else:
            logger.info('%s not in results, no plot could be made.'%(mtype))


    @standalone_figure
    def plot_distance(self,mtype='igrid_search'):
        try:
            # -- necessary information
            (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
            (d,dprob) = self.results['igrid_search']['d']

            ax_d = pl.gca()

            gal = self.info['galpos']
            # -- the plot
            dzoom = dprob>1e-4
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
        except KeyError:
            logger.info('No distance/reddening plotted due to KeyError.')

    @standalone_figure
    def plot_grid_model(self,ptype='prob'):
        """
        Grid of models
        """
        if 'spType' in self.info:
            pl.title(self.info['spType'])
        cutlogg = (self.results['igrid_search']['grid']['logg']<=4.4) & (self.results['igrid_search']['grid']['ci_red']<=0.95)
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        gal = self.info['galpos']

        n = 75000
        region = self.results['igrid_search']['grid']['ci_red']<0.95
        total_prob = 100-(1-self.results['igrid_search']['grid']['ci_red'][cutlogg][-n:])*d_prob_models*100
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


    @standalone_figure
    def plot_MW_side(self):
        im = Image.open(config.get_datafile('images','NRmilkyway.tif'))
        left,bottom,width,height = 0.0,0.0,1.0,1.0
        startm,endm = 183,-177
        startv,endv = -89,91

        xwidth = startm-endm
        ywidth = 90.
        ratio = ywidth/xwidth

        gal = list(self.info['galpos'])
        if gal[0]>180:
            gal[0] = gal[0] - 360.
        # -- boundaries of ESO image
        pl.imshow(im,extent=[startm,endm,startv,endv],origin='lower')
        pl.plot(gal[0],gal[1],'rx',ms=15,mew=2)
        pl.xlim(startm,endm)
        pl.ylim(startv,endv)
        pl.xticks([])
        pl.yticks([])

    @standalone_figure
    def plot_MW_top(self):
        im = Image.open(config.get_datafile('images','topmilkyway.jpg'))
        pl.imshow(im,origin='lower')
        pl.box(on=False)
        pl.xticks([])
        pl.yticks([])
        xlims = pl.xlim()
        ylims = pl.ylim()
        gal = self.info['galpos']
        pl.plot(2800,1720,'ro',ms=10)
        pl.plot([2800,-5000*np.sin(gal[0]/180.*np.pi)+2800],[1720,5000*np.cos(gal[0]/180.*np.pi)+1720],'r-',lw=2)

        # -- necessary information
        if 'igrid_search' in self.results and 'd_mod' in self.results['igrid_search']:
            (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
            (d,dprob) = self.results['igrid_search']['d']
        else:
            d = np.linspace(0,1000,2)
            dprob = np.zeros(len(d))

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

    @standalone_figure
    def plot_finderchart(self,cmap_photometry=pl.cm.spectral,window_size=5.):
        """
        Size is x and y width in arcminutes
        """
        try:
            dec = 'jdedeg' in self.info and self.info['jdedeg'] or None
            ra = 'jradeg' in self.info and self.info['jradeg'] or None
            data,coords,size = mast.get_dss_image(self.ID,ra=ra,dec=dec)
            pl.imshow(data[::-1],extent=[-size[0]/2*60,size[0]/2*60,
                                        -size[1]/2*60,size[1]/2*60],cmap=pl.cm.RdGy_r)#Greys
            pl.xlim(-window_size/2.,+window_size/2.)
            pl.ylim(-window_size/2.,+window_size/2.)
            xlims,ylims = pl.xlim(),pl.ylim()
            keep_this = -self.master['color'] & (self.master['cmeas']>0)
            toplot = self.master[keep_this]
            systems = np.array([system.split('.')[0] for system in toplot['photband']],str)
            set_systems = sorted(list(set(systems)))
            color_cycle = itertools.cycle([cmap_photometry(j) for j in np.linspace(0, 1.0, len(set_systems))])
            for system in set_systems:
                color = next(color_cycle)
                keep = systems==system
                if sum(keep):
                    pl.plot(toplot['_RAJ2000'][keep][0]*60,
                            toplot['_DEJ2000'][keep][0]*60,'x',label=system,
                            mew=2.5,ms=15,alpha=0.5,color=color)
            leg = pl.legend(numpoints=1,prop=dict(size='x-small'),loc='best',fancybox=True)
            leg.get_frame().set_alpha(0.75)
            pl.xlim(xlims)
            pl.ylim(ylims)
            pl.xlabel(r'Right ascension $\alpha$ [arcmin]')
            pl.ylabel(r'Declination $\delta$ [arcmin]')
        except:
            logger.warning('No image found of %s'%(self.ID))
            pass

        if 'pm' in self.info:
            logger.info("Found proper motion info")
            ppm_ra,ppm_de = (self.info['pm']['pmRA'],self.info['pm']['epmRA']),(self.info['pm']['pmDE'],self.info['pm']['epmDE'])
            pl.annotate('',xy=(ppm_ra[0]/50.,ppm_de[0]/50.),
                  xycoords='data',xytext=(0,0),textcoords='data',color='red',
                  arrowprops=dict(facecolor='red', shrink=0.05),
                  horizontalalignment='right', verticalalignment='top')
            pl.annotate('pmRA: %.1f $\pm$ %.1f mas/yr\npmDE: %.1f $\pm$ %.1f mas/yr'%(ppm_ra[0],ppm_ra[1],ppm_de[0],ppm_de[1]),
                        xy=(0.05,0.25),xycoords='axes fraction',color='red')
            if 'igrid_search' in self.results and 'd' in self.results['igrid_search']:
                (d,dprob) = self.results['igrid_search']['d']
                max_distance = d[np.argmax(dprob)]
                e_max_distance = abs(max_distance - d[np.argmin(np.abs(dprob-0.5*max(dprob)))])
            elif 'plx' in self.info and 'v' in self.info['plx'] and 'e' in self.info['plx']:
                plx,eplx = self.info['plx']['v'],self.info['plx']['e']
                dist = 1000./ufloat(plx,eplx)
                max_distance,e_max_distance = dist.nominal_value,dist.std_dev
            else:
                max_distance = 1000.
                e_max_distance = 100.

            tang_velo = 'Tan. vel. at %.0f+/-%.0f pc: '%(max_distance,e_max_distance)

            max_distance = conversions.convert('pc','km',max_distance,e_max_distance)
            ppm_ra = conversions.convert('mas/yr','rad/s',*ppm_ra)
            ppm_de = conversions.convert('mas/yr','rad/s',*ppm_de)
            max_distance = unumpy.uarray(*[max_distance[0],max_distance[1]])
            x = unumpy.uarray(*[ppm_ra[0],ppm_ra[1]])
            y = unumpy.uarray(*[ppm_de[0],ppm_de[1]])
            velocity = max_distance*utan( usqrt(x**2+y**2))


            pl.annotate(tang_velo + '%s km/s'%(velocity),xy=(0.05,0.2),xycoords='axes fraction',color='red')
            if 'Vel' in self.info and 'v' in self.info['Vel']:
                rad_velo = 'Rad. vel.: %.1f'%(self.info['Vel']['v'])
                if 'e' in self.info['Vel']:
                    rad_velo += '+/-%.1f'%(self.info['Vel']['e'])
                pl.annotate(rad_velo+' km/s',xy=(0.05,0.15),xycoords='axes fraction',color='red')


    def make_plots(self):
        """
        Make all available plots
        """
        pl.figure(figsize=(22,12))
        rows,cols = 3,4
        pl.subplots_adjust(left=0.04, bottom=0.07, right=0.97, top=0.96,
                wspace=0.17, hspace=0.24)
        pl.subplot(rows,cols,1);self.plot_grid(ptype='ci_red')
        pl.subplot(rows,cols,2);self.plot_grid(ptype='ebv')
        pl.subplot(rows,cols,3);self.plot_grid(ptype='z')
        pl.subplot(rows,cols,4);self.plot_distance()

        pl.subplot(3,2,3);self.plot_sed(colors=False)
        pl.subplot(3,2,5);self.plot_sed(colors=True)

        pl.subplot(rows,cols,7);self.plot_chi2(colors=False)
        pl.subplot(rows,cols,11);self.plot_chi2(colors=True)

        #pl.subplot(rows,cols,8);self.plot_grid_model(ptype='prob')
        #pl.subplot(rows,cols,12);self.plot_grid_model(ptype='radii')

        pl.figure(figsize=(12,12))
        pl.axes([0,0.0,1.0,0.5]);self.plot_MW_side()
        pl.axes([0,0.5,0.5,0.5]);self.plot_MW_top()
        pl.axes([0.5,0.5,0.5,0.5]);self.plot_finderchart()

    #}

    #{Input and output

    def save_photometry(self,photfile=None):
        """
        Save master photometry to a file.

        :param photfile: name of the photfile. Defaults to C{starname.phot}.
        :type photfile: str
        """
        # -- write to file
        if photfile is not None:
            self.photfile = photfile
        logger.info('Save photometry to file %s'%(self.photfile))
        # -- add some comments
        if self.ID:
            if not 'bibcode' in self.master.dtype.names:
                self.master = crossmatch.add_bibcodes(self.master)
            if not 'comments' in self.master.dtype.names:
                self.master = vizier.quality_check(self.master,self.ID)
        ascii.write_array(self.master,self.photfile,header=True,auto_width=True,use_float='%g',comments=['#'+json.dumps(self.info)])

    def load_photometry(self,photfile=None):
        """
        Load the contents of the photometry file to the master record.

        :param photfile: name of the photfile. Defaults to the value of C{self.photfile}.
        :type photfile: str
        """
        if photfile is not None:
            self.photfile = photfile
        logger.info('Load photometry from file %s'%(self.photfile))
        self.master,comments = ascii.read2recarray(self.photfile,return_comments=True)
        #self.info = json.loads(comments[-3])                                                                                           #### HE COMENTADO ESTO PARA QUE NO SE ATASQUE

        ## to make the builder backwards-compatible with files made with an older version that did not implement a 'phases' column yet
        if 'phase' not in self.master.dtype.names:
            extra_cols = [[0]*len(self.master['meas'])]
            extradtype = [('phase',np.int)]
            names = list(self.master.dtype.names)
            lastnames = []
            if 'bibcode' in names:
                lastnames.append('bibcode')
                names.remove('bibcode')
            if 'comments' in names:
                lastnames.append('comments')
                names.remove('comments')

            lastrecords = self.master[lastnames]
            self.master = numpy_ext.recarr_addcols(self.master[names],extra_cols,extradtype)
            self.master = numpy_ext.recarr_join(self.master,lastrecords)


    def save_fits(self,filename=None,overwrite=True):
        """
        Save content of SED object to a FITS file.

        The .fits file will contain the following extensions if they are present in the object:
           1) data (all previously collected photometric data = content of .phot file)
           2) model_igrid_search (full best model SED table from grid_search)
           3) igrid_search (CI from grid search = actual grid samples)
           4) synflux_igrid_search (integrated synthetic fluxes for best model from grid_search)
           5) model_imc (full best model SED table from monte carlo)
           6) imc (CI from monte carlo = actual monte carlo samples)
           7) synflux_imc (integrated synthetic fluxes for best model from monte carlo)

        :param filename: name of SED FITS file
        :type filename: string
        :param overwrite: overwrite old FITS file if true
        :type overwrite: boolean

        Example usage:

        >>> #mysed.save_fits()
        >>> #mysed.save_fits(filename='myname.fits')
        """
        if filename is None:
            filename = str(os.path.splitext(self.photfile)[0]+'.fits')
        if overwrite:
            if os.path.isfile(filename):
                os.remove(filename)
                logger.info('Old FITS file removed')

        # -- write primary header
        #prim_header = {}
        #for key in self.info:
            #if not (isinstance(self.info[key],float) or isinstance(self.info[key],str)):
                #continue
            #prim_header[key] = self.info[key]
        #fits.write_recarray(np.array([[0]]),filename,header_dict=prim_header,ext=0)

        # -- write master data
        master = self.master.copy()
        fits.write_recarray(master,filename,header_dict=dict(extname='data'))

        # -- write the rest
        for mtype in self.results:#['igrid_search','imc']:
            eff_waves,synflux,photbands = self.results[mtype]['synflux']
            chi2 = self.results[mtype]['chi2']

            results_modeldict = dict(extname='model_'+mtype)
            results_griddict = dict(extname=mtype)
            keys = sorted(self.results[mtype])
            for key in keys:
                if 'CI' in key:
                    for ikey in self.results[mtype][key]:
                        if '_l' not in ikey and '_u' not in ikey and ikey != 'chisq':
                            results_modeldict[ikey] = self.results[mtype][key][ikey]
                        results_griddict[ikey] = self.results[mtype][key][ikey]
                if key=='factor':
                    results_griddict[key] = self.results[mtype][key]

            fits.write_array(list(self.results[mtype]['model']),filename,
                            names=('wave','flux','dered_flux'),
                            units=('AA','erg/s/cm2/AA','erg/s/cm2/AA'),
                            header_dict=results_modeldict)

            index = 1
            while index > 0:
                headerdict = results_modeldict.copy()
                key = 'model{}'.format(index)
                headerdict['extname'] = key+'_'+mtype
                if key in list(self.results[mtype].keys()):
                    fits.write_array(list(self.results[mtype][key]),filename,
                             names=('wave','flux','dered_flux'),
                             units=('AA','erg/s/cm2/AA','erg/s/cm2/AA'),
                             header_dict=headerdict)
                    index += 1
                else:
                    index = 0

            if 'grid' in self.results[mtype]:
                fits.write_recarray(self.results[mtype]['grid'],filename,header_dict=results_griddict)

            results = np.rec.fromarrays([synflux,eff_waves,chi2],dtype=[('synflux','f8'),('mod_eff_wave','f8'),('chi2','f8')])

            fits.write_recarray(results,filename,header_dict=dict(extname='synflux_'+mtype))

        logger.info('Results saved to FITS file: %s'%(filename))


    #def load_fits2(self,filename=None):
        #"""
        #Load a previously made SED FITS file. Only works for SEDs saved with
        #the save_fits function after 14.06.2012.

        #The .fits file can contain the following extensions:
           #1) data (all previously collected photometric data = content of .phot file)
           #2) model_igrid_search (full best model SED table from grid_search)
           #3) igrid_search (CI from grid search = actual grid samples)
           #4) synflux_igrid_search (integrated synthetic fluxes for best model from grid_search)
           #5) model_imc (full best model SED table from monte carlo)
           #6) imc (CI from monte carlo = actual monte carlo samples)
           #7) synflux_imc (integrated synthetic fluxes for best model from monte carlo)

        #:param filename: name of SED FITS file
        #:type filename: string
        #:rtype: bool
        #:return: true if Fits file could be loaded
        #"""
        #if filename is None:
            #filename = os.path.splitext(self.photfile)[0]+'.fits'
        #if not os.path.isfile(filename):
            #logger.warning('No previous results saved to FITS file {:s}'.format(filename))
            #return False
        #ff = pf.open(filename)

        ## -- observed photometry
        #fields = ff['data'].columns.names
        #master = np.rec.fromarrays([ff['data'].data.field(field) for field in fields],names=','.join(fields))
        ## -- remove the whitespace in columns with strings added by fromarrays
        #for i,name in enumerate(master.dtype.names):
            #if master.dtype[i].str.count('S'):
                #for j,element in enumerate(master[name]):
                    #master[name][j] = element.strip()
        #self.master = master

        ## -- add dictionary that will contain the results
        #if not hasattr(self,'results'):
            #self.results = {}

        ## -- grid search and MC results
        #mtypes = [ext.header['extname'] for ext in ff[1:]]
        #mtypes = list(set(mtypes) - set(['data']))

        #for mtype in mtypes:
            #if mtype.startswith('i'):
                #continue
            #else:
                #mtype = mtype.lower().split('_',1) #lstrip('synflux_').lstrip('model_')
                #prefix,mtype = splitted[0],splitted[1]
            #if not mtype in self.results:
                    #self.results[mtype] = {}
            #self.results[mtype][prefix] = np.array(ff[prefix+'_'+mtype].data.field('wave'),dtype='float64'),np.array(ff[prefix+'_'+mtype].data.field('flux'),dtype='float64'),np.array(ff[prefix+'_'+mtype].data.field('dered_flux'),dtype='float64')
            #self.results[mtype]['chi2'] = np.array(ff['synflux_'+mtype].data.field('chi2'),dtype='float64')
            #self.results[mtype]['synflux'] = np.array(ff['synflux_'+mtype].data.field('mod_eff_wave'),dtype='float64'),np.array(ff['synflux_'+mtype].data.field('synflux'),dtype='float64'),self.master['photband']

        #searchmtypes = []
        #for mtype in mtypes:
            #if 'igrid_search' in mtype:
                #searchmtypes.append(mtype)
            #elif 'iminimize' in mtype:
                #searchmtypes.append(mtype)
            #elif 'imc' in mtype:
                #searchmtypes.append(mtype)

        #for mtype in searchmtypes: #['igrid_search','iminimize','imc']:
            #try:
                #fields = ff[mtype].columns.names
                #master = np.rec.fromarrays([ff[mtype].data.field(field) for field in fields],names=','.join(fields))
                #if not mtype in self.results:
                    #self.results[mtype] = {}
                #self.results[mtype]['grid'] = master
                #if 'factor' in ff[mtype].header:
                    #self.results[mtype]['factor'] = np.array([ff[mtype].header['factor']])[0]

                #headerkeys = ff[mtype].header.ascardlist().keys()
                #for key in headerkeys[::-1]:
                    #for badkey in ['xtension','bitpix','naxis','pcount','gcount','tfields','ttype','tform','tunit','factor','extname']:
                        #if key.lower().count(badkey):
                            #headerkeys.remove(key)
                            #continue
                #self.results[mtype]['CI'] = {}
                #for key in headerkeys:
                    ## -- we want to have the same types as the original: numpy.float64 --> np.array([..])[0]
                    #self.results[mtype]['CI'][key.lower()] = np.array([ff[mtype].header[key]])[0]
            #except KeyError:
                #continue

        ##self.results['igrid_search'] = {}
        ##fields = ff['igrid_search'].columns.names
        ##master = np.rec.fromarrays([ff['igrid_search'].data.field(field) for field in fields],names=','.join(fields))
        ##self.results['igrid_search']['grid'] = master
        ##self.results['igrid_search']['factor'] = ff['igrid_search'].header['factor']

        ##self.results['model'] = ff[2].data.field('wave'),ff[2].data.field('flux'),ff[2].data.field('dered_flux')
        ##self.results['chi2'] = ff[4].data.field('chi2')
        ##self.results['synflux'] = ff[4].data.field('mod_eff_wave'),ff[4].data.field('synflux'),ff[1].data.field('photband')

        #ff.close()

        #logger.info('Loaded previous results from FITS file: %s'%(filename))
        #return filename

    def load_fits(self,filename=None):
        """
        Load a previously made SED FITS file. Only works for SEDs saved with
        the save_fits function after 14.06.2012.

        The .fits file can contain the following extensions:
           1) data (all previously collected photometric data = content of .phot file)
           2) model_igrid_search (full best model SED table from grid_search)
           3) igrid_search (CI from grid search = actual grid samples)
           4) synflux_igrid_search (integrated synthetic fluxes for best model from grid_search)
           5) model_imc (full best model SED table from monte carlo)
           6) imc (CI from monte carlo = actual monte carlo samples)
           7) synflux_imc (integrated synthetic fluxes for best model from monte carlo)

        :param filename: name of SED FITS file
        :type filename: string
        :rtype: bool
        :return: true if Fits file could be loaded
        """
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.fits'
        if not os.path.isfile(filename):
            logger.warning('No previous results saved to FITS file {:s}'.format(filename))
            return False
        ff = pf.open(filename)

        # -- observed photometry
        fields = ff['data'].columns.names
        master = np.rec.fromarrays([ff['data'].data.field(field) for field in fields],names=','.join(fields))

        # -- remove the whitespace in columns with strings added by fromarrays
        for i,name in enumerate(master.dtype.names):
            if master.dtype[i].str.count('S'):
                for j,element in enumerate(master[name]):
                    master[name][j] = element.strip()
        self.master = master

        # -- add dictionary that will contain the results
        if not hasattr(self,'results'):
            self.results = {}

        # -- grid search and MC results
        mtypes = [ext.header['extname'] for ext in ff[1:]]
        mtypes = list(set(mtypes) - set(['data','DATA']))

        for mtype in mtypes:
            mtype = mtype.lower()
            if mtype.startswith('i'):
                if not mtype in self.results:
                        self.results[mtype] = {}
                try:
                    fields = ff[mtype].columns.names
                    master = np.rec.fromarrays([ff[mtype].data.field(field) for field in fields],names=','.join(fields))
                    if not mtype in self.results:
                        self.results[mtype] = {}
                    self.results[mtype]['grid'] = master
                    if 'factor' in ff[mtype].header:
                        self.results[mtype]['factor'] = np.array([ff[mtype].header['factor']])[0]

                    headerkeys = list(ff[mtype].header.keys()) #ascardlist().keys()
                    for key in headerkeys[::-1]:
                        for badkey in ['xtension','bitpix','naxis','pcount','gcount','tfields','ttype','tform','tunit','factor','extname']:
                            if key.lower().count(badkey):
                                headerkeys.remove(key)
                                continue
                    self.results[mtype]['CI'] = {}
                    for key in headerkeys:
                        # -- we want to have the same types as the original: numpy.float64 --> np.array([..])[0]
                        self.results[mtype]['CI'][key.lower()] = np.array([ff[mtype].header[key]])[0]
                except KeyError as msg:
                    print(msg)
                    continue
            else:
                splitted = mtype.lower().split('_',1) #lstrip('synflux_').lstrip('model_')
                prefix,mtype = splitted[0],splitted[1]
                if not mtype in self.results:
                    self.results[mtype] = {}
                if 'model' in prefix:
                    self.results[mtype][prefix] = np.array(ff[prefix+'_'+mtype].data.field('wave'),dtype='float64'),np.array(ff[prefix+'_'+mtype].data.field('flux'),dtype='float64'),np.array(ff[prefix+'_'+mtype].data.field('dered_flux'),dtype='float64')
                elif 'synflux' in prefix:
                    self.results[mtype]['chi2'] = np.array(ff[prefix+'_'+mtype].data.field('chi2'),dtype='float64')
                    self.results[mtype][prefix] = np.array(ff[prefix+'_'+mtype].data.field('mod_eff_wave'),dtype='float64'),np.array(ff[prefix+'_'+mtype].data.field('synflux'),dtype='float64'),self.master['photband']

        ff.close()

        logger.info('Loaded previous results from FITS file: %s'%(filename))
        return filename

    def save_hdf5(self, filename=None, update=True):
        """
        Save content of SED object to a HDF5 file. (HDF5 is the successor of FITS files,
        providing a clearer structure of the saved content.)
        This way of saving is more thorough that save_fits(), fx. the CI2D confidence
        intervals are not save to a fits file, but are saved to a hdf5 file.
        Currently the following data is saved to HDF5 file:

        - sed.master (photometry)
        - sed.results (results from all fitting methods)
        - sed.constraints (extra constraints on the fits)

        :param filename: name of SED FITS file
        :type filename: string
        :param update: if True, an existing file will be updated with the current information, if
        False, an existing fill be overwritten
        :type update: bool
        :return: the name of the output HDF5 file.
        :rtype: string
        """

        if filename is None:
            filename = str(os.path.splitext(self.photfile)[0]+'.hdf5')

        data = dict()
        data['master'] = self.master
        data['results'] = self.results
        data['constraints'] = self.constraints

        hdf5.write_dict(data, filename, update=update)

        logger.info('Results saved to HDF5 file: %s'%(filename))
        return filename

    def load_hdf5(self,filename=None):
        """
        Load a previously made SED from HDF5 file.

        :param filename: name of SED FITS file
        :type filename: string
        :return: True if HDF5 file could be loaded
        :rtype: bool
        """
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.hdf5'
        if not os.path.isfile(filename):
            logger.warning('No previous results saved to HFD5 file {:s}'.format(filename))
            return False

        data = hdf5.read2dict(filename)

        self.master = data.get('master', {})
        self.results = data.get('results', {})
        self.constraints = data.get('constraints', {})

        logger.info('Loaded previous results from HDF5 file: %s'%(filename))
        logger.debug('Loaded following datasets from HDF5 file:\n %s'%(list(data.keys())))
        return True

    def save_bibtex(self):
        """
        Convert the bibcodes in a phot file to a bibtex file.

        The first line in the bibtex file contains a \citet command citing
        all photometry.
        """
        filename = os.path.splitext(self.photfile)[0]+'.bib'
        crossmatch.make_bibtex(self.master,filename=filename)

    def save_summary(self,filename=None,CI_limit=None,method='igrid_search',chi2type='ci_red'):
        """
        Save a summary of the results to an ASCII file.
        """
        # -- open the summary file to write the results
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.sum'

        if CI_limit is None:
            CI_limit = self.CI_limit

        # -- gather the results:
        grid_results = self.results[method]['grid']
        start_CI = np.argmin(np.abs(grid_results[chi2type]-CI_limit))
        factor = self.results[method]['factor']
        names = ['factor','chi2_type','ci_limit']
        results = [factor,chi2type,CI_limit*100]
        for name in grid_results.dtype.names:
            lv,cv,uv = grid_results[name][start_CI:].min(),\
                       grid_results[name][-1],\
                       grid_results[name][start_CI:].max()
            names += [name+'_l',name,name+'_u']
            results += [lv,cv,uv]
        # -- write the used photometry to a file
        include_grid = self.master['include']
        photbands = ":".join(self.master[include_grid]['photband'])
        references = ",".join(self.master[include_grid]['bibcode'])
        used_photometry = photometry2str(self.master[include_grid],comment='#')
        used_atmosphere = '#'+model.defaults2str()+'\n'
        used_photbands = '#'+photbands+'\n'
        used_references = '#'+references
        comments = used_photometry+used_atmosphere+used_photbands+used_references

        contents = np.array([results]).T
        contents = np.rec.fromarrays(contents,names=names)
        ascii.write_array(contents,filename,auto_width=True,header=True,
                          comments=comments.split('\n'),mode='a',use_float='%g')

        logger.info('Saved summary to {0}'.format(filename))


    def save_important_info(self,filename=None,CI_limit=None,method='igrid_search',chi2type='ci_red'):
        """
        Save a summary of the results to an ASCII file.
        """
        # -- open the summary file to write the results
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.sum'

        if CI_limit is None:
            CI_limit = self.CI_limit

        # -- gather the results:
        grid_results = self.results[method]['grid']
        start_CI = np.argmin(np.abs(grid_results[chi2type]-CI_limit))
        factor = self.results[method]['factor']
        #names = ['scaling_factor','chi2_type','ci_limit']
        results = [factor,chi2type,CI_limit*100]
        print('Metallicity:')
        print(grid_results['z'][-1])
        wanted_names = ['ebv','logg','teff','z','chisq']
        for name in wanted_names:
            lv,cv,uv = grid_results[name][start_CI:].min(),\
                       grid_results[name][-1],\
                       grid_results[name][start_CI:].max()
            #names += [name+'_l',name,name+'_u']
            results += ["%.3f"%lv,"%.3f"%cv,"%.3f"%uv]

        contents = np.array([results]).T
        contents = np.rec.fromarrays(contents)
        ascii.write_array(contents,filename,auto_width=True,mode='a',use_float='%g')

        logger.info('Saved summary to {0}'.format(filename))


    #}

class BinarySED(SED):

    def __init__(self,ID=None,photfile=None,plx=None,load_fits=True,load_hdf5=True,label='', **kwargs):
        """
        Setup the Binary sed in the same way as a normal SED.
        The masses of both components can be provided, and will then be used in igrid_search,
        iminimize, and while calculating CI and CI2D confidence intervalls
        """
        super(BinarySED, self).__init__(ID=ID,photfile=photfile,plx=plx,label=label,\
                                             load_fits=load_fits,load_hdf5=load_hdf5)

        self.set_constraints(**kwargs)


    def set_constraints(self, **kwargs):
        """
        Add constraints that are used when fitting the Binary SED. Up till now the following
        contraints are supported:
            - masses (in Msol)
            - distance (in Rsol)
        TODO: This function should in the future accept Units.
        """

        if 'masses' in kwargs:
            self.constraints['masses'] = kwargs['masses']
        if 'distance' in kwargs:
            self.constraints['distance'] = kwargs['distance']

    def constraints2str(self):
        """
        Summarizes all constraints in a string.
        """
        res = ""
        for key in list(self.constraints.keys()):
            res += "Using constraint: %s = %s\n"%(key, self.constraints[key])
        res = res[:-1]
        return res

    def igrid_search(self,points=100000,teffrange=None,loggrange=None,ebvrange=None,\
                    zrange=None,rvrange=((3.1,3.1),(3.1,3.1)),vradrange=((0,0),(0,0)),\
                    radrange=(None,None),compare=True,df=None,CI_limit=None,\
                    set_model=True, distance=None,**kwargs):
        """
        Fit fundamental parameters using a (pre-integrated) grid search.

        If called consecutively, the ranges will be set to the CI_limit of previous
        estimations, unless set explicitly.

        If called for the first time, the ranges will be +/- np.inf by defaults,
        unless set explicitly.
        """

        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit

        # -- set defaults limits
        ranges = self.generate_ranges(teffrange=teffrange[0],\
                        loggrange=loggrange[0],ebvrange=ebvrange[0],\
                        zrange=zrange[0],rvrange=rvrange[0],vradrange=vradrange[0],
                        radrange=radrange[0],teff2range=teffrange[1],\
                        logg2range=loggrange[1],ebv2range=ebvrange[1],\
                        z2range=zrange[1],rv2range=rvrange[1],vrad2range=vradrange[1],
                        rad2range=radrange[1])

        # -- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'%\
                   (photometry2str(self.master[include_grid])))

        logger.info('The following constraints are included in the fitting process:\n%s'%\
                   (self.constraints2str()))

        # -- build the grid, run over the grid and calculate the CHI2
        masses = self.constraints.get('masses', None)
        pars = fit.generate_grid_pix(self.master['photband'][include_grid], masses=masses, points=points, **ranges)

        chisqs,scales,escales,lumis = fit.igrid_search_pix(self.master['cmeas'][include_grid],
                             self.master['e_cmeas'][include_grid],
                             self.master['photband'][include_grid],
                             model_func=model.get_itable_pix,constraints=self.constraints,
                             **pars)
        fitres = dict(chisq=chisqs, scale=scales, escale=escales, labs=lumis)

        # -- collect all the results in a record array
        self.collect_results(grid=pars, fitresults=fitres, mtype='igrid_search')

        # -- do the statistics
        self.calculate_statistics(df=df, ranges=ranges, mtype='igrid_search')

        # -- compute the confidence intervals
        ci = self.calculate_confidence_intervals(mtype='igrid_search',chi2_type='red',\
                                                                     CI_limit=CI_limit)
        self.store_confidence_intervals(mtype='igrid_search', **ci)

        # -- remember the best model
        if set_model: self.set_best_model()

    def generate_fit_param(self, start_from='igrid_search', **pars):
        """
        generates a dictionary with parameter information that can be handled by fit.iminimize
        """
        masses = self.constraints.get('masses', None)
        result = super(BinarySED, self).generate_fit_param(start_from=start_from, **pars)

        # -- Couple ebv and rv of both components
        result['ebv2_expr'] = 'ebv'
        result['rv2_expr'] = 'rv'

        if masses != None:
            # -- Couple the radii to the masses
            G, Msol, Rsol = constants.GG_cgs, constants.Msol_cgs, constants.Rsol_cgs
            result['rad_value'] = np.sqrt(G*masses[0]*Msol/10**result['logg_value'])
            result['rad2_value'] = np.sqrt(G*masses[1]*Msol/10**result['logg2_value'])
            result['rad_expr'] = 'sqrt(%0.5f/10**logg) * 165.63560394542122'%(masses[0])
            result['rad2_expr'] = 'sqrt(%0.5f/10**logg2) * 165.63560394542122'%(masses[1])
            result['rad_min'], result['rad_max'] = None, None
            result['rad2_min'], result['rad2_max'] = None, None
            result['rad_vary'], result['rad2_vary'] = False, False
        else:
            result['rad_value'] = self.results[start_from]['CI']['rad']
            result['rad2_value'] = self.results[start_from]['CI']['rad2']

        return result

    def iminimize(self, teff=(None,None), logg=(None,None), ebv=(None,None), z=(None,None),
                  rv=(None,None), vrad=(0,0), teffrange=(None,None), loggrange=(None,None),
                  ebvrange=(None,None), zrange=(None,None), rvrange=(None,None),
                  vradrange=(None,None), radrange=(None,None), compare=True, df=None,
                  distance=None, start_from='igrid_search', points=None, CI_limit=None,
                  calc_ci=True, set_model=True, **kwargs):
        """ Binary minimizer """

        ranges = self.generate_ranges(teffrange=teffrange[0],\
                        loggrange=loggrange[0],ebvrange=ebvrange[0],\
                        zrange=zrange[0],rvrange=rvrange[0],vradrange=vradrange[0],
                        radrange=radrange[0],teff2range=teffrange[1],\
                        logg2range=loggrange[1],ebv2range=ebvrange[1],\
                        z2range=zrange[1],rv2range=rvrange[1],vrad2range=vradrange[1],
                        rad2range=radrange[1])
        pars = self.generate_fit_param(teff=teff[0],logg=logg[0],ebv=ebv[0],z=z[0],\
                            rv=rv[0], vrad=vrad[0],teff2=teff[1],logg2=logg[1],\
                            ebv2=ebv[1],z2=z[1],rv2=rv[1],vrad2=vrad[1], \
                            start_from=start_from, **ranges)

        # -- Print the used data and constraints
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'%\
                   (photometry2str(self.master[include_grid])))
        logger.info('The following constraints are included in the fitting process:\n%s'%\
                   (self.constraints2str()))

        # -- pass all ranges and starting values to the fitter
        kick_list = ['teff', 'logg', 'teff2', 'logg2', 'ebv']
        grid, chisq, nfev, scale, lumis = fit.iminimize(self.master['cmeas'][include_grid],
               self.master['e_cmeas'][include_grid], self.master['photband'][include_grid],
               points=points, kick_list=kick_list, constraints=self.constraints, **pars)

        logger.info('Minimizer Succes with startpoints=%s, chi2=%s, nfev=%s'%(len(chisq), chisq[0], nfev[0]))
        # -- handle the results
        fitres = dict(chisq=chisq, nfev=nfev, scale=scale, labs=lumis)
        self.collect_results(grid=grid, fitresults=fitres, mtype='iminimize', selfact='chisq')

        # -- Do the statistics
        self.calculate_statistics(df=5, ranges=ranges, mtype='iminimize', selfact='chisq')

        # -- Store the confidence intervals
        ci = self._get_imin_ci(mtype='iminimize',**ranges)
        self.store_confidence_intervals(mtype='iminimize', **ci)
        if calc_ci:
            self.calculate_iminimize_CI(mtype='iminimize', CI_limit=CI_limit)

        # -- remember the best model
        if set_model: self.set_best_model(mtype='iminimize')


    def calculate_iminimize_CI(self, mtype='iminimize', CI_limit=0.66, **kwargs):
        """
        Calculate the confidence intervals for each parameter using the lmfit
        calculate confidence interval method.

        The calculated confidence intervals are stored in the results['CI']
        dictionary. If the method fails, or if the asked CI is outside the
        provided ranges, those ranges will be set as CI.

        This method works in the same way as for a single SED, but it adds
        the mass as an extra constraint if the mass of both components is stored.
        """
        masses = self.constraints.get('masses',None)
        super(BinarySED, self).calculate_iminimize_CI(mtype=mtype, CI_limit=CI_limit,\
                                         masses=masses, constraints=self.constraints,\
                                         **kwargs)

    def calculate_iminimize_CI2D(self,xpar, ypar, mtype='iminimize', limits=None, res=10, **kwargs):
        """
        Calculated 2 dimentional confidence intervals for the given parameters,
        using lmfit methods.

        The calculated confidence intervals are stored in the results['CI2D']
        dictionary.

        This method works in the same way as for a single SED, but it adds
        the mass as an extra constraint if the mass of both components is stored.
        """
        masses = self.constraints.get('masses',None)
        super(BinarySED, self).calculate_iminimize_CI2D(xpar, ypar, mtype=mtype,\
                                 limits=limits, res=res, masses=masses, **kwargs)

    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004', **kwargs):
        """
        Get reddenend and unreddened model
        """
        logger.info('Interpolating approximate full SED of best model')

        # -- synthetic flux
        include = self.master['include']
        synflux = np.zeros(len(self.master['photband']))
        keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
        keep = keep & include

        if mtype in ['igrid_search', 'iminimize']:
            scale = self.results[mtype]['CI']['scale']

            # -- get (approximated) reddened and unreddened model
            wave,flux = model.get_table_multiple(teff=(self.results[mtype]['CI']['teff'],self.results[mtype]['CI']['teff2']),
                                    logg=(self.results[mtype]['CI']['logg'],self.results[mtype]['CI']['logg2']),
                                    ebv=(self.results[mtype]['CI']['ebv'],self.results[mtype]['CI']['ebv2']),
                                    radius=(self.results[mtype]['CI']['rad'],self.results[mtype]['CI']['rad2']),
                                    law=law)
            wave_ur,flux_ur = model.get_table_multiple(teff=(self.results[mtype]['CI']['teff'],self.results[mtype]['CI']['teff2']),
                                      logg=(self.results[mtype]['CI']['logg'],self.results[mtype]['CI']['logg2']),
                                      ebv=(0,0),
                                      radius=(self.results[mtype]['CI']['rad'],self.results[mtype]['CI']['rad2']),
                                      law=law)
            # -- get synthetic photometry
            pars = {}
            for key in list(self.results[mtype]['CI'].keys()):
                if not key[-2:] == '_u' and not key[-2:] == '_l':
                    pars[key] = self.results[mtype]['CI'][key]
            synflux_,pars = model.get_itable(photbands=self.master['photband'][keep], **pars)
            flux,flux_ur = flux*scale,flux_ur*scale

            synflux[keep] = synflux_

            synflux[-self.master['color']] *= scale
            chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
            # -- calculate effective wavelengths of the photometric bands via the model
            #   values
            eff_waves = filters.eff_wave(self.master['photband'],model=(wave,flux))
            self.results[mtype]['model'] = wave,flux,flux_ur
            self.results[mtype]['synflux'] = eff_waves,synflux,self.master['photband']
            self.results[mtype]['chi2'] = chi2


class PulsatingSED(SED):

    def __init__(self,ID=None,photfile=None,plx=None,load_fits=True,load_hdf5=True,label='', **kwargs):
        """
        Setup the Binary sed in the same way as a normal SED.
        The masses of both components can be provided, and will then be used in igrid_search,
        iminimize, and while calculating CI and CI2D confidence intervalls
        """
        super(PulsatingSED, self).__init__(ID=ID,photfile=photfile,plx=plx,label=label,\
                                             load_fits=load_fits,load_hdf5=load_hdf5)

        self.set_constraints(**kwargs)


    def set_constraints(self, **kwargs):
        """
        Add constraints that are used when fitting the Pulsating SED. Up till now the following
        contraints are supported (and in fact needed to do a fitting):
            - deltaTeff (in K)
            - deltaLogg (in dex)
        For each 'phase' additional to the zeroth phase, a tuple containing a lower and upper limit
        needs to be added to the deltaTeff and deltaLogg list. These limits have a sign!

        TODO: This function should in the future accept Units.
        """

        #if 'mass' in kwargs:
        #    self.constraints['mass'] = kwargs['mass']
        #if 'distance' in kwargs:
            #self.constraints['distance'] = kwargs['distance']
        if 'deltaTeff' in kwargs:
            deltaTeff = kwargs.get('deltaTeff')
            print(deltaTeff)
            for i in range(len(deltaTeff)):
                self.constraints['delta{}Teff'.format(i+1)] = deltaTeff[i]
        if 'deltaLogg' in kwargs:
            deltaLogg = kwargs.get('deltaLogg')
            print(deltaLogg)
            for i in range(len(deltaTeff)):
                self.constraints['delta{}Logg'.format(i+1)] = deltaLogg[i]

    def constraints2str(self):
        """
        Summarizes all constraints in a string.
        """
        res = ""
        for key in list(self.constraints.keys()):
            res += "Using constraint: %s = %s\n"%(key, self.constraints[key])
        res = res[:-1]
        return res

    def igrid_search(self,points=100000,teffrange=None,loggrange=None,ebvrange=None,\
                    zrange=None,rvrange=(3.1,3.1),vradrange=None,df=None,CI_limit=None,\
                    set_model=True, distance=None,**kwargs):
        """
        Fit fundamental parameters at various phases simultaneously, using a (pre-integrated) grid search.

        Parameter ranges only need to be provided for the first 'phase', but be aware that
        these ranges should encompass the parameters of your star at all phases! The part of the
        parameter space that falls outside the box defined for the first 'phase' is not
        probed at other phases, even though the defined constraints would in principle allow it.

        """

        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit

        ## -- the data sets
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'%\
                   (photometry2str(self.master[include_grid])))

        phase_indices = self.master['phase']
        unique_phases = np.unique(phase_indices)
        logger.info('Data sets with the following phase indices are included in the fitting process:\n%s'%\
                   (str(unique_phases)))

        ## -- set defaults limits
        ranges = self.generate_ranges(teffrange=teffrange,loggrange=loggrange,ebvrange=ebvrange,zrange=zrange,rvrange=rvrange,vradrange=vradrange)

        ## -- the constraints used in the fit
        logger.info('The following constraints are included in the fitting process:\n%s'%\
                   (self.constraints2str()))

        ## -- build the grid, run over the grid and calculate the CHI2
        newkwargs = ranges.copy()
        newkwargs.update(self.constraints)
        pars = fit.generate_grid_pix_pulsating2(self.master['photband'][include_grid],\
                                                points=points, unique_phases=unique_phases,\
                                                **newkwargs)

        chisqs,scales,escales,lumis = fit.igrid_search_pix2(self.master['cmeas'][include_grid],
                             self.master['e_cmeas'][include_grid],\
                             self.master['photband'][include_grid],phases=self.master['phase'][include_grid],\
                             unique_phases=unique_phases,model_func=model.get_itable_pix,\
                             constraints=self.constraints,stat_func=fit.stat_chi2b,\
                             **pars)

        ## -- collect all the results in record arrays, remove the points that fall
        #    out of the model grid (the nan in chisqs), do the statistics and
        #    compute the confidence intervals
        for i in range(len(unique_phases)):
            fitres = dict(chisq=chisqs[:,i], scale=scales[:,i], escale=escales[:,i], labs=lumis[:,i])
            partpars = pars.fromkeys(list(pars.keys()))
            for key in list(pars.keys()):
                partpars[key] = pars[key][:,i]
            self.collect_results(grid=partpars, fitresults=fitres, mtype='igrid_search_{}'.format(unique_phases[i]))

            ## -- do the statistics
            self.calculate_statistics(df=df, ranges=ranges, mtype='igrid_search_{}'.format(unique_phases[i]))

            ## -- compute the confidence intervals
            ci = self.calculate_confidence_intervals(mtype='igrid_search_{}'.format(unique_phases[i]),\
                                                     chi2_type='red',CI_limit=CI_limit)
            self.store_confidence_intervals(mtype='igrid_search_{}'.format(unique_phases[i]), **ci)

            ## -- remember the best model
            if set_model: self.set_best_model(mtype='igrid_search_{}'.format(unique_phases[i]))

        ## -- collect the results of the all-inclusive fit
        index = len(unique_phases)
        fitres = dict(chisq=chisqs[:,index], scale=scales[:,index], escale=escales[:,index], labs=lumis[:,index])
        allpars = pars.fromkeys(list(pars.keys()))
        for key in list(pars.keys()):
            allpars[key] = pars[key][:,0]
        for i in range(len(unique_phases)-1):
            fitres['scale{}'.format(i+1)] = scales[:,index+i+1]
            fitres['escale{}'.format(i+1)] = escales[:,index+i+1]
            fitres['labs{}'.format(i+1)] = lumis[:,index+i+1]
            allpars['teff{}'.format(i+1)] = pars['teff'][:,i+1]
            allpars['logg{}'.format(i+1)] = pars['logg'][:,i+1]
            ranges['teff{}range'.format(i+1)] = self.constraints['delta{}Teff'.format(i+1)]
            ranges['logg{}range'.format(i+1)] = self.constraints['delta{}Logg'.format(i+1)]

        self.collect_results(grid=allpars, fitresults=fitres, mtype='igrid_search_all')

        ## -- do the statistics of the all-inclusive fit
        self.calculate_statistics(df=df, ranges=ranges, mtype='igrid_search_all')

        ## -- compute the confidence intervals of the all-inclusive fit
        ci = self.calculate_confidence_intervals(mtype='igrid_search_all',\
                                                 chi2_type='red',CI_limit=CI_limit)
        self.store_confidence_intervals(mtype='igrid_search_all', **ci)

        ## -- remember the best model
        if set_model: self.set_best_model(mtype='igrid_search_all')

    def generate_fit_param(self, start_from='igrid_search', **pars):
        """
        generates a dictionary with parameter information that can be handled by fit.iminimize
        """
        raise NotImplementedError

    def calculate_iminimize_CI(self, mtype='iminimize', CI_limit=0.66, **kwargs):
        raise NotImplementedError

    def calculate_iminimize_CI2D(self, xpar, ypar, mtype='iminimize', limits=None, res=10, **kwargs):
        raise NotImplementedError

    def iminimize(self, teff=None, logg=None, ebv=None, z=0, rv=3.1, vrad=0, teffrange=None,
                     loggrange=None, ebvrange=None, zrange=None, rvrange=None, vradrange=None,
                     points=None, distance=None, start_from='igrid_search',df=None, CI_limit=None,
                     calc_ci=False, set_model=True, **kwargs):
        raise NotImplementedError

    def imc(self,teffrange=None,loggrange=None,ebvrange=None,zrange=None,start_from='igrid_search',\
                 distribution='uniform',points=None,fitmethod='fmin',disturb=True):
        raise NotImplementedError

    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004',**kwargs):
        """
        Get reddenend and unreddened model
        """
        super(PulsatingSED, self).set_best_model(mtype=mtype,law=law,**kwargs)

        if '_all' in mtype:
            logger.info('Interpolating approximate full SED of best model at additional phases.')

            # -- synthetic flux
            include = self.master['include']
            phases = self.master['phase']
            unique_phases = np.unique(phases)
            synflux = self.results[mtype]['synflux'][1]
            keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
            keep = keep & include

            if ('igrid_search' in mtype) | ('iminimize' in mtype): #mtype in ['igrid_search', 'iminimize']:
                # -- get the metallicity right
                files = model.get_file(z='*')
                if type(files) == str: files = [files] #files needs to be a list!
                metals = np.array([pf.getheader(ff)['Z'] for ff in files])
                metals = metals[np.argmin(np.abs(metals-self.results[mtype]['CI']['z']))]

                for i in range(len(unique_phases)-1):
                    keep = keep & (phases == unique_phases[i+1])
                    scale = self.results[mtype]['CI']['scale{}'.format(i+1)]
                    # -- get (approximated) reddened and unreddened model
                    wave,flux = model.get_table(teff=self.results[mtype]['CI']['teff{}'.format(i+1)],
                                        logg=self.results[mtype]['CI']['logg{}'.format(i+1)],
                                        ebv=self.results[mtype]['CI']['ebv'],
                                        z=metals,
                                        law=law)
                    wave_ur,flux_ur = model.get_table(teff=self.results[mtype]['CI']['teff{}'.format(i+1)],
                                            logg=self.results[mtype]['CI']['logg{}'.format(i+1)],
                                            ebv=0,
                                            z=metals,
                                            law=law)
                    # -- get synthetic photometry
                    synflux_,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff{}'.format(i+1)],
                                   logg=self.results[mtype]['CI']['logg{}'.format(i+1)],
                                   ebv=self.results[mtype]['CI']['ebv'],
                                   z=self.results[mtype]['CI']['z'],
                                   photbands=self.master['photband'][keep])

                    flux,flux_ur = flux*scale,flux_ur*scale

                    synflux[keep] = synflux_

                    #synflux,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
                    #                          logg=self.results[mtype]['CI']['logg'],
                    #                          ebv=self.results[mtype]['CI']['ebv'],
                    #                          photbands=self.master['photband'])
                    synflux[-self.master['color'] & keep] *= scale
                    chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
                    # -- calculate effective wavelengths of the photometric bands via the model
                    #   values
                    eff_waves = filters.eff_wave(self.master['photband'],model=(wave,flux))
                    self.results[mtype]['model{}'.format(i+1)] = wave,flux,flux_ur
                    self.results[mtype]['synflux'] = eff_waves,synflux,self.master['photband']
                    self.results[mtype]['chi2'] = chi2

    @standalone_figure
    def plot_sed(self,colors=False,mtype='igrid_search',plot_redded=True,plot_deredded=False,
            plot_unselected=True,wave_units='AA',flux_units='erg/s/cm2',**kwargs):
        """
        Plot a fitted SED together with the data.

        Example usage:

        First set the SED:

        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()

        Then make the plots:

        >>> p = pl.figure()
        >>> p = pl.subplot(121)
        >>> mysed.plot_sed(colors=False)
        >>> p = pl.subplot(122)
        >>> mysed.plot_sed(colors=True)

        ]]include figure]]ivs_sed_builder_plot_sed_01.png]
        """
        if 'phase' in kwargs:
            uniquephase = kwargs.pop('phase')
            phase = True
        else:
            uniquephase = None
            phase = False
        annotation = kwargs.pop('annotation',True)

        def plot_sed_getcolors(master,color_dict=None):
            myphotbands = [iphotb.split('.')[1] for iphotb in master['photband'][master['color']]]
            if not myphotbands:  # -- If there are no colours none can be returned (added by joris 30-01-2012)
                return [],[],[],None
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

        # -- get the color cycle
        systems = np.array([system.split('.')[0] for system in self.master['photband']],str)
        set_systems = sorted(list(set(systems)))

        if ('absolutesymbolcolor' in kwargs) and (kwargs.pop('absolutesymbolcolor') == True):
            sortedphotsystems,plotcolorvalues = filters.get_plotsymbolcolorinfo()
            selectedcolorinds = np.array([np.where(sortedphotsystems == system)[0][0] for system in set_systems])
            color_cycle = itertools.cycle([pl.cm.spectral(j) for j in plotcolorvalues[selectedcolorinds]])
        else:
            color_cycle = itertools.cycle([pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))])

        # -- for plotting reasons, we translate every color to an integer
        for system in set_systems:
            color = next(color_cycle)
            keep = (systems==system) & (self.master['color']==colors)
            if not plot_unselected:
                keep = keep & self.master['include']
            if phase:
                keep = keep & (self.master['phase'] == uniquephase)
            # -- synthetic:
            if sum(keep) and mtype in self.results and 'synflux' in self.results[mtype]:
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                    y = self.results[mtype]['synflux'][1][keep]
                else:
                    x = self.results[mtype]['synflux'][0][keep]
                    y = self.results[mtype]['synflux'][1][keep]
                    # -- convert to correct units
                    y = conversions.convert('erg/s/cm2/AA',flux_units,y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)
                pl.plot(x,y,'x',ms=10,mew=2,alpha=0.75,color=color,**kwargs)
            # -- include:
            keep = (systems==system) & (self.master['color']==colors) & self.master['include'] & (self.master['phase'] == uniquephase)
            if sum(keep):
                if colors:
                    # -- translate every color to an integer
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    if mtype in self.results and 'synflux' in self.results[mtype]:
                        x = self.results[mtype]['synflux'][0][keep]
                    else:
                        x = self.master['cwave'][keep]
                    y = self.master['cmeas'][keep]
                    e_y = self.master['e_cmeas'][keep]
                    y,e_y = conversions.convert('erg/s/cm2/AA',flux_units,y,e_y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)


                p = pl.errorbar(x,y,yerr=e_y,fmt='o',label=system,
                                capsize=10,ms=7,color=color,**kwargs)

            # -- exclude:
            label = np.any(keep) and '_nolegend_' or system
            keep = (systems==system) & (self.master['color']==colors) & -self.master['include']
            if phase:
                keep = keep & (self.master['phase'] == uniquephase)
            if sum(keep) and plot_unselected:
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results[mtype]['synflux'][0][keep]
                    if np.any(np.isnan(x)):
                        x = self.master['cwave'][keep]
                    y = self.master['cmeas'][keep]
                    e_y = self.master['e_cmeas'][keep]
                    y,e_y = conversions.convert('erg/s/cm2/AA',flux_units,y,e_y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)
                pl.errorbar(x,y,yerr=e_y,fmt='o',label=label,
                            capsize=10,ms=7,mew=2,color=color,mfc='1',mec=color,**kwargs)

        # -- only set logarithmic scale if absolute fluxes are plotted
        #   and only plot the real model then
        if not colors:
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            pl.gca().set_autoscale_on(False)

            # -- the model
            if phase & (uniquephase != 0):
                modelname = 'model{}'.format(uniquephase)
            else:
                modelname = 'model'

            if mtype in self.results and modelname in self.results[mtype]:
                wave,flux,flux_ur = self.results[mtype][modelname]

                flux = conversions.convert('erg/s/cm2/AA',flux_units,flux,wave=(wave,'AA'))
                flux_ur = conversions.convert('erg/s/cm2/AA',flux_units,flux_ur,wave=(wave,'AA'))
                wave = conversions.convert('AA',wave_units,wave)

                if plot_redded:
                    pl.plot(wave,flux,'r-',**kwargs)
                if plot_deredded:
                    pl.plot(wave,flux_ur,'k-',**kwargs)
            pl.ylabel(conversions.unit2texlabel(flux_units,full=True))
            pl.xlabel('wavelength [{0}]'.format(conversions.unit2texlabel(wave_units)))
        else:
            xlabels = list(color_dict.keys())
            xticks = [color_dict[key] for key in xlabels]
            pl.xticks(xticks,xlabels,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Color name')
            pl.xlim(min(xticks)-0.5,max(xticks)+0.5)

        pl.grid()
        if colors:
            leg = pl.legend(loc='best',prop=dict(size='x-small'))
        else:
            leg = pl.legend(loc='upper right',prop=dict(size='x-small'))
            leg.get_frame().set_alpha(0.5)
        loc = (0.45,0.05)
        if mtype in self.results and 'grid' in self.results[mtype] and annotation:
            if phase & (uniquephase != 0):
                exten = uniquephase
            else:
                exten = ''
            teff = self.results[mtype]['grid']['teff{}'.format(exten)][-1]
            logg = self.results[mtype]['grid']['logg{}'.format(exten)][-1]
            ebv = self.results[mtype]['grid']['ebv'][-1]
            met = self.results[mtype]['grid']['z'][-1]
            scale = self.results[mtype]['grid']['scale{}'.format(exten)][-1]
            angdiam = 2*conversions.convert('sr','mas',scale)
            pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\nlogZ=%.3f dex\n$\Theta$=%.3g mas'%(teff,logg,ebv,met,angdiam),
                        loc,xycoords='axes fraction')
        logger.info('Plotted SED as %s'%(colors and 'colors' or 'absolute fluxes'))
        teff = "%d" % teff
        logg = "%.2f" % logg
        ebv = "%.3f" % ebv
        metallicity = "%.3f" % met
        '''f = open('/home/anae/python/SEDfitting/Bastars_info.txt', 'a')
        f.writelines(str(teff)+'\t'+str(logg)+'\t'+str(ebv)+'\n'+'\t'+str(metallicity)+'\n')
        f.closed'''


class Calibrator(SED):
    """
    Convenience class for a photometric standard star or calibrator.
    """
    def __init__(self,ID=None,photfile=None,plx=None,load_fits=False,label='',library='calspec'):
        names,fits_files,phot_files = model.read_calibrator_info(library=library)
        index = names.index(ID)
        # --retrieve fitsfile information
        if library in ['ngsl','calspec']:
            fits_file = pf.open(fits_files[index])
            wave = fits_file[1].data.field('wavelength')
            flux = fits_file[1].data.field('flux')
            fits_file.close()
        elif library=='stelib':
            wave,flux = fits.read_spectrum(fits_files[index])
        # --photfile:
        photfile = phot_files[index]
        super(Calibrator,self).__init__(photfile=photfile,plx=plx,load_fits=load_fits,label=label)
        self.set_model(wave,flux)


class SampleSEDs(object):
    """
    Class representing a list of SEDs.
    """
    def __init__(self,targets,**kwargs):
        """
        Initialize a sample.

        This can be done either with a list of IDs or with a list of SED
        instances. The SED must exist! That is, for each ID, there must be a
        phot or FITS file.
        """
        # -- we don't want to load the FITS files by default because they
        #   can eat up memory
        kwargs.setdefault('load_fits',False)
        self.targets = []
        self.seds = []
        # -- create all SEDs
        for target in targets:
            # -- perhaps we gave an ID: in that case, create the SED instance
            if not isinstance(target,SED):
                mysed = SED(target,**kwargs)
                if not mysed.has_photfile():
                    raise ValueError("No phot file found for {}".format(target))
            # -- perhaps already an SED object: then do nothing
            else:
                mysed = target
            self.seds.append(mysed)
            self.targets.append(mysed.ID)
        logger.info("Loaded {} SEDs".format(len(self.seds)))

    def __iter__(self):
        """
        Allow iteration over the SED instances.
        """
        for sed in self.seds:
            yield sed

    def __len__(self):
        """
        The length of a SampleSEDs instance is the number SED instances.
        """
        return len(self.seds)

    def __getitem__(self,key):
        """
        Implements various ways to get individual seds.

        Allows integer indexing, slicing, indexing with integer and boolean arrays.
        """
        # -- via slicing
        if isinstance(key,slice):
            return SampleSEDs([self[ii] for ii in range(*key.indices(len(self)))])
        # -- via an integer
        elif isinstance(key,int):
            return self.seds[key]
        else:
            # -- try to make the input an array
            try:
                key = np.array(key)
            except:
                raise TypeError("Cannot use instance of type {} for indexing".format(type(key)))
            # -- integer array slicing
            if key.dtype==np.dtype(int):
                return SampleSEDs([self[ii] for ii in key])
            # -- boolean array slicing
            elif key.dtype==np.dtype(bool):
                return SampleSEDs([self[ii] for ii in range(len(key)) if key[ii]])
            # -- that's all I can come up with
            else:
                raise TypeError("Cannot use arrays of type {} for indexing".format(key.dtype))

    def summarize(self):
        # -- collect the names of all the different sources
        sources = {}
        for sed in self.seds:
            # -- collect the different sources
            these_sources = list(set(sed.master['source']))
            for source in these_sources:
                if not source in sources:
                    sources[source] = []
                # -- now for each source, collect the names of the passbands
                keep = sed.master['source']==source
                sources[source] += list(set(sed.master[keep]['photband']))
                sources[source] = sorted(list(set(sources[source])))
        # -- next step: for each source, create a record array where the columns
        #   are the different photbands. Fill in the values for the photbands
        #   for each target when possible.
        summary = []
        source_names = sorted(sources.keys())
        for source in source_names:
            # -- create the record array
            data = np.zeros((len(self.targets),2*len(sources[source])))
            #   but remember to have errors with the photbands
            names = []
            for photband in sources[source]:
                names += [photband,'e_'+photband]
            data = np.rec.fromarrays(data.T,names=names)
            # -- fill in the values
            for i,sed in enumerate(self.seds):
                for photband in sources[source]:
                    keep = (sed.master['source']==source) & (sed.master['photband']==photband)
                    # -- fill in nans for value and error when not present
                    if not sum(keep):
                        data[photband][i] = np.nan
                        data['e_'+photband][i] = np.nan
                    # -- otherwise give the first value you've found
                    else:
                        data[photband][i] = sed.master[keep]['cmeas'][0]
                        data['e_'+photband][i] = sed.master[keep]['e_cmeas'][0]
                        #if not sum(keep): logger.warning('multiple defined photband ({}) and source ({})'.format(photband,source))
            summary.append(data)
        dtypes = [(source,summary[i].dtype) for i,source in enumerate(source_names)]
        output = np.zeros(len(self.targets),dtype=np.dtype(dtypes))
        for i,name in enumerate(output.dtype.names):
            output[name] = summary[i]
        return output

    def get_data(self,source,photband,label=None):
        """
        Get all data on a particular passband from a particular source.

        If label is not None, synthetic flux from a model will be added to the
        columns.
        """
        records = []
        if label is not None:
            synflux_label = [('synflux','f8')]
        else:
            synflux_label = []
        for targetname,sed in zip(self.targets,self.seds):
            # -- for the source, collect the names of the passbands
            keep = (sed.master['source']==source) & (sed.master['photband']==photband)
            # -- is there a model SED for which we want to retrieve synthetic photometry?
            if label is not None:
                if sum(keep)==0:
                    synflux = [0]
                else:
                    synflux = [sed.results[label]['synflux'][1][keep][0]]
            else:
                synflux = []
            # -- if no data on the matter is present, put zero values everywhere
            if sum(keep)==0:
                records.append([targetname]+list(np.zeros(len(sed.master.dtype.names)))+synflux)
            else:
                records.append([targetname]+list(sed.master[keep][0])+synflux)
        # -- make the thing into a record array
        dtypes = np.dtype([('targetname','|S25')] + sed.master.dtype.descr + synflux_label)
        output = np.rec.fromrecords(records,names=dtypes.names)
        output = np.array(output,dtype=dtypes)
        return output

    def get_confidence_interval(self,parameter='teff',mtype='igrid_search'):
        values = np.zeros((len(self),3))
        for i,sed in enumerate(self):
            sed.load_fits()
            values[i,0] = sed.results[mtype]['CI'][parameter+'_l']
            values[i,1] = sed.results[mtype]['CI'][parameter]
            values[i,2] = sed.results[mtype]['CI'][parameter+'_u']
            sed.clear()
        return values




if __name__ == "__main__":
    import sys
    import doctest
    import pprint
    from ivs.aux import loggers

    if not sys.argv[1:]:
        doctest.testmod()
        pl.show()
    else:
        name = " ".join([string for string in sys.argv[1:] if not '=' in string])
        units = [string.split('=')[1] for string in sys.argv[1:] if 'units=' in string]
        if not units:
            units = 'erg/s/cm2/AA'
        else:
            units = units[0]
        logger = loggers.get_basic_logger("")
        mysed = SED(name)
        pprint.PrettyPrinter(indent=4).pprint(mysed.info)
        mysed.get_photometry(units=units)
        mysed.plot_data()
        pl.show()
        answer = input('Keep photometry file %s (y/N)'%(mysed.photfile))
        if not 'y' in answer.lower():
            os.unlink(mysed.photfile)
            logger.info('Removed %s'%(mysed.photfile))
    raise SystemExit

    # -- clean up
    if os.path.isfile('HD180642.fits'):
        os.remove('HD180642.fits')
    raise SystemExit
    # -- PCA analysis
    master['include'] = True
    exclude(master,names=['STROMGREN.HBN-HBW','USNOB1'],wrange=(1.5e4,np.inf))
    do_pca = False
    include_pca = master['include']

    if sum(keep)>2:
        do_pca = True
        logger.info("Start of PCA analysis to find fundamental parameters")
        colors,index = np.unique(master['photband'][include_pca],return_index=True)
        A,grid = fit.get_PCA_grid(colors,ebvrange=(0,0.5),res=10)
        P,T,(means,stds) = fit.get_PCA(A)
        calib = fit.calibrate_PCA(T,grid,function='linear')
        obsT,e_obsT = master['cmeas'][include_pca][index], master['e_cmeas'][include_pca][index]
        pars = fit.get_PCA_parameters(obsT,calib,P,means,stds,e_obsT=e_obsT,mc=mc)
        teff,logg,ebv = pars[0]
        logger.info("PCA result: teff=%.0f, logg=%.2f, E(B-V)=%.3f"%(teff,logg,ebv))

        # -- find angular diameter
        logger.info('Estimation of angular diameter')
        iflux = model.get_itable(teff=teff,logg=logg,ebv=ebv,photbands=master['photband'][include_grid])
        scale_pca = np.median(master['cmeas'][include_grid]/iflux)
        angdiam = 2*np.arctan(np.sqrt(scale_pca))/np.pi*180*3600*1000
        logger.info('Angular diameter = %.3f mas'%(angdiam))

        # -- get model
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
