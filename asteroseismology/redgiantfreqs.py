"""
Functions for Red Giant Seismology.

Contents:
=========

This library contains several functions which allow to identify modes 
in the power-spectra of red giants. For the description of the pressure modes, 
the formulism of Mosser et al. (2011A&A...525L...9M) for the 'Universal 
Oscillation Pattern' is available.

For the description of the frequencies and rotational splitting of the mixed 
dipole modes, the formalism of Mosser et al. (2012A&A...540A.143M) and (2012arXiv1209.3336M)
has been scripted. To apply this method, one needs to know the value of the true 
period spacing. The modulation of the mixing between pressure- and gravity-dipole 
modes is described through a Lorentzian profile. This approach assumes a constant 
Period Spacing for all radial orders. Comparing the modelled frequencies with the 
real frequencies, found in an spectrum, will work well for the several radial
orders but deviations can be found in other regions of the spectrum.
A similar Lorentzian description is valid for rotational splitting.


Contents:
---------
    * purePressureMode: calculates the position of the pure pressure pressure 
    mode of a spherical degree \ell and radial order.

    * universalPattern: creates an array with the frequencies of all pure pressure
    modes with \ell=0,1,2 and 3, whereby the index i selectes all m
    odes of the same spherical degree:
    e.g.: asymptoticRelation[0] = all \ell=0, asymptoticRelation[2] = all \ell=2, 
    asymptoticRelation[-1] gives an array with the radial orders (as float).

    * asymptoticRelation: returns an array with the frequencies of the pure 
    pressure modes of the pure pressure modes 

    * rotational splitting: calculates the expected rotational splitting 
    of a dipole mode based on the degree of mixing
    

Example: KIC 4448777
--------------------
This star is a H-shell burning red giant with rotational splitting.
The powerspectrum can be found here: ...

>>> import ivs.asteroseismology.redgiantfreqs as rg
>>> centralRadial=  226.868   	# [muHz], Frequency of the central radial mode
>>> largeSep	= 17.000		# [muHz], large Frequency separation in powerspectrum
>>> smallSep	= 2.201			# [muHz], small Frequency separation in powerspectrum
>>> 
>>> DeltaPg		= 89.9			# [seconds],	true period spacing of dipole modes
>>> dfMax		= 0.45			# [muHz], normalized rotational splitting (i.e. |freq_{m=0} - freq_{m=+/-1}|
>>> 
>>> qCoupling= 0.15; lambdaParam = 0.5; beParam = 0.08 #Default values

First model the frequencies of the pure pressure modes.

>>> frequenciesUniversalPattern = rg.universalPattern(largeSep, smallSep, centralRadial)

Next, we calculate the frequencies of the mixed dipole modes:

>>> mixedModesPattern = rg.asymptoticRelation(centralRadial,largeSep, DeltaPg, qCoupling)

>>> rotationalSplittingModulation = rg.symmetricRotationalModulation(dfMax, largeSep, frequenciesUniversalPattern, mixedModesPattern, lambdaParam, beParam)




"""

import sys, math, os
import scipy as scipy
import pylab as pl
import numpy as np
import numpy.random as npr
import random as rd
from time import gmtime, strftime
import glob as glob
from scipy import stats
import logging

#-- define formats

format  = ''
datefmt = ''

logging.basicConfig(level=logging.INFO,format=format,datefmt=datefmt)
logger = logging.getLogger('IVS.RG')



def purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, ell):
    """
    Calculates and returns the frequency of a single pure pressure mode for a spherical (\ell=0,1,2, and 3)[Nota Bene: 3 to be implemented].
    
    smallSep can be set to an abitrary value if not calculating \ell=2 modes.

	@param nnn: radial order of the mode with respect to the radial order of the central radial mode.
	@type nnn: float or integer
	@param largeSep: large Frequency separation in powerspectrum, in muHz
	@type largeSep: float
	@param smallSep: small Frequency separation in powerspectrum, in muHz
	@type smallSep: float
	@param centralRadial: Frequency of the central radial mode, in muHz
	@type centralRadial: float	
	@param epsilon: echelle phase shift of the central Radial mode:  epsilon = centralRadial / largeSep - centralRadial // largeSep
	@type epsilon: float
	@param ell:  spherical degree of a mode
	@type ell: float or integer
	@return frequency in muHz
	@rtype float array
    """

    ell = np.float(ell)
    nnn = np.float(nnn)
    smallSep /= largeSep 
    dominantRadialOrderPressure = np.floor(centralRadial / largeSep)
    alphaCorr = 0.015 * largeSep**(-0.32)

    if ell == 0.: constantDEll = 0.0               # for radial modes
    if ell == 1.: constantDEll = -(ell/2.+0.025)   # for dipole modes
    if ell == 2.: constantDEll = smallSep          # for quadrupole modes
    if ell == 3.: constantDEll = 0.0               # for radial modes

    purePressureModeEll = largeSep * ( (dominantRadialOrderPressure + nnn) + epsilon - constantDEll + alphaCorr/2.*(nnn)**2.)

    return purePressureModeEll
    





def universalPattern(largeSep, smallSep, centralRadial, numberOfRadialOrders = 3):

    """ 
    Calculates and returns an array containing frequencies the frequencies of the pure pressure modes 
    of the degrees spherical degrees \ell= 0,1,2 and 3 for a defined range of radial orders.  
    
    When slicing, the index number is equal to the spherical degree \ell, universalPattern[0] = all \ell=0, universalPattern[2] = all \ell=2, 
    The number of the radial order is given as output[-1]. If desired, the array can also be printed to the screen.

	@param largeSep: large Frequency separation in powerspectrum, in muHz
	@type largeSep: float
	@param smallSep: small Frequency separation in powerspectrum, in muHz
	@type smallSep: float
	@param centralRadial: Frequency of the central radial mode, in muHz
	@type centralRadial: float	
	@param numberOfRadialOrders: for how many radial orders above and below the Central Radial Mode should the frequencies be calculated.
	@type numberOfRadialOrders: float or integer

	@return array with frequencies in muHz, [-1] indicates the radial order.
	@rtype float array
    """
    
   
    epsilon = centralRadial / largeSep - centralRadial // largeSep
    radialOrderRange = np.arange(-numberOfRadialOrders,numberOfRadialOrders+1,dtype='float')
    dominantRadialOrderPressure = np.floor(centralRadial/largeSep)
    logger.info('dominantRadialOrderPressure {0}'.format(dominantRadialOrderPressure))
    univPatternPerOrder = []
    
    for nnn in radialOrderRange:
        radial =     purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 0.)
        dipole =     purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 1.)
        quadrupole = purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 2.)
        septupole =  purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 3.)
        
        univPatternPerOrder.append([radial,dipole,quadrupole,septupole,nnn+dominantRadialOrderPressure])
    
    univPatternPerOrder = np.array(univPatternPerOrder, dtype='float')
    
 
    logger.info('radial\t\tdipole\t\tquadrupole\tseptupole\tRadial Order (Pressure)')
    logger.info('--------------------------------------------------------------------------------------')
    for modeValue in univPatternPerOrder: logger.info('{0:.3f}\t\t{1:.3f}\t\t{2:.3f}\t\t{3:.3f}\t\t{4}'.format(modeValue[0],modeValue[1],modeValue[2],modeValue[3],modeValue[4]))
    logger.info('--------------------------------------------------------------------------------------')  

    return univPatternPerOrder.T
    




def asymptoticRelation(centralRadial, largeSep, DeltaPg, qCoupling= 0.15,
                       approximationTreshold=0.001, numberOfRadialOrders = 4):
    """ Calculates and returns the frequencies mixed dipole modes in the powerspectrum of a red-giant star.
    
    This function the modulation of the mixing between pressure- and gravity-dipole and derives the frequencies of the mixed modes from it.

	@param centralRadial: Frequency of the central radial mode, in muHz
	@type centralRadial: float	
	@param largeSep: large Frequency separation in powerspectrum, in muHz
	@type largeSep: float
	@param smallSep: small Frequency separation in powerspectrum, in muHz
	@type smallSep: float
	@param DeltaPg: value of the true period spacing
	@type DeltaPg: float
	@param qCoupling: coupling factor
	@type qCoupling: float
	
	@param approximationTreshold: prefiltering of the solution. if modes seem to be missing, increase this parameter
	@type approximationTreshold:float
	@param numberOfRadialOrders: for how many radial orders above and below the Central Radial Mode should the frequencies be calculated.
	@type numberOfRadialOrders: float or integer	
	@return frequencies of the m=0 component of the mixed dipole modes
	@rtype float array
    """
    
    DeltaPg /= 10.**6.
    epsilon = centralRadial / largeSep - centralRadial // largeSep
    radialOrderRange = np.arange(-numberOfRadialOrders,numberOfRadialOrders+1.,dtype='float')
    dominantRadialOrderPressure = np.floor(centralRadial/largeSep)
    univPatternPerOrder = []; positionMixedModes = []
    rotationalSplitting = []
    logger.info('calculating dipole mixed modes for {0} radial ordres ...'.format(len(radialOrderRange)))


    for nnn in radialOrderRange:
        valueDipole =  purePressureMode(nnn,centralRadial, largeSep, 2., epsilon, 1.)
        freqRange =    np.linspace((centralRadial+largeSep*(nnn-0.25)),  (centralRadial+largeSep*(nnn+1.25)), 100000)
        equalityLine = freqRange/largeSep
        arcTerm = largeSep / scipy.pi * np.arctan( qCoupling * np.tan( scipy.pi / (DeltaPg*freqRange) ) )        
        #print 'nnn',nnn,', pure pressure dipole',valueDipole,min(freqRange),max(freqRange)
        
        # distance to equalityLine
        residualsToEquealityLine = np.array((freqRange/largeSep,(valueDipole+arcTerm)/largeSep-equalityLine)) # Frequenz Abstand der moeglichen Mixed Modes zur _equalityLine_

        # filtering for the close modes, speeds up the process
        filteredResidualsToEquealityLine = residualsToEquealityLine[:,(residualsToEquealityLine[1] < approximationTreshold/1.)]    # nur die, die nahe genug an der Mixed Mode liegen
        
        """ # plotting the method to solve the implicit formula
        pl.figure(num=2)
        pl.plot(freqRange/largeSep,equalityLine,'r', alpha=0.5)
        pl.plot(valueDipole/largeSep,valueDipole/largeSep,'mo')
        pl.plot(freqRange/largeSep,(valueDipole+arcTerm)/largeSep,'b', alpha=0.5)
        pl.plot(residualsToEquealityLine[0],residualsToEquealityLine[1],'k', alpha=0.5)
        pl.plot(filteredResidualsToEquealityLine[0],filteredResidualsToEquealityLine[1])
        #"""

        # find the closest mode to line of equility
        for jjj in np.arange(0,len(filteredResidualsToEquealityLine[0])-1):
            if filteredResidualsToEquealityLine[1,jjj] >= 0. and filteredResidualsToEquealityLine[1,jjj+1] <= 0.:
                    modePosition = scipy.average(filteredResidualsToEquealityLine[0,jjj:jjj+2])*largeSep
                    rotationalSplitting = 0.
                    #modePositionError = average(filteredResidualsToEquealityLine[1,jjj:jjj+2])*largeSep # deviation to 0.00 (i.e. error of assumption in muHz)
                    positionMixedModes.append(modePosition)

    dipoleMixedModes = np.array(positionMixedModes, dtype='float').T

    return dipoleMixedModes





def symmetricRotationalModulation(dfMax, largeSep, frequenciesUniversalPattern, mixedModesPattern, lambdaParam = 0.5, beParam = 0.08):
    """ Calculates and returns rotational splitting
    
    This function the modulation of the mixing between pressure- and gravity-dipole and derives the frequencies of the mixed modes from it.

	@param dfMax: the largest rotational splitting measured for a given star
	@type dfMax: float
	@param largeSep: large Frequency separation in powerspectrum, in muHz
	@type largeSep: float

	@param frequenciesUniversalPattern: output from rg.universalPattern
	@type frequenciesUniversalPattern: array, float

	@param mixedModesPattern: output from rg.asymptoticRelation
	@type mixedModesPattern: array, float

	@param lambdaParam: parameter for fitting, Default: lambdaParam = 0.5
	@type lambdaParam: float

	@param beParam:  parameter for fitting, Default: beParam = 0.08
	@type beParam: float
	
	@return array with the rotational splitting 
	@rtype array, float
    """


    
    radialModes =         frequenciesUniversalPattern[0]
    pureDipolePressure =    frequenciesUniversalPattern[1]
    mixedDipoleMode =        mixedModesPattern
    
    rotationalCentreMode = []; rotationalSplitting = []
    for nnn in np.arange(0,len(radialModes.T)-1):
        pureDipole =    pureDipolePressure[ (pureDipolePressure > radialModes[nnn]) & (pureDipolePressure < radialModes[nnn+1])]
        mixedDipoles=    mixedDipoleMode[ (mixedDipoleMode > radialModes[nnn]) & (mixedDipoleMode < radialModes[nnn+1])] 

        for mode in mixedDipoles:
            modulationProfile = 1.- lambdaParam / (1. + ( (mode - pureDipole[0]) / (beParam*largeSep) )**2. )
            rotationalCentreMode.append(mode) 
            rotationalSplitting.append(modulationProfile)# i.e., frequency of m=0 component, symmetric rotational splitting df. m=0 +/- df

    rotationalCentreMode = np.array(rotationalCentreMode)
    rotationalSplitting = np.array(rotationalSplitting).T*dfMax
    #print shape(rotationalCentreMode),shape(rotationalSplitting)
    rotationalSplitting = np.vstack((rotationalCentreMode,rotationalSplitting))  
    rotationalSplitting = np.array(rotationalSplitting, dtype='float')
    
    return rotationalSplitting
