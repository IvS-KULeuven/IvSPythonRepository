# -*- coding: utf-8 -*-
"""
Calculate the distance to a star with Lutz-Kelker bias
"""
import numpy as np
from numpy import exp,sqrt,sin,pi
import logging

logger = logging.getLogger("IVS.DIST")


def rho(z,z_sun=20.0,hd=31.8,hh=490.,sigma=1.91e-3,f=0.039):
    """
    Stellar galactic density function.
    
    Galactic coordinate z: self-gravitating isothermal disk plus a Gaussian halo
    
    See, e.g., Maiz-Apellaniz, Alfaro and Sota 2007/2008 (Poster)
    
    Other values we found in the literature:
    
    z_sun,hd,sigma,f = 24.7,34.2,1.62e-3,0.058
    
    @param z: galactic coordinate z (parsec)
    @type z: array or float
    @param z_sun: Sun's distance above the Galactic plane (20.0 +/- 2.9 pc)
    @type z_sun: float
    @param hd: disk scale height (31.8+/-1.6 pc)
    @type hd: float
    @param hh: halo half width (490+/-170 pc)
    @type hh: float
    @param sigma: number of stars per square parsec (total surface number density)
    (1.91e-3+/-0.11e-3)
    @type sigma: float
    @param f: fraction of stars in the halo population (0.039+/-0.015)
    @type f: float
    @return: halo+disk,halo and disk density function
    @rtype: 3-tuple
    """
    disk = sigma* (1-f)/(4*hd) * np.cosh(0.5*(z+z_sun)/hd)**-2
    halo = sigma*f/(2*sqrt(2)*hh) * exp( -(0.5*(z+z_sun)/hh)**2)
    return halo+disk,halo,disk

def probability_cd(r,plx,e_plx):
    """
    Compute the probability for an object to be at distance r (pc), given
    its parallax (mas) and error on the parallax (mas) and a constant density
    function.
    
    Unnormalised!
    
    To obtain the probabilty, multiply with a stellar galactic density function.
    
    @param r: distance (pc)
    @type r: float/array
    @param plx: parallax (mas)
    @type plx: float/array
    @param e_plx: error on parallax (mas)
    @type e_plx: float/array
    @return: probability function
    """
    plx /= 1000.
    e_plx /= 1000.
    A = 1.
    prob1 = A*r**2*exp(-0.5*((1-r*plx)/(r*e_plx))**2)
    return prob1

def distprob(r,theta,plx,**kwargs):
    """
    Compute the probability for an object to be located at a distance r (pc),
    given its parallax and galactic lattitude.
    
    theta in degrees
    plx is a tuple!
    
    returns (unnormalised) probability density function.
    """
    z = r*sin(theta/180.*pi)
    pcd = probability_cd(r,plx[0],plx[1])
    galaxy,halo,disk = rho(z,**kwargs)
    dpr = pcd*galaxy
    logger.info("Calculate Lutz-Kelker bias for THETA=%.3fdeg"%(theta))
    return dpr