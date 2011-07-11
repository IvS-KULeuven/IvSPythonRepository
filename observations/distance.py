# -*- coding: utf-8 -*-
"""
Calculate the distance to a star with Lutz-Kelker bias

Section 1. Bias introduced by the shape of the Galaxy
=====================================================

In this section, we compare the distance calculation via the inverse of the
parallax with a method that takes into account the distribution of stars in the
Milky Way. The latter can be important when the error on the parallax is larger
than 25%: then, it becomes more likely to have observed a star located further
in the galaxy but which appears erronously closer, just because there are many
more stars in the volume Delta_d (as the volume is much larger). This depends
on the the location of the star in the Milky Way (in the Halo, the number of
stars is low anyway) and thus also on the shape of the Milky Way. Here, the
Galaxy is approximated as the sum of a self-gravitating isothermal disk and a
Gaussian Halo (L{rho}). For more details, see, e.g., Maiz-Apellaniz, Alfaro and
Sota 2007/2008 (Poster).

To correctly compute the distance to a star given it's parallax, we need to
know it's position wrt the galactic plane. So, we need information on the
coordinates of the star:

>>> from ivs.catalogs import sesame

We take three examples: Vega, which has a well known parallax, HIP14, which has a
mediocre value, and HIP15, which has a badly determined value. When retrieving
the information from SIMBAD, make sure we have Van Leeuwen's (2007) parallax,
and that we have the galactic coordinates of the target.

>>> vega = sesame.search('vega',fix=True)
>>> hip14= sesame.search('hip14',fix=True)
>>> hip15= sesame.search('hip15',fix=True)

We have the following parameters::

    vega : 130.23 +/- 0.36 (0.28%)
    HIP14: 4.86 +/- 0.67 (13.79%)
    HIP15: 1.91 +/- 1.14 (59.69%)

We can compute the distance probability functions via L{distprob}:

>>> r1 = np.linspace(7,8,1000)
>>> r2 = np.linspace(0,500,1000)
>>> r3 = np.linspace(0,1.5e4,1000)
>>> prob1 = distprob(r1,vega['galpos'][1],(vega['plx']['v'],vega['plx']['e']))
>>> prob2 = distprob(r2,hip14['galpos'][1],(hip14['plx']['v'],hip14['plx']['e']))
>>> prob3 = distprob(r3,hip15['galpos'][1],(hip15['plx']['v'],hip15['plx']['e']))

It is useful to compare with the Gaussian approximation. First we need some
extra modules to compute with uncertainties, and to evaluate the Gaussian
function.

>>> from ivs.sigproc import evaluate
>>> from ivs.units.uncertainties import ufloat

Then we can approximate the distance to the star (in pc) as the inverse of the
parallax (in arcsec).

>>> vega_ = 1./ufloat((vega['plx']['v']/1000.,vega['plx']['e']/1000.))
>>> hip14_ = 1./ufloat((hip14['plx']['v']/1000.,hip14['plx']['e']/1000.))
>>> hip15_ = 1./ufloat((hip15['plx']['v']/1000.,hip15['plx']['e']/1000.))
>>> prob1_ = evaluate.gauss(r1,[1.,vega_.nominal_value,vega_.std_dev()])
>>> prob2_ = evaluate.gauss(r2,[1.,hip14_.nominal_value,hip14_.std_dev()])
>>> prob3_ = evaluate.gauss(r3,[1.,hip15_.nominal_value,hip15_.std_dev()])

We summarize everything in a plot: up to a relative error of ~25%, the Gaussian
approximation works pretty well. From then on, even the peak of the distribution
shows a bias (the Lutz-Kelker bias), but also have tails can appear.

>>> p = pl.figure()
>>> p = pl.subplot(221)
>>> p = pl.plot(r1,prob1/prob1.max(),'k-',lw=2,label='Vega (LK)')
>>> p = pl.plot(r1,prob1_,'r--',lw=2,label='Vega (GA)')
>>> p = pl.legend()
>>> p = pl.xlabel('Distance (pc)');p = pl.ylabel('Unnormalised probability')
>>> p = pl.xlim(7.55,7.8)
>>> p = pl.subplot(222)
>>> p = pl.plot(r2,prob2/prob2.max(),'k-',lw=2,label='HIP14 (LK)')
>>> p = pl.plot(r2,prob2_,'r--',lw=2,label='HIP14 (GA)')
>>> p = pl.xlabel('Distance (pc)');p = pl.ylabel('Unnormalised probability')
>>> p = pl.legend()
>>> p = pl.subplot(212)
>>> p = pl.plot(r3,prob3/prob3.max(),'k-',lw=2,label='HIP15 (LK)')
>>> p = pl.plot(r3,prob3_,'r--',lw=2,label='HIP15 (GA)')
>>> p = pl.xlabel('Distance (pc)');p = pl.ylabel('Unnormalised probability')
>>> p = pl.legend()

]include figure]]ivs_observations_distance_LK.png]

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

def probability_cd(r,plx):
    """
    Compute the probability for an object to be at distance r (pc), given
    its parallax (mas) and error on the parallax (mas) and a constant density
    function.
    
    Unnormalised!
    
    To obtain the probabilty, multiply with a stellar galactic density function.
    
    @param r: distance (pc)
    @type r: float/array
    @param plx: parallax (mas) and error
    @type plx: tuple of float/array
    @return: probability function
    @rtype: float/array
    """
    plx,e_plx = plx[0]/1000,plx[1]/1000.
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
    
    @param r: distance (pc)
    @type r: float/array
    @param theta: galactic latitude (deg)
    @type theta: float
    @param plx: parallax (mas) and error
    @type plx: tuple of float/array
    @return: probability function
    @rtype: float/array
    """
    z = r*sin(theta/180.*pi)
    pcd = probability_cd(r,plx)
    galaxy,halo,disk = rho(z,**kwargs)
    dpr = pcd*galaxy
    logger.info("Calculate Lutz-Kelker bias for THETA=%.3fdeg"%(theta))
    return dpr


if __name__ == "__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()