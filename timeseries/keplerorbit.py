# -*- coding: utf-8 -*-
"""
Fit and evaluate binary (Keplerian) orbits.

Be careful with the inclination angle: if you only want an image projected onto
the sky, it is not case-sensitive. Only for the radial velocities this becomes
important (the primary and secondary can be mixed).

Example usage: compute the absolute orbits of the two components of
mu Cassiopeia, and compute the relative orbit of the second component.

See ps file on http://www.chara.gsu.edu/~gudehus/binary.html and Drummond et al.,
1995 (the latter's orbital parameters are used).

Neccesary imports:

>>> from ivs.units.constants import *

Orbital elements, Euler angles and masses of the system

>>> P = 21.753*365 # period in days
>>> e = 0.561 # ecc
>>> omega,Omega,i = 332.7, 47.3, 106.8 # degrees
>>> T0 = 1975.738*365 # year in days
>>> plx = 133.2 # mas
>>> M1,M2 = 0.719,0.168 # Msun
>>> times = np.linspace(T0,T0+0.95*P,100)

Derive semi-major axis of absolute orbits (semi-major axis of relative orbit is a+a2)

>>> fA = M2**3/ (M1+M2)**2 * Msol
>>> a1 = ((P*24*3600)**2 * GG * fA / (4*np.pi**2))**(1/3.)
>>> a1 = a1/au
>>> a2 = a1*M1/M2

Convert the angles to radians

>>> omega,Omega,i = omega/180*np.pi,Omega/180*np.pi,i/180*np.pi

Compute the orbits of the primary and secondary in the plane of the sky

>>> x1,y1,z1 = orbit_on_sky(times,(P,e,a1,T0,omega,Omega,i),component='prim',distance=(1000./133.,'pc'))
>>> x2,y2,z2 = orbit_on_sky(times,(P,e,a2,T0,omega,Omega,i),component='sec',distance=(1000./133.,'pc'))

Plot it:

>>> import pylab as pl
>>> colors = times/365.
>>> p = pl.title('mu Cassiopeia')
>>> p = pl.scatter(x1,y1,c=colors,edgecolors='none',cmap=pl.cm.spectral,label='_nolegend_')
>>> p = pl.scatter(x2,y2,c=colors,edgecolors='none',cmap=pl.cm.spectral,label='_nolegend_')
>>> p = pl.scatter(x2-x1,y2-y1,c=colors,edgecolors='none',cmap=pl.cm.spectral,label='_nolegend_')
>>> p = pl.plot(x1,y1,'b-',label='Primary')
>>> p = pl.plot(x2,y2,'g-',label='Secondary')
>>> p = pl.plot(x2-x1,y2-y1,'r-',label='Secondary relative')
>>> p = pl.colorbar()
>>> p = pl.legend(loc='lower right')
>>> p = pl.grid()
>>> p = pl.xlabel(r'Angular size (arcsec) --> N')
>>> p = pl.ylabel(r'Angular size (arcsec) --> E')

]include figure]]ivs_binary_keplerorbit_muCas.png]

Example usage 2: compute the absolute orbits of the two components of
Capella, and compute the relative orbit of the primary component.

See Hummel et al. 1994 (?).

Neccesary imports:

>>> from ivs.units.constants import *

Orbital elements, Euler angles and masses of the system

>>> P = 104.022 # period in days
>>> e = 0.000 #ecc
>>> i = 137.18 # degree
>>> T0 = 2447528.45 # year
>>> omega = 0 # degrees
>>> Omega = 40.8 # degrees
>>> plx = 76.20 # mas
>>> M1 = 2.69 #Msun
>>> M2 = 2.56 #Msun
>>> RV0 = 0.
>>> times = np.linspace(T0,T0+0.95*P,100)

Derive semi-major axis of absolute orbits (semi-major axis of relative orbit is a+a2)

>>> fA = M2**3/ (M1+M2)**2 * Msol
>>> a1 = ((P*24*3600)**2 * GG * fA / (4*np.pi**2))**(1/3.)
>>> a1 = a1/au
>>> a2 = a1*M1/M2

Convert the angles to radians

>>> omega,Omega,i = omega/180*np.pi,Omega/180*np.pi,i/180*np.pi

Compute the orbits of the primary and secondary in the plane of the sky

>>> x1,y1,z1 = orbit_on_sky(times,(P,e,a1,T0,omega,Omega,i),component='prim',distance=(1000./plx,'pc'))
>>> x2,y2,z2 = orbit_on_sky(times,(P,e,a2,T0,omega,Omega,i),component='sec',distance=(1000./plx,'pc'))

Plot it:

>>> import pylab as pl
>>> colors = times/365.
>>> p = pl.figure()
>>> p = pl.title('Capella')
>>> p = pl.scatter(x1,y1,c=colors,edgecolors='none',cmap=pl.cm.spectral,label='_nolegend_')
>>> p = pl.scatter(x2,y2,c=colors,edgecolors='none',cmap=pl.cm.spectral,label='_nolegend_')
>>> p = pl.scatter(x1-x2,y1-y2,c=colors,edgecolors='none',cmap=pl.cm.spectral,label='_nolegend_')
>>> p = pl.plot(x1,y1,'b-',label='Primary')
>>> p = pl.plot(x2,y2,'g-',label='Secondary')
>>> p = pl.plot(x1-x2,y1-y2,'r-',label='Secondary relative')
>>> p = pl.colorbar()
>>> p = pl.legend(loc='lower right')
>>> p = pl.grid()
>>> p = pl.xlabel(r'Angular size DEC (arcsec) --> N')
>>> p = pl.ylabel(r'Angular size RA (arcsec) --> E')

]include figure]]ivs_binary_keplerorbit_capella.png]

"""
import logging
import numpy as np
from scipy import optimize
from ivs.units import conversions
from ivs.units.constants import *
from ivs.coordinates import vectors

logger = logging.getLogger('IVS.KEPLER')

#{ Radial velocities

def radial_velocity(parameters,times=None,theta=None,itermax=8):
    """
    Evaluate Radial velocities due to kepler orbit.
    
    These parameters define the Keplerian orbit if you give times points (C{times}, days):
        1. Period of the system (days)
        2. time of periastron passage T0 (not x0!) (HJD)
        3. eccentricity
        4. omega, the longitude of periastron gives the orientation
           of the orbit within its own plane (radians)
        5. the semiamplitude of the velocity curve (km/s)
        6. systemic velocity RV0 (RV of centre of mass of system) (km/s)
    
    These parameters define the Keplerian orbit if you give angles (C{theta}, radians):
        1. Period of the system (days)
        2. eccentricity
        3. a, the semi-major axis of the absolute orbit of the component (au)
        4. omega, the longitude of periastron gives the orientation of the
        orbit within in its own plane (radians)
        5. inclination of the orbit (radians)
        6. systemic velocity RV0 (RV of centre of mass of system) (km/s)
    
    The periastron passage T0 can be derived via x0 by calculating
    
    T0 = x0/(2pi*Freq) + times[0]
    
    See e.g. p 41,42 of Hilditch, 'An Introduction To Close Binary Stars'
    
    @parameter parameters: parameters of Keplerian orbit (dependent on input)
    @type parameters: iterable
    @parameter times: observation times (days)
    @type times: numpy array
    @parameter theta: orbital angles (radians)
    @type theta: numpy array
    @param itermax: number of iterations to find true anomaly
    @type itermax: integer
    @return: fitted radial velocities (km/s)
    @rtype: ndarray
    """
    
    #-- in the case time points are given:
    if times is not None:
        #-- expand parameters and calculate x0
        P,T0,e,omega,K,RV0 = parameters
        freq = 1./P
        x0 = T0*2*np.pi*freq
        #-- calculate true anomaly
        E,true_an = true_anomaly(times*2*np.pi*freq-x0,e,itermax=itermax)
        #-- evaluate Keplerian radial velocity orbit
        RVfit = RV0 + K*(e*np.cos(omega) + np.cos(true_an+omega))
    
    elif theta is not None:
        P,e,a,omega,i,RV0 = parameters
        P = conversions.convert('d','s',P)
        a = conversions.convert('au','m',a)
        K = 2*np.pi*a*np.sin(i)/ (P*np.sqrt(1.-e**2))
        RVfit = RV0 + K*(e*np.cos(omega) + np.cos(theta+omega))/1000.
        
    return RVfit

def orbit_in_plane(times,parameters,component='primary',coordinate_frame='polar'):
    """
    Construct an orbit in the orbital plane.
    
    Give times in days
    
    Parameters contains:
        1. period (days)
        2. eccentricity
        3. semi-major axis (au)
        4. time of periastron passage T0 (not x0!) (HJD)
    
    Return r (m) and theta (radians)
    
    @param times: times of observations (days)
    @type times: array
    @param parameters: list of parameters (P,e,a,T0)
    @type parameters: list
    @param component: component to calculate the orbit of. If it's the secondary,
    the angles will be shifted by 180 degrees
    @type component: string, one of 'primary' or 'secondary'
    @param coordinate_frame: type of coordinates
    @type coordinate_frame: str, ('polar' or 'cartesian')
    @return: coord1,coord2
    @rtype: 2xarray
    """
    P,e,a,T0 = parameters
    P = conversions.convert('d','s',P)
    a = conversions.convert('au','m',a)
    T0 = conversions.convert('d','s',T0)
    times = conversions.convert('d','s',times)
    
    n = 2*np.pi/P
    ma = n*(times-T0)
    E,theta = true_anomaly(ma,e)
    r = a*(1-e*np.cos(E))
    PR = r*np.sin(theta)
    #PR[E>0] *= -1
    #theta[E>0] *= -1
    
    #-- correct angles if secondary component is calculated
    if 'sec' in component.lower():
        theta += np.pi
    
    if coordinate_frame=='polar':
        return r,theta
    elif coordinate_frame=='cartesian':
        return r*np.cos(theta),r*np.sin(theta)

def velocity_in_plane(times,parameters,component='primary',coordinate_frame='polar'):
    """
    Calculate the velocity in the orbital plane.
    
    @param times: times of observations (days)
    @type times: array
    @param parameters: list of parameters (P,e,a,T0)
    @type parameters: list
    @param component: component to calculate the orbit of. If it's the secondary,
    the angles will be shifted by 180 degrees
    @type component: string, one of 'primary' or 'secondary'
    @param coordinate_frame: type of coordinates
    @type coordinate_frame: str, ('polar' or 'cartesian')
    @return: coord1,coord2
    @rtype: 2xarray
    """
    #-- calculate the orbit in the plane
    r,theta = orbit_in_plane(times,parameters,component=component,coordinate_frame='polar')
    P,e,a,T0 = parameters
    P = conversions.convert('d','s',P)
    a = conversions.convert('au','m',a)
    
    #-- compute rdot and thetadot
    l = r*(1+e*np.cos(theta))
    L = 2*np.pi*a**2/P*np.sqrt(1-e**2)
    rdot = L/l*e*np.sin(theta)
    thetadot = L/r**2
    
    #-- convert to the right coordinate frame
    if coordinate_frame=='polar':
        return rdot,thetadot
    elif coordinate_frame=='cartesian':
        vx,vy,vz = vectors.spher2cart((r,theta,np.pi/2.),(rdot,r*thetadot,0))
        return vx,vy


def project_orbit(r,theta,parameters):
    """
    Project an orbit onto the plane of the sky.
    
    Parameters contains the Euler angles:
        1. omega: the longitude of periastron gives the orientation of the
        orbit within in its own plane (radians)
        2. Omega: PA of ascending node (radians)
        3. i: inclination (radians), i=pi/2 is edge on
    
    Returns x,y (orbit in plane of the sky) and z
    
    See Hilditch p41 for a sketch of the coordinates. The difference with this
    is that all angles are inverted. In this approach, North is in the positive
    X direction, East is in the negative Y direction.
    """
    omega,Omega,i = parameters
    x = r*(np.cos(Omega)*np.cos(theta+omega) - np.sin(+Omega)*np.sin(theta+omega)*np.cos(i))
    y = r*(np.sin(Omega)*np.cos(theta+omega) + np.cos(+Omega)*np.sin(theta+omega)*np.cos(i))
    z = r*(np.sin(theta+omega)*np.sin(i))
    return x,y,z
    
    

def orbit_on_sky(times,parameters,distance=None,component='primary'):
    """
    Construct an orbit projected on the sky.
    
    Parameters contains:
        1. period (days)
        2. eccentricity
        3. semi-major axis (au)
        4. time of periastron passage T0 (not x0!) (HJD)
        5. omega: the longitude of periastron gives the orientation of the
        orbit within in its own plane (radians)
        6. Omega: PA of ascending node (radians)
        7. i: inclination (radians), i=pi/2 is edge on
    
    You can give an extra parameter 'distance' as a tuple (value,'unit'). This
    will be used to convert the distances to angular scale (arcsec).
    
    Else, this function returns the distances in AU.
    
    See Hilditch p41 for a sketch of the coordinates. The difference with this
    is that all angles are inverted. In this approach, North is in the positive
    X direction, East is in the negative Y direction.
    """
    #-- divide parameters in 'in-plane' parameters and euler angles
    pars_in_plane,euler_angles = parameters[:4],parameters[4:]
    #-- construct the orbit in the orbital plane
    r,theta = orbit_in_plane(times,pars_in_plane,component=component)
    #-- and project in onto the sky according to the euler angles
    x,y,z = project_orbit(r,theta,euler_angles)
    
    #-- if necessary, convert the true distance to angular scale
    if distance is not None:
        d = conversions.convert(distance[1],'m',distance[0])
        x = 2*np.arctan(x/(2*d))/np.pi*180*3600
        y = 2*np.arctan(y/(2*d))/np.pi*180*3600
        z = 2*np.arctan(z/(2*d))/np.pi*180*3600
        return x,y,z
    else:
        return x/au,y/au,z/au
    
    


def true_anomaly(M,e,itermax=8):
    """
    Calculation of true and eccentric anomaly in Kepler orbits.
    
    M is the phase of the star, e is the eccentricity
    
    See p.39 of Hilditch, 'An Introduction To Close Binary Stars'
    
    @parameter M: phase
    @type M: float
    @parameter e: eccentricity
    @type e: float
    @keyword itermax: maximum number of iterations
    @type itermax: integer
    @return: eccentric anomaly (E), true anomaly (theta)
    @rtype: float,float
    """
    #-- initial value
    Fn = M + e*np.sin(M) + e**2/2.*np.sin(2*M)
    #-- iterative solving of the transcendent Kepler's equation
    for i in range(itermax):
        F = Fn
        Mn = F-e*np.sin(F)
        Fn = F+(M-Mn)/(1.-e*np.cos(F))
        keep = F!=0 #-- take care of zerodivision
        if hasattr(F,'__iter__'):
            if np.all(abs((Fn-F)[keep]/F[keep])<0.00001):
                break
        elif (abs((Fn-F)/F)<0.00001):
            break
    #-- relationship between true anomaly (theta) and eccentric
    #   anomalie (Fn)
    true_an = 2.*np.arctan(np.sqrt((1.+e)/(1.-e))*np.tan(Fn/2.))
    return Fn,true_an




def calculate_phase(T,e,omega,pshift=0):
    """
    Compute orbital phase from true anomaly T
    
    @parameter T: true anomaly
    @type T: float
    @parameter omega: argument of periastron (radians)
    @type omega: float
    @parameter e: eccentricity
    @type e: float
    @parameter pshift: phase shift
    @type pshift: float
    @return: phase of superior conjunction, phase of periastron passage
    @rtype: float,float
    """
    E = 2.0*np.arctan(np.sqrt((1-e)/(1+e)) * np.tan(T/2.0))
    M = E - e*np.sin(E)
    return (M+omega)/(2.0*np.pi) - 0.25 + pshift

    
    
def calculate_critical_phases(omega,e,pshift=0):
    """
    Compute phase of superior conjunction and periastron passage.
    
    Example usage:
    >>> omega = np.pi/4.0
    >>> e = 0.3
    >>> print calculate_critical_phases(omega,e)
    (-0.125, -0.057644612788576133, -0.42054512757020118, -0.19235538721142384, 0.17054512757020118)
    
    @parameter omega: argument of periastron (radians)
    @type omega: float
    @parameter e: eccentricity
    @type e: float
    @parameter pshift: phase shift
    @type pshift: float
    @return: phase of superior conjunction, phase of periastron passage
    @rtype: float,float
    """
    #-- Phase of periastron passage
    Phi_omega = (omega - np.pi/2.0)/(2.0*np.pi) + pshift
    #-- Phase of inferior/superior conjunction
    Phi_conj  = calculate_phase(np.pi/2.0-omega,e,omega,pshift)
    Phi_inf   = calculate_phase(3.0*np.pi/2.0-omega,e,omega,pshift)
    Phi_asc   = calculate_phase(-omega,e,omega,pshift)
    Phi_desc  = calculate_phase(np.pi-omega,e,omega,pshift)
    return Phi_omega,Phi_conj,Phi_inf,Phi_asc,Phi_desc


def eclipse_separation(e,omega):
    """
    Calculate the eclipse separation between primary and secondary in a light curve.
    
    Minimum separation at omega=pi
    Maximum spearation at omega=0
    
    @parameter e: eccentricity
    @type e: float
    @parameter omega: argument of periastron (radians)
    @type omega: float
    @return: separation in phase units (0.5 is half)
    @rtype: float
    """
    radians = np.pi+2*np.arctan(e*np.cos(omega)/np.sqrt(1-e**2)) + 2*e*np.cos(omega)*np.sqrt(1-e**2)/(1-e**2*np.sin(omega)**2)
    return radians/(2*np.pi)


def omega_from_eclipse_separation(separation,e):
    """
    Caculate the argument of periastron from the eclipse separation and eccentricity.
    
    separation in phase units.
    
    @parameter separation: separation in phase units (0.5 is half)
    @type separation: float
    @parameter e: eccentricity
    @type e: float
    @return: omega (longitude of periastron) in radians
    @rtype: float
    """
    minsep_omega = np.pi
    maxsep_omega = 0.
    minsep = eclipse_separation(e,minsep_omega)
    maxsep = eclipse_separation(e,maxsep_omega)
    if separation<minsep or maxsep<separation:
        logger.warning('Phase separation must be between %.3g and %.3g when e=%.3g'%(minsep,maxsep,e))
        return np.nan
    else:
        omega = optimize.bisect(lambda x:separation-eclipse_separation(e,x),maxsep_omega,minsep_omega)
        return omega
    
#}

#{ Kepler's laws

def third_law(M=None,a=None,P=None):
    """
    Kepler's third law.
    
    Give two quantities, derived the third.
    
    M = total mass system (solar units)
    a = semi-major axis (au)
    P = period (d)
    
    >>> print third_law(M=1.,a=1.)
    365.256891359 
    >>> print third_law(a=1.,P=365.25)
    1.00003773538
    >>> print third_law(M=1.,P=365.25)
    0.999987421856
    """
    if a is not None:
        a *= au
    if P is not None:
        P *= (24*3600.)
    if M is not None:
        M *= Msol
    
    if M is None:
        return 4*np.pi**2*a**3/P**2/GG/Msol
    if a is None:
        return (GG*M*P**2/(4*np.pi**2))**(1./3.)/au
    if P is None:
        return np.sqrt(4*np.pi**2*a**3/(GG*M))/(24*3600.)
    

if __name__=="__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()
    