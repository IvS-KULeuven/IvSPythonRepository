# -*- coding: utf-8 -*-
"""
Fit and evaluate binary (Keplerian) orbits.

Example usage: compute the absolute orbits of the two components of
mu Cassiopeia, and compute the relative orbit of the second component.

See ps file on http://www.chara.gsu.edu/~gudehus/binary.html and Drummond et al.,
1995 (the latters orbital parameters are used).

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
>>> p = pl.show()


"""
import numpy as np
from ivs.units import conversions
from ivs.units.constants import *

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

def orbit_in_plane(times,parameters,component='primary'):
    """
    Construct an orbit in the orbital plane.
    
    Parameters contains:
        1. period (days)
        2. eccentricity
        3. semi-major axis (au)
        4. time of periastron passage T0 (not x0!) (HJD)
    
    Return r (m) and theta (radians)
    
    @param component: component to calculate the orbit of. If it's the secondary,
    the angles will be shifted by 180 degrees
    @type component: string, one of 'primary' or 'secondary'
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
    
    return r,theta




def project_orbit(r,theta,parameters):
    """
    Project an orbit on the plane of the sky.
    
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
        if all(abs((Fn-F)/F)<0.00001):
            break
    #-- relationship between true anomaly (theta) and eccentric
    #   anomalie (Fn)
    true_an = 2.*np.arctan(np.sqrt((1.+e)/(1.-e))*np.tan(Fn/2.))
    return Fn,true_an

#}

if __name__=="__main__":
    import doctest
    doctest.testmod()