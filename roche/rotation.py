"""
Roche models of fast and differential rotation

Section 1. Non-rotating spherically symmetric single star
=========================================================

The most trivial usage of this module would be to construct a non-rotating,
spherically symmetric single star. All quantities should be constant over the
stellar surface, except of course the projected ones. We compute this via the
zero-rotation limit of the fast rotating model:

First set some parameters: the rotation rate, effective temperature at the pole,
polar radius, mass and viewing angle.

>>> omega = 0.00       # ratio to critical velocity
>>> T_pole = 5500.     # K
>>> r_pole = 1.        # solar radii
>>> M = 1.             # solar mass
>>> view_angle = pi/2  # radians

Construct a coordinate grid with a resolution of 50 points, covering the whole
star. Unravel, so that we have access to 1D coordinate arrays

>>> theta,phi = local.get_grid(50,full=True)
>>> thetas,phis = np.ravel(theta),np.ravel(phi)

Calculate the shape of the stellar surface and the local gravity for every point
in the grid. This is a 3D vector in Cartesian coordinates, so reshape all
components to match the grid shape. As a reference, also explicitly calculate
the polar surface gravity, which is the z-component of the gravity vector.

>>> radius = np.array([get_fastrot_roche_radius(itheta,r_pole,omega) for itheta in thetas]).reshape(theta.shape)
>>> grav_local = np.array([fastrot_roche_surface_gravity(iradius,itheta,iphi,r_pole,omega,M) for iradius,itheta,iphi in zip(radius.ravel(),thetas,phis)]).T
>>> grav_local = np.array([i.reshape(theta.shape) for i in grav_local])
>>> g_pole = fastrot_roche_surface_gravity(r_pole,0,0,r_pole,omega,M)[-1]

Calculate the value of the surface gravity in SI-units:

>>> grav = vectors.norm(grav_local)

Compute the size of all surface elements, and the angle with respect to the
radius vector (should be all zero). The normal to the surface is the negative of
the surface gravity vector, which is directed inwards:

>>> areas_local,cos_gamma = local.surface_elements((radius,(theta,phi)),-grav_local)

Now we can calculate the local effective temperature (assuming radiative
atmosphere in Von Zeipel law), and bolometric intensities:

>>> teff_local = local.temperature(vectors.norm(grav_local),g_pole,T_pole,beta=1.)
>>> ints_local = local.intensity(teff_local,grav,photband='OPEN.BOL')

Define the line-of-sight vector:

>>> line_of_sight = np.zeros_like(grav_local)
>>> line_of_sight[0] = -sin(view_angle)
>>> line_of_sight[2] = -cos(view_angle)

Compute the angles between surface normals and line-of-sight. We do this in 1D:

>>> gravity = np.array([igrav.ravel() for igrav in grav_local])
>>> line_of_sight = np.array([ilos.ravel() for ilos in line_of_sight])
>>> angles = vectors.angle(-gravity,line_of_sight)
>>> mus = cos(angles)

Now compute the bolometric projected intensity, taking care of limbdarkening
effects:

>>> intens_proj = local.intensity(teff_local.ravel(),grav.ravel(),mus,photband='OPEN.BOL').reshape(theta.shape)

The total and projected luminosity can then be calculated the following way:

>>> print pi*(ints_local*areas_local*constants.Rsol_cgs**2).sum()/constants.Lsol_cgs
0.992471247895
>>> print pi*np.nansum(intens_proj*areas_local*constants.Rsol_cgs**2)/constants.Lsol_cgs
0.360380373413

Now make some plots showing the local quantities:

>>> quantities = areas_local,np.log10(grav*100),teff_local,ints_local,intens_proj,angles/pi*180,radius
>>> names = 'Area','log g', 'Teff', 'Flux', 'Proj. flux', 'Angle'
>>> p = pl.figure()
>>> rows,cols = 2,3    
>>> for i,(quantity,name) in enumerate(zip(quantities,names)):
...    p = pl.subplot(rows,cols,i+1)
...    p = pl.title(name)
...    q = quantity.ravel()
...    vmin,vmax = q[-np.isnan(q)].min(),q[-np.isnan(q)].max()
...    if vmin==vmax: vmin,vmax = q[-np.isnan(q)].mean()-0.01*q[-np.isnan(q)].mean(),q[-np.isnan(q)].mean()+0.01*q[-np.isnan(q)].mean()
...    p = pl.scatter(phis/pi*180,thetas/pi*180,c=q,edgecolors='none',vmin=vmin,vmax=vmax)
...    p = pl.colorbar()
...    p = pl.xlim(0,360)
...    p = pl.ylim(180,0)

]include figure]]ivs_binary_rochepotential_nonrotstar.png]

Section 2. Fast and uniformly rotating star
===========================================

Repeating the same example as above, but now with parameters

>>> omega = 0.98

Gives the figure below:

]include figure]]ivs_binary_rochepotential_fastrotstar.png]

Changing the viewing angle to 80 degrees

>>> view_angle = 80/180.*pi

gives: Note that the maximum projected intensity is higher than in the previous
example: this is because we view directly at higher latitudes now, which are
hotter than the equator. The least flux is coming from the equatorial zone.

]include figure]]ivs_binary_rochepotential_fastrotstar_80incl.png]

To generate an image of the star, it is wise to convert the coordinates to
Cartesian coordinates:

>>> x,y,z = vectors.spher2cart_coord(radius.ravel(),phis,thetas)
>>> p = pl.figure()
>>> p = pl.subplot(121,aspect='equal')
>>> p = pl.title('top view')
>>> p = pl.scatter(x,y,c=teff_local.ravel(),edgecolors='none')
>>> p = pl.subplot(122,aspect='equal')
>>> p = pl.title('side view')
>>> p = pl.scatter(y,z,c=teff_local.ravel(),edgecolors='none')

The movie below is obtained by calculating the projected intensity in the line
of sight

>>> proj_intens,mus = local.projected_intensity(teff_local,grav_local,areas_local,np.array([0,-sin(view_angle),-cos(view_angle)]),photband='OPEN.BOL')

for different lines of sight, and then for each line of sight computing the
Cartesian coordinates, and finally rotating in Y-Z plane.

>>> y,z = vectors.rotate(y,z,-(pi/2-view_angle))

]include figure]]ivs_binary_rochepotential_fastrotstar_shape.png]

]include figure]]ivs_binary_rochepotential_fastrot.gif]

Section 3. Differentially rotating star
=======================================

The differential rotation rate has to follow the stereotypical law

Omega(theta) = Omega_pole + b * sin(theta)

The critical velocity can now easily be  exceeded at the equator, if the rotation
rate increases towards the pole.

First set some parameters: the rotation rate of the equator and pole, effective
temperature at the pole, polar radius, mass and viewing angle.

>>> omega_eq = 0.1       # ratio to critical rotation rate (equatorial)
>>> omega_pl = 0.9       # ratio to critical rotation rate (polar)
>>> T_pole = 5500.       # K
>>> r_pole = 1.0         # solar radii
>>> M = 1.               # solar mass
>>> view_angle = pi/2    # radians

We now do very similar stuff as in Section 1, except for the different Roche
potential. (We can skip making the grid now)

>>> radius = np.array([get_diffrot_roche_radius(itheta,r_pole,M,omega_eq,omega_pl) for itheta in thetas]).reshape(theta.shape)
>>> grav_local = np.array([diffrot_roche_surface_gravity(iradius,itheta,iphi,r_pole,M,omega_eq,omega_pl) for iradius,itheta,iphi in zip(radius.ravel(),thetas,phis)]).T
>>> grav_local = np.array([i.reshape(theta.shape) for i in grav_local])
>>> g_pole = diffrot_roche_surface_gravity(r_pole,0,0,r_pole,M,omega_eq,omega_pl)[-1]

Calculate the value of the surface gravity in SI-units:

>>> grav = vectors.norm(grav_local)

Compute the size of all surface elements, and the angle with respect to the
radius vector (should be all zero). The normal to the surface is the negative of
the surface gravity vector, which is directed inwards:

>>> areas_local,cos_gamma = local.surface_elements((radius,(theta,phi)),-grav_local)

Now we can calculate the local effective temperature (assuming radiative
atmosphere in Von Zeipel law), and bolometric intensities:

>>> teff_local = local.temperature(vectors.norm(grav_local),g_pole,T_pole,beta=1.)
>>> ints_local = local.intensity(teff_local,grav,photband='OPEN.BOL')

Compute the angles between surface normals and line-of-sight. We do this in 1D:

>>> gravity = np.array([igrav.ravel() for igrav in grav_local])
>>> line_of_sight = np.array([ilos.ravel() for ilos in line_of_sight])
>>> angles = vectors.angle(-gravity,line_of_sight)
>>> mus = cos(angles)

Now compute the bolometric projected intensity, taking care of limbdarkening
effects:

>>> intens_proj = local.intensity(teff_local.ravel(),grav.ravel(),mus,photband='OPEN.BOL').reshape(theta.shape)

Now make some plots showing the local quantities:

>>> quantities = areas_local,np.log10(grav*100),teff_local,ints_local,intens_proj,angles/pi*180,radius
>>> names = 'Area','log g', 'Teff', 'Flux', 'Proj. flux', 'Angle'
>>> p = pl.figure()
>>> rows,cols = 2,3    
>>> for i,(quantity,name) in enumerate(zip(quantities,names)):
...    p = pl.subplot(rows,cols,i+1)
...    p = pl.title(name)
...    q = quantity.ravel()
...    vmin,vmax = q[-np.isnan(q)].min(),q[-np.isnan(q)].max()
...    if vmin==vmax: vmin,vmax = q[-np.isnan(q)].mean()-0.01*q[-np.isnan(q)].mean(),q[-np.isnan(q)].mean()+0.01*q[-np.isnan(q)].mean()
...    p = pl.scatter(phis/pi*180,thetas/pi*180,c=q,edgecolors='none',vmin=vmin,vmax=vmax)
...    p = pl.colorbar()
...    p = pl.xlim(0,360)
...    p = pl.ylim(180,0)

]include figure]]ivs_binary_rochepotential_diffrotstar.png]

To generate an image of the star, it is wise to convert the coordinates to
Cartesian coordinates:

>>> x,y,z = vectors.spher2cart_coord(radius.ravel(),phis,thetas)
>>> p = pl.figure()
>>> p = pl.subplot(121,aspect='equal')
>>> p = pl.title('top view')
>>> p = pl.scatter(x,y,c=teff_local.ravel(),edgecolors='none')
>>> p = pl.subplot(122,aspect='equal')
>>> p = pl.title('side view')
>>> p = pl.scatter(y,z,c=teff_local.ravel(),edgecolors='none')

]include figure]]ivs_binary_rochepotential_diffrotstar_shape.png]
"""
import numpy as np
from numpy import pi,cos,sin,sqrt,nan
from scipy.optimize import newton

from ivs.roche import local
from ivs.coordinates import vectors
from ivs.units import constants
from ivs.units import conversions


#{ Fast rotation Roche potential

def fastrot_roche_surface_gravity(r,theta,phi,r_pole,omega,M,norm=False):
    """
    Calculate components of the local surface gravity of the fast rotating
    Roche model.
    
    Input units are solar units.
    Output units are SI.
    Omega is fraction of critical velocity.
    
    See Cranmer & Owocki, Apj (1995)
    
    @param r: radius of the surface element to calculate the surface gravity
    @type r: float/ndarray
    @param theta: colatitude of surface element
    @type theta: float/ndarray
    @param phi: longitude of surface element
    @type phi: float/ndarray
    @param r_pole: polar radius of model
    @type r_pole: float
    @param omega: fraction of critical rotation velocity
    @type omega: float
    @param M: mass of the star
    @type M: float
    @param norm: compute magnitude of surface gravity (True) or vector (False)
    @type norm: boolean
    @return: surface gravity magnitude or vector
    @rtype: 3Xfloat/3Xndarray
    """
    GG = constants.GG_sol
    #-- calculate r-component of local gravity
    x = r/r_pole
    grav_r = GG*M/r_pole**2 * (-1./x**2 + 8./27.*x*omega**2*sin(theta)**2)
    #-- calculate theta-component of local gravity
    grav_th = GG*M/r_pole**2 * (8./27.*x*omega**2*sin(theta)*cos(theta))
    
    grav = np.array([grav_r,grav_th])
    #-- now we transform to spherical coordinates
    grav = vectors.spher2cart( (r,phi,theta),(grav[0],0.,grav[1]) )
    grav = grav*constants.Rsol
    if norm:
        return vectors.norm(grav)
    else:
        return grav

def get_fastrot_roche_radius(theta,r_pole,omega):
    """
    Calculate Roche radius for a fast rotating star.
    
    @param theta: angle from rotation axis
    @type theta: float
    @param r_pole: polar radius in solar units
    @type r_pole: float
    @param omega: angular velocity (in units of the critical angular velocity)
    @omega_: float
    @return: radius at angle theta in solar units
    @rtype: float
    """
    #-- calculate surface
    Rstar = 3*r_pole/(omega*sin(theta)) * cos((pi + np.arccos(omega*sin(theta)))/3.)
    #-- solve singularities
    if np.isinf(Rstar) or sin(theta)<1e-10:
        Rstar = r_pole
    return Rstar
    
def critical_angular_velocity(M,R_pole,units='Hz'):
    """
    Compute the critical angular velocity (Hz).
    
    Definition taken from Cranmer and Owocki, 1995 and equal to
    
    Omega_crit = sqrt( 8GM / 27Rp**3 )
    
    Example usage (includes conversion to period in days):
    
    >>> Omega = critical_angular_velocity(1.,1.)
    >>> P = 2*pi/Omega
    >>> P = conversions.convert('s','d',P)
    >>> print 'Critical rotation period of the Sun: %.3f days'%(P)
    Critical rotation period of the Sun: 0.213 days
        
    @param M: mass (solar masses)
    @type M: float
    @param R_pole: polar radius (solar radii)
    @type R_pole: float
    @param units: if you wish other units than Hz, you can give them here
    @type units: string understandable by L{<units.conversions.convert>}
    @return: critical angular velocity in Hz
    @rtype: float
    """
    M = M*constants.Msol
    R = R_pole*constants.Rsol
    omega_crit = np.sqrt( 8*constants.GG*M / (27.*R**3))
    if units.lower()!='hz':
        omega_crit = conversions.convert('Hz',units,omega_crit)
    return omega_crit

def critical_velocity(M,R_pole,units='km/s',definition=1):
    """
    Compute the critical velocity (km/s)
    
    Definition 1 from Cranmer and Owocki, 1995:
    
    v_c = 2 pi R_eq(omega_c) * omega_c
    
    Definition 2 from Townsend 2004:
    
    v_c = sqrt ( 2GM/3Rp )
    
    which both amount to the same value:
    
    >>> critical_velocity(1.,1.,definition=1)
    356.71131858379499
    >>> critical_velocity(1.,1.,definition=2)
    356.71131858379488
    
    @param M: mass (solar masses)
    @type M: float
    @param R_pole: polar radius (solar radii)
    @type R_pole: float
    @param units: if you wish other units than Hz, you can give them here
    @type units: string understandable by L{units.conversions.convert}
    @return: critical velocity in km/s
    @rtype: float
    """
    if definition==1:
        omega_crit = critical_angular_velocity(M,R_pole)
        P = 2*pi/omega_crit
        R_eq = get_fastrot_roche_radius(pi/2.,R_pole,1.)*constants.Rsol
        veq = 2*pi*R_eq/P
    elif definition==2:
        veq = np.sqrt( 2*constants.GG * M*constants.Msol / (3*R_pole*constants.Rsol))
    veq = conversions.convert('m/s',units,veq)
        
    return veq


#}
#{ Differential rotation Roche potential

def diffrot_roche_potential(r,theta,r_pole,M,omega_eq,omega_pole):
    """
    Definition of Roche potential due to differentially rotating star
    
    We first solve the cubic equation
    
    M{re/rp = 1 + f (x^2 + x + 1)/(6x^2)}
        
    where M{f = re^3 Omega_e^2 / (G M)}
    and   M{x = Omega_e / Omega_p}
    
    This transforms to solving
    
    M{re^3 + b * re + c = 0}
        
    where M{b = -1 / (aXrp) and c = 1/(aX),}
    and M{a = Omega_e^2/(GM) and X = (x^2 + x + 1)/(6x^2)}
    
    @param r: radius of the surface element to calculate the surface gravity
    @type r: float/ndarray
    @param theta: colatitude of surface element
    @type theta: float/ndarray
    @param r_pole: polar radius of model
    @type r_pole: float
    @param M: mass of the star
    @type M: float
    @param omega_eq: fraction of critical rotation velocity at equator
    @type omega_eq: float
    @param omega_pole: fraction of critical rotation velocity at pole
    @type omega_pole: float
    @return: roche potential value
    @rtype: float/ndarray
    """
    GG = constants.GG_sol
    
    Omega_crit = sqrt(8*GG*M/ (27*r_pole**3))
    omega_eq = omega_eq*Omega_crit
    omega_pole = omega_pole*Omega_crit
    x = omega_eq / omega_pole
    
    #-- find R_equator solving a cubic equation:
    a = omega_eq**2/(GG*M)
    X = (x**2+x+1)/(6*x**2)
    b,c = -1./(a*X*r_pole), +1./(a*X)
    m,k = 27*c+0*1j,-3*b+0*1j
    n = m**2-4*k**3 + 0*1j
    om1 = -0.5 + 0.5*sqrt(3)*1j
    om2 = -0.5 - 0.5*sqrt(3)*1j
    c1 = (0.5*(m+sqrt(n)))**(1./3.)
    c2 = (0.5*(m-sqrt(n)))**(1./3.)
    x1 = -1./3. * ( c1 + c2 )
    x2 = -1./3. * ( om2*c1 + om1*c2 )
    x3 = -1./3. * ( om1*c1 + om2*c2 )
    re = x2.real
    
    #   ratio of centrifugal to gravitational force at the equator
    f = re**3 * omega_eq**2 / (GG*M)
    #   ratio Re/Rp
    rat = 1 + (f*(x**2+x+1))/(6.*x**2)
    #   some coefficients for easy evaluation
    alpha = f*(x-1)**2/(6*x**2)*(1/rat)**7
    beta  = f*(x-1)   /(2*x**2)*(1/rat)**5
    gamma = f         /(2*x**2)*(1/rat)**3
    #   implicit equation for the surface
    sinth = sin(theta)
    y = r/r_pole
    surf = alpha*y**7*sinth**6 + beta*y**5*sinth**4 + gamma*y**3*sinth**2 - y +1
    return surf

def diffrot_roche_surface_gravity(r,theta,phi,r_pole,M,omega_eq,omega_pole,norm=False):
    """
    Surface gravity from differentially rotation Roche potential.
    
    Magnitude is OK, please carefully check direction.
    
    @param r: radius of the surface element to calculate the surface gravity
    @type r: float/ndarray
    @param theta: colatitude of surface element
    @type theta: float/ndarray
    @param phi: longitude of surface element
    @type phi: float/ndarray
    @param r_pole: polar radius of model
    @type r_pole: float
    @param M: mass of the star
    @type M: float
    @param omega_eq: fraction of critical rotation velocity at equator
    @type omega_eq: float
    @param omega_pole: fraction of critical rotation velocity at pole
    @type omega_pole: float
    @param norm: compute magnitude of surface gravity (True) or vector (False)
    @type norm: boolean
    @return: surface gravity magnitude or vector
    @rtype: 3Xfloat/3Xndarray
    """
    GG = constants.GG_sol
    Omega_crit = sqrt(8*GG*M/ (27*r_pole**3))
    omega_eq = omega_eq*Omega_crit
    omega_pole = omega_pole*Omega_crit
    x = omega_eq / omega_pole
    
    #-- find R_equator solving a cubic equation:
    a = omega_eq**2/(GG*M)
    X = (x**2+x+1)/(6*x**2)
    b,c = -1./(a*X*r_pole), +1./(a*X)
    m,k = 27*c+0*1j,-3*b+0*1j
    n = m**2-4*k**3 + 0*1j
    om1 = -0.5 + 0.5*sqrt(3)*1j
    om2 = -0.5 - 0.5*sqrt(3)*1j
    c1 = (0.5*(m+sqrt(n)))**(1./3.)
    c2 = (0.5*(m-sqrt(n)))**(1./3.)
    x1 = -1./3. * ( c1 + c2 )
    x2 = -1./3. * ( om2*c1 + om1*c2 )
    x3 = -1./3. * ( om1*c1 + om2*c2 )
    re = x2.real
    
    #   ratio of centrifugal to gravitational force at the equator
    f = re**3 * omega_eq**2 / (GG*M)
    #   ratio Re/Rp
    rat = 1 + (f*(x**2+x+1))/(6.*x**2)
    #   some coefficients for easy evaluation
    alpha = f*(x-1)**2/(6*x**2)*(1/rat)**7
    beta  = f*(x-1)   /(2*x**2)*(1/rat)**5
    gamma = f         /(2*x**2)*(1/rat)**3
    #   implicit equation for the surface
    sinth = sin(theta)
    y = r/r_pole
    grav_th = (6*alpha*y**7*sinth**5 + 4*beta*y**5*sinth**3 + 2*gamma*y**3*sinth)*cos(theta)
    grav_r = 7*alpha/r_pole*y**6*sinth**6 + 5*beta/r_pole*y**4*sinth**4 + 3*gamma/r_pole*y**2*sinth**2 - 1./r_pole
    
    fr = 6*alpha*y**7*sinth**4 + 4*beta*y**5*sinth**2 + 2*gamma*y**3
    magn = GG*M/r**2 * sqrt(cos(theta)**2 + (1-fr)**2 *sin(theta)**2)
    magn_fake = np.sqrt(grav_r**2+(grav_th/r)**2)
    grav_r,grav_th = grav_r/magn_fake*magn,(grav_th/r)/magn_fake*magn
    
    grav = np.array([grav_r*constants.Rsol,grav_th*constants.Rsol])
    #-- now we transform to spherical coordinates
    grav = vectors.spher2cart( (r,phi,theta),(grav[0],0.,grav[1]) )
    
    if norm:
        return vectors.norm(grav)
    else:
        return grav


def get_diffrot_roche_radius(theta,r_pole,M,omega_eq,omega_pole):
    """
    Calculate Roche radius for a differentially rotating star.
    
    @param theta: angle from rotation axis
    @type theta: float
    @param r_pole: polar radius in solar units
    @type r_pole: float
    @param M: mass in solar units
    @type M: float
    @param omega_eq: equatorial angular velocity (in units of the critical angular velocity)
    @omega_eq: float
    @param omega_pole: polar angular velocity (in units of the critical angular velocity)
    @omega_pole: float
    @return: radius at angle theta in solar units
    @rtype: float
    """
    try:
        r = newton(diffrot_roche_potential,r_pole,args=(theta,r_pole,M,omega_eq,omega_pole))
    except RuntimeError:
        r = np.nan    
    return r

def diffrot_law(omega_eq,omega_pole,theta):
    """
    Evaluate a differential rotation law of the form Omega = b1+b2*omega**2
    
    The relative differential rotation rate is the ratio of the rotational shear
    to the equatorial velocity
    
    alpha = omega_eq - omega_pole / omega_eq
    
    The units of angular velocity you put in, you get out (i.e. in terms of
    the critical angular velocity or not).
    
    @param omega_eq: equatorial angular velocity
    @type omega_eq: float
    @param omega_pole: polar angular velocity
    @type omega_pole: float
    @param theta: colatitude (0 at the pole) in radians
    @type theta: float (radians)
    @return: the angular velocity at colatitude theta
    @rtype: float
    """
    b1 = omega_eq
    b2 = omega_pole - omega_eq
    omega_theta = b1 + b2 * cos(theta)**2
    return omega_theta


def diffrot_velocity(coordinates,omega_eq,omega_pole,R_pole,M):
    """
    Calculate the velocity vector of every surface element.
    
    @param coordinates: polar coordinates of stellar surface (phi,theta,radius)
    make sure the radius is in SI units!
    @type coordinates: 3xN array
    @param omega_eq: equatorial angular velocity (as a fraction of the critical one)
    @type omega_eq: float
    @param omega_pole: polar angular velocity (as a fraction of the critical one)
    @type omega_pole: float
    @param R_pole: polar radius in solar radii
    @type R_pole: float
    @param M: stellar mass in solar mass
    @type M: float
    @return: velocity vectors of all surface elements
    @rtype 3xN array
    """
    #-- angular velocity of every surface element
    Omega_crit = critical_angular_velocity(M,R_pole)
    phi,theta,radius = coordinates
    omega_local = diffrot_law(omega_eq,omega_pole,theta)*Omega_crit
    #-- direction of local angular velocity in Cartesian coordinates (directed in upwards z)
    omega_local_vec = np.array([np.zeros_like(omega_local),np.zeros_like(omega_local),omega_local]).T
    
    x,y,z = vectors.spher2cart_coord(radius,phi,theta)
    surface_element = np.array([x,y,z]).T

    velo_local = np.array([np.cross(ielement,iomega_local_vec) for ielement,iomega_local_vec in zip(surface_element,omega_local_vec)]).T
    return velo_local

#}


if __name__=="__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()
