"""
Compute stellar shapes via Roche potential and gradients.

The premisse of the Roche potential is that the stellar mass can be represented
by a point source.

Currently, three types of stellar shapes are implemented:
    1. asynchronously rotating eccentric binary
    2. fast rotating star
    3. differentially (fast) rotating star

This module can be used to calculate the following information:
    1. the distorted shape of the star due to a Roche potential
    2. the local surface gravity assuming Von Zeipel gravity darkening
    3. the local effective temperature
    4. the local velocity vector due to rotation
    5. the radial velocity
    
This information can then be used to calculate:
    1. the total distorted surface area
    2. the total distorted luminosity
    3. the projected intensity in the line of sight
    4. the mean radial velocity (e.g. Rossiter effect in binaries)
    5. a synthetic spectrum starting from a library of spectra
    6. a synthetic light curve from the projected intensities

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

>>> theta,phi = get_grid(50,full=True)
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

>>> areas_local,cos_gamma = surface_elements((radius,theta,phi),-grav_local)

Now we can calculate the local effective temperature (assuming radiative
atmosphere in Von Zeipel law), and bolometric intensities:

>>> teff_local = local_temperature(vectors.norm(grav_local),g_pole,T_pole,beta=1.)
>>> ints_local = local_intensity(teff_local,grav,photband='OPEN.BOL')

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

>>> intens_proj = local_intensity(teff_local.ravel(),grav.ravel(),mus,photband='OPEN.BOL').reshape(theta.shape)

The total and projected luminosity can then be calculated the following way:

>>> print pi*(ints_local*areas_local*constants.Rsol_cgs**2).sum()/constants.Lsol_cgs
1.00299508526
>>> print pi*np.nansum(intens_proj*areas_local*constants.Rsol_cgs**2)/constants.Lsol_cgs
0.361744218355

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

gives:

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

]include figure]]ivs_binary_rochepotential_fastrotstar_shape.png]

"""
import pylab as pl
import numpy as np
from numpy import pi,cos,sin,sqrt,nan
from scipy.optimize import newton
from scipy.spatial import KDTree
import time
from ivs.binary import keplerorbit
from ivs.units import constants
from ivs.units import vectors
from ivs.sed import model as sed_model
from ivs.sed import limbdark
from ivs.misc import loggers
from ivs.io import ascii

logger = loggers.get_basic_logger("")

#{ Eccentric asynchronous binary Roche equation in spherical coordinates

def binary_roche_potential(r,theta,phi,Phi,q,d,F):
    """
    Unitless eccentric asynchronous Roche potential in spherical coordinates.
    
    See Wilson, 1979.
    
    The  synchronicity parameter F is 1 for synchronised circular orbits. For
    pseudo-synchronous eccentrical orbits, it is equal to (Hut, 1981)
    
    F = sqrt( (1+e)/ (1-e)^3)
    
    Periastron is reached when d = 1-e.
        
    @param r: radius of Roche volume at potential Phi (in units of semi-major axis)
    @type r: float
    @param theta: colatitude (0 at the pole, pi/2 at the equator)
    @type theta: float
    @param phi: longitude (0 in direction of COM)
    @type phi: float
    @param Phi: Roche potential value (unitless)
    @type Phi: float
    @param q: mass ratio
    @type q: float
    @param d: separation (in units of semi-major axis)
    @type d: float
    @param F: synchronicity parameter
    @type F: float
    @return: residu between Phi and roche potential
    @rtype: float
    """
    lam,nu = cos(phi)*sin(theta),cos(theta)
    term1 = 1. / r
    term2 = q * ( 1./sqrt(d**2 - 2*lam*d*r + r**2) - lam*r/d**2)
    term3 = 0.5 * F**2 * (q+1) * r**2 * (1-nu**2)
    return (Phi - (term1 + term2 + term3))

def binary_roche_potential_gradient(x,y,z,q,d,F,norm=False):
    """
    Gradient of eccenctric asynchronous Roche potential in cartesian coordinates.
    
    See Phoebe scientific reference, http://phoebe.fiz.uni-lj.si/docs/phoebe_science.ps.gz
    
    x,y,z,d in real units! (otherwise you have to scale it yourself)
    
    @param x: x-axis
    @type x: float'
    @param y: y-axis
    @type y: float'
    @param z: z-axis
    @type z: float'
    @param q: mass ratio
    @type q: float
    @param d: separation (in units of semi-major axis)
    @type d: float
    @param F: synchronicity parameter
    @type F: float
    @param norm: flag to return magnitude (True) or vector form (False)
    @type norm: boolean
    @return: Roche potential gradient
    @rtype: ndarray or float
    """
    r = sqrt(x**2 + y**2 + z**2)
    r_= sqrt((d-x)**2 + y**2 + z**2)
    dOmega_dx = - x / r**3 + q * (d-x) / r_**3 + F**2 * (1+q)*x - q/d**2
    dOmega_dy = - y / r**3 - q * y     / r_**3 + F**2 * (1+q)*y
    dOmega_dz = - z / r**3 - q * z     / r_**3
    
    dOmega = np.array([dOmega_dx,dOmega_dy,dOmega_dz])
    if norm:
        return vectors.norm(dOmega)
    else:
        return dOmega


def binary_roche_surface_gravity(x,y,z,d,omega,M1,M2,a=1.,norm=False):
    """
    Calculate surface gravity in an eccentric asynchronous binary roche potential.
    """
    q = M2/M1
    x_com = q*d/(1+q)
    
    r = np.array([x,y,z])
    d_cf = np.array([d-x_com,0,0])
    d = np.array([d,0,0])
    h = d - r
    
    term1 = - G*M1/vectors.norm(r)**3*r
    term2 = - G*M2/vectors.norm(h)**3*h
    term3 = - omega**2 * d_cf
    
    #print term1,term2,term3
    #print omega**2,d_cf,x_com,d
    #print 'r',r
    #print 'd',d
    #print 'x_com',x_com
    #print 'd_cf',d_cf
    #print 'h',h
    #print '-->gravity components (logg, cgs dex)',np.log10(vectors.norm(term1)*100),np.log10(vectors.norm(term2)*100),np.log10(vectors.norm(term3)*100)
    #sys.exit()
    g_pole = term1 + term2 + term3
    if norm:
        return vectors.norm(g_pole)
    else:
        return g_pole


def get_binary_roche_radius(theta,phi,Phi,q,d,F,r_pole=None):
    """
    Calculate the eccentric asynchronous binary Roche radius in spherical coordinates.
    
    This is done via the Newton-Raphson secant method. If r_pole is not given
    as a starting value, it will be calculated here (slowing down the function).
    
    If no radius can be calculated for the given coordinates, 'nan' is returned.
    
    @param theta: colatitude (0 at the pole, pi/2 at the equator)
    @type theta: float
    @param phi: longitude (0 in direction of COM)
    @type phi: float
    @param Phi: Roche potential value (unitless)
    @type Phi: float
    @param q: mass ratio
    @type q: float
    @param d: separation (in units of semi-major axis)
    @type d: float
    @param F: synchronicity parameter
    @type F: float
    @param r_pole: polar radius (serves as starting value for NR method)
    @type r_pole: float
    @return r: radius of Roche volume at potential Phi (in units of semi-major axis)
    @rtype r: float
    """
    if r_pole is None:
        r_pole = newton(binary_roche_potential,1e-5,args=(0,0,Phi,q,ds.min(),F))
    try:
        r = newton(binary_roche_potential,r_pole,args=(theta,phi,Phi,q,d,F))
        if r<0 or r>d:
            r = nan
    except RuntimeError:
        r = nan    
    return r

#}

#{ Fast rotation Roche potential

def fastrot_roche_surface_gravity(r,theta,phi,r_pole,omega,M,norm=False):
    """
    Calculate components of the local surface gravity of the fast rotating
    Roche model.
    
    Input units are solar units.
    Output units are SI.
    Omega is fraction of critical velocity.
    
    See Cranmer & Owocki, Apj (1995)
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
    

#}
#{ Differential rotation Roche potential

def diffrot_roche_potential(r,theta,r_pole,M,omega_eq,omega_pole):
    """
    Definition of Roche potential due to differentially rotating star
    
    We first solve the cubic equation
    
        re/rp = 1 + f (x^2 + x + 1)/(6x^2)
        
    where f = re^3 Omega_e^2 / (G M)
    and   x = Omega_e / Omega_p
    
    This transforms to solving
    
        re^3 + b * re + c = 0
        
    where b = -1 / (aXrp) and c = 1/(aX),
    and a = Omega_e^2/(GM) and X = (x^2 + x + 1)/(6x^2)
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
    
    Magnitude is OK, direction doesn't seem to be.
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
    
    grav = np.array([grav_r*constatns.Rsol,grav_th*constants.Rsol])
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
        r = nan    
    return r
#}

#{ General functions

def surface_elements((r,theta,phi),(surfnormal_x,surfnormal_y,surfnormal_z)):
    """
    Compute surface area of elements in a grid.
    
    theta,phi must be generated like mgrid(theta_range,phi_range)
    
    usually, the surfnormals are acquired via differentiation of a gravity potential,
    and is then equal to the *negative* of the local surface gravity.
    """
    #-- compute the grid size at each location
    dtheta = theta[1:]-theta[:-1]
    dtheta = np.vstack([dtheta,dtheta[-1]])
    
    dphi = phi[:,1:]-phi[:,:-1]
    dphi = np.column_stack([dphi,dphi[:,-1]])
    #-- compute the angle between the surface normal and the radius vector
    x,y,z = vectors.spher2cart_coord(r,phi,theta)
    
    a = np.array([x,y,z])
    b = np.array([surfnormal_x,surfnormal_y,surfnormal_z])
    
    cos_gamma = vectors.cos_angle(a,b)
    
    return r**2 * sin(theta) * dtheta * dphi / cos_gamma, cos_gamma

def local_temperature(surface_gravity,g_pole,T_pole,beta=1.):
    """
    Calculate local temperature.
    
    beta is gravity darkening parameter.
    """
    Grav = abs(surface_gravity/g_pole)**beta
    Teff = Grav**0.25 * T_pole
    return Teff

def local_intensity(teff,grav,mu=None,photband='OPEN.BOL'):
    """
    Calculate local intensity.
    
    beta is gravity darkening parameter.
    """
    if mu is None:
        mu = np.ones_like(teff)
    if (teff<3500).any() or np.isnan(teff).any():
        print 'WARNING: point outside of grid, minimum temperature is 3500K'
        teff = np.where((teff<3500) | np.isnan(teff),3500,teff)
    if (grav<0.01).any() or np.isnan(grav).any():
        print 'WARNING: point outside of grid, minimum gravity is 0 dex'
        grav = np.where((np.log10(grav*100)<0.) | np.isnan(grav),0.01,grav)
    intens = np.array([limbdark.get_itable(teff=iteff,logg=np.log10(igrav*100),absolute=True,mu=imu,photbands=[photband])[0] for iteff,igrav,imu in zip(teff.ravel(),grav.ravel(),mu.ravel())])
    return intens.reshape(teff.shape)
    

def projected_intensity(teff,gravity,areas,line_of_sight,photband='OPEN.BOL'):
    """
    Compute projected intensity in the line of sight.
    
    gravity is vector directed inwards in the star
    line of sight is vector.
    """
    ones = np.ones_like(gravity[0])
    losx = line_of_sight[0]*ones
    losy = line_of_sight[1]*ones
    losz = line_of_sight[2]*ones
    angles = vectors.angle(-gravity,np.array([losx,losy,losz]))
    mus = cos(angles)
    grav_ = vectors.norm(gravity)
    #-- intensity is less if we look at the limb
    intens = local_intensity(teff,grav_,mus,photband=photband)
    #-- intensity is less if the surface element area is small or on the limb
    return intens*areas*mus,mus


def reflection_effect(primary,secondary,theta,phi,A1=1.,A2=1.,max_iter=1):
    #-- reflection effect
    #--------------------
    reflection_iter = 0
    while (reflection_iter<max_iter):
        R1,R2 = np.ones(len(primary['teff'])/4),np.ones(len(primary['teff'])/4)
        for i in xrange(len(R1)):
            if i%100==0: print i,len(R1)
            #-- radiation from secondary onto primary
            s12 = np.array([ primary['x'][i]-secondary['x'],
                             primary['y'][i]-secondary['y'],
                             primary['z'][i]-secondary['z']])
            psi2 = vectors.angle(+s12,-np.array([secondary['gravx'],secondary['gravy'],secondary['gravz']]))
            psi1 = vectors.angle(-s12,-np.array([primary['gravx'][i:i+1],primary['gravy'][i:i+1],primary['gravz'][i:i+1]]))
            s12 = s12.ravel()
            
            keep = (psi2<pi/2.) & (psi1<pi/2.)
            Lambda_2_bol = np.array([limbdark.get_itable(teff=iteff,logg=np.log10(igrav*100),theta=ipsi2,photbands=['OPEN.BOL'],absolute=False)[0]\
                            for iteff,igrav,ipsi2 in zip(secondary['teff'][keep],secondary['grav'][keep],psi2[keep])])
    
            s = vectors.norm(s12[keep])
            J_21_entrant = A1 * R2[i] * np.sum(secondary['flux'][keep] * cos(psi1[keep])*cos(psi2[keep])*Lambda_2_bol*secondary['areas'][keep]/s**2)
            if np.isnan(J_21_entrant).any(): raise ValueError
            J_1_tot      = primary['flux'][i]
            R1_new = 1+ J_21_entrant/J_1_tot
        
            #-- radiation from primary onto secondary
            s12 = np.array([primary['x']-secondary['x'][i],
                            primary['y']-secondary['y'][i],
                            primary['z']-secondary['z'][i]])
            psi2 = vectors.angle(+s12,-np.array([primary['gravx'],primary['gravy'],primary['gravz']]))
            psi1 = vectors.angle(-s12,-np.array([secondary['gravx'][i:i+1],secondary['gravy'][i:i+1],secondary['gravz'][i:i+1]]))
            s12 = s12.ravel()
            
            keep = (psi2<pi/2.) & (psi1<pi/2.)
            Lambda_2_bol = np.array([limbdark.get_itable(teff=iteff,logg=np.log10(igrav*100),theta=ipsi2,photbands=['OPEN.BOL'],absolute=False)[0]\
                    for iteff,igrav,ipsi2 in zip(primary['teff'][keep],primary['grav'][keep],psi2[keep])])
        
            s = vectors.norm(s12)
            J_12_entrant = A2 * R1[i] * np.nansum(primary['flux'][keep] * cos(psi1[keep])*cos(psi2[keep])*Lambda_2_bol*primary['areas'][keep]/s**2)
            J_2_tot      = secondary['flux'][i]
            R2_new = 1+ J_12_entrant/J_2_tot
                
            R1[i] = R1_new
            R2[i] = R2_new
        
        #================ START DEBUGGING PLOTS ===================
        #pl.figure()
        #pl.subplot(423);pl.title('old primary')
        #pl.scatter(phis_,thetas_,c=teff_,edgecolors='none',cmap=pl.cm.spectral)
        #pl.colorbar()
        #pl.subplot(424);pl.title('old secondary')
        #pl.scatter(phis_,thetas_,c=teff2_,edgecolors='none',cmap=pl.cm.spectral)
        #pl.colorbar()
        #================   END DEBUGGING PLOTS ===================
        
        #-- adapt the teff only (and when) the increase is more than 1% (=>1.005**0.25=1.01)
        break_out = True
        trash,trash2,R1,R2 = stitch_grid(theta,phi,R1.reshape(theta.shape),R2.reshape(theta.shape))
        del trash,trash2
        R1 = R1.ravel()
        R2 = R2.ravel()
        
        if (R1[-np.isnan(R1)]>1.05).any():
            print "Significant reflection effect on primary (max %.3f%%)"%((R1.max()**0.25-1)*100)
            primary['teff']*= R1**0.25
            primary['flux'] = local_intensity(primary['teff'],primary['grav'],np.ones_like(primary['teff']),photbands=['OPEN.BOL'])
            break_out = False
        else:
            print 'Maximum reflection effect on primary: %.3f%%'%((R1.max()**0.25-1)*100)
        
        if (R2[-np.isnan(R2)]>1.05).any():
            print "Significant reflection effect on secondary (max %.3g%%)"%((R2.max()**0.25-1)*100)
            secondary['teff']*= R2**0.25
            secondary['flux'] = local_intensity(secondary['teff'],secondary['grav'],np.ones_like(secondary['teff']),photbands=['OPEN.BOL'])
            break_out = False
        else:
            print 'Maximum reflection effect on secondary: %.3g%%'%((R1.max()**0.25-1)*100)
        
        if break_out:
            break
            
        
        reflection_iter += 1
        
        #================ START DEBUGGING PLOTS ===================
        #pl.subplot(411)
        #pl.scatter(x,y,c=teff_,edgecolors='none',cmap=pl.cm.spectral,vmin=T_pole,vmax=T_pole2)
        #pl.scatter(x2,y2,c=teff2_,edgecolors='none',cmap=pl.cm.spectral,vmin=T_pole,vmax=T_pole2)
        #pl.colorbar()
        
        #pl.subplot(425);pl.title('new primary')
        #pl.scatter(phis_,thetas_,c=teff_,edgecolors='none',cmap=pl.cm.spectral)
        #pl.colorbar()
        #pl.subplot(426);pl.title('new secondary')
        #pl.scatter(phis_,thetas_,c=teff2_,edgecolors='none',cmap=pl.cm.spectral)
        #pl.colorbar()
        
        #pl.subplot(427)
        #pl.scatter(phis_,thetas_,c=R1,edgecolors='none',cmap=pl.cm.spectral)
        #pl.colorbar()
        #pl.subplot(428)
        #pl.scatter(phis_,thetas_,c=R2,edgecolors='none',cmap=pl.cm.spectral)
        #pl.colorbar()
        #================   END DEBUGGING PLOTS ===================
    return primary,secondary

def get_grid(*args,**kwargs):
    """
    Construct a coordinate grid
    """
    gtype = kwargs.get('gtype','spherical')
    full = kwargs.get('full',False)
    
    if 'spher' in gtype.lower():
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution
        
        dtheta = pi/2/res1 # step size coordinate 1
        dphi = pi/res2     # step size coordinate 2
        
        #-- full grid or only one quadrant
        if full:
            theta0,thetan = dtheta/2, pi-dtheta/2
            phi0,phin = 0,2*pi-dphi
        else:
            theta0,thetan = dtheta/2,pi/2-dtheta/2
            phi0,phin = 0,pi-dphi
            
        theta,phi = np.mgrid[theta0:thetan:res1*1j,phi0:phin:res2*1j]
        return theta,phi

def stitch_grid(theta,phi,*quant,**kwargs):
    """
    Stitch a grid together that was originally defined on 1 quadrant.
    
    We add the three other quandrants.
    """
    seamless = kwargs.get('seamless',False)
    vtype = kwargs.get('vtype',['scalar' for i in quant])
    ravel = kwargs.get('ravel',False)
    
    #-- basic coordinates
    alltheta = np.vstack([np.hstack([theta,theta]),np.hstack([theta+pi/2,theta+pi/2])])
    allphi   = np.vstack([np.hstack([phi,phi+pi]),np.hstack([phi,phi+pi])])
    
    #-- all other quantities
    allquan = []
    for i,iquant in enumerate(quant):
        #-- if they are scalar values, they do not change direction
        top1,top2 = iquant,iquant[:,::-1]
        bot1,bot2 = iquant[::-1],iquant[::-1][:,::-1]
        if vtype[i]=='scalar':
            allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([bot1,bot2])]))
        #-- vector components do change direction
        elif vtype[i]=='x':
            allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([bot1,bot2])]))
        elif vtype[i]=='y':
            allquan.append(np.vstack([np.hstack([top1,-top2]),np.hstack([bot1,-bot2])]))
        elif vtype[i]=='z':
            allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([-bot1,-bot2])]))
        elif vtype[i]=='vx':
            allquan.append(np.vstack([np.hstack([top1,-top2]),np.hstack([bot1,-bot2])]))
        elif vtype[i]=='vy':
            allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([bot1,bot2])]))
        elif vtype[i]=='vz':
            allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([-bot1,-bot2])]))
    
    out = [alltheta,allphi]+allquan
    
    #-- for plotting reasons, remove the vertical seam and the the bottom hole.
    if seamless:
        #-- vertical seam
        out = [np.column_stack([i,i[:,0]]) for i in out]
        out[1][:,-1] += 2*pi
        #-- fill bottom hole
        out = [np.vstack([i,i[0]]) for i in out]
        out[0][-1] += pi
    
    #-- ravel to 1d arrays if asked for
    if ravel:
        out = [i.ravel() for i in out]
    
    return out

#}

def diffrot_test():
    res = 20
    theta,phi = np.mgrid[0:pi/2:res*1j,0:pi:res*1j]
    thetas,phis = np.ravel(theta),np.ravel(phi)
    thetas = np.linspace(0,pi/2,100)
    phis = np.zeros_like(thetas)
    dtheta = thetas.ptp() / (res-1)
    dphi = phis.ptp() / (res-1)
    r_pole = 1.
    M = 1.
    omega_eq = 0.00#1
    omega_pole = 0.00#1
    omegas = [0.01,0.1,0.2,0.5,0.7,0.9,0.99,1.0,2.]
    pl.figure()
    rows,cols = 2,4
    color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(omegas))]
    for i in range(1,rows*cols+1):
        pl.subplot(rows,cols,i)
        pl.gca().set_color_cycle(color_cycle)
    
    phi = pi/4
    for omega in omegas:
        radius = np.array([get_diffrot_roche_radius(itheta,r_pole,M,omega,omega) for itheta in thetas])
        gravit = np.array([diffrot_roche_surface_gravity(iradius,itheta,phi,r_pole,M,omega,omega) for iradius,itheta in zip(radius,thetas)])
        radius2= np.array([get_fastrot_roche_radius(itheta,r_pole,omega) for itheta in thetas])
        gravit2= np.array([fastrot_roche_surface_gravity(iradius,itheta,phi,r_pole,omega,M) for iradius,itheta in zip(radius2,thetas)])
        pl.subplot(rows,cols,1)
        pl.plot(radius*sin(thetas),radius*cos(thetas),'-',label='%.2f'%(omega))
        pl.gca()._get_lines.count -= 1
        pl.plot(radius2*sin(thetas),radius2*cos(thetas),'--',lw=2)
        pl.subplot(rows,cols,2)
        pl.scatter(radius2*sin(thetas),radius2*cos(thetas),c=np.log10(vectors.norm(gravit2.T)*100),edgecolors='none',cmap=pl.cm.spectral,vmin=3,vmax=4.45)
        pl.subplot(rows,cols,3)
        pl.scatter(radius*sin(thetas),radius*cos(thetas),c=np.log10(vectors.norm(gravit.T)*100),edgecolors='none',cmap=pl.cm.spectral,vmin=3,vmax=4.45)
        pl.subplot(rows,cols,4)
        pl.plot(gravit.T[0],gravit.T[1],'-')
        pl.gca()._get_lines.count -= 1
        pl.plot(gravit2.T[0],gravit2.T[1],'--',lw=2)
    
    for omega in omegas:
        radius = np.array([get_diffrot_roche_radius(itheta,r_pole,M,omega,1.5*omega) for itheta in thetas])
        gravit = np.array([diffrot_roche_surface_gravity(iradius,itheta,phi,r_pole,M,omega,1.5*omega) for iradius,itheta in zip(radius,thetas)])
        radius2= np.array([get_fastrot_roche_radius(itheta,r_pole,omega) for itheta in thetas])
        gravit2= np.array([fastrot_roche_surface_gravity(iradius,itheta,phi,r_pole,omega,M) for iradius,itheta in zip(radius2,thetas)])
        pl.subplot(rows,cols,5)
        pl.plot(radius*sin(thetas),radius*cos(thetas),'-',label='%.2f'%(omega))
        pl.gca()._get_lines.count -= 1
        pl.plot(radius2*sin(thetas),radius2*cos(thetas),'--',lw=2)
        pl.subplot(rows,cols,6)
        pl.scatter(radius2*sin(thetas),radius2*cos(thetas),c=np.log10(vectors.norm(gravit2.T)*100),edgecolors='none',cmap=pl.cm.spectral,vmin=3,vmax=4.45)
        pl.subplot(rows,cols,7)
        pl.scatter(radius*sin(thetas),radius*cos(thetas),c=np.log10(vectors.norm(gravit.T)*100),edgecolors='none',cmap=pl.cm.spectral,vmin=3,vmax=4.45)
        pl.subplot(rows,cols,8)
        pl.plot(gravit.T[0],gravit.T[1],'-')
        pl.gca()._get_lines.count -= 1
        pl.plot(gravit2.T[0],gravit2.T[1],'--',lw=2)
    
    
    pl.subplot(rows,cols,1);pl.ylim(pl.xlim());pl.legend()
    pl.subplot(rows,cols,2);pl.ylim(pl.xlim());pl.colorbar()
    pl.subplot(rows,cols,3);pl.ylim(pl.xlim());pl.colorbar()
    pl.subplot(rows,cols,5);pl.ylim(pl.xlim());pl.legend()
    pl.subplot(rows,cols,6);pl.ylim(pl.xlim());pl.colorbar()
    pl.subplot(rows,cols,7);pl.ylim(pl.xlim());pl.colorbar()
    pl.subplot(rows,cols,4)#;pl.xlim(-300,0);pl.ylim(0,50)
    pl.subplot(rows,cols,8)#;pl.xlim(-300,0);pl.ylim(0,50)
    

def full_test_diff_rot():
    #-- full calculation of one example
    omega = 1.0
    T_pole = 15000.
    r_pole = 2.
    M = 1.
    view_angle = 90./180*pi
    
    g_pole = fastrot_roche_surface_gravity(r_pole,0,0,r_pole,omega,M)[-1]
    theta,phi = get_grid(50)
    thetas,phis = np.ravel(theta),np.ravel(phi)
    
    radius = np.array([get_fastrot_roche_radius(itheta,r_pole,omega) for itheta in thetas]).reshape(theta.shape)
    grav_local = np.array([fastrot_roche_surface_gravity(iradius,itheta,iphi,r_pole,omega,M) for iradius,itheta,iphi in zip(radius.ravel(),thetas,phis)]).T
    #radius = np.array([get_diffrot_roche_radius(itheta,r_pole,M,omega,2*omega) for itheta in thetas]).reshape(theta.shape)
    #grav_local = np.array([diffrot_roche_surface_gravity(iradius,itheta,iphi,r_pole,M,omega,2*omega) for iradius,itheta,iphi in zip(radius.ravel(),thetas,phis)]).T
    grav_local = np.array([i.reshape(theta.shape) for i in grav_local])
    grav = vectors.norm(grav_local)
    
    #-- derive all local quantities
    areas_local,cos_gamma = surface_elements((radius,theta,phi),-grav_local)
    teff_local = local_temperature(vectors.norm(grav_local),g_pole,T_pole,beta=1.)
    ints_local = local_intensity(teff_local,grav,np.ones_like(cos_gamma),photbands=['OPEN.BOL'])
    line_of_sight = np.zeros_like(grav_local)
    line_of_sight[0] = sin(view_angle)
    line_of_sight[2] = cos(view_angle)
    
    #-- stitch the grid!
    thetas_,phis_,radius_,gravx_,gravy_,gravz_,losx_,losy_,losz_,grav_,areas_,teff_,intens_ = stitch_grid(theta,phi,radius,
                            grav_local[0],grav_local[1],grav_local[2],
                            line_of_sight[0],line_of_sight[1],line_of_sight[2],
                            grav,areas_local,teff_local,ints_local,seamless=False,ravel=True,
                            vtype=['scalar','x','y','z','scalar','scalar','scalar','scalar','scalar','scalar','scalar'])
    
    angles_ = vectors.angle(np.array([-gravx_,-gravy_,-gravz_]),np.array([losx_,losy_,losz_]))
    mus_ = cos(angles_)
    intens_proj_ = local_intensity(teff_,grav_,mus_,photbands=['OPEN.BOL'])
    print "Luminosity",pi*(intens_*areas_*constants.Rsol_cgs**2).sum()/constants.Lsol_cgs
    print "Projected",pi*np.nansum(intens_proj_*areas_*constants.Rsol_cgs**2)/constants.Lsol_cgs
    
    pl.figure()
    rows,cols = 4,3
    
    pl.subplot(rows,cols,1);pl.title('area')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=areas_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,2);pl.title('grav')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=np.log10(grav_*100),edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,3);pl.title('teff')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=teff_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,4);pl.title('intens')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=intens_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,5);pl.title('intens_proj')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=intens_proj_*areas_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,6);pl.title('angle')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=angles_/pi*180,edgecolors='none',vmin=0,vmax=90)
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,7);pl.title('radius')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=radius_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,10);pl.title('gravx')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=-gravx_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,11);pl.title('gravy')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=-gravy_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    pl.subplot(rows,cols,12);pl.title('gravz')
    pl.scatter(phis_/pi*180,thetas_/pi*180,c=-gravz_,edgecolors='none')
    pl.colorbar();pl.xlim(0,360);pl.ylim(180,0)
    
    
    x,y,z = vectors.spher2cart_coord(radius_,phis_,thetas_)
    output = np.array([x,y,z,gravx_,gravy_,gravz_,grav_,areas_,teff_]).T
    ascii.write_array(output,'test/fastrot_shape.dat')
    
    pl.show()


def binary_light_curve_synthesis(**parameters):
    
    #-- some parameters are optional
    #   file output parameters
    name = parameters.setdefault('name','mybinary')
    direc = parameters.setdefault('direc','')
    #   calculation details
    res = parameters.setdefault('gres',20)                   # resolution of the grid
    tres= parameters.setdefault('tres',125)                  # resolution of the phase diagram
    passband = parameters.setdefault('passband','JOHNSON.V') # photometric passband
    max_iter_reflection = parameters.setdefault('ref_iter',1) # maximum number of iterations of reflection effect
    #   orbital parameters
    gamma = parameters.setdefault('gamma',0.)            # systemic velocity [km/s]
    incl = parameters.setdefault('incl',90.)             # system inclination angle [deg]
    q = parameters.setdefault('q',1.)                    # mass ratio (M2/M1)
    e = parameters.setdefault('e',0.)                    # eccentricity
    omega = parameters.setdefault('omega',0.)            # argument of periastron
    F = parameters.setdefault('F1',sqrt((1+e)/(1-e)**3)) # synchronicity parameter primary
    F2= parameters.setdefault('F2',sqrt((1+e)/(1-e)**3)) # synchronicity parameter secondary
    A1= parameters.setdefault('A1',1.)                   # albedo primary
    A2= parameters.setdefault('A2',1.)                   # albedo secondary
    
    #-- others are mandatory:
    T_pole = parameters['Tpole1']       # Primary Polar temperature   [K]
    T_pole2 = parameters['Tpole2']      # Secondary Polar temperature [K]
    P = parameters['P']                 # Period [days]
    asini = parameters['asini']         # total semi-major axis*sini  [AU]
    Phi  = parameters['Phi1']           # Gravitational potential of primary [-]
    Phi2 = parameters['Phi2']           # Gravitational potential of secondary [-]
    
    #-- derive parameters needed in the calculation
    a = asini/sin(incl/180.*pi)
    a1 = a / (1.+1./q)
    a2 = a / (1.+   q)
    q2 = 1./q
    view_angle = incl/180.*pi
    #   calculate the total mass (in solar mass units) from Kepler's third law:
    a_ = conversions.convert('au','m',a)
    P_ = conversions.convert('d','s',P)
    M  = (4*pi**2 * a_**3 / P_**2 / constants.GG ) / constants.Msol
    M1 = M / (1.+   q)
    M2 = M / (1.+1./q)
    #   calculate location of first Lagrangian point L1 and L2
    mu = M2 / (M1+M2)
    z = (mu/3.)**(1./3.)
    L1 = z- 1/3.*z**2 - 1/9.*z**3 + 58./81*z**4
    L2 = z+ 1/3.*z**2 - 1/9.*z**3 + 50./81*z**4
    parameters['M1'] = M1
    parameters['M2'] = M2
    parameters['L1'] = L1*a*constants.au/constants.Rsol
    parameters['L2'] = L2*a*constants.au/constants.Rsol
    
    #-- calculate the Keplerian orbits of the primary and secondary
    times = np.linspace(0.5*P,1.5*P,tres)
    r1,theta1 = keplerorbit.orbit_in_plane(times,[P,e,a1,gamma],component='primary')
    r2,theta2 = keplerorbit.orbit_in_plane(times,[P,e,a2,gamma],component='secondary')
    RV1 = keplerorbit.radial_velocity([P,e,a1,pi/2,view_angle,gamma],theta=theta1,itermax=8)
    RV2 = keplerorbit.radial_velocity([P,e,a2,pi/2,view_angle,gamma],theta=theta2,itermax=8)
    
    r1 = r1/constants.au
    r2 = r2/constants.au
    #   put them in cartesian coordinates
    x1o,y1o = r1*cos(theta1)/a,r1*sin(theta1)/a
    x2o,y2o = r2*cos(theta2)/a,r2*sin(theta2)/a
    #   and rotate towards the line of sight in the XZ plane
    rot_i = -(pi/2 - view_angle)
    x1o_,z1o = vectors.rotate(x1o,np.zeros_like(y1o),rot_i)
    x2o_,z2o = vectors.rotate(x2o,np.zeros_like(y2o),rot_i)
    
    #-- calculate separation at all phases
    ds = np.sqrt( (x1o-x2o)**2 + (y1o-y2o)**2)
    
    #-- calculate the polar radius and the radius towards L1 at minimum separation
    r_pole = newton(binary_roche_potential,1e-5,args=(0,0,Phi,q,ds.min(),F))
    r_pole2 = newton(binary_roche_potential,1e-5,args=(0,0,Phi2,q2,ds.min(),F2))
    parameters['Rp1'] = r_pole*a*constants.au/constants.Rsol
    parameters['Rp2'] = r_pole2*a*constants.au/constants.Rsol
    
    #-- calculate the critical Phis and Roche radius
    phi1_crit = binary_roche_potential(L1,0,0,Phi,q,ds.min(),F)
    phi2_crit = binary_roche_potential(L1,0,0,Phi2,q2,ds.min(),F2)
    R_roche1 = a * 0.49 * q2**(2./3.) / (0.6*q2**(2./3.) + np.log(1+q2**(1./3.)))
    R_roche2 = a * 0.49 * q**(2./3.) / (0.6*q**(2./3.) + np.log(1+q**(1./3.)))
    
    #-- timescales
    tdyn1 = np.sqrt(2*(r_pole*constants.au)**3/(constants.GG*M1*constants.Msol))
    tdyn2 = np.sqrt(2*(r_pole2*constants.au)**3/(constants.GG*M2*constants.Msol))
    lumt1 = (M1<2.) and M1**4 or (1.5*M1**3.5) # from mass luminosity relation
    lumt2 = (M2<2.) and M2**4 or (1.5*M2**3.5) # from mass luminosity relation
    tthr1 = constants.GG * (M1*constants.Msol)**2 / (r_pole*constants.au*lumt1*constants.Lsol)
    tthr2 = constants.GG * (M2*constants.Msol)**2 / (r_pole2*constants.au*lumt2*constants.Lsol)
    tnuc1 = 7e9 * M1 / lumt1
    tnuc2 = 7e9 * M2 / lumt2
    print '================================================='
    print "Period P                  =    %.3g d"%(P)
    print "Inclination angle i          = %.3g deg"%(incl)
    print "Systemic velocity gamma      = %.3g km/s"%(gamma)
    print "Eccentricity e               = %.3g"%(e)
    print "Mass ratio q                 = %.3g"%(q)
    print "Mass primary M1              = %.3g Msun"%(M1)
    print "Mass secondary M2            = %.3g Msun"%(M2)
    print "Polar Radius primary R1      = %.3g Rsun (Roche = %.3g)"%(r_pole*a*constants.au/constants.Rsol,R_roche1*constants.au/constants.Rsol)
    print "Polar Radius secondary R2    = %.3g Rsun (Roche = %.3g)"%(r_pole2*a*constants.au/constants.Rsol,R_roche2*constants.au/constants.Rsol)
    print "Center-of-mass               = %.3g Rsun"%(q/(1+q)*a*constants.au/constants.Rsol)
    print 'Lagrange L1                  = %.3g Rsun'%(L1*a*constants.au/constants.Rsol)
    print 'Lagrange L2                  = %.3g Rsun'%(L2*a*constants.au/constants.Rsol)
    print 'Phi primary                  = %.3g (critical = %.3g)'%(Phi,phi1_crit)
    print 'Phi secondary                = %.3g (critical = %.3g)'%(Phi2,phi2_crit)
    print 'Luminosity primary L1 (MS)   = %.3g Lsun'%(lumt1)
    print 'Luminosity secondary L2 (MS) = %.3g Lsun'%(lumt2)
    print 'System semi-major axis a     = %.3g Rsun'%(a*constants.au/constants.Rsol)
    print '------'
    print 'TDYN primary                 = %.3g hr'%(tdyn1/3600.)
    print 'TDYN secondary               = %.3g hr'%(tdyn2/3600.)
    print 'TNUC primary                 = %.3g yr'%(tnuc1)
    print 'TNUC secondary               = %.3g yr'%(tnuc2)
    print 'TTHERM primary               = %.3g yr'%(tthr1/(24*3600*365))
    print 'TTHERM secondary             = %.3g yr'%(tthr2/(24*3600*365))
    print '================================================='
    
    #-- construct the grid to calculate stellar shapes
    theta,phi = get_grid(res)
    thetas,phis = np.ravel(theta),np.ravel(phi)
    
    light_curve = np.zeros_like(times)
    RV1_corr = np.zeros_like(times)
    RV2_corr = np.zeros_like(times)
    to_SI = a*constants.au
    to_CGS = a*constants.au*100.
    
    if os.path.isfile(os.path.join(direc,'%s.fits'%(name))):
        os.remove(os.path.join(direc,'%s.fits'%(name)))
        print "Removed existing file"
    
    ext_dict = {}
    for di,d in enumerate(ds):
        print "STEP",di,
        
        #-- this is the angular velocity due to rotation and orbit
        #   you get the rotation period of the star via 2pi/omega_rot (in sec)
        omega_rot = F * 2*pi/P_ * 1/d**2 * sqrt( (1+e)*(1-e))
        omega_rot_vec = np.array([0.,0.,-omega_rot])
        
        if e>0 or di==0:
            #-- compute the star's radius and surface gravity
            out = [[get_binary_roche_radius(itheta,iphi,Phi=Phi,q=q,d=d,F=F,r_pole=r_pole),
                    get_binary_roche_radius(itheta,iphi,Phi=Phi2,q=q2,d=d,F=F2,r_pole=r_pole2)] for itheta,iphi in zip(thetas,phis)]
            rprim,rsec = np.array(out).T
            
            #-- for the primary
            #------------------
            radius  = rprim.reshape(theta.shape)
            this_r_pole = radius[0][0]
            x,y,z = vectors.spher2cart_coord(radius,phi,theta)
            g_pole = binary_roche_surface_gravity(0,0,this_r_pole*to_SI,d*to_SI,omega_rot,M1*constants.Msol,M2*constants.Msol,norm=True)
            Gamma_pole = binary_roche_potential_gradient(0,0,this_r_pole,q,d,F,norm=True)
            zeta = g_pole / Gamma_pole
            dOmega = binary_roche_potential_gradient(x,y,z,q,d,F,norm=False)
            grav_local = dOmega*zeta
            
            #-- here we can compute local quantities: surface gravity, area,
            #   effective temperature, flux and velocity
            grav_local = np.array([i.reshape(theta.shape) for i in grav_local])
            grav = vectors.norm(grav_local)
            areas_local,cos_gamma = surface_elements((radius,theta,phi),-grav_local)
            teff_local = local_temperature(grav,g_pole,T_pole,beta=1.)
            ints_local = local_intensity(teff_local,grav,np.ones_like(cos_gamma),photband='OPEN.BOL')
            velo_local = np.cross(np.array([x,y,z]).T*to_SI,omega_rot_vec).T
            
            #-- here we can compute the global quantities: total surface area
            #   and luminosity
            lumi_prim = 4*pi*(ints_local*areas_local*to_CGS**2).sum()/constants.Lsol_cgs
            area_prim = 4*areas_local.sum()*to_CGS**2/(4*pi*constants.Rsol_cgs**2)
            print ''
            print '  Polar Radius primary   = %.3g Rsun'%(this_r_pole*a*constants.au/constants.Rsol)
            print "  Polar logg primary     = %.3g dex"%(np.log10(g_pole*100))
            print "  Luminosity primary     = %.3g Lsun"%(lumi_prim)
            print "  Surface area primary   = %.3g Asun"%(area_prim)
            print "  Mean Temp primary      = %.3g K"%(np.average(teff_local,weights=areas_local))
            ext_dict['Rp1'] = this_r_pole*a*constants.au/constants.Rsol
            ext_dict['loggp1'] = np.log10(g_pole*100)
            ext_dict['LUMI1'] = lumi_prim
            ext_dict['SURF1'] = area_prim
                                    
            #-- for the secondary
            #--------------------
            radius2 = rsec.reshape(theta.shape)
            this_r_pole2 = radius2[0][0]
            x2,y2,z2 = vectors.spher2cart_coord(radius2,phi,theta)
            g_pole2 = binary_roche_surface_gravity(0,0,this_r_pole2*to_SI,d*to_SI,omega_rot,M2*constants.Msol,M1*constants.Msol,norm=True)
            Gamma_pole2 = binary_roche_potential_gradient(0,0,this_r_pole2,q2,d,F2,norm=True)
            zeta2 = g_pole2 / Gamma_pole2
            dOmega2 = binary_roche_potential_gradient(x2,y2,z2,q2,d,F2,norm=False)
            grav_local2 = dOmega2*zeta2
            
            #-- here we can compute local quantities: : surface gravity, area,
            #   effective temperature, flux and velocity  
            grav_local2 = np.array([i.reshape(theta.shape) for i in grav_local2])
            grav2 = vectors.norm(grav_local2)
            areas_local2,cos_gamma2 = surface_elements((radius2,theta,phi),-grav_local2)
            teff_local2 = local_temperature(grav2,g_pole2,T_pole2,beta=1.)
            ints_local2 = local_intensity(teff_local2,grav2,np.ones_like(cos_gamma2),photband='OPEN.BOL')
            velo_local2 = np.cross(np.array([x2,y2,z2]).T*to_SI,omega_rot_vec).T
            
            #-- here we can compute the global quantities: total surface area
            #   and luminosity
            lumi_sec = 4*pi*(ints_local2*areas_local2*to_CGS**2).sum()/constants.Lsol_cgs
            area_sec = 4*areas_local2.sum()*to_CGS**2/(4*pi*constants.Rsol_cgs**2)
            print '  Polar Radius secondary = %.3g Rsun'%(this_r_pole2*a*constants.au/constants.Rsol)
            print "  Polar logg secondary   = %.3g dex"%(np.log10(g_pole2*100))
            print "  Luminosity secondary   = %.3g Lsun"%(lumi_sec)
            print "  Surface area secondary = %.3g Asun"%(area_sec)
            print "  Mean Temp secondary    = %.3g K"%(np.average(teff_local2,weights=areas_local2))
            ext_dict['Rp2'] = this_r_pole2*a*constants.au/constants.Rsol
            ext_dict['loggp2'] = np.log10(g_pole2*100)
            ext_dict['LUMI2'] = lumi_sec
            ext_dict['SURF2'] = area_sec
            
            #================ START DEBUGGING PLOTS ===================
            #plot_quantities(phi,theta,np.log10(grav2*100.),areas_local2,np.arccos(cos_gamma2)/pi*180,teff_local2,ints_local2,
            #           names=['grav','area','angle','teff','ints'],rows=2,cols=3)
            #pl.show()
            #================   END DEBUGGING PLOTS ===================
                    
            #-- stitch the grid!
            theta_,phi_,radius,gravx,gravy,gravz,grav,areas,teff,ints,vx,vy,vz = \
                         stitch_grid(theta,phi,radius,grav_local[0],grav_local[1],grav_local[2],
                                    grav,areas_local,teff_local,ints_local,velo_local[0],velo_local[1],velo_local[2],seamless=False,
                                    vtype=['scalar','x','y','z','scalar','scalar','scalar','scalar','vx','vy','vz'])
            #-- stitch the grid!
            theta2_,phi2_,radius2,gravx2,gravy2,gravz2,grav2,areas2,teff2,ints2,vx2,vy2,vz2 = \
                         stitch_grid(theta,phi,radius2,grav_local2[0],grav_local2[1],grav_local2[2],
                                    grav2,areas_local2,teff_local2,ints_local2,velo_local2[0],velo_local2[1],velo_local2[2],seamless=False,
                                    vtype=['scalar','x','y','z','scalar','scalar','scalar','scalar','vx','vy','vz'])
            
            #-- vectors and coordinates in original frame
            x_of,y_of,z_of = vectors.spher2cart_coord(radius.ravel(),phi_.ravel(),theta_.ravel())
            x2_of,y2_of,z2_of = vectors.spher2cart_coord(radius2.ravel(),phi2_.ravel(),theta2_.ravel())
            x2_of = -x2_of            
            #-- store information on primary and secondary in a record array
            primary = np.rec.fromarrays([theta_.ravel(),phi_.ravel(),radius.ravel(),
                                         x_of,y_of,z_of,
                                         vx.ravel(),vy.ravel(),vz.ravel(),
                                         gravx.ravel(),gravy.ravel(),gravz.ravel(),grav.ravel(),
                                         areas.ravel(),teff.ravel(),ints.ravel()],
                                  names=['theta','phi','r',
                                         'x','y','z',
                                         'vx','vy','vz',
                                         'gravx','gravy','gravz','grav',
                                         'areas','teff','flux'])
            
            secondary = np.rec.fromarrays([theta2_.ravel(),phi2_.ravel(),radius2.ravel(),
                                         x2_of,y2_of,z2_of,
                                         vx2.ravel(),-vy2.ravel(),vz2.ravel(),
                                         -gravx2.ravel(),gravy2.ravel(),gravz2.ravel(),grav2.ravel(),
                                         areas2.ravel(),teff2.ravel(),ints2.ravel()],
                                  names=['theta','phi','r',
                                         'x','y','z',
                                         'vx','vy','vz',
                                         'gravx','gravy','gravz','grav',
                                         'areas','teff','flux'])
                        
            
            primary,secondary = reflection_effect(primary,secondary,theta,phi,A1=A1,A2=A2,max_iter=max_iter_reflection)

            
        #-- now compute the integrated intensity in the line of sight:
        #-------------------------------------------------------------
        #-- move the components along the orbit in the XY plane
        rot_theta = np.arctan2(y1o[di],x1o[di])
        x,y   = vectors.rotate(primary['x'],primary['y'],rot_theta,x0=x1o[di],y0=y1o[di])
        x2,y2 = vectors.rotate(secondary['x'],secondary['y'],rot_theta,x0=x2o[di],y0=y2o[di])
        gravx_,gravy_ = vectors.rotate(primary['gravx'],primary['gravy'],rot_theta)
        gravx2_,gravy2_ = vectors.rotate(secondary['gravx'],secondary['gravy'],rot_theta)
        vx_,vy_ = vectors.rotate(primary['vx'],primary['vy'],rot_theta)
        vx2_,vy2_ = vectors.rotate(secondary['vx'],secondary['vy'],rot_theta)
        
        #-- rotate towards the line of sight in the XZ plane
        rot_i = -(pi/2 - view_angle)
        x,z   = vectors.rotate(x,primary['z'],rot_i)
        x2,z2 = vectors.rotate(x2,secondary['z'],rot_i)
        gravx_,gravz_ = vectors.rotate(gravx_,primary['gravz'],rot_i)
        gravx2_,gravz2_ = vectors.rotate(gravx2_,secondary['gravz'],rot_i)
        vx_,vz_ = vectors.rotate(vx_,primary['vz'],rot_i)
        vx2_,vz2_ = vectors.rotate(vx2_,secondary['vz'],rot_i)            
        
        #-- we define the line of sight to be always x=1,y=0,z=0, and calculate
        #   and project the intensities in that direction
        line_of_sight = np.array([1.,0,0])
        proj_intens_prim,mus_prim = projected_intensity(primary['teff'],np.array([gravx_,gravy_,gravz_]),
                          primary['areas'],line_of_sight,photband=passband)
        proj_intens_secn,mus_secn = projected_intensity(secondary['teff'],np.array([gravx2_,gravy2_,gravz2_]),
                          secondary['areas'],line_of_sight,photband=passband)
        
        #-- project the coordinates onto the plane of the sky: this is simply the y,z coordinates
        visible_prim = -np.isnan(proj_intens_prim)
        visible_secn = -np.isnan(proj_intens_secn)
        #   only keep the meshpoints (and their areas) faced towards us.
        x_prim = x[visible_prim]
        y_prim = y[visible_prim]
        z_prim = z[visible_prim]
        a_prim = (primary['areas']*mus_prim)[visible_prim]
        x_secn = x2[visible_secn]
        y_secn = y2[visible_secn]
        z_secn = z2[visible_secn]
        a_secn = (secondary['areas']*mus_secn)[visible_secn]
                
        #-- the total intensity is simply the sum of the projected intensities
        #   over all visible meshpoints. To calculate the visibility, we
        #   we collect the Y-Z coordinates in one array for easy matching in
        #   the KDTree
        coords_prim = np.column_stack([y_prim,z_prim])
        coords_secn = np.column_stack([y_secn,z_secn])
        #   We need to know which star is in front. It is the one with the
        #   largest x coordinate
        if x_secn.min()<x_prim.min():
            coords_front = coords_prim
            coords_back  = coords_secn
            intens_front = proj_intens_prim[visible_prim]
            intens_back  = proj_intens_secn[visible_secn]
            areas_front = a_prim
            areas_back = a_secn
            velo_front = -vx_[visible_prim] + RV1[di]*1000
            velo_back  = -vx2_[visible_secn] + RV2[di]*1000
            front_cmap,back_cmap = pl.cm.hot,pl.cm.cool
            print 'Primary in front,',
            front_component = 1
        else:
            coords_front = coords_secn
            coords_back  = coords_prim
            intens_front = proj_intens_secn[visible_secn]
            intens_back  = proj_intens_prim[visible_prim]
            areas_front = a_secn
            areas_back = a_prim
            velo_front = -vx2_[visible_secn] + RV2[di]*1000
            velo_back  = -vx_[visible_prim] + RV1[di]*1000
            front_cmap,back_cmap = pl.cm.cool,pl.cm.hot
            front_component = 2
            print 'Secondary in front,',

        #   now find the coordinates of the front component closest to the
        #   the coordinates of the back component
        tree = KDTree(coords_front)
        distance,order = tree.query(coords_back)
        #   meshpoints of the back component inside an eclipse have a
        #   nearest neighbouring point in the (projected) front component
        #   which is closer than sqrt(area) of the surface element connected
        #   to that neighbouring point on the front component
        in_eclipse = distance < np.sqrt(areas_front[order])
        if np.sum(in_eclipse)>0:
            print 'during eclipse',
        else:
            print 'outside eclipse',
        
        #-- so now we can easily compute the total intensity as the sum of
        #   all visible meshpoints:
        total_intensity = intens_front.sum() + intens_back[-in_eclipse].sum()
        light_curve[di] = total_intensity
        print "---> Total intensity:",total_intensity
        
        intens_back[in_eclipse] = 0
        velo_back[in_eclipse] = 0
        
        #-- now calculate the *real* observed radial velocity and projected intensity
        proj_intens_prim[-visible_prim] = 0.
        proj_intens_secn[-visible_secn] = 0.
        proj_velo_prim = -vx_
        proj_velo_prim[-visible_prim] = 0.
        proj_velo_secn = -vy_
        proj_velo_secn[-visible_secn] = 0.
        if front_component==1:
            proj_intens_secn[visible_secn] = np.where(in_eclipse,np.nan,proj_intens_secn[visible_secn])
            RV1_corr[di] = np.average(velo_front/1000.,weights=intens_front)
            RV2_corr[di] = np.average(velo_back[-in_eclipse]/1000.,weights=intens_back[-in_eclipse])
            proj_velo_secn[visible_secn] = np.where(in_eclipse,np.nan,proj_velo_secn[visible_secn])
        else:
            proj_intens_prim[visible_prim] = np.where(in_eclipse,np.nan,proj_intens_prim[visible_prim])
            RV2_corr[di] = np.average(velo_front/1000.,weights=intens_front)
            RV1_corr[di] = np.average(velo_back[-in_eclipse]/1000.,weights=intens_back[-in_eclipse])
            proj_velo_prim[visible_prim] = np.where(in_eclipse,np.nan,proj_velo_prim[visible_prim])
        
        ext1_dict = parameters.copy()
        ext1_dict['extname'] = 'orbit'
        fits.write_array([times,light_curve,RV1,RV2,RV1_corr,RV2_corr,x1o_,y1o,z1o,x2o_,y2o,z2o],os.path.join(direc,'%s.fits'%(name)),
                        names=('time','flux','RV1','RV2','RVC1','RVC2','x1','y1','z1','x2','y2','z2'),
                        units=('d','n/a','a','a','a','a','a','a'),
                        header_dict=ext1_dict,ext=1)
        
        ext_dict['time'] = times[di]
        ext_dict['flux'] = total_intensity
        ext_dict['extname'] = 'PHS%05d'%(di)
        output = [x ,y ,z ,gravx_ ,gravy_ ,gravz_ ,primary['grav'],vx_,vy_,vz_,primary['areas']  ,primary['teff']  ,primary['flux']  ,proj_intens_prim,proj_velo_prim]
        output+= [x2,y2,z2,gravx2_,gravy2_,gravz2_,secondary['grav'],vx2_,vy2_,vz2_,secondary['areas'],secondary['teff'],secondary['flux'],proj_intens_secn,proj_velo_secn]
        fits.write_array(output,os.path.join(direc,'%s.fits'%(name)),
                       names=['x1','y1','z1','gravx1','gravy1','gravz1','grav1','vx1','vy1','vz1','areas1','teff1','flux1','projflux1','projv1',
                              'x2','y2','z2','gravx2','gravy2','gravz2','grav2','vx2','vy2','vz2','areas2','teff2','flux2','projflux2','projv2'],
                       header_dict=ext_dict)
        
        
        #================ START DEBUGGING PLOTS ===================
        pl.figure(figsize=(16,11))
        pl.subplot(221,aspect='equal');pl.title('line of sight intensity')
        pl.scatter(coords_front[:,0],coords_front[:,1],c=intens_front/areas_front,edgecolors='none',cmap=front_cmap)
        pl.scatter(coords_back[:,0],coords_back[:,1],c=intens_back/areas_back,edgecolors='none',cmap=back_cmap)
        pl.xlim(-1.,1.)
        pl.ylim(-0.5,0.5)
        pl.xlabel('X [semi-major axis]')
        pl.ylabel('Z [semi-major axis]')
        
        pl.subplot(222,aspect='equal');pl.title('line of sight velocity')
        vmin = min((vx_.min()/1000.+RV1.min()),(vx2_.min()/1000.+RV2.min()))
        vmax = max((vx_.max()/1000.+RV2.max()),(vx2_.max()/1000.+RV2.max()))
        pl.scatter(coords_front[:,0],coords_front[:,1],c=velo_front/1000.,edgecolors='none',cmap=pl.cm.RdBu_r,vmin=vmin,vmax=vmax)
        pl.scatter(coords_back[:,0],coords_back[:,1],c=velo_back/1000.,edgecolors='none',cmap=pl.cm.RdBu_r,vmin=vmin,vmax=vmax)
        cbar = pl.colorbar()
        cbar.set_label('Radial velocity [km/s]')
        pl.xlim(-1.,1.)
        pl.ylim(-0.5,0.5)
        pl.xlabel('X [semi-major axis]')
        pl.ylabel('Z [semi-major axis]')
        
        pl.subplot(223,aspect='equal');pl.title('Top view')
        try:
            pl.scatter(x[visible_prim],coords_front[:,0],c=intens_front/areas_front,edgecolors='none',cmap=front_cmap)
            pl.scatter(x2[visible_secn],coords_back[:,0],c=intens_back/areas_back,edgecolors='none',cmap=back_cmap)
        except:
            pl.scatter(x2[visible_secn],coords_front[:,0],c=intens_front,edgecolors='none',cmap=front_cmap)
            pl.scatter(x[visible_prim],coords_back[:,0],c=intens_back,edgecolors='none',cmap=back_cmap)
        pl.xlim(-1.,1.)
        pl.ylim(-1.,1.)
        pl.xlabel('X [semi-major axis]')
        pl.ylabel('Y [semi-major axis]')
        pl.subplot(224);pl.title('light curve and RV curve')
        ref_light = (lumi_prim+lumi_sec*constants.Lsol_cgs)/3
        pl.plot(times[:di+1],light_curve[:di+1],'k-',label=passband)
        pl.plot(times[di],light_curve[di],'ko',ms=10)
        pl.xlim(times.min(),times.max())
        pl.ylim(0,2*light_curve[0])
        pl.ylabel('Flux [erg/s/cm2/A/sr]')
        pl.legend(loc='lower left',prop=dict(size='small'))
        pl.twinx(pl.gca())
        
        pl.plot(times[:di+1],RV1_corr[:di+1],'b-',lw=2,label='Numeric 1')
        pl.plot(times[:di+1],RV1[:di+1]-gamma,'g--',lw=2,label='Kepler 1')
        pl.plot(times[di],RV1[di]-gamma,'go',ms=10)
        pl.plot(times[di],RV1_corr[di],'bo',ms=10)
        
        pl.plot(times[:di+1],RV2_corr[:di+1],'r-',lw=2,label='Numeric 2')
        pl.plot(times[:di+1],RV2[:di+1]-gamma,'c--',lw=2,label='Kepler 2')
        pl.plot(times[di],RV2[di]-gamma,'cs',ms=10)
        pl.plot(times[di],RV2_corr[di],'rs',ms=10)
        pl.ylim(1.2*min(min(RV1-gamma),min(RV2-gamma)),+1.2*max(max(RV1-gamma),max(RV2-gamma)))
        pl.xlim(times.min(),times.max())
        pl.ylabel('Radial velocity [km/s]')
        pl.xlabel('Time [d]')
        pl.legend(loc='upper left',prop=dict(size='small'))
        pl.savefig(os.path.join(direc,'%s_los_%04d'%(name,di)))
        pl.close()
        
        pl.figure(figsize=(8,2.0/2.4*8))
        ax = pl.axes([0,0,1,1])
        ax.set_aspect('equal')
        ax.set_axis_bgcolor('k')
        pl.xticks([])
        pl.yticks([])
        if front_component==1:
            sfront = 4
            sback = 16
        else:
            sfront = 16
            sback = 4
        pl.scatter(coords_front[:,0],coords_front[:,1],s=sfront,c=intens_front/areas_front,edgecolors='none',cmap=pl.cm.gray)
        pl.scatter(coords_back[:,0],coords_back[:,1],s=sback,c=intens_back/areas_back,edgecolors='none',cmap=pl.cm.gray)
        pl.xlim(-1.2,1.2)
        pl.ylim(-1.0,1.0)
        pl.savefig(os.path.join(direc,'%s_image_%04d'%(name,di)),facecolor='k')
        pl.close()
        #================   END DEBUGGING PLOTS ===================


if __name__=="__main__":
    import doctest
    doctest.testmod()
    pl.show()
    
    import sys
    #binary_test()
    #full_test_diff_rot()
    #full_test_binary()
    