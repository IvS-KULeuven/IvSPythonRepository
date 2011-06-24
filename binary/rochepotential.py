"""
Compute stellar shapes via Roche potential and gradients.

The premisse of the Roche potential is that the stellar mass can be represented
by a point source.

Currently, three types of stellar shapes are implemented:
    1. asynchronously rotating eccentric binary (L{get_binary_roche_radius})
    2. fast rotating star (L{get_fastrot_roche_radius})
    3. differentially (fast) rotating star (L{get_diffrot_roche_surface_gravity})

This module can be used to calculate the following information:
    1.  the distorted shape of the star due to a Roche potential
    2.  the local surface gravity assuming Von Zeipel gravity darkening
    3.  the L{local effective temperature<local_temperature>}
    4.  the local velocity vector due to rotation
    5.  the radial velocity
    6.  the total distorted surface area
    7.  the total distorted luminosity
    8.  the L{projected intensity<projected_intensity>} in the line of sight
    9.  the mean radial velocity (e.g. Rossiter effect in binaries)
    10. a L{synthetic spectrum<spectral_synthesis>} starting from a library of spectra
    11. a L{synthetic light curve<binary_light_curve_synthesis>} from the projected intensities

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

>>> proj_intens,mus = projected_intensity(teff_local,grav_local,areas_local,np.array([0,-sin(view_angle),-cos(view_angle)]),photband='OPEN.BOL')

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

>>> areas_local,cos_gamma = surface_elements((radius,theta,phi),-grav_local)

Now we can calculate the local effective temperature (assuming radiative
atmosphere in Von Zeipel law), and bolometric intensities:

>>> teff_local = local_temperature(vectors.norm(grav_local),g_pole,T_pole,beta=1.)
>>> ints_local = local_intensity(teff_local,grav,photband='OPEN.BOL')

Compute the angles between surface normals and line-of-sight. We do this in 1D:

>>> gravity = np.array([igrav.ravel() for igrav in grav_local])
>>> line_of_sight = np.array([ilos.ravel() for ilos in line_of_sight])
>>> angles = vectors.angle(-gravity,line_of_sight)
>>> mus = cos(angles)

Now compute the bolometric projected intensity, taking care of limbdarkening
effects:

>>> intens_proj = local_intensity(teff_local.ravel(),grav.ravel(),mus,photband='OPEN.BOL').reshape(theta.shape)

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


Section 4. Binary star
======================

As an example we take the binary star SX Aurigae. This is the minimum set of
parameters.

>>> P = 1.2100802        # period in days
>>> e = 0.               # eccentricity
>>> q = 0.54369          # mass ratio M2/M1
>>> asini = 11.9*constants.Rsol/constants.au # semi-major axis
>>> i = 81.27            # inclination angle
>>> F = 1.               # synchronicity parameter
>>> d = 1.               # separation in semi-major axis units
>>> Phi1 = 3.            # gravitational potential primary
>>> Phi2 = 5.05          # gravitational potential secondary
>>> T_pole1 = 25000.     # polar temperature primary
>>> T_pole2 = 18850.     # polar temperature secondary
>>> M1 = 10.3            # primary mass
>>> M2 = 5.6             # secondary mass

Note that we actually do not need the masses, since we can derive them from
the third kepler law with the other parameters. We do so for convenience.

Everything is calculated in units of semi-major axis. Thus, to convert to SI
units we need the following conversion factor:

>>> a = asini * sin(i/180.*pi)
>>> to_SI = a*constants.au

We need some constants for the calculations: the polar radii and angular velocity

>>> r_pole1 = newton(binary_roche_potential,1e-5,args=(0,0,Phi1,  q,d,F))
>>> r_pole2 = newton(binary_roche_potential,1e-5,args=(0,0,Phi2,1/q,d,F))
>>> omega_rot = F * 2*pi/(P*24*3600) * 1/d**2 * sqrt( (1+e)*(1-e))
>>> omega_rot_vec = np.array([0.,0.,-omega_rot])

Derive the shape of the two stars

>>> radius1 = np.array([get_binary_roche_radius(itheta,iphi,Phi=Phi1,q=  q,d=d,F=F,r_pole=r_pole1) for itheta,iphi in zip(thetas,phis)]).reshape(theta.shape)
>>> radius2 = np.array([get_binary_roche_radius(itheta,iphi,Phi=Phi2,q=1/q,d=d,F=F,r_pole=r_pole2) for itheta,iphi in zip(thetas,phis)]).reshape(theta.shape)

We focus on the primary, then repeat everything for the secondary: The local
surface gravity can only be calculated if we have Cartesian coordinates.

>>> x1,y1,z1 = vectors.spher2cart_coord(radius1,phi,theta)
>>> g_pole1 = binary_roche_surface_gravity(0,0,r_pole1*to_SI,d*to_SI,omega_rot,M1*constants.Msol,M2*constants.Msol,norm=True)
>>> Gamma_pole1 = binary_roche_potential_gradient(0,0,r_pole1,q,d,F,norm=True)
>>> zeta1 = g_pole1 / Gamma_pole1
>>> dOmega1 = binary_roche_potential_gradient(x1,y1,z1,q,d,F,norm=False)
>>> grav_local1 = dOmega1*zeta1
>>> grav_local1 = np.array([i.reshape(theta.shape) for i in grav_local1])
>>> grav1 = vectors.norm(grav_local1)

Now we can, as before, calculate the other local quantities:

>>> areas_local1,cos_gamma1 = surface_elements((radius1,theta,phi),-grav_local1)
>>> teff_local1 = local_temperature(grav1,g_pole1,T_pole1,beta=1.)
>>> ints_local1 = local_intensity(teff_local1,grav1,np.ones_like(cos_gamma1),photband='OPEN.BOL')
>>> velo_local1 = np.cross(np.array([x1,y1,z1]).T*to_SI,omega_rot_vec).T


Now make some plots showing the local quantities:

>>> quantities = areas_local1,np.log10(grav1*100),teff_local1,ints_local1
>>> names = 'Area','log g', 'Teff', 'Flux'
>>> p = pl.figure()
>>> rows,cols = 2,2    
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

]include figure]]ivs_binary_rochepotential_binarystar.png]

>>> p = pl.figure()
>>> p = pl.subplot(121,aspect='equal')
>>> p = pl.title('top view')
>>> p = pl.scatter(x1,y1,c=teff_local1.ravel(),edgecolors='none')
>>> p = pl.subplot(122,aspect='equal')
>>> p = pl.title('side view')
>>> p = pl.scatter(x1,z1,c=teff_local1.ravel(),edgecolors='none')

]include figure]]ivs_binary_rochepotential_binarystar_shape.png]

Section 5. Synthesizing spectra
===============================

In this case study, we calculate the projected radial velocity of a uniformly
rotating star, but in the limit of no deformation due to the rotation. We do this,
so that we can easily compare rotational broadening calculated with the ROTIN
package, and numerically. So let us build a star with these parameters:

>>> omega = 0.0        # ratio to critical velocity
>>> T_pole = 13000.    # K
>>> r_pole = 2.0       # solar radii
>>> M = 1.5            # solar mass
>>> view_angle = pi/2  # radians
>>> theta,phi = get_grid(20,100,full=True,gtype='spher')
>>> thetas,phis = np.ravel(theta),np.ravel(phi)

Then calculate the shape of this star

>>> radius = np.array([get_fastrot_roche_radius(itheta,r_pole,omega) for itheta in thetas]).reshape(theta.shape)
>>> grav_local = np.array([fastrot_roche_surface_gravity(iradius,itheta,iphi,r_pole,omega,M) for iradius,itheta,iphi in zip(radius.ravel(),thetas,phis)]).T
>>> grav_local = np.array([i.reshape(theta.shape) for i in grav_local])
>>> g_pole = fastrot_roche_surface_gravity(r_pole,0,0,r_pole,omega,M)[-1]
>>> grav = vectors.norm(grav_local)

and the local quantities

>>> areas_local,cos_gamma = surface_elements((radius,theta,phi),-grav_local)
>>> teff_local = local_temperature(vectors.norm(grav_local),g_pole,T_pole,beta=1.)
>>> ints_local = local_intensity(teff_local,grav,photband='OPEN.BOL')
>>> x,y,z = vectors.spher2cart_coord(radius.ravel(),phis,thetas)

Assume, with a shape of a non-rotating star, that we have a velocity component
on the surface of the star:

>>> myomega = 0.5
>>> velo_local = diffrot_velocity((phi,theta,radius*constants.Rsol),myomega,myomega,r_pole,M)

Collect all the necessary information in one record array.

>>> starlist = [x,y,z] + [i.ravel() for i in velo_local] + [i.ravel() for i in grav_local] + [teff_local.ravel(),areas_local.ravel()]
>>> starnames = ['x','y','z','vx','vy','vz','gravx','gravy','gravz','teff','areas']
>>> star = np.rec.array(starlist,names=starnames)

Project the star in some line-of-sight. The velocity component in the X-direction
is the radial velocity.

>>> view_angle = pi/2 # edge on
>>> mystar = project(star,view_long=(0,0,0),view_lat=(view_angle,0,0),photband='OPEN.BOL',only_visible=True,plot_sort=True)

We can calculate the synthetic spectra for all surface elements between 7055 and
7075 angstrom:

>>> loggs = np.log10(np.sqrt(mystar['gravx']**2 + mystar['gravy']**2 + mystar['gravz']**2)*100)
>>> iterator = zip(mystar['teff'],loggs,mystar['vx']/1000.)
>>> wave_spec = np.linspace(7055,7075,750)
>>> spectra = np.array([spectra_model.get_table(teff=iteff,logg=ilogg,vrad=ivrad,wave=wave_spec)[1] for iteff,ilogg,ivrad in iterator])

The total observed spectrum is then simply the weighted sum with the local projected
intensities:

>>> average_spectrum = np.average(spectra,weights=mystar['projflux'],axis=0)

We compare with the ROTIN package, which assumes a linear limbdarkening law

>>> original = spectra_model.get_table(teff=mystar['teff'][0],logg=loggs[0],wave=wave_spec)[1]
>>> rotin1 = spectra_model.get_table(teff=mystar['teff'][0],logg=loggs[0],vrot=vmax/1000.*sin(view_angle),wave=wave_spec,fwhm=0,epsilon=0.6,stepr=-1)[1]
    
And a plot can be made via

>>> colors = mystar['eyeflux']/mystar['eyeflux'].max()
>>> p = pl.figure(figsize=(7,7))
>>> ax = pl.axes([0,0.5,1,0.5])
>>> p = ax.set_aspect('equal')
>>> p = ax.set_axis_bgcolor('k')
>>> p = pl.box(on=False)
>>> p = pl.xticks([]);p = pl.yticks([])
>>> p = pl.scatter(mystar['y'],mystar['z'],c=colors,edgecolors='none',cmap=pl.cm.gray,vmin=0,vmax=1)
>>> p = pl.scatter(mystar['y']+1.02*mystar['y'].ptp(),mystar['z'],c=mystar['vx'],edgecolors='none',vmin=vmin,vmax=vmax,cmap=pl.cm.RdBu)
>>> ax = pl.axes([0.13,0.1,0.78,0.4])
>>> p = pl.plot(wave_spec,average_spectrum,'k-',lw=2,label='Numerical')
>>> p = pl.plot(wave_spec,original,'r-',label='No rotation')
>>> p = pl.plot(wave_spec,rotin1,'b-',label='Rotin v=%.1f km/s'%(vmax/1000.),lw=2)
>>> p = pl.ylim(0.95,1.0)
>>> p = pl.legend(loc='lower right',prop=dict(size='small'))
>>> p = pl.xlabel('Wavelength [angstrom]',color='1')
>>> p = pl.ylabel('Normalised flux',color='1')

]include figure]]ivs_binary_rochepotential_norot_spectrum_b.png]

]include figure]]ivs_binary_rochepotential_norot_spectrum_a.png]


"""
import logging
import os
import pylab as pl
import numpy as np
from numpy import pi,cos,sin,sqrt,nan
from scipy.optimize import newton
from scipy.spatial import KDTree
try:
    from scipy.spatial import Delaunay
except ImportError:
    print 'No scipy >= 0.9, falling back on scipy < 0.9: no Delaunay triangulations available'
import time
from ivs.binary import keplerorbit
from ivs.units import constants
from ivs.units import vectors
from ivs.units import conversions
from ivs.sed import model as sed_model
from ivs.sed import limbdark
from ivs.spectra import model as spectra_model
from ivs.misc import loggers
from ivs.io import ascii
from ivs.io import fits

logger = logging.getLogger("BIN.ROCHE")

#{ Eccentric asynchronous binary Roche potential in spherical coordinates

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
    
    term1 = - constants.GG*M1/vectors.norm(r)**3*r
    term2 = - constants.GG*M2/vectors.norm(h)**3*h
    term3 = - omega**2 * d_cf
    
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
    
def critical_angular_velocity(M,R_pole,units='Hz'):
    """
    Compute the critical angular velocity (Hz).
        
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

def critical_velocity(M,R_pole,units='km/s'):
    """
    Compute the critical velocity (km/s)
    
    @param M: mass (solar masses)
    @type M: float
    @param R_pole: polar radius (solar radii)
    @type R_pole: float
    @param units: if you wish other units than Hz, you can give them here
    @type units: string understandable by L{units.conversions.convert}
    @return: critical velocity in km/s
    @rtype: float
    """
    omega_crit = critical_angular_velocity(M,R_pole)
    P = 1./omega_crit
    R_eq = get_fastrot_roche_radius(pi/2,R_pole,omega_crit)*constants.Rsol
    veq = 2*pi*R_eq/P
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
        r = nan    
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
    @type coordinates: 3xNxM array
    @param omega_eq: equatorial angular velocity (as a fraction of the critical one)
    @type omega_eq: float
    @param omega_pole: polar angular velocity (as a fraction of the critical one)
    @type omega_pole: float
    @param R_pole: polar radius in solar radii
    @type R_pole: float
    @param M: stellar mass in solar mass
    @type M: float
    @return: velocity vectors of all surface elements
    @rtype 3xNxM array
    """
    #-- angular velocity of every surface element
    Omega_crit = critical_angular_velocity(M,R_pole)
    phi,theta,radius = coordinates
    omega_local = diffrot_law(omega_eq,omega_pole,theta)*Omega_crit
    #-- direction of local angular velocity in Cartesian coordinates (directed in upwards z)
    omega_local_vec = np.array([np.zeros_like(omega_local.ravel()),np.zeros_like(omega_local.ravel()),omega_local.ravel()]).T
    
    x,y,z = vectors.spher2cart_coord(radius,phi,theta)
    surface_element = np.array([x.ravel(),y.ravel(),z.ravel()]).T

    velo_local = np.array([np.cross(ielement,iomega_local_vec) for ielement,iomega_local_vec in zip(surface_element,omega_local_vec)]).T.ravel()
    return velo_local.reshape((3,theta.shape[0],theta.shape[1]))

#}

#{ Derivation of local quantities

def surface_elements((r,mygrid),(surfnormal_x,surfnormal_y,surfnormal_z),gtype='spher'):
    """
    Compute surface area of elements in a grid.
    
    theta,phi must be generated like mgrid(theta_range,phi_range)
    
    usually, the surfnormals are acquired via differentiation of a gravity potential,
    and is then equal to the *negative* of the local surface gravity.
    """
    theta,phi = mygrid[:2]
    if gtype=='spher':
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
    
    elif gtype=='delaunay':
        #-- compute the angle between the surface normal and the radius vector
        x,y,z = vectors.spher2cart_coord(r,phi,theta)
        a = np.array([x,y,z])
        b = np.array([surfnormal_x,surfnormal_y,surfnormal_z])        
        cos_gamma = vectors.cos_angle(a,b)
        
        delaunay_grid = mygrid[2]

        sizes = np.zeros(len(delaunay_grid.convex_hull))
        
        for i,indices in enumerate(delaunay_grid.convex_hull):
            a = sqrt((x[indices[0]]-x[indices[1]])**2 + (y[indices[0]]-y[indices[1]])**2 + (z[indices[0]]-z[indices[1]])**2)
            b = sqrt((x[indices[0]]-x[indices[2]])**2 + (y[indices[0]]-y[indices[2]])**2 + (z[indices[0]]-z[indices[2]])**2)
            c = sqrt((x[indices[1]]-x[indices[2]])**2 + (y[indices[1]]-y[indices[2]])**2 + (z[indices[1]]-z[indices[2]])**2)
            s = 0.5*(a+b+c)
            sizes[i] = sqrt( s*(s-a)*(s-b)*(s-c))
        
        return sizes, cos_gamma

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
    intens = local_intensity(teff,grav_,mu=mus,photband=photband)
    #-- intensity is less if the surface element area is small (and the area
    #   needs to be projected!!)
    return intens*areas*mus,mus


def project(star,view_long=(0,0,0),view_lat=(pi/2,0,0),photband='OPEN.BOL',
            only_visible=False,plot_sort=False):
    """
    Project and transform coordinates and vectors to align with the line-of-sight.
    
    Parameter C{star} should be a record array containing fields 'teff','gravx',
    'gravy','gravz','areas','vx','vy','vz'
    
    and either you suply ('r','theta','phi') or ('x','y','z')
    
    The XY direction is then the line-of-sight, and the YZ plane is the plane
    of the sky.
    
    An extra column 'projflux' and 'eyeflux' will be added. Projected flux
    takes care of limb darkening, and projected surface area. Eye flux only
    takes care of limbdarkening, and should only be used for plotting reasons.
    
    view_long[0] of 0 means looking in the XY line, pi means looking in the YX line.
    view_lat[0] of pi/2 means edge on, 0 or pi is pole-on.
    
    This function updates all Cartesian coordinates present in the star, but not
    the polar coordinates! The projected fluxes are added as a field 'projflux'
    to the returned record array.
    
    If you set 'only_visible' to True, only the information on the visible parts
    of the star will be contained.
    
    If you set 'plot_sort' to True, the arrays will be returned in a sorted order,
    where the areas at the back come first. This is especially handy for plotting.
    
    @parameters star: record array containing all necessary information on the
    star
    @type star: numpy record array
    @parameter view_long: longitude viewing angle (radians) and coordinate zeropoint
    @type view_long: tuple floats  (radians,x,y)
    @parameter view_lat: inclination viewing angle (radians) and coordinate zeropoint
    @type view_lat: tuple floats (radians,x,z)
    @parameter photband: photometric passband
    @type photband: string
    @parameter only_visible: flag to return only information on visible surface elements
    @type only_visible: boolean
    @parameter plot_sort: flag to sort the surface elements from back to front
    @type plot_sort: boolean
    """
    myshape = star['gravx'].shape
    gravx,gravy,gravz = np.array([star['gravx'].ravel(),star['gravy'].ravel(),star['gravz'].ravel()])
    areas = star['areas'].ravel()
    teff = star['teff'].ravel()
    vx,vy,vz = star['vx'].ravel(),star['vy'].ravel(),star['vz'].ravel()
    #-- if 'x' is not in the star's record array, we assume the polar coordinates
    #   are in there and convert them to Cartesian coordinates
    if not 'x' in star.dtype.names:
        x,y,z = vectors.spher2cart_coord(star['r'].ravel(),star['phi'].ravel(),star['theta'].ravel())
    else:
        x,y,z = star['x'].ravel(),star['y'].ravel(),star['z'].ravel(),
    #-- first we rotate in the XY plane (only for surface coordinates is the
    #   coordinate zeropoint important, the rest are vectors!):
    if view_long[0]!=0:
        x,y = vectors.rotate(x,y,view_long[0],x0=view_long[1],y0=view_long[2])
        gravx,gravy = vectors.rotate(gravx,gravy,view_long[0])
        vx,vy = vectors.rotate(vx,vy,view_long[0])
    
    #-- then we rotate in the YZ plane:
    if view_lat[0]!=pi/2:
        rot_i = -(pi/2 - view_lat[0])
        x,z = vectors.rotate(x,z,rot_i)
        gravx,gravz = vectors.rotate(gravx,gravz,rot_i)
        vx,vz = vectors.rotate(vx,vz,rot_i)
        
    #-- ... and project the fluxes in the line of sight, which is now in the XY
    #   direction:
    view_vector = np.array([1.,0,0])#np.array([-sin(pi/2),0,-cos(pi/2)])
    grav_local = np.array([gravx,gravy,gravz])
    proj_flux,mus = projected_intensity(teff,grav_local,areas,view_vector,photband=photband)
    
    #-- we now construct a copy of the star record array with the changed
    #   coordinates
    new_star = star.copy()
    new_star['gravx'],new_star['gravy'],new_star['gravz'] = gravx,gravy,gravz
    new_star['vx'],new_star['vy'],new_star['vz'] = vx,vy,vz
    if 'x' in star.dtype.names:
        new_star['x'],new_star['y'],new_star['z'] = x,y,z
    else:
        new_star = pl.mlab.rec_append_fields(new_star,'x',x)
        new_star = pl.mlab.rec_append_fields(new_star,'y',x)
        new_star = pl.mlab.rec_append_fields(new_star,'z',x)
    new_star = pl.mlab.rec_append_fields(new_star,'projflux',proj_flux)
    new_star = pl.mlab.rec_append_fields(new_star,'eyeflux',proj_flux/areas/mus)
    
    #-- clip visible areas and sort in plotting order if necessary
    if only_visible:
        new_star = new_star[-np.isnan(new_star['projflux'])]
    if plot_sort:
        new_star = new_star[np.argsort(new_star['x'])]
    
    return new_star




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

#}

#{ Gridding functions

def get_grid(*args,**kwargs):
    """
    Construct a coordinate grid
    
    If you give two resolution, the first is for theta, the second for phi
    """
    gtype = kwargs.get('gtype','spherical')
    full = kwargs.get('full',False)
    
    if 'spher' in gtype.lower():
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution
        
        dtheta = pi/res1 # step size coordinate 1
        dphi = pi/res2     # step size coordinate 2
        
        #-- full grid or only one quadrant
        if full:
            theta0,thetan = dtheta, pi-dtheta
            phi0,phin = 0,2*pi-2*dphi
        else:
            theta0,thetan = dtheta/4,pi/2-dtheta/4
            phi0,phin = 0,pi-dphi
            
        theta,phi = np.mgrid[theta0:thetan:res1*1j,phi0:phin:res2*1j]
        return theta,phi
    
    if 'cil' in gtype.lower():
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution
        
        dsintheta = 2./res1 # step size sin(theta)
        dphi = pi/res2     # step size coordinate 2
        
        #-- full grid or only one quadrant
        if full:
            sintheta0,sinthetan = -1+dsintheta, 1.-dsintheta
            phi0,phin = 0,2*pi-dphi
        else:
            sintheta0,sinthetan = dsintheta,1-dsintheta
            phi0,phin = 0,pi-dphi
            
        sintheta,phi = np.mgrid[sintheta0:sinthetan:res1*1j,phi0:phin:res2*1j]
        return np.arccos(sintheta),phi
    
    if 'delaunay' in gtype.lower():
        #-- resolution doesn't matter that much anymore 
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution
        
        theta = np.random.uniform(low=0,high=pi,size=res1*res2)
        phi = np.random.uniform(low=0,high=2*pi,size=res2*res1)
        
        x = sin(theta)*sin(phi)
        y = sin(theta)*cos(phi)
        z = cos(theta)
        
        points = np.array([x,y,z]).T
        grid = Delaunay(points)
        centers = np.zeros((len(grid.convex_hull),3))
        for i,indices in enumerate(grid.convex_hull):
            centers[i] = [x[indices].sum()/3,y[indices].sum()/3,z[indices].sum()/3]
        theta,phi = np.arccos(centers[:,2]),np.arctan2(centers[:,1],centers[:,0])
        return theta,phi, grid
        
        
        
        
        
    

def stitch_grid(theta,phi,*quant,**kwargs):
    """
    Stitch a grid together that was originally defined on 1 quadrant.
    
    We add the three other quandrants.
    """
    seamless = kwargs.get('seamless',False)
    vtype = kwargs.get('vtype',['scalar' for i in quant])
    gtype = kwargs.get('gtype','spher')
    ravel = kwargs.get('ravel',False)
    
    if gtype == 'spher':
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
    else:
        out = [theta,phi] + list(quant)
        
    return out

#}
def spectral_synthesis(*stars,**kwargs):
    """
    Generate a synthetic spectrum of one or more stars.
    
    If you give more than one star, you get more than one synthesized spectrum
    back. To compute the total spectrum, just add them up.
    
    WARNING: the spectra are scaled with the value of the projected intensity.
    This is usually calculated within a specific photometric passband. The total
    spectrum will thus be dependent on the photometric passband, unless you
    specified 'OPEN.BOL' (which is the bolometric open filter). If you want to
    know the RV for a specific line, you might want to look into calculating the
    projected intensities at specific wavelength intervals.
    
    @param stars: any number of (projected) star record arrays
    @type stars: tuple of star record arrays
    @keyword wave: wavelength template to compute the synthetic spectrum on
    @type wave: numpy array
    @return: tuple of synthetic spectra, scaled to total intensity
    @rtype: tuple of numpy arrays
    """
    wave = kwargs.get('wave',np.linspace(7050,7080,1000))
    #-- these are the log surface gravity values in cgs [dex]
    spectra = []
    get_spectrum = spectra_model.get_table
    for mystar in stars:
        loggs = np.log10(np.sqrt(mystar['gravx']**2 + mystar['gravy']**2 + mystar['gravz']**2)*100)
        iterator = zip(mystar['teff'],loggs,mystar['vx']/1000.)
        print "Temperature range:",mystar['teff'].min(),mystar['teff'].max()
        print "logg range:",loggs.min(),loggs.max()
        #-- retrieve the synthetic spectra for all surface elements, interpolated
        #   on the supplied wavelength grid, and taking care of radial velocity shifts
        spectra_elements = np.array([get_spectrum(teff=iteff,logg=ilogg,vrad=ivrad,wave=wave)[1] for iteff,ilogg,ivrad in iterator])
        #-- and add them up according to the flux projected in the line of sight
        average_spectrum = np.sum(spectra_elements*mystar['projflux'],axis=0)
        spectra.append(average_spectrum)
    return spectra


def binary_light_curve_synthesis(**parameters):
    """
    Generate a synthetic light curve of a binary system.
    
    @keyword name: name of the binary system, used for output
    @type name: string
    @keyword direc: directory to put output files
    @type direc: string (if empty string, then current directory will be taken,
    if None, then not output will be written)
    @parameter gres: grid resolution for primary and secondary. If integer, the
    resolution of both components and both longitude and latitude will be the
    same. A 2-tuple will be interpreted as (latitude resolution, longitude resolution)
    @type gres: integer, 2-tuple or 4-tuple
    @parameter tres: number of phase steps to comptue the light curve on
    @type tres: integer
    """
    #-- some parameters are optional
    #   file output parameters
    name = parameters.setdefault('name','mybinary')
    direc = parameters.setdefault('direc','')
    #   calculation details
    res = parameters.setdefault('gres',20)                    # resolution of the grid
    gtype = parameters.setdefault('gtype','spher')
    tres= parameters.setdefault('tres',125)                   # resolution of the phase diagram
    photband = parameters.setdefault('photband','JOHNSON.V')  # photometric passband
    max_iter_reflection = parameters.setdefault('ref_iter',1) # maximum number of iterations of reflection effect
    #   orbital parameters
    gamma = parameters.setdefault('gamma',0.)            # systemic velocity [km/s]
    incl = parameters.setdefault('incl',90.)             # system inclination angle [deg]
    q = parameters.setdefault('q',1.)                    # mass ratio (M2/M1)
    e = parameters.setdefault('e',0.)                    # eccentricity
    omega = parameters.setdefault('omega',0.)            # argument of periastron
    #   component parameters
    F = parameters.setdefault('F1',sqrt((1+e)/(1-e)**3)) # synchronicity parameter primary
    F2= parameters.setdefault('F2',sqrt((1+e)/(1-e)**3)) # synchronicity parameter secondary
    A1= parameters.setdefault('A1',1.)                   # albedo primary
    A2= parameters.setdefault('A2',1.)                   # albedo secondary
    beta1 = parameters.setdefault('beta1',1.0)           # gravity darkening primary
    beta2 = parameters.setdefault('beta2',1.0)           # gravity darkening secondary
    
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
    
    logger.info('=================================================')
    logger.info('GENERAL SYSTEM AND COMPONENT PROPERTIES')
    logger.info('=================================================')
    logger.info("Period P                     = %.3g d"%(P))
    logger.info("Inclination angle i          = %.3g deg"%(incl))
    logger.info("Systemic velocity gamma      = %.3g km/s"%(gamma))
    logger.info("Eccentricity e               = %.3g"%(e))
    logger.info("Mass ratio q                 = %.3g"%(q))
    logger.info("Mass primary M1              = %.3g Msun"%(M1))
    logger.info("Mass secondary M2            = %.3g Msun"%(M2))
    logger.info("Polar Radius primary R1      = %.3g Rsun (Roche = %.3g)"%(r_pole*a*constants.au/constants.Rsol,R_roche1*constants.au/constants.Rsol))
    logger.info("Polar Radius secondary R2    = %.3g Rsun (Roche = %.3g)"%(r_pole2*a*constants.au/constants.Rsol,R_roche2*constants.au/constants.Rsol))
    logger.info("Center-of-mass               = %.3g Rsun"%(q/(1+q)*a*constants.au/constants.Rsol))
    logger.info('Lagrange L1                  = %.3g Rsun'%(L1*a*constants.au/constants.Rsol))
    logger.info('Lagrange L2                  = %.3g Rsun'%(L2*a*constants.au/constants.Rsol))
    logger.info('Phi primary                  = %.3g (critical = %.3g)'%(Phi,phi1_crit))
    logger.info('Phi secondary                = %.3g (critical = %.3g)'%(Phi2,phi2_crit))
    logger.info('Luminosity primary L1 (MS)   = %.3g Lsun'%(lumt1))
    logger.info('Luminosity secondary L2 (MS) = %.3g Lsun'%(lumt2))
    logger.info('System semi-major axis a     = %.3g Rsun'%(a*constants.au/constants.Rsol))
    logger.info('------')
    logger.info('TDYN primary                 = %.3g hr'%(tdyn1/3600.))
    logger.info('TDYN secondary               = %.3g hr'%(tdyn2/3600.))
    logger.info('TNUC primary                 = %.3g yr'%(tnuc1))
    logger.info('TNUC secondary               = %.3g yr'%(tnuc2))
    logger.info('TTHERM primary               = %.3g yr'%(tthr1/(24*3600*365)))
    logger.info('TTHERM secondary             = %.3g yr'%(tthr2/(24*3600*365)))
    logger.info('=================================================')
    
    #-- construct the grid to calculate stellar shapes
    if hasattr(res,'__iter__'):
        mygrid = get_grid(res[0],res[1],gtype=gtype)
        theta,phi = mygrid[:2]
    else:
        mygrid = get_grid(res,gtype=gtype)
        theta,phi = mygrid[:2]
    thetas,phis = np.ravel(theta),np.ravel(phi)
    
    light_curve = np.zeros_like(times)
    RV1_corr = np.zeros_like(times)
    RV2_corr = np.zeros_like(times)
    to_SI = a*constants.au
    to_CGS = a*constants.au*100.
    
    fitsfile = os.path.join(direc,'%s.fits'%(name))
    if direc is not None and os.path.isfile(fitsfile):
        os.remove(fitsfile)
        logger.warning("Removed existing file %s"%(fitsfile))
    
    ext_dict = {}
    for di,d in enumerate(ds):
        report = "STEP %04d"%(di)
        
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
            this_r_pole = radius.ravel()[0]
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
            areas_local,cos_gamma = surface_elements((radius,mygrid),-grav_local,gtype=gtype)
            teff_local = local_temperature(grav,g_pole,T_pole,beta=beta1)
            ints_local = local_intensity(teff_local,grav,np.ones_like(cos_gamma),photband='OPEN.BOL')
            velo_local = np.cross(np.array([x,y,z]).T*to_SI,omega_rot_vec).T
            
            #-- here we can compute the global quantities: total surface area
            #   and luminosity
            lumi_prim = 4*pi*(ints_local*areas_local*to_CGS**2).sum()/constants.Lsol_cgs
            area_prim = 4*areas_local.sum()*to_CGS**2/(4*pi*constants.Rsol_cgs**2)
            logger.info('----PRIMARY DERIVED PROPERTIES')
            logger.info('Polar Radius primary   = %.3g Rsun'%(this_r_pole*a*constants.au/constants.Rsol))
            logger.info("Polar logg primary     = %.3g dex"%(np.log10(g_pole*100)))
            logger.info("Luminosity primary     = %.3g Lsun"%(lumi_prim))
            logger.info("Surface area primary   = %.3g Asun"%(area_prim))
            logger.info("Mean Temp primary      = %.3g K"%(np.average(teff_local,weights=areas_local)))
            ext_dict['Rp1'] = this_r_pole*a*constants.au/constants.Rsol
            ext_dict['loggp1'] = np.log10(g_pole*100)
            ext_dict['LUMI1'] = lumi_prim
            ext_dict['SURF1'] = area_prim
                                    
            #-- for the secondary
            #--------------------
            radius2 = rsec.reshape(theta.shape)
            this_r_pole2 = radius2.ravel()[0]
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
            areas_local2,cos_gamma2 = surface_elements((radius2,mygrid),-grav_local2,gtype=gtype)
            teff_local2 = local_temperature(grav2,g_pole2,T_pole2,beta=beta2)
            ints_local2 = local_intensity(teff_local2,grav2,np.ones_like(cos_gamma2),photband='OPEN.BOL')
            velo_local2 = np.cross(np.array([x2,y2,z2]).T*to_SI,omega_rot_vec).T
            
            #-- here we can compute the global quantities: total surface area
            #   and luminosity
            lumi_sec = 4*pi*(ints_local2*areas_local2*to_CGS**2).sum()/constants.Lsol_cgs
            area_sec = 4*areas_local2.sum()*to_CGS**2/(4*pi*constants.Rsol_cgs**2)
            logger.info('----SECONDARY DERIVED PROPERTIES')
            logger.info('Polar Radius secondary = %.3g Rsun'%(this_r_pole2*a*constants.au/constants.Rsol))
            logger.info("Polar logg secondary   = %.3g dex"%(np.log10(g_pole2*100)))
            logger.info("Luminosity secondary   = %.3g Lsun"%(lumi_sec))
            logger.info("Surface area secondary = %.3g Asun"%(area_sec))
            logger.info("Mean Temp secondary    = %.3g K"%(np.average(teff_local2,weights=areas_local2)))
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
                                    grav,areas_local,teff_local,ints_local,velo_local[0],velo_local[1],velo_local[2],
                                    seamless=False,gtype=gtype,
                                    vtype=['scalar','x','y','z','scalar','scalar','scalar','scalar','vx','vy','vz'])
            #-- stitch the grid!
            theta2_,phi2_,radius2,gravx2,gravy2,gravz2,grav2,areas2,teff2,ints2,vx2,vy2,vz2 = \
                         stitch_grid(theta,phi,radius2,grav_local2[0],grav_local2[1],grav_local2[2],
                                    grav2,areas_local2,teff_local2,ints_local2,velo_local2[0],velo_local2[1],velo_local2[2],
                                    seamless=False,gtype=gtype,
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
            
            #-- take care of the reflection effect
            primary,secondary = reflection_effect(primary,secondary,theta,phi,
                                       A1=A1,A2=A2,max_iter=max_iter_reflection)

        #-- now compute the integrated intensity in the line of sight:
        #-------------------------------------------------------------
        rot_theta = np.arctan2(y1o[di],x1o[di])
        prim = project(primary,view_long=(rot_theta,x1o[di],y1o[di]),
                       view_lat=(view_angle,0,0),photband=photband,
                       only_visible=True,plot_sort=True)
        secn = project(secondary,view_long=(rot_theta,x2o[di],y2o[di]),
                       view_lat=(view_angle,0,0),photband=photband,
                       only_visible=True,plot_sort=True)
        prim['vx'] = -prim['vx'] + RV1[di]*1000.
        secn['vx'] = -secn['vx'] + RV2[di]*1000.
        
        if direc is not None:
            fits.write_recarray(prim,os.path.join(direc,'primary.fits'))
            fits.write_recarray(secn,os.path.join(direc,'secondary.fits'))
        
        #-- the total intensity is simply the sum of the projected intensities
        #   over all visible meshpoints. To calculate the visibility, we
        #   we collect the Y-Z coordinates in one array for easy matching in
        #   the KDTree
        #   We need to know which star is in front. It is the one with the
        #   largest x coordinate
        if secn['x'].min()<prim['x'].min():
            front,back = prim,secn
            front_component = 1
            report += ' Primary in front'
        else:
            front,back = secn,prim
            front_component = 2
            report += ' Secondary in front'
        coords_front = np.column_stack([front['y'],front['z']])
        coords_back = np.column_stack([back['y'],back['z']])    

        #   now find the coordinates of the front component closest to the
        #   the coordinates of the back component
        tree = KDTree(coords_front)
        distance,order = tree.query(coords_back)
        #   meshpoints of the back component inside an eclipse have a
        #   nearest neighbouring point in the (projected) front component
        #   which is closer than sqrt(area) of the surface element connected
        #   to that neighbouring point on the front component
        in_eclipse = distance < np.sqrt(front['areas'][order])
        if np.sum(in_eclipse)>0:
            report += ' during eclipse'
        else:
            report += ' outside eclipse'
        
        #-- so now we can easily compute the total intensity as the sum of
        #   all visible meshpoints:
        total_intensity = front['projflux'].sum() + back['projflux'][-in_eclipse].sum()
        light_curve[di] = total_intensity
        report += "---> Total intensity: %g "%(total_intensity)
        
        if di==0:
            ylim_lc = (0.95*min(prim['projflux'].sum(),secn['projflux'].sum()),1.2*(prim['projflux'].sum()+secn['projflux'].sum()))
        back['projflux'][in_eclipse] = 0
        back['eyeflux'][in_eclipse] = 0
        back['vx'][in_eclipse] = 0
        back['vy'][in_eclipse] = 0
        back['vz'][in_eclipse] = 0
        
        #-- now calculate the *real* observed radial velocity and projected intensity
        if front_component==1:
            RV1_corr[di] = np.average(front['vx']/1000.,weights=front['projflux'])
            RV2_corr[di] = np.average(back['vx'][-in_eclipse]/1000.,weights=back['projflux'][-in_eclipse])
            front_cmap = pl.cm.hot
            back_cmap = pl.cm.cool_r
        else:
            RV2_corr[di] = np.average(front['vx']/1000.,weights=front['projflux'])
            RV1_corr[di] = np.average(back['vx'][-in_eclipse]/1000.,weights=back['projflux'][-in_eclipse])        
            front_cmap = pl.cm.cool_r
            back_cmap = pl.cm.hot
        report += 'RV1=%.3f, RV2=%.3f'%(RV1_corr[di],RV2_corr[di])
        logger.info(report)
        
        #================ START DEBUGGING PLOTS ===================
        if direc is not None:
            #--   first calculate the size of the picture, and the color scales
            if di==0:
                size_x = 1.2*(max(prim['y'].ptp(),secn['y'].ptp())/2. + max(ds))
                size_y = 1.2*(max(prim['z'].ptp(),secn['z'].ptp())/2. + max(ds) * cos(view_angle))
                vmin_image,vmax_image = 0,max([front['eyeflux'].max(),back['eyeflux'].max()])        
                size_top = 1.2*(max(prim['x'].ptp(),secn['x'].ptp())/2. + max(ds))
                
            
            pl.figure(figsize=(16,11))
            pl.subplot(221,aspect='equal');pl.title('line of sight intensity')
            pl.scatter(back['y'],back['z'],c=back['eyeflux'],edgecolors='none',cmap=back_cmap)
            pl.scatter(front['y'],front['z'],c=front['eyeflux'],edgecolors='none',cmap=front_cmap)
            pl.xlim(-size_x,size_x)
            pl.ylim(-size_y,size_y)
            pl.xlabel('X [semi-major axis]')
            pl.ylabel('Z [semi-major axis]')
            
            #-- line-of-sight velocity of the system
            pl.subplot(222,aspect='equal');pl.title('line of sight velocity')
            if di==0:
                vmin_rv = min((prim['vx'].min()/1000.+RV1.min()),(prim['vx'].min()/1000.+RV2.min()))
                vmax_rv = max((secn['vx'].max()/1000.+RV2.max()),(secn['vx'].max()/1000.+RV2.max()))
            pl.scatter(front['y'],front['z'],c=front['vx']/1000.,edgecolors='none',cmap=pl.cm.RdBu_r,vmin=vmin_rv,vmax=vmax_rv)
            pl.scatter(back['y'],back['z'],c=back['vx']/1000.,edgecolors='none',cmap=pl.cm.RdBu_r,vmin=vmin_rv,vmax=vmax_rv)
            cbar = pl.colorbar()
            cbar.set_label('Radial velocity [km/s]')
            pl.xlim(-size_x,size_x)
            pl.ylim(-size_y,size_y)
            pl.xlabel('X [semi-major axis]')
            pl.ylabel('Z [semi-major axis]')
            
            #-- top view of the system
            pl.subplot(223,aspect='equal');pl.title('Top view')
            pl.scatter(prim['x'],prim['y'],c=prim['eyeflux'],edgecolors='none',cmap=pl.cm.hot)
            pl.scatter(secn['x'],secn['y'],c=secn['eyeflux'],edgecolors='none',cmap=pl.cm.cool_r)
            pl.xlim(-size_top,size_top)
            pl.ylim(-size_top,size_top)
            pl.xlabel('X [semi-major axis]')
            pl.ylabel('Y [semi-major axis]')
            
            #-- light curve and radial velocity curve
            pl.subplot(224);pl.title('light curve and RV curve')
            pl.plot(times[:di+1],2.5*np.log10(light_curve[:di+1]),'k-',label=photband)
            pl.plot(times[di],2.5*np.log10(light_curve[di]),'ko',ms=10)
            pl.xlim(times.min(),times.max())
            pl.ylim(2.5*np.log10(ylim_lc[0]),2.5*np.log10(ylim_lc[1]))
            pl.ylabel('Flux [erg/s/cm2/A/sr]')
            pl.legend(loc='lower left',prop=dict(size='small'))
            
            pl.twinx(pl.gca())
            #   primary radial velocity (numerical, kepler and current)
            pl.plot(times[:di+1],RV1_corr[:di+1],'b-',lw=2,label='Numeric 1')
            pl.plot(times[:di+1],RV1[:di+1]-gamma,'g--',lw=2,label='Kepler 1')
            pl.plot(times[di],RV1[di]-gamma,'go',ms=10)
            pl.plot(times[di],RV1_corr[di],'bo',ms=10)
            #   secondary radial velocity (numerical, kepler and current)
            pl.plot(times[:di+1],RV2_corr[:di+1],'r-',lw=2,label='Numeric 2')
            pl.plot(times[:di+1],RV2[:di+1]-gamma,'c--',lw=2,label='Kepler 2')
            pl.plot(times[di],RV2[di]-gamma,'cs',ms=10)
            pl.plot(times[di],RV2_corr[di],'rs',ms=10)
            pl.ylim(1.2*min(min(RV1-gamma),min(RV2-gamma)),+1.2*max(max(RV1-gamma),max(RV2-gamma)))
            pl.xlim(times.min(),times.max())
            pl.ylabel('Radial velocity [km/s]')
            pl.xlabel('Time [d]')
            pl.legend(loc='upper left',prop=dict(size='small'))
            pl.savefig(os.path.join(direc,'%s_los_%04d'%(name,di)),facecolor='0.75')
            pl.close()
            
            #-- REAL IMAGE picture
            pl.figure(figsize=(7,size_y/size_x*7))
            ax = pl.axes([0,0,1,1])
            ax.set_aspect('equal')
            ax.set_axis_bgcolor('k')
            pl.xticks([]);pl.yticks([])
            if front_component==1: sfront,sback = 16,9
            else:                  sfront,sback = 9,16
            pl.scatter(back['y'][-in_eclipse],back['z'][-in_eclipse],s=sback,c=back['eyeflux'][-in_eclipse],edgecolors='none',cmap=pl.cm.gray,vmin=vmin_image,vmax=vmax_image)
            pl.scatter(front['y'],front['z'],s=sfront,c=front['eyeflux'],edgecolors='none',cmap=pl.cm.gray,vmin=vmin_image,vmax=vmax_image)
            pl.xlim(-size_x,size_x);pl.ylim(-size_y,size_y)
            pl.savefig(os.path.join(direc,'%s_image_%04d'%(name,di)),facecolor='k')
            pl.close()
            #================   END DEBUGGING PLOTS ===================
        
    return times, light_curve, RV1_corr, RV2_corr

if __name__=="__main__":
    import doctest
    doctest.testmod()
    pl.show()
    
    logger = loggers.get_basic_logger("")
    import sys
    