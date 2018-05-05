"""
Calculate the stellar surface displacements due to pulsations.


"""
from numpy import sqrt,pi,sin,cos,exp,log10,linspace,mgrid,polyval,zeros_like
from scipy.special import legendre
from scipy.misc import factorial

#{ Helper functions

def legendre_(l,m,x):
    """
    Legendre polynomial.

    Check equation (3) from Townsend, 2002:

    >>> ls,x = [0,1,2,3,4,5],cos(linspace(0,pi,100))
    >>> check = 0
    >>> for l in ls:
    ...     for m in range(-l,l+1,1):
    ...         Ppos = legendre_(l,m,x)
    ...         Pneg = legendre_(l,-m,x)
    ...         mycheck = Pneg,(-1)**m * factorial(l-m)/factorial(l+m) * Ppos
    ...         check += sum(abs(mycheck[0]-mycheck[1])>1e-10)
    >>> print check
    0
    """
    m_ = abs(m)
    legendre_poly = legendre(l)
    deriv_legpoly_ = legendre_poly.deriv(m=m_)
    deriv_legpoly = polyval(deriv_legpoly_,x)
    P_l_m = (-1)**m_ * (1-x**2)**(m_/2.) * deriv_legpoly
    if m<0:
        P_l_m = (-1)**m_ * factorial(l-m_)/factorial(l+m_) * P_l_m
    return P_l_m

def sph_harm(theta,phi,l=2,m=1):
    """
    Spherical harmonic according to Townsend, 2002.

    This function is memoized: once a spherical harmonic is computed, the
    result is stored in memory

    >>> theta,phi = mgrid[0:pi:20j,0:2*pi:40j]
    >>> Ylm20 = sph_harm(theta,phi,2,0)
    >>> Ylm21 = sph_harm(theta,phi,2,1)
    >>> Ylm22 = sph_harm(theta,phi,2,2)
    >>> Ylm2_2 = sph_harm(theta,phi,2,-2)
    >>> p = figure()
    >>> p = gcf().canvas.set_window_title('test of function <sph_harm>')
    >>> p = subplot(411);p = title('l=2,m=0');p = imshow(Ylm20.real,cmap=cm.RdBu)
    >>> p = subplot(412);p = title('l=2,m=1');p = imshow(Ylm21.real,cmap=cm.RdBu)
    >>> p = subplot(413);p = title('l=2,m=2');p = imshow(Ylm22.real,cmap=cm.RdBu)
    >>> p = subplot(414);p = title('l=2,m=-2');p = imshow(Ylm2_2.real,cmap=cm.RdBu)

    """
    factor = (-1)**m * sqrt( (2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m))
    Plm = legendre_(l,m,cos(theta))
    return factor*Plm*exp(1j*m*phi)

def dsph_harm_dtheta(theta,phi,l=2,m=1):
    """
    Derivative of spherical harmonic wrt colatitude.

    Using Y_l^m(theta,phi).

    Equation::

        sin(theta)*dY/dtheta = (l*J_{l+1}^m * Y_{l+1}^m - (l+1)*J_l^m * Y_{l-1,m})

    E.g.: Phd thesis of Joris De Ridder
    """
    if abs(m)>=l:
        Y = 0.
    else:
        factor = 1./sin(theta)
        term1 = l     * norm_J(l+1,m) * sph_harm(theta,phi,l+1,m)
        term2 = (l+1) * norm_J(l,m)   * sph_harm(theta,phi,l-1,m)
        Y = factor * (term1 - term2)
    return Y/sin(theta)

def dsph_harm_dphi(theta,phi,l=2,m=1):
    """
    Derivative of spherical harmonic wrt longitude.

    Using Y_l^m(theta,phi).

    Equation::

        dY/dphi = i*m*Y
    """
    return 1j*m*sph_harm(theta,phi,l,m)


def norm_J(l,m):
    """
    Normalisation factor
    """
    if abs(m)<l:
        J = sqrt( float((l**2-m**2))/(4*l**2-1))
    else:
        J = 0
    return J

def norm_atlp1(l,m,Omega,k):
    """
    Omega is actually spin parameter (Omega_rot/omega_freq)
    """
    return Omega * (l-abs(m)+1)/(l+1) * 2./(2*l+1) * (1-l*k)

def norm_atlm1(l,m,Omega,k):
    """
    Omega is actually spin parameter (Omega_rot/omega_freq)
    """
    return Omega * (l+abs(m))/l * 2./(2*l+1) * (1 + (l+1)*k)

#}

#{ Displacement fields

def radial(theta,phi,l,m,t):
    """
    Radial displacement, see Zima 2006.

    t in phase units
    """
    return sph_harm(theta,phi,l,m) * exp(1j*t)

def colatitudinal(theta,phi,l,m,t,Omega,k):
    term1 = k * dsph_harm_dtheta(theta,phi,l,m) * exp(1j*t)
    term2 = norm_atlp1(l,m,Omega,k) / sin(theta) * dsph_harm_dphi(theta,phi,l+1,m)  * exp(1j*t + pi/2)
    term3 = norm_atlm1(l,m,Omega,k) / sin(theta) * dsph_harm_dphi(theta,phi,l-1,m)  * exp(1j*t - pi/2)
    return term1 + term2 + term3

def longitudinal(theta,phi,l,m,t,Omega,k):
    term1 = k /sin(theta) * dsph_harm_dphi(theta,phi,l,m)*exp(1j*t)
    term2 = -norm_atlp1(l,m,Omega,k) * dsph_harm_dtheta(theta,phi,l+1,m)*exp(1j*t+pi/2)
    term3 = -norm_atlm1(l,m,Omega,k) * dsph_harm_dtheta(theta,phi,l-1,m)*exp(1j*t-pi/2)
    return term1 + term2 + term3

def surface(theta,phi,l,m,t,Omega=0.1,k=1.,asl=0.2,radius=1.):
    ksi_r = asl*sqrt(4*pi)*radial(theta,phi,l,m,t)
    if l>0:
        ksi_theta = asl*sqrt(4*pi)*colatitudinal(theta,phi,l,m,t,Omega,k)
        ksi_phi = asl*sqrt(4*pi)*longitudinal(theta,phi,l,m,t,Omega,k)
    else:
        ksi_theta = zeros_like(theta)
        ksi_phi = zeros_like(phi)

    return (radius+ksi_r.real),\
           (theta + ksi_theta.real),\
           (phi + ksi_phi.real)

def observables(theta,phi,teff,logg,l,m,t,Omega=0.1,k=1.,asl=0.2,radius=1.,delta_T=0.01+0j,delta_g=0.01+0.5j):
    rad_part = radial(theta,phi,l,m,t)
    ksi_r = asl*sqrt(4*pi)*rad_part
    if l>0:
        ksi_theta = asl*sqrt(4*pi)*colatitudinal(theta,phi,l,m,t,Omega,k)
        ksi_phi = asl*sqrt(4*pi)*longitudinal(theta,phi,l,m,t,Omega,k)
    else:
        ksi_theta = zeros_like(theta)
        ksi_phi = zeros_like(phi)
    gravity = 10**(logg-2)
    return (radius+ksi_r.real),\
           (theta + ksi_theta.real),\
           (phi + ksi_phi.real),\
           (teff + (delta_T*rad_part*teff).real),\
           log10(gravity+(delta_g*rad_part*gravity).real)+2

#}

def mainx():
    # Test function: legendre_
    ls,x = [0,1,2,3,4,5],cos(linspace(0,pi,100))
    check = 0
    for l in ls:
        for m in range(-l,l+1,1):
            Ppos = legendre_(l,m,x)
            Pneg = legendre_(l,-m,x)
            mycheck = Pneg,(-1)**m * factorial(l-m)/factorial(l+m) * Ppos
            check += sum(abs(mycheck[0]-mycheck[1])>1e-10)
    print(check) # should yield "0"

    # Test function: sph_harm
    import matplotlib.pyplot as plt
    theta,phi = mgrid[0:pi:20j,0:2*pi:40j]
    Ylm20 = sph_harm(theta,phi,2,0)
    Ylm21 = sph_harm(theta,phi,2,1)
    Ylm22 = sph_harm(theta,phi,2,2)
    Ylm2_2 = sph_harm(theta,phi,2,-2)

    fig,axs = plt.subplots(1,4,figsize=(16,4))
    fig.suptitle('Test of function <sph_harm>', fontsize=20)
    axs[0].set_title('l=2,m=0');axs[0].imshow(Ylm20.real,cmap='jet')
    axs[1].set_title('l=2,m=1');axs[1].imshow(Ylm21.real,cmap='jet')
    axs[2].set_title('l=2,m=2');axs[2].imshow(Ylm22.real,cmap='jet')
    axs[3].set_title('l=2,m=-2');axs[3].imshow(Ylm2_2.real,cmap='jet')
    plt.tight_layout()

    # Test functions: dsph_harm_dtheta, norm_J (called by dsph_harm_dtheta), and dsph_harm_dphi
    spherical_harmonic_derivative_wrt_colatitude = dsph_harm_dtheta(theta,phi,l=2,m=1)
    spherical_harmonic_derivative_wrt_longitude = dsph_harm_dphi(theta,phi,l=2,m=1)

    fig,axs = plt.subplots(2,1,figsize=(12,8))
    axs[0].set_title('Test of function <dsph_harm_dtheta>', fontsize=20)
    axs[0].imshow(spherical_harmonic_derivative_wrt_colatitude.real,cmap='jet')
    axs[1].set_title('Test of function <dsph_harm_dphi>', fontsize=20)
    axs[1].imshow(spherical_harmonic_derivative_wrt_longitude.real,cmap='jet')
    plt.tight_layout()

    # Test functions: norm_atlp1,norm_atlm1
    l,m = 2,1
    Omega,k = 0.1,1
    print(norm_atlp1(l,m,Omega,k))
    print(norm_atlm1(l,m,Omega,k))

    # Test function: radial
    times = linspace(0,2*pi,100)
    t = times[50]
    radial_displacement = radial(theta,phi,l,m,t)

    plt.figure()
    plt.title('Test of function <radial>', fontsize=20)
    plt.imshow(radial_displacement.real,cmap='jet')
    plt.tight_layout()

    # Test function: colatitudinal,longitudinal
    fig,axs = plt.subplots(2,1,figsize=(12,8))
    axs[0].set_title('Test of function <colatitudinal>', fontsize=20)
    axs[0].imshow(colatitudinal(theta,phi,l,m,t,Omega,k).real,cmap='jet')
    axs[1].set_title('Test of function <longitudinal>', fontsize=20)
    axs[1].imshow(longitudinal(theta,phi,l,m,t,Omega,k).real,cmap='jet')
    plt.tight_layout()
    # plt.show()

    # Test function: surface
    test_surface = surface(theta,phi,l,m,t,Omega=0.1,k=1.,asl=0.2,radius=1.)
    print(type(test_surface)) # returns a tuple

    # Test function: observables
    teff,logg = 5777,4.438 # approximate solar values
    test_observables = observables(theta,phi,teff,logg,l,m,t,Omega=0.1,k=1.,asl=0.2,radius=1.,delta_T=0.01+0j,delta_g=0.01+0.5j)
    print(type(test_observables)) # returns a tuple

if __name__ == '__main__': mainx()

#=========THE FOLLOWING EXAMPLE IS BROKEN (IMPORTED MODULES RAISE ERRORS) | USE AT YOUR OWN RISK==============

# if __name__=="__main__":
#     import numpy as np
#     from ivs.roche import local
#     from ivs.coordinates import vectors
#     from enthought.mayavi import mlab
#     from divers import multimedia
#     from scipy.spatial import Delaunay
#     theta,phi = local.get_grid(50,25,gtype='triangular')
#     r = np.ones_like(theta)
#     x,y,z = vectors.spher2cart_coord(r,phi,theta)
#     points = np.array([x,y,z]).T
#     grid = Delaunay(points)
#     #keep = phi>pi
#     #theta,phi = theta[keep],phi[keep]
#     l,m = 2,2
#     asl = 0.01
#
#     for k in [0,1.,2.]:
#         for l in range(1,5):
#             for m in range(0,l+1,1):
#                 mlab.figure(size=(1000,800))
#                 mlab.gcf().scene.disable_render = True
#                 if l==0 or l==1:
#                     asl = 0.1
#                 else:
#                     asl = 0.01
#                 old_center=None
#                 for i,t in enumerate(np.linspace(0,2*pi,100)):
#                     print(k,l,m,i)
#                     r,th,ph = surface(theta,phi,l,m,t,asl=asl,k=k)
#                     center,size,normal = local.surface_normals(r,ph,th,grid,gtype='triangular')
#
#                     if i==0:
#                         colors = r
#                         r_c,phi_c,theta_c = vectors.cart2spher_coord(*center.T)
#                         colors_ = r_c
#
#
#                     mlab.clf()
#                     mlab.points3d(center.T[0],center.T[1],center.T[2],colors_,scale_factor=0.05,scale_mode='none',colormap='RdBu',vmin=colors_.min(),vmax=colors_.max())
#                     #mlab.quiver3d(center.T[0],center.T[1],center.T[2],normal.T[0],normal.T[1],normal.T[2],colormap='spectral',scale_mode='none')
#                     mlab.colorbar()
#
#                     if i>=1:
#                         vx,vy,vz = center.T[0]-old_center.T[0],\
#                                    center.T[1]-old_center.T[1],\
#                                    center.T[2]-old_center.T[2]
#                         v = np.sqrt(vx**2+vy**2+vz**2)
#                         mlab.quiver3d(center.T[0],center.T[1],center.T[2],\
#                                       vx,vy,vz,scalars=v,colormap='spectral',scale_mode='scalar')
#                         #mlab.show()
#                     old_center = center.copy()
#                     mlab.view(distance=5,azimuth=-90,elevation=90)
#                     mlab.savefig('pulsation_lm%d%d_k%03d_%03d.png'%(l,m,k,i))
#                 mlab.close()
#                 multimedia.make_movie('pulsation_lm%d%d_k%03d_*.png'%(l,m,k),output='pulsation_lm%d%d_k%03d.avi'%(l,m,k))
