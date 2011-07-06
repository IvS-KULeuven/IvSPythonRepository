"""
Calculate the stellar surface displacements due to pulsations.


"""

from numpy import sqrt,pi,sin,cos,exp
from scipy.special import lpmv,legendre
from scipy.misc.common import factorial
from scipy.integrate import dblquad,quad

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
    deriv_legpoly = np.polyval(deriv_legpoly_,x)
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
    Derivative of spherical harmonic to colatitude.
    
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
    Derivative of spherical harmonic to longitude.
    
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
    return asl * sph_harm(theta,phi,l,m) * exp(1j*t)

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

def surface(theta,phi,l,m,t,Omega=0.1,k=1.,radius=1.,asl=1):
    ksi_r = asl*sqrt(4*pi)*radial(theta,phi,l,m,t)
    ksi_theta = asl*sqrt(4*pi)*colatitudinal(theta,phi,l,m,t,Omega,k)
    ksi_phi = asl*sqrt(4*pi)*longitudinal(theta,phi,l,m,t,Omega,k)
    
    return radius+ksi_r, theta + ksi_theta, phi + ksi_phi
    
#}

if __name__=="__main__":
    from ivs.roche import local
    from ivs.coordinates import vectors
    from enthought.mayavi import mlab
    theta,phi,grid = get_grid(20,40,gtype='delaunay')
    l,m = 0,0
    
    mlab.figure(size=(1000,800))
    mlab.gcf().scene.disable_render = True
    
    for i,t in enumerate(np.linspace(0,1,10)):
        r,th,ph = surface(theta,phi,l,m,t)
        x,y,z = vectors.spher2cart_coords(r,th,ph)
        mlab.clf()
        mlab.points3d(x,y,z,scale_fator=0.05,scale_mode='none')
        mlab.savefig('pulsation_%03d.png'%(i))
    mlab.show()
        