"""
Vector and coordinate identities and transformations

General conventions:

    - theta is the colatitude, i.e. the angle from the Z-axis
    - phi is the longitude


Transform cartesian coordinates to spherical coordinates and back

>>> x,y,z = np.random.uniform(low=-1,high=1,size=(3,10))
>>> r,phi,theta = cart2spher_coord(x,y,z)
>>> x_,y_,z_ = spher2cart_coord(r,phi,theta)
>>> np.allclose(x,x_),np.allclose(y,y_),np.allclose(z,z_)
(True, True, True)

Transform cartesian vectors to spherical vectors and back

>>> x0,y0,z0 = np.random.uniform(low=-1,high=1,size=(3,10))
>>> x1,y1,z1 = np.random.uniform(low=-1,high=1,size=(3,10))
>>> r0,phi0,theta0 = cart2spher_coord(x0,y0,z0)
>>> r1,phi1,theta1 = cart2spher((x0,y0,z0),(x1,y1,z1))
>>> x1_,y1_,z1_ = spher2cart((r0,phi0,theta0),(r1,phi1,theta1))
>>> np.allclose(x1,x1_),np.allclose(y1,y1_),np.allclose(z1,z1_)
(True, True, True)

"""

import numpy as np
from numpy import cos,sin,pi,sqrt


#{ Coordinate transformations

def spher2cart_coord(r,phi,theta):
    """
    Spherical to Cartesian coordinates
    """
    x = r*cos(phi)*sin(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(theta)
    return x,y,z

def cart2spher_coord(x,y,z):
    """
    theta colatitude
    phi longitude
    """
    rho = np.sqrt(x**2+y**2+z**2)
    phi = np.arctan2(y,x)
    theta = np.arctan2(np.sqrt(x**2+y**2),z)
    return rho,phi,theta

def rotate(x,y,theta,x0=0.,y0=0.):
    """
    Rotate counter clockwise.
    """
    xnew = cos(theta)*x - sin(theta)*y - x0
    ynew = sin(theta)*x + cos(theta)*y - y0
    return xnew,ynew
#}

#{ Coordinate vector transformations

def spher2cart(spherical_coord, a_spherical_coord):
    """
    theta is angle from z-axis (colatitude)
    phi is longitude

    E.g. http://www.engin.brown.edu/courses/en3/Notes/Vector_Web2/Vectors6a/Vectors6a.htm

    >>> np.random.seed(111)
    >>> r,phi,theta = np.random.uniform(low=-1,high=1,size=(3,2))
    >>> a_r,a_phi,a_theta = np.random.uniform(low=-1,high=1,size=(3,2))
    >>> a_x,a_y,a_z = spher2cart((r,phi,theta),(a_r,a_phi,a_theta))
    """
    (r,phi,theta) = spherical_coord
    (a_r,a_phi,a_theta) = a_spherical_coord
    ax = sin(theta)*cos(phi)*a_r + cos(theta)*cos(phi)*a_theta - sin(phi)*a_phi
    ay = sin(theta)*sin(phi)*a_r + cos(theta)*sin(phi)*a_theta + cos(phi)*a_phi
    az = cos(theta)         *a_r - sin(theta)         *a_theta
    return ax,ay,az


def cart2spher(cartesian_coordinats_0, cartesian_coordinats_1):
    """
    theta is angle from z-axis (colatitude)
    phi is longitude

    return r,phi,theta
    """
    (x0,y0,z0) = cartesian_coordinats_0
    (x1,y1,z1) = cartesian_coordinats_1
    r,phi,theta = cart2spher_coord(x0,y0,z0)
    ar     = sin(theta)*cos(phi)*x1 + sin(theta)*sin(phi)*y1 + cos(theta)*z1
    atheta = cos(theta)*cos(phi)*x1 + cos(theta)*sin(phi)*y1 - sin(theta)  *z1
    aphi   = -sin(phi)          *x1 + cos(phi)           *y1
    return ar,aphi,atheta

#}

#{ Vector identities

def norm(vec):
    """
    Euclidic norm of a vector (or of a grid of vectors)

    Input vectors should be numpy arrays.
    """
    return sqrt((vec**2).sum(axis=0))


def angle(vec1,vec2):
    """
    Compute angle between two vectors (or between two grids of vectors).

    Input vectors should be numpy arrays.
    """
    return np.arccos( (vec1*vec2).sum(axis=0) / (norm(vec1)*norm(vec2)))


def cos_angle(vec1,vec2):
    """
    Compute cosine of angle between two vectors (or between two grids of vectors).

    Input vectors should be numpy arrays.
    """
    return (vec1*vec2).sum(axis=0) / (norm(vec1)*norm(vec2))

#}


if __name__=="__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
