"""
Vector and coordinate identities and transformations
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
    
    
def rotate(x,y,theta,x0=0.,y0=0.):
    """
    Rotate counter clockwise.
    """
    xnew = cos(theta)*x - sin(theta)*y - x0
    ynew = sin(theta)*x + cos(theta)*y - y0
    return xnew,ynew
#}

#{ Coordinate vector transformations

def spher2cart((r,phi,theta),(a_r,a_phi,a_theta)):
    """
    theta is angle from z-axis (colatitude)
    phi is longitude
    """
    transfo = np.matrix([[sin(theta)*cos(phi),  cos(theta)*cos(phi), -sin(phi)],
                         [sin(theta)*sin(phi),  cos(theta)*sin(phi),  cos(phi)],
                         [cos(theta)         , -sin(theta)         ,  0       ]])
    vector = np.matrix([a_r,a_theta,a_phi]).T
    return np.asarray((transfo*vector).T)[0]

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
    Compute angle between two vectors (or between two grids of vectors).
    
    Input vectors should be numpy arrays.
    """
    return (vec1*vec2).sum(axis=0) / (norm(vec1)*norm(vec2))

#}