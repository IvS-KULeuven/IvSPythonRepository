
"""
covellipse package
Author: Joris De Ridder

The covellipse package allows to draw contours of a 2D covariance matrix.

Example:

Import the required packages:

>>> from numpy import *
>>> from numpy.random import multivariate_normal
>>> from matplotlib import pylab as p
>>> from covellipse import sigmaEllipse

Make some fake data:

>>> covMatrix = array([[2.0, -0.75*sqrt(4)*sqrt(2)],[-0.75*sqrt(4)*sqrt(2), 4.0]])
>>> m = multivariate_normal([0,0],covMatrix, 3000)

Compute the 2-sigma ellipse, centered around (0,0):

>>> x,y=sigmaEllipse(0,0,covMatrix,2)

Make a plot of the data as well as the ellipse:

>>> p.scatter(m[:,0], m[:,1], s=5)
>>> p.plot(x,y,linewidth=2) 
"""


import numpy as np
from math import atan,sqrt,pi


def ellipse(xc, yc, a, b, phi):

    """
    Returns x and y values of a complete ellipse
    
    @param xc: x-coordinate of the center of the ellipse
    @type xc: float
    @param yc: y-coordinate of the center of the ellipse
    @type yc: float
    @param a: half the long axis of the ellipse
    @type a: float
    @param b: half the small axis of the ellipse
    @type b: float
    @param phi: angle between long axis and the x-axis (radians)
    @type phi: float    
    @return: x,y: coordinate arrays of ellipse 
    @rtype: tuple of 2 arrays (dtype: float)
    """
    
    t = np.linspace(0,2*pi,100)
    x = xc + a * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi)
    y = yc + a * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi)
    
    return x,y
    



def sigmaEllipse(xc, yc, covMatrix, nSigma):

    """
    Return x and y value of the |nSigma|-sigma ellipse
    Given a covariance matrix, and the center of the ellipse (xc,yc)
    
    @param xc: x-coordinate of the center of the ellipse
    @type xc: float
    @param yc: y-coordinate of the center of the ellipse
    @type yc: float
    @param covMatrix: 2D covariance matrix
    @type covMatrix: 2x2 float nd-array
    @param nSigma: the nSigma-contour of the ellipse will be returned   
    @type nSigma: float
    @return: (x,y): the x and y values of the ellipse 
    @rtype: tuple of two float ndarrays
    """

    a = covMatrix[0,0]
    b = covMatrix [0,1]
    c = b
    d = covMatrix[1,1]
    
    eigenvalue1 = 0.5*(a+d+np.sqrt((a-d)**2+4*b*c))
    eigenvalue2 = 0.5*(a+d-np.sqrt((a-d)**2+4*b*c))
    
    semimajor = sqrt(eigenvalue1) * nSigma
    semiminor = sqrt(eigenvalue2) * nSigma
    
    if b == 0.0:
        eigenvector1 = np.array([1, 0])
        eigenvector2 = np.array([0, 1])
        phi = 0.0
    else:
        eigenvector1 = np.array([b, eigenvalue1-a])
        eigenvector2 = np.array([b, eigenvalue2-a])
        phi = atan((eigenvalue1-a)/b)
           
    x,y = ellipse(xc, yc, semimajor, semiminor, phi)
    
    return x,y

    
    
