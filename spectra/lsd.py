"""
Compute Least-Square-Deconvolutions of spectra.

Section 1. Single stars
=======================

Section 1.1 Without Tikhonov regularization
-------------------------------------------

Generate some spectra of a single star:

>>> velos,V,S,masks = __generate_test_spectra(1,binary=False,noise=0.01)

Compute the LSD profiles in certain radial velocity range:

>>> rvs = np.linspace(-10,10,100)
>>> Z,cc = lsd(velos,V,S,rvs,masks)

Plot both the input spectrum and the LSD and CCF profiles:

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(velos,V[:,0],'k-',label='Observation')
>>> p = pl.vlines(masks[0][0],1,1-np.array(masks[0][1]),color='r',label='Mask')
>>> p = pl.legend(loc='best')
>>> p = pl.subplot(122)
>>> p = pl.plot(rvs,Z[0][0],'k-',label='LSD')
>>> p = pl.plot(rvs,cc[0][0],'r-',label='CCF')
>>> p = pl.legend(loc='best')

]]include figure]]ivs_spectra_lsd01.png]

Section 1.2 With Tikhonov regularization
----------------------------------------

Do the same as above, but for more noise and different Tikhonov regularizations.

>>> velos,V,S,masks = __generate_test_spectra(1,binary=False,noise=0.1)

>>> rvs = np.linspace(-10,10,100)
>>> lambdas = [0.,0.1,0.5,1.0]
>>> output = [lsd(velos,V,S,rvs,masks,Lambda=lam) for lam in lambdas]

We plot both the input spectrum and the LSD and CCF profiles:

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(velos,V[:,0],'k-',label='Observation')
>>> p = pl.vlines(masks[0][0],1,1-np.array(masks[0][1]),color='r',label='Mask')
>>> p = pl.legend(loc='best')
>>> p = pl.subplot(122)
>>> for lam,(Z,cc) in zip(lambdas,output):
...     p = pl.plot(rvs,Z[0][0],'-',label='$\Lambda$=%.1f'%(lam),lw=2)
>>> p = pl.legend(loc='best')

]]include figure]]ivs_spectra_lsd02.png]

Section 2. Binary stars
=======================

Generate some spectra of a binary star:

>>> velos,V,S,masks = __generate_test_spectra(1,binary=True,noise=0.01)

Compute the LSD profiles in certain radial velocity range, first using only
one line mask, then both:

>>> rvs = np.linspace(-10,10,100)
>>> Z1,cc1 = lsd(velos,V,S,rvs,masks[:1])
>>> Z2,cc2 = lsd(velos,V,S,rvs,masks)

Plot both the spectrum and the LSD and CCF profiles. Note that the CCF in the
binary case is exactly the same as in the single star case (see implementation).
First, we plot the LSD profile under the assumption of a single star. Second,
we plot the LSD profile when taking binarity into account.

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(velos,V[:,0],'k-',label='Observation')
>>> p = pl.vlines(masks[0][0],1,1-np.array(masks[0][1]),color='r',label='Mask 1',lw=2)
>>> p = pl.legend(loc='lower right')
>>> p = pl.subplot(122)
>>> p = pl.plot(rvs,Z1[0][0],'k-',label='LSD 1')
>>> p = pl.plot(rvs,cc1[0][0],'r-',label='CCF 1')
>>> p = pl.legend(loc='best')

]]include figure]]ivs_spectra_lsd03.png]

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(velos,V[:,0],'k-',label='Observation')
>>> p = pl.vlines(masks[0][0],1,1-np.array(masks[0][1]),color='r',label='Mask 1',lw=2)
>>> p = pl.vlines(masks[1][0],1,1-np.array(masks[1][1]),color='b',label='Mask 2',lw=2)
>>> p = pl.legend(loc='lower right')
>>> p = pl.subplot(122)
>>> p = pl.plot(rvs,Z2[0][0],'r-',label='LSD 1',lw=2)
>>> p = pl.plot(rvs,Z2[0][1],'b-',label='LSD 2',lw=2)
>>> p = pl.legend(loc='best')

]]include figure]]ivs_spectra_lsd04.png]

"""
import pylab as pl
import numpy as np
import numpy.linalg as la
from ivs.sigproc import evaluate
import itertools

def lsd(velos,V,S,rvs,masks,Lambda=0.):
    """
    Compute LSD profiles and cross correlation functions.
    
    Possibility to include Tikhonov regularization to clean up the profiles,
    when setting C{Lambda>0}.
    
    Possibility to include multiprofile LSD. Parameter C{masks} should be a list
    of C{(centers,weights)} (if you give only one mask, give
    C{masks=[(centers,weights)]}.
    
    See Donati, 1997 for the original paper and Kochukhov, 2010 for extensions.
    
    @parameter velos: velocity vector of observations
    @type velos: array of length N_spec
    @parameter V: observation array
    @type V: N_obs x N_spec array
    @parameter S: weights of individual pixels
    @type S: array of length N_spec
    @parameter rvs: radial velocity vector to compute the profile on
    @type rvs: array of length N_rv
    @parameter Lambda: Tikhonov regularization parameter
    @parameter masks: list of tuples (center velocities, weights)
    @type masks: list (length N_mask) of tuples of 1D arrays
    @type Lambda: float
    @return: LSD profile, CCF of shape (N_obs x (N_rv.N_mask))
    @rtype: 2D array, 2D array
    """
    #-- some global parameters
    m,n = len(rvs),len(velos)
    Nspec = V.shape[1]
    Nmask = len(masks)
    V = np.matrix(V)-1
    
    #-- weights of the individual pixels
    S = np.matrix(np.diag(S))
    #-- line masks (yes, this can be vectorized but I'm too lazy for the moment)
    M = np.matrix(np.zeros((n,m*len(masks))))
    for N,(line_centers,weights) in enumerate(masks):
        for l,lc in enumerate(line_centers):
            for i in range(n):
                for j in range(m-1):
                    vi = velos[i]-lc
                    if not (rvs[j]<vi<rvs[j+1]): continue
                    M[i,j+  N*m] = weights[l]*(rvs[j+1]-vi)/(rvs[j+1]-rvs[j])
                    M[i,j+1+N*m] = weights[l]*(vi    -rvs[j])/(rvs[j+1]-rvs[j])
    M = np.matrix(M)
    #-- regularization parameter
    if Lambda:
        R = np.matrix(np.zeros((m*Nmask,m*Nmask)))
        for i in range(1,m-1):
            R[i,i] = 2
            R[i-1,i] = -1
            R[i+1,i] = -1
        R[0,0] = 1
        R[1,0] = -1
        R[-1,-1] = 1
        R[-2,-1] = -1
    #-- compute the LSD    
    X = M.T*(S**2)
    XM = X*M
    if Lambda:
        XM = XM+Lambda*R
    cc = X*V # this is in fact the cross correlation profile
    #-- XM is of shape (mxm), cc is of shape (mxNspec)
    #-- we can solve this system quickly ourselves or call numpy. I find the
    #   latter more elegant, but it might be slower.
    #Z = la.inv(XM)*cc
    Z,res,rank,s = la.lstsq(XM,cc)
    #-- retrieve LSD profile and cross-correlation function
    Z = np.array(Z.T)
    cc = np.array(cc.T)
    #-- split up the profiles
    Z_ = []
    C_ = []
    for i in range(len(Z)):
        Z_.append([])
        C_.append([])
        for N in range(Nmask):
            Z_[-1].append(Z[i][N*m:(N+1)*m])
            C_[-1].append(cc[i][N*m:(N+1)*m])
    #-- that's it!
    return Z_,C_

def __generate_test_spectra(Nspec,binary=False,noise=0.01):
    spec_length = 1000 # n
    velo_length = 100 # m
    obs_velo1 = [-3.14,-3.5,-4,3.25][:Nspec]
    obs_velo2 = [+3.14,+3.5,+4,-3.25][:Nspec]
    m,n = velo_length,spec_length
    #-- for the observations
    velos = np.linspace(-30,30,spec_length)
    line_centers1 = [-11,3,9]
    weights1 = [0.5,0.1,0.3]
    line_centers2 = [-5,2,10]
    weights2 = [0.3,0.4,0.1]
    
    #-- weights
    S = np.ones(spec_length)
    #-- profiles
    V = [np.zeros(spec_length) for i in range(Nspec)]
    masks = [(line_centers1,weights1)]
    for i,prof_velo in enumerate(obs_velo1):
        for line_center,weight in zip(line_centers1,weights1):
            V[i] += evaluate.gauss(velos,[weight,line_center-prof_velo,1.])
        V[i] += np.random.normal(size=n,scale=noise)
    if binary:
        masks.append((line_centers2,weights2))
        for i,prof_velo in enumerate(obs_velo2):
            for line_center,weight in zip(line_centers2,weights2):
                V[i] += evaluate.gauss(velos,[weight,line_center-prof_velo,1.])            
    V = 1-np.array(V)
    V = np.matrix(V).T
    
    return velos,V,S,masks

if __name__=="__main__":
    import doctest
    doctest.testmod()
    pl.show()
