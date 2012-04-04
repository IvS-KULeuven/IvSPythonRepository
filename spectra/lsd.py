import pylab as pl
import numpy as np
import numpy.linalg as la
from ivs.sigproc import evaluate
import itertools

def lsd(velos,V,S,rvs,masks,times=None,Lambda=0.):
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
    V = np.matrix(S)
    
    #-- weights of the individual pixels
    S = np.matrix(np.diag(S))
    #-- line masks
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
        #assert(R.shape==(m*Nmask,m*Nmask))
    #-- check some shapes
    #assert(V.shape==(n,Nspec))
    #assert(M.shape==(n,m*Nmask))
    #assert(S.shape==(n,n))
    #-- compute the LSD    
    S2 = S**2
    X = M.T*S2
    XM = X*M
    if Lambda:
        XM = XM+Lambda*R
    cc = np.dot(X,V)
    unc = la.inv(XM)
    Z = np.dot(unc,cc)
    #-- check the shape
    #assert(Z.shape==(m*Nmask,Nspec))
    #-- retrieve LSD profile and cross-correlation function
    Z = np.array(Z.T)
    cc = np.array(cc.T)
    #-- plotting example
    #colors = itertools.cycle([pl.cm.spectral(j) for  j in np.linspace(0,1,Nspec)])
    #for i in range(len(Z)):
        #color = colors.next()
        #for N in range(Nmask):
            #print i,N*m,(N+1)*m
            #pl.plot(rvs,Z[i][N*m:(N+1)*m],'o-',color=color)
    #-- that's it!
    return Z,cc

#def generate_spectra(Nspec):
    #spec_length = 1000 # n
    #velo_length = 100 # m
    #obs_velo1 = [-3.14,-3.5,-4,3.25][:Nspec]
    #obs_velo2 = [+3.14,+3.5,+4,-3.25][:Nspec]
    #m,n = velo_length,spec_length
    ##-- for the observations
    #velos = np.linspace(-30,30,spec_length)
    #line_centers1 = [-11,3,9]
    #weights1 = [0.5,0.1,0.3]
    #line_centers2 = [-5,2,10]
    #weights2 = [0.3,0.4,0.1]
    
    ##-- weights
    #S = np.ones(spec_length)
    ##-- profiles
    #pl.figure()
    #V = [np.ones(spec_length) for i in range(Nspec)]
    #for i,prof_velo in enumerate(obs_velo1):
        #for line_center,weight in zip(line_centers1,weights1):
            #V[i] += evaluate.gauss(velos,[weight,line_center-prof_velo,1.])
    #for i,prof_velo in enumerate(obs_velo2):
        #for line_center,weight in zip(line_centers2,weights2):
            #V[i] += evaluate.gauss(velos,[weight,line_center-prof_velo,1.])
        #V[i] += np.random.normal(size=n,scale=0.1)
        #pl.plot(velos,1-V[i]+i/10.,'k-')
    #V = 1-np.array(V)
    #V = np.matrix(V).T
    
    #return velos,V,S,[(line_centers1,weights1),(line_centers2,weights2)]

#velos,V,S,masks = generate_spectra(4)    
#rvs = np.linspace(-10,10,100)

#Z,cc = lsd(velos,V,S,rvs,masks,times=None,Lambda=1.)

#pl.show()