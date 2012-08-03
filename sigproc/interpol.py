"""
Non-standard interpolation methods.
"""
import numpy as np
from scipy import ndimage
import pyfits as pf
import time
import itertools

def __df_dx(oldx,oldy,index,sharp=False):
    """
    Preliminary estimation of df/dx
    """
    xm1,fm1 = oldx[index-1],oldy[index-1]
    x  ,f   = oldx[index]  ,oldy[index]
    xp1,fp1 = oldx[index+1],oldy[index+1]
    
    if not sharp:
        df_dx = 1./(xp1-xm1) * ( fp1*(x-xm1)/(xp1-x) - fm1*(xp1-x)/(x-xm1)) + \
            f*(xp1-2*x+xm1) / ( (x-xm1)*(xp1-x))
    else:
        df_dx = min( (fp1-f)/(xp1-x), (f-fm1)/(x-xm1)  )
    return df_dx

#-- interpolation by polynomial of degree 3
def __P1(x,x0,x1): return  (x-x1)**2 * (2*x-3*x0+x1) / (x1-x0)**3
def __P2(x,x0,x1): return -(x-x0)**2 * (2*x-3*x1+x0) / (x1-x0)**3
def __P3(x,x0,x1): return  (x-x0)    * (x-x1)**2     / (x1-x0)**2
def __P4(x,x0,x1): return  (x-x0)**2 * (x-x1)        / (x1-x0)**2

def local_interpolation(newx,oldx,oldy,full_output=False):
    """
    A local interpolation method by a polynomial of degree 3.
    
    After Marc-Antoine Dupret, 2002 (extended version of PhD thesis).
    
    Cannot extrapolate!
    
    >>> np.random.seed(1114)
    >>> oldx = np.sort(np.random.uniform(size=10))
    >>> oldy = oldx**2#np.random.uniform(size=10)
    >>> oldy[5:] = -oldx[5:]**2
    >>> newx = np.linspace(oldx.min(),oldx.max(),1000)
    >>> newy,disconts = local_interpolation(newx,oldx,oldy,full_output=True)
    
    >>> sharpy = newy.copy()
    >>> sharpy[disconts] = np.nan
    >>> smoothy = newy.copy()
    >>> smoothy[-disconts] = np.nan
    >>> p = pl.figure()
    >>> p = pl.plot(oldx,oldy,'ks-',ms=10)
    >>> p = pl.plot(newx,smoothy,'go-',lw=2,ms=2,mec='g')
    >>> p = pl.plot(newx,sharpy,'ro-',lw=2,ms=2,mec='r')
    
    @param newx: new x-array to interpolate on
    @type newx: ndarray
    @param oldx: old x-array to interpolate from
    @type oldx: ndarray
    @param oldy: old y-array to interpolate from
    @type oldy: ndarray
    @param full_output: also return points where discontinuities were detected
    @type full_output: boolean
    @return: interpolated array(, discontinuities)
    @rtype: ndarray(,ndarray)
    """
    #-- extend axis to be able to interpolate last point
    lastx = 2*oldx[-1]-oldx[-2]
    lasty = (oldy[-1]-oldy[-2])/(oldx[-1],oldx[-2])*(oldx[-1]-lastx) + oldy[-1]
    oldy = np.hstack([oldy,lasty])
    oldx = np.hstack([oldx,lastx])
    #-- prepare new y array
    newy = np.zeros_like(newx)
    #-- keep information on where sharp features occur, if asked
    if full_output:
        disconts = np.zeros(len(newx),bool)
        #disconts = np.zeros(len(newx))
    
    for i,x in enumerate(newx):
        index = oldx.searchsorted(x)
        #if index>=(len(oldx)-1): continue
        x0,f0 = oldx[index-1],oldy[index-1]
        x1,f1 = oldx[index],oldy[index]
        x2,f2 = oldx[index-2],oldy[index-2]
        
        #-- check sharpness of feature
        sharpness = 1./((f1-f0)/(x1-x0)*(x1-x2)/(f1-f2))
        sharp = (0.2<=sharpness<=0.5) 
        if full_output:
            disconts[i] = sharp#ness#sharp
        #-- preliminary estimation of df/dx
        dfdx0 = __df_dx(oldx,oldy,index-1,sharp=sharp)
        dfdx1 = __df_dx(oldx,oldy,index,sharp=sharp)
        
        
        #-- interpolation by polynomial of degree 3
        P1 =  (x-x1)**2 * (2*x-3*x0+x1) / (x1-x0)**3
        P2 = -(x-x0)**2 * (2*x-3*x1+x0) / (x1-x0)**3
        P3 =  (x-x0)    * (x-x1)**2     / (x1-x0)**2
        P4 =  (x-x0)**2 * (x-x1)        / (x1-x0)**2
        
        #-- interpolating polynomial
        Px = f0*P1 + f1*P2 + dfdx0*P3 + dfdx1*P4
        newy[i] = Px
    if full_output:
        return newy,disconts
    else:
        return newy


def local_interpolation_ND(newx,oldx,oldy):
    #-- oldx should be a list of [x,y,z,...]
    #   oldy should be array of N(x) x N(y) x N(z) x ...
    points = np.zeros((len(args),3))
    for _x in args:
        index = oldx.searchsorted(_x)
        x0,f0 = oldx[index-1],oldy[index-1]
        x1,f1 = oldx[index],oldy[index]
        x2,f2 = oldx[index-2],oldy[index-2]
    
    polynomials = []
    

def create_pixeltypegrid(grid_pars,grid_data):
    """
    Creates pixelgrid and arrays of axis values.
    
    Starting from:
        * grid_pars: 2D numpy array, 1 column per parameter, unlimited number of cols
        * grid_data: 2D numpy array, 1 column per variable, data corresponding to the rows in grid_pars
    
    The grid should be rectangular and complete, i.e. every combination of the unique values in the 
    parameter columns should exist. If not, a nan value will be inserted.
    
    @param grid_pars: Npar x Ngrid array of parameters
    @type grid_pars: array
    @param grid_data: Ndata x Ngrid array of data
    @type grid_data:array
    @return: axis values and pixelgrid
    @rtype: array, array
    """

    uniques = [np.unique(column, return_inverse=True) for column in grid_pars]
    #[0] are the unique values, [1] the indices for these to recreate the original array

    axis_values = [uniques_[0] for uniques_ in uniques]
    unique_val_indices = [uniques_[1] for uniques_ in uniques]
    
    data_dim = np.shape(grid_data)[0]

    par_dims   = [len(uv[0]) for uv in uniques]

    par_dims.append(data_dim)
    pixelgrid = np.ones(par_dims)
    
    # We put np.inf as default value. If we get an inf, that means we tried to access
    # a region of the pixelgrid that is not populated by the data table
    pixelgrid[pixelgrid==1] = np.inf
    
    # now populate the multiDgrid
    indices = [uv[1] for uv in uniques]
    pixelgrid[indices] = grid_data.T
    return axis_values, pixelgrid

def interpolate(p, axis_values, pixelgrid):
    """
    Interpolates in a grid prepared by create_pixeltypegrid().
    
    p is an array of parameter arrays
    
    @param p: Npar x Ninterpolate array
    @type p: array
    @return: Ndata x Ninterpolate array
    @rtype: array
    """
    # convert requested parameter combination into a coordinate
    p_ = [np.searchsorted(av_,val) for av_, val in zip(axis_values,p)]
    lowervals_stepsize = np.array([[av_[p__-1], av_[p__]-av_[p__-1]] \
                         for av_, p__ in zip(axis_values,p_)])
    p_coord = (p-lowervals_stepsize[:,0])/lowervals_stepsize[:,1] + np.array(p_)-1


    # interpolate
    return np.array([ndimage.map_coordinates(pixelgrid[...,i],p_coord, order=1, prefilter=False) \
                for i in range(np.shape(pixelgrid)[-1])])




if __name__=='__main__':
    import pylab as pl
    from doctest import testmod
    testmod()
    pl.show()
    