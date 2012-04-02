"""
Non-standard interpolation methods.
"""
import numpy as np

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

def local_interpolation(newx,oldx,oldy,full_output=False):
    """
    A local interpolation method by a polynomial of degree 3.
    
    Cannot extrapolate!
    
    >>> np.random.seed(1114)
    >>> oldx = np.sort(np.random.uniform(size=10))
    >>> oldy = np.random.uniform(size=10)
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
        disconts = np.zeros_like(newx,bool)
    
    for i,x in enumerate(newx):
        index = oldx.searchsorted(x)
        #if index>=(len(oldx)-1): continue
        x0,f0 = oldx[index-1],oldy[index-1]
        x1,f1 = oldx[index],oldy[index]
        x2,f2 = oldx[index-2],oldy[index-2]
        
        #-- check sharpness of feature
        sharpness = 1./((f1-f0)/(x1-x0)*(x1-x2)/(f1-f2))
        sharp = -(0.2<=sharpness<=0.5) 
        if full_output:
            disconts[i] = sharp
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
        

if __name__=='__main__':
    import pylab as pl
    from doctest import testmod
    testmod()
    pl.show()
    