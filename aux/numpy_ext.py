# -*- coding: utf-8 -*-
"""
Operations on numpy arrays not present in the standard package.
"""
import numpy as np
import pylab as pl

#{ Normal arrays
def unique_arr(a,axis=0,return_index=False):
    """
    Distil unique rows/cols from an array
    
    C{axis=0}: unique rows
    C{axis=1}: unique cols
    
    Example with rows:
    >>> c = np.sort(np.random.normal(size=(3,5)),axis=0)
    >>> d = np.r_[c,c,c]
    >>> du = np.sort(unique_arr(d,axis=0),axis=0)
    >>> np.all(du==c)
    True
    
    Example with columns:
    >>> c = np.sort(np.random.normal(size=(3,5)),axis=1)
    >>> d = np.hstack([c,c])
    >>> du = np.sort(unique_arr(d,axis=1),axis=1)
    >>> np.all(du==c)
    True
    
    @param a: array to remove duplicate entries from
    @type a: numpy array
    @param axis: keep unique elements in rows (axis=0) or columns (axis=1)
    @type axis: integer
    @param return_index: return index array with unique elements
    @type return_index: bool
    @return: unique array(,indices)
    @rtype: numpy array(, numpy array)
    """
    
    if axis==0:
        a = np.ascontiguousarray(a)
        a_,index = np.unique1d(a.view([('',a.dtype)]*a.shape[1]),return_index=True)
        a = a_.view(a.dtype).reshape(-1,a.shape[1])
    elif axis==1:
        a = np.transpose(a)
        a = np.ascontiguousarray(a)
        a_,index = np.unique1d(a.view([('',a.dtype)]*a.shape[1]),return_index=True)
        a = a_.view(a.dtype).reshape(-1,a.shape[1])
        a = np.transpose(a)
    if return_index:
        return a,index
    else:
        return a

def sort_order(a,order):
    """
    Sort an array along several axes
    
    >>> a = np.array([[ 0.,  1.],\
                   [ 1.,  1.],\
                   [ 1.,  0.],\
                   [ 0.,  0.],\
                   [ 0.,  3.],\
                   [ 1.,  0.],\
                   [ 1.,  3.],\
                   [ 1.,  2.],\
                   [ 1.,  3.],\
                   [ 0.,  4.]])
    >>> sort_order(a,[0,1])
    array([[ 0.,  0.],
           [ 0.,  1.],
           [ 0.,  3.],
           [ 0.,  4.],
           [ 1.,  0.],
           [ 1.,  0.],
           [ 1.,  1.],
           [ 1.,  2.],
           [ 1.,  3.],
           [ 1.,  3.]])
    
    @param a: numpy array to sort
    @type a: numpy array
    @param order: order of the columns to sort
    @type order: list of integers (indices)
    @return: sorted array
    @rtype: numpy array
    """
    a = np.ascontiguousarray(a.copy())+0.
    a_ = np.sort(a.view([('',a.dtype)]*a.shape[1]), order=['f%d'%(i) for i in order], axis=0)
    a = a_.view(a.dtype).reshape(-1,a.shape[1])
    return a


def match_arrays(a,b):
    """
    Return closest-match indices from b in a.
    
    Example usage:
        >>> a = np.random.uniform(size=10)
        >>> b = a[np.array(np.random.uniform(size=3)*10,int)]
        >>> ind = match_arrays(a,b)
        >>> all(a[ind] == b)
        True
    """
    sa = np.argsort(a)
    a_ = a[sa]
    a_ = a_[:-1] + np.diff(a_)/2.
    closest_index = np.searchsorted(a_,b)
    indices = sa[closest_index]
    return indices

def stdw(data,weights=None): 
    """
    Calculates the weighted (sample) standard deviation of a list of numbers.  
     
    @type data: list 
    @param data: input data, must be a two dimensional list in format [value, weight] 
    @rtype: float 
    @return: weighted standard deviation 
    """ 
    if len(data)==1:
        return 0
    listMean=np.average(data,weights=weights) 
    sum=0 
    wSum=0 
    wNonZero=0 
    for el,w in zip(data,weights): 
        if w>0.0: 
            sum=sum+float((el-listMean)/w)*float((el-listMean)/w) 
            wSum=wSum+float(1.0/w)*float(1.0/w) 
             
    nFactor=float(len(data))/float(len(data)-1) 
    stdev=np.sqrt(nFactor*(sum/wSum))
    
    return stdev    
    
    
def deriv(x,y):
    """
    3 point Lagrangian differentiation.
    
    Returns z = dy/dx
    
    Example usage:
    
    >>> X = np.array([ 0.1, 0.3, 0.4, 0.7, 0.9])
    >>> Y = np.array([ 1.2, 2.3, 3.2, 4.4, 6.6])
    >>> deriv(X,Y)
    array([  3.16666667,   7.83333333,   7.75      ,   8.2       ,  13.8       ])
    """
    #-- body derivation
    x0_x1 = np.roll(x,1)-x
    x1_x2 = x-np.roll(x,-1)
    x0_x2 = np.roll(x,1)-np.roll(x,-1)
    derivee = np.roll(y,1)*x1_x2/(x0_x1*x0_x2)
    derivee = derivee + y*(1/x1_x2-1/x0_x1)
    derivee = derivee - np.roll(y,-1)*x0_x1/(x0_x2*x1_x2)
    #-- edges
    derivee[0]=y[0]*(1./x0_x1[1]+1./x0_x2[1])
    derivee[0]=derivee[0]-y[1]*x0_x2[1]/(x0_x1[1]*x1_x2[1])
    derivee[0]=derivee[0]+y[2]*x0_x1[1]/(x0_x2[1]*x1_x2[1])
    nm3=len(x)-3
    nm2=len(x)-2
    nm1=len(x)-1
    derivee[nm1]=-y[nm3]*x1_x2[nm2]/(x0_x1[nm2]*x0_x2[nm2])
    derivee[nm1]=derivee[nm1]+y[nm2]*x0_x2[nm2]/(x0_x1[nm2]*x1_x2[nm2])
    derivee[nm1]=derivee[nm1]-y[nm1]*(1./x0_x2[nm2]+1./x1_x2[nm2])
    
    return derivee


#def deriv_noise(x,y):
    #"""
    #Compact finite difference scheme on non-uniform mesh.
    
    #See Gamet, Ducros, Nicoud et al., 1999.
    #"""
    #alpha = 
    #beta = 
    #hi_1 = 
    #hi0 = 
    #hi1 = 
    #hi_1*

#}
#{ Record arrays

def recarr(x,mydtype):
    """
    Convert an array to a record array.
    dtype = [('name1',int),('name2',float)]
    
    >>> x = np.array([[1,1],[2,2]])
    >>> y = recarr(x,[('ones',int),('twos',int)])
    >>> print y['ones']
    [1 1]
    
    @param x: array to convert
    @type x: numpy array
    @param mydtype: dtype of record array
    @type mydtype: list of tuples ('column name',type)
    @return: convert record array
    @rtype: numpy record array
    """
    y = [tuple([col[i] for col in x]) for i in range(len(x[0]))]
    y = np.array(y,dtype=mydtype)
    return y

def recarr_addrows(x,rows):
    """
    Add rows to a record array
    
    >>> x = np.array([[1,1],[2,2]])
    >>> y = recarr(x,[('ones',int),('twos',int)])
    >>> z = recarr_addrows(y,[[1,2]])
    >>> print z['ones']
    [1 1 1]
    
    @param x: original record array
    @type x: numpy record array
    @param rows: list of lists/tuples or 2D-numpy array
    @type rows: list or ndarray
    @return: extended record array
    @rtype: numpy record array
    """
    rows = [tuple(row) for row in rows]
    x = np.core.records.fromrecords(x.tolist()+rows,dtype=x.dtype)
    return x

def recarr_addcols(x,cols,dtypes_ext):
    """
    Add columns to a record array
    
    >>> x = np.array([[1,1],[2,2]])
    >>> y = recarr(x,[('ones',int),('twos',int)])
    >>> z = recarr_addcols(y,[[3,3],[4,4]],[('threes',int),('fours',int)])
    >>> print z['fours']
    [4 4]
    
    @param x: original record array
    @type x: numpy record array
    @param cols: list of lists/tuples or 2D-numpy array
    @type cols: list or ndarray
    @param dtypes_ext: dtype of extra columns
    @type dtypes_ext: list of tuples ('column name',type)
    @return: extended record array
    @rtype: numpy record array
    """
    names = list(x.dtype.names)
    dtypes = [(name,x.dtype[names.index(name)].str) for name in names]
    dtypes += dtypes_ext
    rows = []
    for i in range(len(x)):
        rows.append(list(x[i]) + [col[i] for col in cols])
    x = np.core.records.fromrecords(rows,dtype=dtypes)
    return x

def recarr_join(arr1,arr2):
    """
    Join to record arrays column wise.
    """
    arr1 = arr1.copy()
    for field in arr2.dtype.names:
        arr1 = pl.mlab.rec_append_fields(arr1,field,arr2[field])
    return arr1

#}

if __name__=="__main__":
    import doctest
    doctest.testmod()