# -*- coding: utf-8 -*-
"""
Operations on numpy arrays not present in the standard package.
"""
import numpy as np

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
#}

if __name__=="__main__":
    import doctest
    doctest.testmod()