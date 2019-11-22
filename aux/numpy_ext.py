# -*- coding: utf-8 -*-
"""
Operations on numpy arrays not present in the standard package.
"""
import numpy as np
import pylab as pl
from scipy import spatial
import itertools

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
        a_,index = np.unique(a.view([('',a.dtype)]*a.shape[1]),return_index=True)
        a = a_.view(a.dtype).reshape(-1,a.shape[1])
    elif axis==1:
        a = np.transpose(a)
        a = np.ascontiguousarray(a)
        a_,index = np.unique(a.view([('',a.dtype)]*a.shape[1]),return_index=True)
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

def argmax2D(a):
    """
    Calculate the argmax of a 2D array.

    Example usage:

    >>> output = np.zeros(100)
    >>> for i in xrange(100):
    ...     a = np.random.normal(size=(5,6))
    ...     x,y = argmax2D(a)
    ...     output[i] = (a[x,y] == a.max())
    >>> sum(output)
    100.0

    @param a: array to calculate the position of the maximum of
    @type a: numpy 2D array
    @rtype: (index,index)
    @return: x,y coordinates
    """
    index = np.argmax(a)
    y,x = index%a.shape[1],index//a.shape[1]
    return x,y

def argmin2D(a):
    """
    Calculate the argmin of a 2D array.

    Example usage:

    >>> output = np.zeros(100)
    >>> for i in xrange(100):
    ...     a = np.random.normal(size=(5,6))
    ...     x,y = argmin2D(a)
    ...     output[i] = (a[x,y] == a.min())
    >>> sum(output)
    100.0

    @param a: array to calculate the position of the minimum of
    @type a: numpy 2D array
    @rtype: (index,index)
    @return: x,y coordinates
    """
    index = np.argmin(a)
    y,x = index%a.shape[1],index//a.shape[1]
    return x,y


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

def random_rectangular_grid(gridpoints,size):
    """
    Generate random points in a non-convex continuous rectangular grid.

    This routine will subdivide the grid in smaller cubicles and then populate
    those with a number of points proportional to the relative size of the
    cubicle. Because of this, C{size} will never be the true size of the
    sample returned, but will always be a little bit higher or lower.

    B{Warning: This function requiresevery square to be populated by at least
    one point. If this is not the case, it will be forced that every square
    has a point.}

    >>> xs = np.array([1,2,0,1,2,3,0,1,2,3,1,2.])
    >>> ys = np.array([4,4,3,3,3,3,2,2,2,2,1,1.])
    >>> gridpoints = np.column_stack([xs,ys])
    >>> sample = random_rectangular_grid(gridpoints,10000)

    >>> p = pl.figure()
    >>> p = pl.plot(xs,ys,'ro')
    >>> p = pl.plot(sample[:,0],sample[:,1],'ko',ms=2)

    ]include figure]]ivs_aux_numpy_ext_grid.png]

    Gridpoints should be Ngridpoints x Ncols
    """
    #-- make the KDTree, axis values and axis indices
    tree = spatial.KDTree(gridpoints)
    axis_values = [np.sort(np.unique(col)) for col in gridpoints.T]
    axis_indices = [np.arange(1,len(axis)) for axis in axis_values]

    #-- we need to find the total volume and then we need to sample that volume
    #   uniformly. We can only include grid cubes that have corners that are all
    #   present in the grid
    hypercube_sides  = [abs(np.diff(axis)) for axis in axis_values]
    hypercube_volume = 0.
    random_ranges_sizes = []
    for combination,limit_indices in zip(itertools.product(*hypercube_sides),itertools.product(*axis_indices)):
        limit_indices = np.array(limit_indices)
        #-- create all corners of the particular grid cube, they need to be in the
        #   KDTree!
        for corner in itertools.product(*list(zip(limit_indices-1,limit_indices))):
            dist = tree.query([axis_values[i][j] for i,j in enumerate(corner)])[0]
            if dist!=0:
                break
        #-- if all the corners of the cube are in the grid, compute the volume and
        #   remember the lower and upper limits to generate random uniform points.
        else:
            low = [axis_values[i][j-1] for i,j in enumerate(limit_indices)]
            high = [axis_values[i][j] for i,j in enumerate(limit_indices)]
            volume = np.product(combination)
            hypercube_volume+= volume
            random_ranges_sizes.append([low,high,volume])
    #-- now we're ready to actually generate the grid

    sample = np.vstack([np.random.uniform(low=low,high=high,size=(int(max(vol/hypercube_volume*size,1)),len(low))) for low,high,vol in random_ranges_sizes])

    return sample

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
        rows.append(tuple(list(x[i]) + [col[i] for col in cols]))
    x = np.core.records.fromrecords(rows,dtype=dtypes)
    return x

def recarr_join(arr1,arr2):
    """
    Join to record arrays column wise.
    """
    arr1 = arr1.copy()
    for field in arr2.dtype.names:
        arr1 = np.lib.recfunctions.rec_append_fields(arr1,field,arr2[field])
    return arr1

#}

#{ Line intersections

def find_intersections(A, B):
    """
    Find intersection of two 2D lines.

    Example usage:

    >>> pcoeff1 = [1,2,-2.]
    >>> pcoeff2 = [-2.12,2.45,3.7321]
    >>> x = np.linspace(-3,3,7)
    >>> A = np.column_stack([x,np.polyval(pcoeff1,x)])
    >>> B = np.column_stack([x,np.polyval(pcoeff2,x)])

    We find two intersections:

    >>> xs,ys = find_intersections(A,B)
    >>> print xs,ys
    [-1.22039755  1.34367003] [-2.77960245  2.71835017]

    >>> p = pl.figure()
    >>> p = pl.plot(x,A[:,1],'bo-')
    >>> p = pl.plot(x,B[:,1],'go-')
    >>> p = pl.plot(xs,ys,'rs',ms=5)

    Returns empty arrays when there are no intersections.

    ]include figure]]ivs_aux_numpy_ext_intersect.png]
    """
    # min, max and all for arrays
    amin = lambda x1, x2: np.where(x1<x2, x1, x2)
    amax = lambda x1, x2: np.where(x1>x2, x1, x2)
    aall = lambda abools: np.dstack(abools).all(axis=2)
    slope = lambda line: (lambda d: d[:,1]/d[:,0])(np.diff(line, axis=0))

    x11, x21 = np.meshgrid(A[:-1, 0], B[:-1, 0])
    x12, x22 = np.meshgrid(A[1:, 0], B[1:, 0])
    y11, y21 = np.meshgrid(A[:-1, 1], B[:-1, 1])
    y12, y22 = np.meshgrid(A[1:, 1], B[1:, 1])

    m1, m2 = np.meshgrid(slope(A), slope(B))
    m1inv, m2inv = 1/m1, 1/m2

    yi = (m1*(x21-x11-m2inv*y21) + y11)/(1 - m1*m2inv)
    xi = (yi - y21)*m2inv + x21

    xconds = (amin(x11, x12) < xi, xi <= amax(x11, x12),
              amin(x21, x22) < xi, xi <= amax(x21, x22) )
    yconds = (amin(y11, y12) < yi, yi <= amax(y11, y12),
              amin(y21, y22) < yi, yi <= amax(y21, y22) )

    return xi[aall(xconds)], yi[aall(yconds)]

#}

if __name__=="__main__":
    import doctest
    doctest.testmod()
    pl.show()
