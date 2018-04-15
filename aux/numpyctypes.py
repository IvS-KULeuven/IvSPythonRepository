"""
Module to convert a numpy array to a ctypes struct.
This struct can then be passed to a native C library.

Author: Joris De Ridder
"""

import numpy as np
import ctypes as C



ctypesDict = {'d' : C.c_double,
              'b' : C.c_char,
              'h' : C.c_short,
              'i' : C.c_int,
              'l' : C.c_long,
              'q' : C.c_longlong,
              'B' : C.c_ubyte,
              'H' : C.c_ushort,
              'I' : C.c_uint,
              'L' : C.c_ulong,
              'Q' : C.c_ulonglong}






def c_ndarray(a, dtype = None, ndim = None, shape = None, requirements = None):

    """
    Returns a ctypes structure of the array 'a'
    containing the arrays info (data, shape, strides, ndim).
    A check is made to ensure that the array has the specified dtype
    and requirements.

    Example:

    >>> myArray = np.arange(10.0)
    >>> myCstruct = c_ndarray(myArray, dtype=np.double, ndim = 3, shape = (4,3,2),
    ...                       requirements = ['c_contiguous'])

    @param a: the numpy array to be converted
    @type a: ndarray
    @param dtype: the required dtype of the array, convert if it doesn't match
    @type dtype: numpy dtype
    @param ndim: the required number of axes of the array,
                 complain if it doesn't match
    @type ndim: integer
    @param shape: required shape of the array, complain if it doesn't match
    @type shape: tuple
    @param requirements: "ensurearray", "aligned", "fortran", "f_contiguous",
                         or "c_contiguous". Convert if it doesn't match.
    @type requirements: list
    @return: ctypes structure with the fields:
                - data: pointer to the data : the type is determined with the
                        dtype of the array, and with ctypesDict.
                - shape: pointer to long array : size of each of the dimensions
                - strides: pointer to long array : strides in elements (not bytes)
    @rtype: ctypes structure

    """

    if not requirements:

        # Also allow derived classes of ndarray

        array = np.asanyarray(a, dtype=dtype)

    else:

        # Convert requirements to captial letter codes:
        # (ensurearray' -> 'E';  'aligned' -> 'A'
        #  'fortran', 'f_contiguous', 'f' -> 'F'
        #  'contiguous', 'c_contiguous', 'c' -> 'C')

        requirements = [x[0].upper() for x in requirements]
        subok = (0 if 'E' in requirements else 1)

        # Make from 'a' an ndarray with the specified dtype, but don't copy the
        # data (yet). This also ensures that the .flags attribute is present.

        array = np.array(a, dtype=dtype, copy=False, subok=subok)

        # See if copying all data is really necessary.
        # Note: 'A' = (A)ny = only (F) it is was already (F)

        copychar = 'A'
        if 'F' in requirements:
            copychar = 'F'
        elif 'C' in requirements:
            copychar = 'C'

        for req in requirements:
            if not array.flags[req]:
                array = array.copy(copychar)
                break


    # If required, check the number of axes and the shape of the array

    if ndim is not None:
        if array.ndim != ndim:
            raise TypeError("Array has wrong number of axes")

    if shape is not None:
        if array.shape != shape:
            raise TypeError("Array has wrong shape")

    # Define a class that serves as interface of an ndarray to ctypes.
    # Part of the type depends on the array's dtype.

    class ndarrayInterfaceToCtypes(C.Structure):
        pass

    typechar = array.dtype.char

    if typechar in ctypesDict:
        ndarrayInterfaceToCtypes._fields_ =                       \
                  [("data", C.POINTER(ctypesDict[typechar])),
                   ("shape" , C.POINTER(C.c_long)),
                   ("strides", C.POINTER(C.c_long))]
    else:
        raise TypeError("dtype of input ndarray not supported")

    # Instantiate the interface class and attach the ndarray's internal info.
    # Ctypes does automatic conversion between (c_long * #) arrays and POINTER(c_long).

    ndarrayInterface = ndarrayInterfaceToCtypes()
    ndarrayInterface.data = array.ctypes.data_as(C.POINTER(ctypesDict[typechar]))
    ndarrayInterface.shape = (C.c_long * array.ndim)(*array.shape)
    ndarrayInterface.strides = (C.c_long * array.ndim)(*array.strides)
    for n in range(array.ndim):
        ndarrayInterface.strides[n] /= array.dtype.itemsize

    return ndarrayInterface


