# -*- coding: utf-8 -*-
"""
Read and write ASCII files.
"""
import gzip
import logging
import os

import numpy as np

logger = logging.getLogger("IO.ASCII")

#{ Input
def read2list(filename,**kwargs):
    """
    Load an ASCII file to list of lists.
    
    The comments and data go to two different lists.
    
    Also opens gzipped files.
    
    @param filename: name of file with the data
    @type filename: string
    @keyword commentchar: character(s) denoting comment rules
    @type commentchar: list of str
    @keyword splitchar: character seperating entries in a row (default: whitespace)
    @type splitchar: str or None
    @return: list of lists (data rows)
             list of lists (comments lines without commentchar),
    @rtype: (list,list)
    """
    commentchar = kwargs.get('commentchar',['#'])
    splitchar = kwargs.get('splitchar',None)
    
    if os.path.splitext(filename)[1] == '.gz':
        ff = gzip.open(filename)
    else:
        ff = open(filename)
        
    data = []  # data
    comm = []  # comments
    
    while 1:  # might call read several times for a file
        line = ff.readline()
        if not line: break  # end of file
        
        #-- when reading a comment line
        if line[0] in commentchar:
            comm.append(line[1:].strip())
            continue # treat next line
        
        #-- when reading data, split the line
        data.append(line.strip().split(splitchar))
    ff.close()
    
    #-- report that the file has been read
    logger.debug('Data file %s read'%(filename))
    
    #-- and return the contents
    return data,comm

def read2array(filename,**kwargs):
    """
    Load ASCII file to a numpy array.
    
    For a list of extra keyword arguments, see C{<read2list>}.
    
    If you want to return a list of the columns instead of rows, just do
    
    C{>>> col1,col2,col3 = ascii.read2array(myfile).T}
    
    @param filename: name of file with the data
    @type filename: string
    @keyword dtype: type of numpy array (default: float)
    @type dtype: numpy dtype
    @keyword return_comments: flag to return comments (default: False)
    @type return_comments: bool
    @return: data array (, list of comments)
    @rtype: ndarray (, list)
    """
    dtype = kwargs.get('dtype',np.float)
    return_comments = kwargs.get('return_comments',False)
    data,comm = read2list(filename,**kwargs)
    
    data = np.array(data,dtype=dtype)
    return return_comments and (data,comm) or data

def read2recarray(filename,**kwargs):
    """
    Load ASCII file to a numpy record array.
    
    For a list of extra keyword arguments, see C{<read2list>}.
    
    FI dtypes is None, we have some room to automatically detect the contents
    of the columns. This is not implemented yet.
    
    the keyword 'dtype' should be equal to a list of tuples, e.g.
    
    C{dtype = [('col1','a10'),('col2','>f4'),..]}
    
    @param filename: name of file with the data
    @type filename: string
    @keyword dtype: dtypes of record array 
    @type dtype: list of tuples
    @keyword return_comments: flag to return comments (default: False)
    @type return_comments: bool
    @return: data array (, list of comments)
    @rtype: ndarray (, list)
    """
    dtype = kwargs.get('dtype',None)
    return_comments = kwargs.get('return_comments',False)

    #-- first read in as a normal array
    data,comm = read2list(filename,**kwargs)
    data = np.array(data,dtype=np.str).T
    
    #-- if dtypes is None, we have some room to automatically detect the contents
    #   of the columns. This is not implemented yet.
    #-- cast all columns to the specified type
    dtype = np.dtype(dtype)
    data = [np.cast[dtype[i]](data[i]) for i in range(len(data))]
    #-- and build the record array
    data = np.rec.array(data, dtype=dtype)
    return return_comments and (data,comm) or data
#}

#{ Output

def write_array(data, filename, **kwargs):
    """
    Save a numpy array to an ASCII file.
    
    Add comments via keyword comments (a list of strings denoting every comment
    line). By default, the comment lines will be preceded by the C{commentchar}.
    If you want to override this behaviour, set C{commentchar=''}.
    
    @keyword header: optional header for column names
    @type header: list of str
    @keyword comments: comment lines
    @type comments: list of str
    @keyword commentchar: comment character
    @type commentchar: str
    @keyword sep: separator for the columns and header names
    @type sep: str
    @keyword axis0: string denoting the orientation of the matrix. If you gave
    a list of columns, set C{axis0='cols'}, otherwise C{axis='rows'} (default).
    @type axis0: str, one of C{cols}, C{rows}.
    @keyword mode: file mode (a for appending, w for (over)writing...)
    @type mode: char (one of 'a','w'...)
    """
    header = kwargs.get('header',[])
    comments = kwargs.get('comments',None)
    commentchar = kwargs.get('commentchar','#')
    sep = kwargs.get('sep',' ')
    axis0 = kwargs.get('axis0','rows')
    mode = kwargs.get('model','w')
    
    #-- switch to rows first if a list of columns is given
    if axis0.lower()!='rows':
        data = data.T
    
    ff = open(filename,mode)
    if header:
        ff.write('#'+sep.join(header)+'\n')
    
    #-- write to a file
    for row in data:
        ff.write(sep.join(['%g'%(col) for col in row])+'\n')
    ff.close()
        
            
        
    
    
    
    

#}