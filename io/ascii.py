# -*- coding: utf-8 -*-
"""
Read and write ASCII files.
"""
import gzip
import logging
import os
import re
import numpy as np

logger = logging.getLogger("IO.ASCII")
logger.addHandler(logging.NullHandler())

#{ General Input
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
    @keyword skip_empty: skip empty lines
    @type skip_empty: bool
    @keyword skip_lines: skip nr of lines (including comment and empty lines)
    @type skip_lines: integer
    @return: list of lists (data rows)
             list of lists (comments lines without commentchar),
    @rtype: (list,list)
    """
    commentchar = kwargs.get('commentchar',['#'])
    splitchar = kwargs.get('splitchar',None)
    skip_empty = kwargs.get('skip_empty',True)
    skip_lines = kwargs.get('skip_lines',0)
    
    if os.path.splitext(filename)[1] == '.gz':
        ff = gzip.open(filename)
    else:
        ff = open(filename)
        
    data = []  # data
    comm = []  # comments
    
    line_nr = -1
    
    #-- fixwidth split or character split?
    if splitchar is None or isinstance(splitchar,str):
        fixwidth = False
    else:
        fixwidth = True
    
    while 1:  # might call read several times for a file
        line = ff.readline()
        if not line: break  # end of file
        line_nr += 1
        if line_nr<skip_lines:        
            continue
        
        #-- strip return character from line
        if skip_empty and line.isspace():
            continue # empty line
        
        #-- remove return characters
        line = line.replace('\n','')
        #-- when reading a comment line
        if line[0] in commentchar:
            comm.append(line[1:])
            continue # treat next line
        
        #-- when reading data, split the line
        if not fixwidth:
            data.append(line.split(splitchar))
        else:
            data.append(fw2python(line,splitchar))
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
    
    IF dtypes is None, we have some room to automatically detect the contents
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
    splitchar = kwargs.get('splitchar',None)
    
    #-- first read in as a normal array
    data,comm = read2list(filename,**kwargs)
    
    #-- if dtypes is None, we have some room to automatically detect the contents
    #   of the columns. This is not fully implemented yet, and works only
    #   if the second-to-last and last columns of the comments denote the
    #   name and dtype, respectively
    if dtype is None:
        data = np.array(data,dtype=str).T
        header = comm[-2].replace('|',' ').split()
        types = comm[-1].replace('|','').split()
        dtype = [(head,typ) for head,typ in zip(header,types)]
        dtype = np.dtype(dtype)
    elif isinstance(dtype,list):
        data = np.array(data,dtype=str).T
        dtype = np.dtype(dtype)
    #-- if dtype is a list, assume it is a list of fixed width stuff.
    elif isinstance(splitchar,list):
        types,lengths = fws2info(splitchar)
        dtype = []
        names = range(300)
        for i,(fmt,lng) in enumerate(zip(types,lengths)):
            if fmt.__name__=='str':
                dtype.append((str(names[i]),(fmt,lng)))
            else:
                dtype.append((str(names[i]),fmt))
        dtype = np.dtype(dtype)
    
    #-- cast all columns to the specified type
    data = [np.cast[dtype[i]](data[i]) for i in range(len(data))]
        
    #-- and build the record array
    data = np.rec.array(data, dtype=dtype)
    return return_comments and (data,comm) or data
#}

#{ General Output

def write_array(data, filename, **kwargs):
    """
    Save a numpy array to an ASCII file.
    
    Add comments via keyword comments (a list of strings denoting every comment
    line). By default, the comment lines will be preceded by the C{commentchar}.
    If you want to override this behaviour, set C{commentchar=''}.
    
    If you give a record array, you can simply set C{header} to C{True} to write
    the header, instead of specifying a list of strings.
    
    @keyword header: optional header for column names
    @type header: list of str (or boolean for record arrays)
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
    @keyword auto_width: automatically determine the width of the columns
    @type auto_width: bool
    @keyword formats: formats to use to write each column
    @type formats: list of string formatters
    """
    header = kwargs.get('header',None)
    comments = kwargs.get('comments',None)
    commentchar = kwargs.get('commentchar','#')
    sep = kwargs.get('sep',' ')
    axis0 = kwargs.get('axis0','rows')
    mode = kwargs.get('mode','w')
    auto_width = kwargs.get('auto_width',False)
    formats = kwargs.get('formats',None)
    # use '%g' or '%f' or '%e' for writing floats automatically from record arrays with auto width
    use_float = kwargs.get('use_float','%f') 
    
    #-- switch to rows first if a list of columns is given
    if not isinstance(data,np.ndarray):
        data = np.array(data)
    if not 'row' in axis0.lower():
        data = data.T
    
    if formats is None:
        try:
            formats = [('S' in str(data[col].dtype) and '%s' or use_float) for col in data.dtype.names]
        except TypeError:
            formats = [('S' in str(col.dtype) and '%s' or '%s') for col in data.T]
    #-- determine width of columns: also take the header label into account
    col_widths = []
    #-- for record arrays
    if auto_width is True and header==True:
        for fmt,head in zip(formats,data.dtype.names):
            col_widths.append(max([len('%s'%(fmt)%(el)) for el in data[head]]+[len(head)]))
    #-- for normal arrays and specified header
    elif auto_width is True and header is not None:
        for i,head in enumerate(header):
            col_widths.append(max([len('%s'%(formats[i])%(el)) for el in data[:,i]]+[len(head)]))
    #-- for normal arrays without header
    elif auto_width is True and header is not None:
        for i in range(data.shape[1]):
            col_widths.append(max([len('%s'%(formats[i])%(el)) for el in data[:,i]]))
    
    if header is True:
        col_fmts = [str(data.dtype[i]) for i in range(len(data.dtype))]
        header = data.dtype.names
    else:
        col_fmts = None
    
    ff = open(filename,mode)
    if comments is not None:
        ff.write('\n'.join(comments)+'\n')
    
    #-- WRITE HEADER
    #-- when header is desired and automatic width
    if header is not None and col_widths:
        ff.write('#'+sep.join(['%%%s%ss'%(('s' in fmt and '-' or ''),cw)%(head) for fmt,head,cw in zip(formats,header,col_widths)])+'\n')
    #-- when header is desired
    elif header is not None:
        ff.write('#'+sep.join(header)+'\n')
    
    #-- WRITE COLUMN FORMATS
    if col_fmts is not None and col_widths:
        ff.write('#'+sep.join(['%%%s%ss'%(('s' in fmt and '-' or ''),cw)%(colfmt) for fmt,colfmt,cw in zip(formats,col_fmts,col_widths)])+'\n')
    elif col_fmts is not None:
        ff.write('#'+sep.join(['%%%ss'%('s' in fmt and '-' or '')%(colfmt) for fmt,colfmt in zip(formats,col_fmts)])+'\n')
    
    #-- WRITE DATA
    #-- with automatic width
    if col_widths:
        for row in data:
            ff.write(' '+sep.join(['%%%s%s%s'%(('s' in fmt and '-' or ''),cw,fmt[1:])%(col) for col,cw,fmt in zip(row,col_widths,formats)])+'\n')
    #-- without automatic width
    else:
        for row in data:
            ff.write(sep.join(['%s'%(col) for col in row])+'\n')
    ff.close()
    
    

#}

#{ Helper functions

def fw2info(fixwidth):
    """
    Convert fix width information to python information
  
    >>> fw2info('F13.4')
    ([float], [13])
    >>> fw2info('I2')
    ([int], [2])
    >>> fw2info('6I2')
    ([int, int, int, int, int, int], [2, 2, 2, 2, 2, 2])
    >>> fw2info('5F12.3')
    ([float, float, float, float, float], [12, 12, 12, 12, 12])
    """
    translator = {}
    translator['I'] = int
    translator['F'] = float
    translator['A'] = str
    
    numbers = re.compile('[\d]*[\d\.*]+',re.IGNORECASE)
    formats = re.compile('[a-z]+',re.IGNORECASE)
  
    numbers = numbers.findall(fixwidth)
    formats = formats.findall(fixwidth)
    if len(numbers)==1:
        numbers = ['1'] + numbers
    widths = int(numbers[0])*[int(nr.split('.')[0]) for nr in numbers[-1:]]
    formats= int(numbers[0])*formats
    formats=[translator[fmt] for fmt in formats]
    return formats,widths

def fws2info(fixwidths):
    info = [fw2info(fw) for fw in fixwidths]
    types = []
    for iinfo in info: types += iinfo[0]
    length = []
    for iinfo in info: length += iinfo[1]
    length = [0]+[sum(length[:i+1]) for i in range(len(length))]
    return types,length
    

def fw2python(line,fixwidths,missing_value=0):
    """
    Fixwidth info to python objects
    """
    types,length = fws2info(fixwidths)
    processed_line = []
    for i,itype in enumerate(types):
        number = line[length[i]:length[i+1]]
        if number.isspace():
            processed_line.append(itype(missing_value))
        else:
            processed_line.append(itype(line[length[i]:length[i+1]]))
    return processed_line
#}

#{ Source specific

def read_harps(filename):
    """
    Return spectrum and header from a HARPS normalised file.
    
    @parameter filename: name of the file
    @type filename: str
    @return: wavelengths,flux,header
    @rtype: array,array,dict
    """
    data,comments = read2array(filename,return_comments=True)
    wave,flux = data[:,0],data[:,1]
    header = {'instrument':'harps'}
    for line in comments:
        line = line.split()
        if line and line[0]=='HJD':
            header['HJD'] =float(line[-1])
    return wave,flux,header

def read_mad(infile,add_dp=False,**kwargs):
    """
    Reads .phot from MADROT or .nad from MAD.
    
    MAD is the nonadiabatic pulsation code from Marc-Antoine Dupret.
    
    MADROT is the extension for slowly rotating stars in the traditional
    approximation.
    
    This function serves a generic read function, it reads *.phot and *.nad
    files, both from MAD and MADROT.
    
    Returns list of dicts
    
    @param add_dp: add period spacings information
    @type add_dp: boolean
    @return: star info, blocks
    """
    try:
        if os.path.splitext(infile)[1]=='.phot':
            star_info,blocks = read_photrot(infile,**kwargs)
        elif os.path.splitext(infile)[1]=='.nad':
            star_info,blocks = read_nad(infile,**kwargs)
        else:
            logger.error('Cannot interpret file "%s"'%(infile))
    except:
        logger.error('Error encountered in file %s'%(infile))
        raise
    #-- add period spacing
    if add_dp:
        add_period_spacings(star_info,blocks)
    return star_info,blocks


def read_photrot(infile,teff=None,logg=None):
    """
    Reads .phot output file of MADROT.
    
    Returns list of dictionaries
    
    For quick reading, you can already give a range on teff and logg, and skip
    the rest of file.
        
    if there is a NAD file with the same name (but extension .nad instead of .phot),
    we read in the header of this file too for information on the star::
        
        #   M/M_sun =  1.40      Teff    =   6510.2    Log (L/L_sun) =  0.5168                              
        #   Log g   =  4.2748    R/R_sun =  1.4283     age (y)       =  2.3477E+08                          
        #   X       =  0.70      Z       =  0.020      alphaconv     =  1.80     overshoot =  0.40 
    
   if the filename is given as C{M1.4_0.020_0.40_0.70_001-0097.phot}, we add
   information on mass, Z, overshoot, X and evolution number to star_info,
   except for those quantities already derived from a NAD file.
   
   One dictionary representes information about one pulsation mode.
   
   @param teff: range in effective temperature for file selection (value, sigma)
   @type teff: tuple
   @param logg: range in surface gravity for file selection (value, sigma)
   @type logg: tuple
   """
    dict_keys = ('l', 'm', 'par', 'n', 'f_cd_com', 'f_cd_obs', 'ltilde', 'fT', 
                 'psiT', 'fg', 'psig', 'K', 'omega_puls', 'omega_rot')
    with open(infile,'r') as f:
        line = f.readline()
        #first line gives teff, logg, XYZ, remember this as global information
        #about the star
        logTeff, logg_, Fe_H = line.split()
        logTeff, logg_, Fe_H = float(logTeff), float(logg_), float(Fe_H)
        star_info = dict(teff=10**logTeff,logg=logg_,Fe_H=Fe_H)
        
        #-- read the NAD file if it exists:
        nadfile =  os.path.splitext(infile)[0]+'.nad'
        if os.path.isfile(nadfile):
            star_info.update(read_nad(nadfile,only_header=True)[0])
        else:
            pass
            # not implemented yet.
            
        #second line gives number of blocks
        line = f.readline()
        nr_b = int(line)
        blocks = []
        
        #-- check for fundamental parameters
        valid_funds = True
        if teff is not None:
            if not ((teff[0]-teff[1]) <= star_info['teff'] <= (teff[0]+teff[1])):
                valid_funds = False
        if logg is not None:
            if not ((logg[0]-logg[1]) <= star_info['logg'] <= (logg[0]+logg[1])):
                valid_funds = False
        print star_info['teff'],teff,valid_funds
        #loop over all blocks and store them
        pos = 0
        fail = False
        nrfails = 0
        while 1 and valid_funds:
            line = f.readline()
            if not line: break
            # skip empty lines
            if len(line.strip())==0: continue
            elif pos==0:
                # wait for the line which gives the column info
                if line.strip()[0]=='l':
                    pos = 2
                    continue
            elif pos==2:
                #we are now at the line which gives values for l,m,par,...
                #create dictionary and store values
                d = {}
                for key_, val_ in zip(dict_keys,line.split()):
                    #-- sometimes values can be merged together
                    try:
                        d[key_]=float(val_)
                    except ValueError:
                        fail = True
                        nrfails += 1
                #add spin parameter if not failed
                if not fail:
                    d['s'] = 2*d['omega_rot']/d['omega_puls']
                
                for key_ in ('l', 'm', 'par', 'n'): d[key_]=int(d[key_])
                #create the lists to store B
                d['B']=[]
                d['l_B']=[]
                pos=3
            elif pos==3:
                #reading the 50 next lines which give vector B and associated l-values
                l_B_, B_ = line.split()
                d['B'].append(float(B_))
                d['l_B'].append(int(l_B_))
                if len(d['B'])==50:
                    d['B'] = array(d['B'])
                    d['l_B'] = array(d['l_B'])
                    if not fail: blocks.append(d)
                    pos = 0
                    continue
        
    logger.debug("Read %d modes from MADROT %s (%d modes failed)"%(len(blocks),infile,nrfails))
    
    return star_info,blocks
        

def read_phot(photfile,teff=None,logg=None):
    """
    Reads .phot output file of MAD.
    
    Returns list of dictionaries
    
    One dictionary representes information about one pulsation mode.
    """
    dict_keys = ('l', 'fT','psiT', 'fg', 'psig')
    
    blocks = []
    
    with open(photfile,'r') as f:
        fc = f.readlines()
        #first line gives teff, logg, XYZ, remember this as global information
        #about the star
        logTeff, logg = fc[0].split()[:2]
        logTeff, logg = float(logTeff), float(logg)
        star_info = dict(teff=10**logTeff,logg=logg)
        
        #second line gives number of blocks
        nr_b = int(fc[1])
        blocks = []
        
        #loop over all blocks and store them
        pos = 0
        fail = False
        nrfails = 0
        for line in fc[2:]:
            entries = line.strip().split()
            block = {}
            for i,key in enumerate(dict_keys):
                block[key] = entries[i]
            block['m'] = 0 # only one m-value
            block['n'] = 0 # no information on nodes
            block['f_cd_com'] = 0 # no information on frequency
            block['omega_rot'] = 0 # no rotation
            blocks.append(block)
        
    logger.debug("Read %d modes from MAD %s"%(len(blocks),photfile))
    return star_info,blocks
    
def read_nad(nadfile,only_header=False):
    """
    Reads .nad output file of MAD.
    
    Returns list of dictionaries
    One dictionary representes information about one pulsation mode.
    """
    
    dict_keys = {'l':'l',
                 'm':'m',
                 'n':'n',
                 'freq_com (c/d)':'f_cd_com',
                 'freq_in (c/d)':'f_cd_obs',
                 'freq (c/d)': 'f_cd_obs', # for non rotating nad
                 'freq (micHz)': 'f_hz_obs', # for non rotating nad
                 'ltilde':'l-tilde',
                 'l-tilde':'l-tilde',
                 'Im(omega)':'Im(omega)',
                 'Im(sigma) (micHz)':'Im(sigma)_hz',
                 'Re(omega)':'Re(omega)',
                 'Re_ad(omega)':'Re_ad(omega)',
                 'f_T':'fT',
                 'psi_T':'psiT',
                 'f_g':'fg',
                 'psi_g':'psig',
                 'K':'K',
                 'lifetime (d)':'lifetime',
                 'Q (d)':'Q'}  # for non rotating nad
    star_info = {}
    
    blocks = []
    myheader = None
    lines = 0
    
    with open(nadfile,'r') as ff:
    
        while 1:
            lines += 1
            line = ff.readline()
            if not line: break            # end of file, quit
            if line.isspace(): continue   # blank line, continue
            
            line = line.strip()
            if line[0]=='#' and lines<=3:  # general header containing star info
                contents = [i for i in line[1:].split()]
                entry_name = None
                entry_value = None
                for i in range(len(contents)):
                    if not entry_name:
                        entry_name = contents[i]
                        continue
                    if entry_name and contents[i]=='=':
                        entry_value = 0.
                        continue
                    if entry_value is not None:
                        entry_value = float(contents[i])
                        star_info[entry_name.lower().replace(' ','')] = entry_value
                        entry_name = None
                        entry_value = None
                    else:
                        entry_name += contents[i]
            elif only_header and lines>3: # we were only interested in the header
                break
            elif line[0]=='#':
                myheader = [head.strip() for head in line[1:].split('  ') if head]
            else:
                line = line.split()
                thismode = {}
                for i,key in enumerate(myheader):
                    thismode[dict_keys[key]] = float(line[i])
                #-- for nonrotating nad files, m is not available
                if not 'm' in myheader:
                    thismode['m'] = 0
                #-- for rotating AND nonrotating nad files, omega_rot is not available
                thismode['omega_rot'] = 0
                blocks.append(thismode)
    
    logger.debug("Read %d modes from MAD (nad) %s"%(len(blocks),nadfile))
    return star_info,blocks




#}