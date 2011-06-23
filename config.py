# -*- coding: utf-8 -*-
"""
Global configuration of the IvS package.
"""
import os
import sys
import glob as glob_module

#-- You can add directories here, but be sure that the relative paths within
#   those directories are correct!
data_dirs = [os.getenv('ivsdata'),'/STER/pieterd/IVSDATA/', '/STER/kristofs/IVSdata']
ivs_dirs = dict(coralie='/STER/coralie/',
                hermes='/STER/mercator/hermes/')
         


def get_datafile(relative_path,basename):
    """
    Reconstruct the location of a data file and check whether it exists.
    
    If the file exists, it will return the absolute path of the filename.
    
    If the file does not exist, it will raise an IOError.
    
    @param relative_path: relative path starting from main data directory tree
    @type relative_path: str
    @param basename: filename
    @type basename: str
    @return: absolute path to the file
    @rtype: str
    """
    for data_dir in data_dirs:
        if data_dir is None: continue
        
        filename = os.path.join(data_dir,relative_path,basename)
        
        if os.path.isfile(filename):
            break
    else:
        str_data_dirs = ", ".join([idir for idir in data_dirs if idir is not None])
        relative_file = os.path.join(relative_path,basename)
        raise IOError, "File %s not found in any of the specified data directories %s"%(relative_file,str_data_dirs)
    
    return filename

def glob(relative_path,arg='*'):
    """
    Glob the files in a relative path.
    
    @param relative_path: relative path starting from main data directory tree
    @type relative_path: str
    @param arg: argument to use in glob
    @type arg: str
    @return: sorted list of files
    @rtype: list of str
    """
    files = []
    for data_dir in data_dirs:
        if data_dir is None: continue
        files += glob_module.glob(os.path.join(data_dir,relative_path,arg))
    files.sort()
    return files
    
    

if __name__=="__main__":    
    to_install = ['spectra/pyrotin4']
    main_dir = os.getcwd()
    if len(sys.argv)>1:
        if sys.argv[1]=='compile':
            for name in to_install:
                direc,pname = os.path.dirname(name),os.path.basename(name)
                os.path.chdir(direc)
                if not os.path.isfile(name+'.so'):
                    os.system('f2py -c %s.f -m %s'%(pname,pname))
                os.path.chdir(main_dir)
       
