# -*- coding: utf-8 -*-
"""
Global configuration of the IvS package.

Usage: $ python config.py compile
"""

import os
import sys
import glob as glob_module

#-- You can add directories here, but be sure that the relative paths within
#   those directories are correct!
data_dirs = [os.getenv('ivsdata'),'/STER/pieterd/IVSDATA/', '/STER/kristofs/IVSdata','/STER/jorisv/IVSDATA/',
             '/STER/kenneth/Python_repository/','/home/ben/public_html/opacities','/STER/michelh/IVSDATA/']

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
        raise IOError("File %s not found in any of the specified data directories %s"%(relative_file,str_data_dirs))

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
    import subprocess
    from ivs.aux import loggers
    import shutil
    import time
    import glob
    logger = loggers.get_basic_logger()

    to_install = ['spectra/pyrotin4',
                  'timeseries/deeming','timeseries/eebls','timeseries/multih',
                  'timeseries/pyclean','timeseries/pyKEP','timeseries/pyscargle',
                  'timeseries/pyGLS','timeseries/pyfasper_single','timeseries/pyfasper',
                  'timeseries/pyscargle_single','timeseries/pydft','sigproc/pyfinterpol']
    main_dir = os.getcwd()
    if len(sys.argv)>1:
        if len(sys.argv)>=3:
            compiler = sys.argv[2]
        else:
            compiler = 'gfortran'
        if sys.argv[1]=='compile':
            answer = 'y'
            for name in to_install:
                #-- break up fortran filepath in file and directory name
                direc,pname = os.path.dirname(name),os.path.basename(name)
                if not os.path.isfile(name+'.so'):
                    #-- build the command to compile the program and log it to
                    #   the user
                    cmd = 'f2py --fcompiler=%s -c %s.f -m %s'%(compiler,os.path.join(direc,pname),pname)
                    logger.info('Compiling %s: %s'%(pname.upper(),cmd))
                    if answer!='Y':
                        answer = input('Continue? [Y/n] ')
                        if answer.lower()=='n':
                            continue
                    #-- call the compiling command
                    p = subprocess.check_output(cmd,shell=True)#,stdout=devnull)
                    #-- check if compilation went fine
                    # find compiled file name
                    compiled_filename = list(filter(os.path.isfile, glob.glob('./'+pname+'*.so')))
                    # if it exists move it to to appropriate dir
                    # and change the compiled name e.g. deeming.cpython-36m-x86_64-linux-gnu.so
                    # to e.g. deeming.so
                    if compiled_filename:
                        shutil.move(compiled_filename[0],name+'.so')
                        logger.info('... succeeded')
                    else:
                        logger.error('FAILED')
                else:
                    logger.info('%s already compiled'%(pname.upper()))
