* Clone the git repository to create a local copy:
    
    git clone /home/pieterd/python/ivs/ivs .

* Don't forget to add the repository to your pythonpath. If it is located in the
directory DIRECTORY, then, e.g. in your bash profile, add:
    
    export PYTHONPATH="$PYTHONPATH:DIRECTORY"


* The IvS Python repository contains mostly python routines. Some of the time-critical
functions, however, are written in fortran. To compile them you can run

    $ python config.py compile

If you want to specify your own fortran compiler you can do so with

    $ python config.py compile f77

Note: sometimes the compilation process fails. If so, try to compile spectra/pyrotin4.f
manually, and then retry the automatic compilation.

* In the config file you may also change the paths where the data catalogs (variable: data_dir) 
can be found, if you are not using the default locations. If you clone the repository
on your laptop and you are away from the institute, you can copy the IVS data directories
to your local disk and just add that directory to "data_dir". You don't need to remove
the other directories, if they don't exist, they will simply be skipped.


* Make sure that your python path points to the ivs folder, so that you can simply import
using, for example:

    >>> from ivs.statistics import linearregression

  If your python repository is in, e.g., ~/python/ivs/ivs , you can put in your .bash_profile:
  
    PYTHONPATH="/Users/Joris/Development/Python/ivs/"
    export PYTHONPATH


* To generate the documentation, copy the file

    /home/pieterd/python/ivs/makedoc.py
    
to a directory of your choice (referenced as /home/user/mydir/ henceforth),
but preferable *not* in the ivs root directory (e.g. one directory above it).
Create a directory "doc" in that directory:
    
    $ mkdir /home/user/mydir/doc

And copy the directory "/home/pieterd/python/ivs/doc/images" to that directory:
    
    $ cp -r /home/pieterd/python/ivs/doc/images /home/user/mydir/doc/
    
Now run the makedoc.py script in the right directory:
    
    $ cd /home/user/mydir
    $ python makedoc.py
    
Open "/home/user/mydir/doc/index.html" in your favorite browser and start browsing!
Whenever you change something yourself in your local branch or you pull changes
from someone elses, you can re-run the makedoc.py script.






Encountered errors and their solutions:
=======================================

1. Q: When I run "python config.py compile", I get the following error: 
numpy.distutils.fcompiler.CompilerNotFound: gnu95: f90 nor f77
A: Install gfortran.

