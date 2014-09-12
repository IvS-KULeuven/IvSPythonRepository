* Clone the git repository to create a local copy. First 'cd' to the directory where you keep all your 
  python scripts, e.g. ~/python:
  
    $ cd python/

  If you're using a Desktop computer of the IvS, the easiest way to clone is:
  
    git clone /home/joris/python/ivs/ivs/proj.git ivs
  
  which will clone the entiry repository into the ~/python/ivs folder. It's important to choose 'ivs'
  as folder name, because all importing in python is done like
  
    from ivs.statistics import linearregression
    
  That is, it will assume an 'ivs' folder. 
  If you're using your own laptop, you can clone using:

    $ git clone YOURNAME@copernicus.ster.kuleuven.be:/home/joris/python/ivs/ivs/proj.git ivs
    
  Again, this will clone all repository files in the ~/python/ivs folder. Be aware, however, that only
  the python scripts and the documentation are being cloned, not the (numerous and sometimes huge) datafiles
  that come along with it, containing, for example, limbdarkening coefficients.
    
* Updating your own clone of the IvS python repository to the most recent version can be done with

    $ cd ivs
    $ git pull

* The IvS Python repository contains mostly python routines. Some of the time-critical
functions, however, are written in fortran. To compile them you can run

    $ python config.py compile

If you want to specify your own fortran compiler you can do so with

    $ python config.py compile f77

Note: sometimes the compilation process fails. If so, try to compile spectra/pyrotin4.f
manually, and then retry the automatic compilation:
    
    $ cd spectra/
    $ f2py --fcompiler=gfortran -c pyrotin4.f -m pyrotin4
    $ cd ../
    $ python config.py compile

* In the config file you may also change the paths where the data catalogs (variable: data_dir) 
can be found, if you are not using the default locations. If you clone the repository
on your laptop and you are away from the institute, you can copy the IVS data directories
to your local disk and just add that directory to "data_dir". You don't need to remove
the other directories, if they don't exist, they will simply be skipped.


* Make sure that your python path points to the ivs folder, so that you can simply import
using, for example:

    >>> from ivs.statistics import linearregression

  If your python repository is in, e.g., ~/python/ivs/ , you can put in your .bash_profile:
  
    export PYTHONPATH=/home/YOURNAME/python:$PYTHONPATH

Warning: don't put ~/python/ivs in your Python path, but just ~/python.


* To generate the documentation, simply run the script

    $ python makedoc.py

  in the repository's root folder. This assumes that 'epydoc' is available which is 
  already installed on all IvS computers. On your own laptop, you can get it from 
  http://epydoc.sourceforge.net.
 
Open "/doc/html/index.html" in your favorite browser and start browsing!
Whenever you change something yourself in your local branch or you pull changes
from someone else, you can re-run the makedoc.py script.


* Questions about the IvS Python repository can be asked through the small mailing list
    python@ster.kuleuven.be.
  
* Happy computing!





Encountered errors and their solutions:
=======================================

1. Q: When I run "python config.py compile", I get the following error: 
numpy.distutils.fcompiler.CompilerNotFound: gnu95: f90 nor f77
A: Install gfortran.

