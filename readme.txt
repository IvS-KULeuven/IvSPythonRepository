
* The IvS Python repository contains mostly python routines. Some of the time-critical
functions, however, are written in fortran. To compile them you can run

    $ python config.py compile

If you want to specify your own fortran compiler you can do so with

    $ python config.py compile f77


* In the config file you may also change the paths where the data catalogs (variable: data_dir) 
can be found, if you are not using the default locations.


* Make sure that your python path points to the ivs folder, so that you can simply import
using, for example:

    >>> from ivs.statistics import linearregression

  If your python repository is in, e.g., ~/python/ivs/ivs , you can put in your .bash_profile:
  
    PYTHONPATH="/Users/Joris/Development/Python/ivs/"
    export PYTHONPATH
