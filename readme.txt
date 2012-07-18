
* The IvS Python repository contains mostly python routines. Some of the time-critical
functions, however, are written in fortran. To compile them you can run

$ python config.py compile

* In the config file you may also change the paths where the data catalogs (variable: data_dir) 
can be found, if you are not using the default locations.

* Make sure that your python path points to the ivs folder, so that you can simply import
using, for example:

>>> from ivs.statistics import linearregression

