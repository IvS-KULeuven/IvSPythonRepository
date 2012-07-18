
The IvS Python repository contains mostly python routines. Some of the time-critical
functions, however, are written in fortran. To compile them you can run

$ python config.py compile

Make sure that your python path points to the ivs folder, so that you can simply import
using, for example:

>>> from ivs.statistics import linearregression

