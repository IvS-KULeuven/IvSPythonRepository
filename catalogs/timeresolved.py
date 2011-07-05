# -*- coding: utf-8 -*-
"""
Retrieve epoch photometry from the internet.

Author: Joris De Ridder & Pieter Degroote

Error messages are written to the logger "timeresolved".
"""

from __future__ import with_statement
import httplib
import logging
import os
import pyfits
import numpy as np
from ivs.misc import loggers
from ivs.catalogs import sesame
from ivs.io import fits
from ivs import config
        
logger = logging.getLogger("timeresolved")
logger.addHandler(loggers.NullHandler)







if __name__=="__main__":
    import doctest
    doctest.testmod()
    