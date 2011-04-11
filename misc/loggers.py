# -*- coding: utf-8 -*-
"""
Set up and configure basic logging facilities.

Typcical use: import the module, and call

>>> logger = get_basic_logger()

in your script. All relevant messages will get logged to the terminal.
"""
import logging

class NullHandler(logging.Handler):
    """
    The NullHandler is part of the standard logging module only from Python 2.7 on.
    """
    level = 100
    def emit(self, record):
        pass

def get_basic_logger(name="",clevel='INFO',
                             flevel='DEBUG',filename=None,filemode='w'):
    """
    Return a basic logger via a log file and/or terminal.
    
    Example 1: log only to the console, accepting levels "INFO" and above
    
    >>> logger = loggers.get_basic_logger()
    
    Example 2: log only to the console, accepting levels "DEBUG" and above
    
    >>> logger = loggers.get_basic_logger(clevel='DEBUG')
    
    Example 3: log only to a file, accepting levels "DEBUG" and above
    
    >>> logger = loggers.get_basic_logger(clevel=None,filename='mylog.log')
    
    Example 4: log only to a file, accepting levels "INFO" and above
    
    >>> logger = loggers.get_basic_logger(clevel=None,flevel='INFO',filename='mylog.log')
    
    Example 5: log to the terminal (INFO and above) and file (DEBUG and above)
    
    >>> logger = loggers.get_basic_logger(filename='mylog.log')
    
    @param name: name of the logger
    @type name: str
    """
    #-- define formats
    format  = '%(asctime)s %(name)s %(levelname)-8s %(message)s'
    datefmt = '%a, %d %b %Y %H:%M'
    
    if clevel: clevel = logging.__dict__[clevel.upper()]
    if flevel: flevel = logging.__dict__[flevel.upper()]
    
    #-- set up basic configuration.
    #   The basicConfig sets up one default logger. If you give a filename, it's
    #   a FileHandler, otherwise a StreamHandler.
    #-- If we want console and filename, first set up a basic FileHandler, then
    #   add terminal StreamHandler
    if filename is not None:
        logging.basicConfig(level=min(flevel,clevel),
                            format=format,datefmt=datefmt,
                            filename=filename,filemode=filemode)
    if filename is not None and clevel:    
        # define a Handler which writes INFO messages or higher to the sys.stderr
        ch = logging.StreamHandler()
        ch.setLevel(clevel)
        # set a format which is simpler for console use
        formatter = logging.Formatter(format,datefmt)
        # tell the handler to use this format
        ch.setFormatter(formatter)
        logging.getLogger(name).addHandler(ch)
    #-- If we only want a console:
    else:
        logging.basicConfig(level=clevel,format=format,datefmt=datefmt,filename=filename,filemode=filemode)        
    return logging.getLogger(name)
