"""
Retrieve Hipparcos epoch photometry from the internet

Author: Joris De Ridder

Error messages are written to the logger "timeresolved".
"""

from __future__ import with_statement
import numpy as np
import httplib
import logging


# Setup the logger.
# Add at least one handler to avoid the message "No handlers could be found" 
# on the console. The NullHandler is part of the standard logging module only 
# from Python 2.7 on.

class NullHandler(logging.Handler):
    def emit(self, record):
        pass
        
logger = logging.getLogger("timeresolved")
nullHandler = NullHandler()
logger.addHandler(nullHandler)



def getHipData(hipnr, outputFileName=None):

    """
    Retrieve Hipparcos epoch photometry from the ESA website.
    Output are ndarrays with the time series. If a filename is
    specified, the time series together with some header information 
    is written to the file.

    Example:
    
    >>> time, magnitude, error, flag = getHipData(23124)
    >>> time = time[flag <= 2]
    >>> magnitude = magnitude[flag <= 2]
    >>> errorbar = errorbar[flag <= 2]
    >>> # Do the same but also write it to a file.
    >>> time, magnitude, error, flag = getHipData(23124, "myfile.txt")
    
    
    @param hipnr: the hipparcos number of the star. 
                  E.g. 1234 or "1234"
    @type hipnr: integer or string
    @param outputFileName: optional. If given, the retrieved Hipparcos time series 
                           will be written to the specified file.                           
    @type outputFileName: string
    @return: 4 ndarrays: time (days), magnitudes, errorbars, flags
    @rtype: tuple
    
    """

    server = "www.rssd.esa.int"
    webpage = "/hipparcos_scripts/HIPcatalogueSearch.pl?hipepId="

    # Connect to the website, en retrieve the wanted webpage
    
    conn = httplib.HTTPConnection(server)
    conn.request("GET", webpage + str(hipnr))
    response = conn.getresponse()
    if response.reason != "OK":
        logger.error("Data retrieval for HIP%s not possible. Reason: %s" % (str(hipnr), response.reason))
        return array([])
    else:
        logger.info("Data retrieval for HIP%s: OK" % str(hipnr))
        
    contents = response.read()
    conn.close()
    
    # Parse the webpage, to remove the html codes (line starts with <").
    # Put a "#" in front of the header information, and format nicely.
    # Write to the output file.
    
    if outputFileName is not None:
        outputFile = open(outputFileName, 'w')

    time = []
    magnitude = []
    errorbar = []
    flag = []
    
    for line in contents.split('\n'):
        if line == "": continue
        if not line.startswith("<"):
            if line[0].isdigit():
                items = line.split('|')
                time += [float(items[0])]
                magnitude += [float(items[1])]
                errorbar += [float(items[2])]
                flag += [float(items[3])]
                if outputFileName is not None:
                    line = line.replace("|", " ").replace("\r", "")
            else:
                if outputFileName is not None:
                    line = "# " + line
                    
            if outputFileName is not None:
                outputFile.write(line + "\n")
                
    if outputFileName is not None:
        outputFile.close()
    
    return np.array(time, dtype=np.double),        \
           np.array(magnitude, dtype=np.double),   \
           np.array(errorbar, dtype=np.double),    \
           np.array(flag, dtype=np.int)
    