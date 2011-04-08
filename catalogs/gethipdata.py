# -*- coding: utf-8 -*-
"""
Retrieve Hipparcos epoch photometry from the internet

Author: Joris De Ridder

Error messages are written to the logger "gethipdata".
"""

from __future__ import with_statement
import httplib
import logging
import numpy as np


# Setup the logger.
# Add at least one handler to avoid the message "No handlers could be found" 
# on the console. The NullHandler is part of the standard logging module only 
# from Python 2.7 on.

class NullHandler(logging.Handler):
    def emit(self, record):
        pass
        
logger = logging.getLogger("gethipdata")
nullHandler = NullHandler()
logger.addHandler(nullHandler)



def getHipData(hipnr, outputFileName=None):

    """
    Retrieve Hipparcos epoch photometry from the ESA website.
    The time series together with some header information is stored in a record
    array and dictionary, and optionally written in a specified file.
    
    The time points are given in barycentric Julian Date and are corrected
    for the offset of 2440000.0 in the original files, but B{only in the output
    record array}. The output files display the B{original contents}.
    
    The header dictionary is of style::
    
        {'HH14': ('A', 'Annex flag (light curves)'),
         ...}
    
    For more information:
    C{http://www.rssd.esa.int/SA-general/Projects/Hipparcos/CATALOGUE_VOL1/sect2_05.ps.gz}
    
    Example 1: using the output file
    
    >>> data,header = getHipData(23124, "myfile.txt")
    >>> from numpy import loadtxt
    >>> data = loadtxt("myfile.txt")
    >>> time,magnitude,errorbar,flag = data[:,0],data[:,1],data[:,2],data[:,3]
    >>> time = time[flag <= 2.]
    >>> magnitude = magnitude[flag <= 2]
    >>> errorbar = errorbar[flag <= 2]
    
    Example 2: using the output record array
    
    >>> data,header = getHipData(23124)
    >>> data = data[ data['q_flag'<=2.] ]

    
    @param hipnr: the hipparcos number of the star. 
                  E.g. 1234 or "1234"
    @type hipnr: integer or string
    @param outputFileName: the name of the file that will be created
                           to save the Hipparcos time series
    @type outputFileName: string
    @return: record array with fields C{time}, C{mag}, C{e_mag}, C{q_mag} and
    a dictionary containing the header information
    @rtype: rec array, dict
    
    """

    server = "www.rssd.esa.int"
    webpage = "/hipparcos_scripts/HIPcatalogueSearch.pl?hipepId="

    # Connect to the website, en retrieve the wanted webpage
    
    conn = httplib.HTTPConnection(server)
    conn.request("GET", webpage + str(hipnr))
    response = conn.getresponse()
    if response.reason != "OK":
        logger.error("Data retrieval for HIP%s not possible. Reason: %s" % (str(hipnr), response.reason))
        return
    else:
        logger.info("Data retrieval for HIP%s: OK" % str(hipnr))
        
    contents = response.read()
    conn.close()
    
    # Parse the webpage, to remove the html codes (line starts with <").
    # Put a "#" in front of the header information, and format nicely.
    # Write to the output file if asked for.
    data = []
    header = {}
    
    if outputFileName:
        outputFile = open(outputFileName,'w')
    
    for line in contents.split('\n'):
        if line == "": continue
        if not line.startswith("<"):
            line = line.replace("|", " ").replace("\r", "")
            #-- this is the header
            if not line[0].isdigit():
                sline = line.split(':')
                #-- only keep header entries of the style "key: value information"
                #   in the dictionary
                if len(sline)==2:
                    key,info = sline
                    info = info.split()
                    header[key] = (info[0]," ".join(info[1:]))
                if outputFileName:
                    line = "# " + line
            #-- this is the real contents
            else:
                data.append(tuple(line.split()))
            if outputFileName:
                outputFile.write(line + "\n")
    if outputFileName:
        outputFile.close()
    
    #-- now make a record array
    #   we choose the header name to be in the VizieR style
    dtypes = [('time','>f4'),('mag','>f4'),('e_mag','>f4'),('q_mag','>f4')]
    data = np.rec.array(data,dtype=dtypes)
    
    #-- fix the time offset
    data['time'] += 2440000.0
    return data,header
    