# -*- coding: utf-8 -*-
"""
Retrieve Hipparcos epoch photometry from the internet

Author: Joris De Ridder & Pieter Degroote

Error messages are written to the logger "timeresolved".
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
        
logger = logging.getLogger("timeresolved")
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
       
    For more information:
    C{http://www.rssd.esa.int/SA-general/Projects/Hipparcos/CATALOGUE_VOL1/sect2_05.ps.gz}
    
    Example:
    
    >>> data,header = getHipData(1234)
    >>> data = data[data['q_mag'] <= 2]         # keep only the good points
    
    To write the retrieved data to a file:
    
    >>> data, header = getHipData(1234 , "myfile.txt")
    
   To store the different columns in separate arrays:
    
    >>> data, header = getHipData(1234)
    >>> time = data['time']
    >>> magnitude = data['mag']
    >>> errorbar = data['e_mag']
    >>> qualityflag = data['q_mag']
    
    
    @param hipnr: the hipparcos number of the star. 
                  E.g. 1234 or "1234"
    @type hipnr: integer or string
    @param outputFileName: the name of the file that will be created
                           to save the Hipparcos time series
    @type outputFileName: string
    @return: record array with fields time, mag, e_mag (errorbar), 
             q_mag (quality flag), and a dictionary containing the 
             header information. The header dictionary is of style
             {'HH14': ('A', 'Annex flag (light curves)'), ...}
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
            
            # This is the header
            
            if not line[0].isdigit():
                sline = line.split(':')
                
                # Only keep header entries of the style 
                # "key: value information" in the dictionary
                
                if len(sline)==2:
                    key,info = sline
                    info = info.split()
                    header[key] = (info[0]," ".join(info[1:]))
                if outputFileName:
                    line = "# " + line
                    
            # This is the real contents
            
            else:
                data.append(tuple(line.split()))
            if outputFileName:
                outputFile.write(line + "\n")
    if outputFileName:
        outputFile.close()
    
    # Make a record array.
    # Choose the header names to be in the VizieR style.
    
    dtypes = [('time','>f4'),('mag','>f4'),('e_mag','>f4'),('q_mag','i')]
    data = np.rec.array(data,dtype=dtypes)
    
    # Fix the time offset
    
    data['time'] += 2440000.0
    
    return data,header




if __name__=="__main__":
    import doctest
    doctest.testmod()
    