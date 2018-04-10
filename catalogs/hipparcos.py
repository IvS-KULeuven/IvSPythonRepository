# -*- coding: utf-8 -*-
"""
Retrieve Hipparcos epoch/intermediate photometry from the internet.

Author: Joris De Ridder & Pieter Degroote
"""


import http.client
import logging
import numpy as np
from ivs.aux import loggers
from ivs.catalogs import sesame

logger = logging.getLogger("CAT.HIP")
logger.addHandler(loggers.NullHandler)

def getHipData(ID,dtype='ep',outputFileName=None):

    """
    Retrieve Hipparcos epoch/intermediate photometry from the ESA website.

    The time series together with some header information is stored in a record
    array and dictionary, and optionally written in a specified file.

    The time points are given in barycentric Julian Date and are corrected
    for the offset of 2440000.0 in the original files, but B{only in the output
    record array}. The output files display the B{original contents}.

    For epoch photometry, set C{dtype='ep'}.
    For intermediate date, set C{dtype='i'}.

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

    In the case of intermediate data products:
        - orbit: orbit number
        - source: source of abscissa (F=FAST, f=rejected FAST, N=NDAC,n=NDAC rejected)
        - d_acosd: partial derivative wrt alpha cos(delta)
        - d_d: partial derivative wrt delta
        - d_pi: partial derivative wrt parallax
        - d_mua: partial derivative wrt proper motion alpha cos(delta)
        - d_mud: partial derivative wrt proper motion delta

    @param ID: identification of the star: if you give an integer or string that
    can be converted to an integer, it is assumed to be the hipparcos number of
    the star.  E.g. 1234 or "1234". If it is not an integer, the star will
    be resolved via sesame to get the HIP number if possible
    @type ID: integer or string
    @param dtype: data type (epoch ('ep') photometry or intermediate ('i') data)
    @type dtype: string (one of ('ep','i'))
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
    webpage = "/hipparcos_scripts/HIPcatalogueSearch.pl?hip%sId="%(dtype)

    # Resolve the name if necessary (i.e., if it's not a HIP number). If the
    # star has no HIP number, log an error and return None
    try:
        hipnr = int(ID)
    except ValueError:
        info = sesame.search(ID,db='S')
        IDs = [alias for alias in info['alias'] if 'HIP' in alias]
        if len(IDs)!=1:
            logger.error("Data retrieval for %s not possible. Reason: no HIP number resolved" % (ID))
            return
        hipnr = IDs[0].split(' ')[1]

    # Connect to the website, en retrieve the wanted webpage

    conn = http.client.HTTPConnection(server)
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
            line = line.replace("\r", "")

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
                data.append(line.split('|'))
                #-- correct for empty fields
                data[-1] = tuple([(entry.replace(' ','')=='' and np.nan or entry) for entry in data[-1]])
            if outputFileName:
                outputFile.write(line + "\n")
    if outputFileName:
        outputFile.close()

    # Make a record array.
    # Choose the header names to be in the VizieR style.

    if dtype=='ep':
        dtypes = [('time','f8'),('mag','f8'),('e_mag','f8'),('q_mag','i')]
    elif dtype=='i':
        dtypes = [('orbit','i'),('source','a1'),
                  ('d_acosd','f8'),('d_d','f8'),('d_pi','f8'),
                  ('d_mua','f8'),('d_mud','f8'),
                  ('abs_res','f8'),('abs_std','f8'),('cor','f8')]
    data = np.rec.array(data,dtype=dtypes)

    # Fix the time offset
    if dtype=='ep':
        data['time'] += 2440000.0


    return data,header





def getEpochData(ID,outputFileName=None):
    """
    Convenience function to retrieve epoch data.
    """
    return getHipData(ID,dtype='ep',outputFileName=outputFileName)





def getIntermediateData(ID,outputFileName=None):
    """
    Convenience function to retrieve intermediate data.
    """
    return getHipData(ID,dtype='i',outputFileName=outputFileName)


