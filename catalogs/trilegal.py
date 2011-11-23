# -*- coding: utf-8 -*-
"""
Web query of TRILEGAL populsation synthesis code.
Author: Joris De Ridder
"""


from time import sleep
from urllib import urlretrieve
from mechanize import Browser, urlopen

def trilegal(outputFileName,
             useGalacticCoordinates = True,
             eqRightAscension = 0, eqDeclination = 0, 
             galLongitude = 0, galLatitude = 90, 
             fieldArea = 1,
             passband = 4, magnitudeLimit = 26, magnitudeResolution = 0.1,
             IMFtype = 3,
             includeBinaries = True, binaryFraction = 0.3, lowerBinaryMassRatio = 0.7, upperBinaryMassRatio = 1.0):
    
    """
    Query the web interface of the TRILEGAL population synthesis code.
    
    The TRILEGAL webform is automatically filled and submitted. The computations are done locally
    on Girardi's computer. As soon as they are finished, the script retrieves the data file.
    
    Example:
    
    >>> trilegal("output.txt", useGalacticCoordinates=True, galLongitude=3, galLatitude=14, fieldArea=1, magnitudeLimit=7)
    
    @param outputFileName: name of file wherein trilegal output will be saved
    @type outputFileName: string
    @param galCoordinates: if True: use galactic coordinates and ignore equatorial coordinates, 
                           if False: use equatorial coordinates and ignore galactic coordinates
    @type galCoordinates: boolean                       
    @param eqRightAscension: equatorial right ascension (alpha) in degrees
    @type eqRightAscension: integer
    @param eqDeclination: equatorial declination (delta) in degrees
    @type eqDeclination: integer
    @param galLongitude: galactic longitude (l) in degrees
    @type galLongitude: integer
    @param galLatitude: galactic latitude (b) in degrees
    @type galLatitude: integer
    @param fieldArea: total field area in square degrees (max. 10 deg^2)
    @type fieldArea: integer
    @param passband: U,B,V,R,I,J,H,K = 1,2,3,4,5,6,7,8 for magnitude limit
    @type passband: integer
    @param magnitudeLimit: magnitude limit in specified passband
    @type magnitudeLimit: float
    @param magnitudeResolution: Distance modulus resolution of Galaxy components (mag)
    @type magnitudeResolution: float
    @param IMFtype: type of Initial Mass Function of single stars (1 = Salpeter with cutoff at 0.01, Msun, 2 = Chabrier exponential, 3 = Chabrier lognormal, 4 = Kroupa corrected for binaries, 5 = Kroupa not corrected for binaries)
    @type IMFtype: integer
    @param includeBinaries: include binaries in the population (True or False)
    @type includeBinaries: boolean
    @param binaryFraction: fraction of binaries 
    @type binaryFraction: float
    @param lowerBinaryMassRatio: lower limit of binary mass fraction
    @type lowerBinaryMassRatio: float
    @param upperBinaryMassRatio: upper limit of binary mass fraction
    @type upperBinaryMassRatio: float
    @return None. A file is retrieved
    
    """
    
    # The latest Trilegal web version
    
    trilegalURL = "http://stev.oapd.inaf.it/cgi-bin/trilegal_1.4"
    
    # Get the web form
    
    print("Opening TRILEGAL web interface")
    
    myBrowser = Browser()
    myBrowser.open(trilegalURL)
    myBrowser.select_form(nr=0)    # there is only one form...
    
    # Fill in the form. To know how the different fields in the form are
    # named, we used
    # >>> request = mechanize.Request(trilegalURL)
    # >>> response = mechanize.urlopen(request)
    # >>> forms = mechanize.ParseResponse(response, backwards_compat=False)
    # >>> print forms[0]
    
    print("Filling TRILEGAL web form")
    
    myBrowser["gal_coord"]    = [str(int(useGalacticCoordinates)+1)]   # 1 or 2
    myBrowser["eq_alpha"]     = str(eqRightAscension)
    myBrowser["eq_delta"]     = str(eqDeclination)
    myBrowser["gc_l"]         = str(galLongitude)
    myBrowser["gc_b"]         = str(galLatitude)
    myBrowser["field"]        = str(fieldArea)
    myBrowser["icm_lim"]      = str(passband) 
    myBrowser["mag_lim"]      = str(magnitudeLimit)
    myBrowser["mag_res"]      = str(magnitudeResolution)
    myBrowser["binary_kind"]  = [str(int(includeBinaries))]
    myBrowser["binary_frac"]  = str(binaryFraction)
    myBrowser["binary_mrinf"] = str(lowerBinaryMassRatio)
    myBrowser["binary_mrsup"] = str(upperBinaryMassRatio)
     
    # Submit the completed form
    
    print("Submitting completed TRILEGAL web form")
    
    nextWebPage = myBrowser.submit()
    
    # Trilegal is now computing the result. Click on the special "Refresh" 
    # button until the webpage says that the computations are finished.
    
    print ("Waiting until TRILEGAL computations are finished")
    
    myBrowser.select_form(nr=0)                   # one form on the "be patient" web page 
    message = "Your job was finished"
    while (message not in nextWebPage.read()):
        nextWebPage = urlopen(myBrowser.click())  # click on the Refresh button
        myBrowser.select_form(nr=0)               # select form again, so that we can make a click again
        sleep(5)                                  # to not overload the website with refresh requests
        
    # Get the url of the outputfile, and retrieve it. This can take a while.
    
    print("Retrieving TRILEGAL output file")
    
    outputLink = myBrowser.links(url_regex="lgirardi/tmp/output").next()
    urlretrieve(outputLink.absolute_url, outputFileName)
    myBrowser.close()
    
    # Save the parameters in an info file
    
    parameterInfo = """
    useGalacticCoordinates {0}
    eqRightAscension {1}
    eqDeclination {2}
    galLongitude {3}
    galLatitude {4}
    fieldArea {5}
    passband {6}
    magnitudeLimit {7}
    magnitudeResolution {8}
    IMFtype {9}
    includeBinaries {10}
    binaryFraction {11}
    lowerBinaryMassRatio {12}
    upperBinaryMassRatio {13}
    """.format(useGalacticCoordinates,
               eqRightAscension,
               eqDeclination,
               galLongitude, 
               galLatitude, 
               fieldArea,
               passband, 
               magnitudeLimit, 
               magnitudeResolution,
               IMFtype,
               includeBinaries, 
               binaryFraction, 
               lowerBinaryMassRatio, 
               upperBinaryMassRatio)
                  
    infoFileName = "info_" + outputFileName
    with open(infoFileName, 'w') as infoFile:
        infoFile.write(parameterInfo)
