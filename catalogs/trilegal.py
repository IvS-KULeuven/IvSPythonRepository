# -*- coding: utf-8 -*-
"""
Web query of TRILEGAL populsation synthesis code.
Author: Joris De Ridder
"""


from time import sleep
from urllib import urlretrieve
from mechanize import Browser, urlopen

def trilegal(outputFileName,
             longitude = 0, latitude = 0,
             coordinateType = "galactic",
             fieldArea = 1,
             passband = 4, magnitudeLimit = 26, magnitudeResolution = 0.1,
             IMFtype = 3,
             includeBinaries = True, binaryFraction = 0.3, lowerBinaryMassRatio = 0.7, upperBinaryMassRatio = 1.0,
             useThinDisc = False,
             useThickDisc = False,
             useBulge = True):
    
    """
    Query the web interface of the TRILEGAL population synthesis code.
    
    The TRILEGAL webform is automatically filled and submitted. The computations are done locally
    on Girardi's computer. As soon as they are finished, the script retrieves the data file.
    
    Example:
    
    >>> trilegal("output.txt", longitude=3, latitude=14, coordinateType="galactic", fieldArea=1, magnitudeLimit=7, useThinDisc=True)
    
    @param outputFileName: name of file wherein trilegal output will be saved
    @type outputFileName: string
    @param longitude: galactic longitude (degrees) or right ascension (hours)
    @type longitude: integer
    @param latitude:  galactic latitude (degrees) or declination (degrees)
    @type latitude: integer
    @param coordinateType: either "galactic", or "equatorial" 
    @type coordinateType: string
    @param fieldArea: total field area in square degrees (max. 10 deg^2)
    @type fieldArea: float
    @param passband: U,B,V,R,I,J,H,K = 1,2,3,4,5,6,7,8 for magnitude limit
    @type passband: integer
    @param magnitudeLimit: magnitude limit in specified passband
    @type magnitudeLimit: float
    @param magnitudeResolution: Distance modulus resolution of Galaxy components (mag)
    @type magnitudeResolution: float
    @param IMFtype: type of Initial Mass Function of single stars 
                    1 = Salpeter with cutoff at 0.01, Msun, 
                    2 = Chabrier exponential, 
                    3 = Chabrier lognormal, 
                    4 = Kroupa corrected for binaries, 
                    5 = Kroupa not corrected for binaries
    @type IMFtype: integer
    @param includeBinaries: include binaries in the population (True or False)
    @type includeBinaries: boolean
    @param binaryFraction: fraction of binaries 
    @type binaryFraction: float
    @param lowerBinaryMassRatio: lower limit of binary mass fraction
    @type lowerBinaryMassRatio: float
    @param upperBinaryMassRatio: upper limit of binary mass fraction
    @type upperBinaryMassRatio: float
    @param useThinDisk: if True use squared hyperbolic secant along z, if False don't include
    @type useThinDisk: boolean
    @param useThickDisk: if True use squared hyperbolic secant along z, if False don't include
    @type useThickDisk: boolean
    @param useBulge: if True use triaxal bulge, if False don't include
    @type useBulge: boolean
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
    
    if coordinateType == "galactic":
        myBrowser["gal_coord"] = ["1"]
        myBrowser["gc_l"]      = str(longitude)
        myBrowser["gc_b"]      = str(latitude)
    else:
        myBrowser["gal_coord"] = ["2"]
        myBrowser["eq_alpha"]  = str(longitude)
        myBrowser["eq_delta"]  = str(latitude)
    
    myBrowser["field"]        = str(fieldArea)
    myBrowser["icm_lim"]      = str(passband) 
    myBrowser["mag_lim"]      = str(magnitudeLimit)
    myBrowser["mag_res"]      = str(magnitudeResolution)
    myBrowser["binary_kind"]  = [str(int(includeBinaries))]
    myBrowser["binary_frac"]  = str(binaryFraction)
    myBrowser["binary_mrinf"] = str(lowerBinaryMassRatio)
    myBrowser["binary_mrsup"] = str(upperBinaryMassRatio)

    if useThinDisc:
        myBrowser["thindisk_kind"] = ["3"]
    else:
        myBrowser["thindisk_kind"] = ["0"]
        
    if useThickDisc:
        myBrowser["thickdisk_kind"] = ["3"]
    else:
        myBrowser["thickdisk_kind"] = ["0"]
        
    if useBulge:
        myBrowser["bulge_kind"] = ["2"]
    else:
        myBrowser["bulge_kind"] = ["0"]
     
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
    coordinateType {0}
    longitude {1}
    latitude {2}
    fieldArea {3}
    passband {4}
    magnitudeLimit {5}
    magnitudeResolution {6}
    IMFtype {7}
    includeBinaries {8}
    binaryFraction {9}
    lowerBinaryMassRatio {10}
    upperBinaryMassRatio {11}
    """.format(coordinateType,
               longitude,
               latitude,
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
