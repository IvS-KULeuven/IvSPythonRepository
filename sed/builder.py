# -*- coding: utf-8 -*-
"""
SED builder program.

To construct an SED, use the SED class. The functions defined in this module
are mainly convenience functions specifically for that class, but can be used
outside of the SED class if you know what you're doing.

Section 1. Retrieving and plotting photometry of a target
=========================================================

>>> mysed = SED('HD180642')
>>> mysed.get_photometry()
>>> mysed.plot_data()

and call Pylab's C{show} function to show to the screen:

]]include figure]]ivs_sed_builder_example_photometry.png]

You can give a B{search radius} to C{get_photometry} via the keyword C{radius}.
The default value is 10 arcseconds for stars dimmer than 6th magnitude, and 60
arcseconds for brighter stars. The best value of course depends on the density
of the field.

>>> mysed.get_photometry(radius=5.)

If your star's name is not recognised by any catalog, you can give coordinates
to look for photometry. In that case, the ID of the star given in the C{SED}
command will not be used to search photometry (only to save the phot file):

>>> mysed.get_photometry(ra=289.31167983,dec=1.05941685)

Note that C{ra} and C{dec} are given in B{degrees}.

You best B{switch on the logger} (see L{ivs.aux.loggers.get_basic_logger}) to see the progress:
sometimes, access to catalogs can take a long time (the GATOR sources are
typically slow). If one of the C{gator}, C{vizier} or C{gcpd} is impossibly slow,
you can B{include/exclude these sources} via the keywords C{include} or C{exclude},
which take a list of strings (choose from C{gator}, C{vizier} and/or C{gcpd}).

>>> mysed.get_photometry(exclude=['gator'])

The results will be written to the file B{HD180642.phot}. An example content is::

  #  meas    e_meas flag unit photband          source                        _r   _RAJ2000    _DEJ2000   cwave       cmeas     e_cmeas cunit       color include
  #float64   float64 |S20 |S30 |S30              |S50                     float64    float64     float64 float64     float64     float64 |S50         bool    bool
    7.823      0.02 nan  mag  WISE.W3           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05  123337 4.83862e-17 8.91306e-19 erg/s/cm2/A     0       1
    7.744     0.179 nan  mag  WISE.W4           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05  222532 4.06562e-18 6.70278e-19 erg/s/cm2/A     0       1
     7.77     0.023 nan  mag  WISE.W1           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05 33791.9   6.378e-15  1.3511e-16 erg/s/cm2/A     0       1
    7.803      0.02 nan  mag  WISE.W2           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05   46293 1.82691e-15 3.36529e-17 erg/s/cm2/A     0       1
    8.505     0.016 nan  mag  TYCHO2.BT         I/259/tyc2                 0.042   7.17e-06    1.15e-06  4204.4 2.76882e-12 4.08029e-14 erg/s/cm2/A     0       1
    8.296     0.013 nan  mag  TYCHO2.VT         I/259/tyc2                 0.042   7.17e-06    1.15e-06 5321.86 1.93604e-12 2.31811e-14 erg/s/cm2/A     0       1
     8.27       nan nan  mag  JOHNSON.V         II/168/ubvmeans             0.01    1.7e-07    3.15e-06 5504.67 1.80578e-12 1.80578e-13 erg/s/cm2/A     0       1
     0.22       nan nan  mag  JOHNSON.B-V       II/168/ubvmeans             0.01    1.7e-07    3.15e-06     nan     1.40749    0.140749 flux_ratio      1       0
    -0.66       nan nan  mag  JOHNSON.U-B       II/168/ubvmeans             0.01    1.7e-07    3.15e-06     nan     1.21491    0.121491 flux_ratio      1       0
     8.49       nan nan  mag  JOHNSON.B         II/168/ubvmeans             0.01    1.7e-07    3.15e-06 4448.06 2.54162e-12 2.54162e-13 erg/s/cm2/A     0       1
     7.83       nan nan  mag  JOHNSON.U         II/168/ubvmeans             0.01    1.7e-07    3.15e-06 3641.75 3.08783e-12 3.08783e-13 erg/s/cm2/A     0       1
    2.601       nan nan  mag  STROMGREN.HBN-HBW J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan     1.66181    0.166181 flux_ratio      1       0
   -0.043       nan nan  mag  STROMGREN.M1      J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan    0.961281   0.0961281 flux_ratio      1       0
    8.221       nan nan  mag  STROMGREN.Y       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685 5477.32 1.88222e-12 1.88222e-13 erg/s/cm2/A     0       1
    0.009       nan nan  mag  STROMGREN.C1      J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan     0.93125    0.093125 flux_ratio      1       0
    0.238       nan nan  mag  STROMGREN.B-Y     J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan     1.28058    0.128058 flux_ratio      1       0
    8.459       nan nan  mag  STROMGREN.B       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685  4671.2 2.41033e-12 2.41033e-13 erg/s/cm2/A     0       1
    8.654       nan nan  mag  STROMGREN.V       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685 4108.07 2.96712e-12 2.96712e-13 erg/s/cm2/A     0       1
    8.858       nan nan  mag  STROMGREN.U       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685 3462.92 3.40141e-12 3.40141e-13 erg/s/cm2/A     0       1
     7.82      0.01 nan  mag  JOHNSON.J         J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06 12487.8 2.36496e-13 2.36496e-15 erg/s/cm2/A     0       1
     7.79      0.01 nan  mag  JOHNSON.K         J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06 21951.2 3.24868e-14 3.24868e-16 erg/s/cm2/A     0       1
     7.83      0.04 nan  mag  JOHNSON.H         J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06 16464.4 8.64659e-14 3.18552e-15 erg/s/cm2/A     0       1
   8.3451    0.0065 nan  mag  HIPPARCOS.HP      I/239/hip_main             0.036   2.17e-06    1.15e-06 5275.11 1.91003e-12 1.91003e-14 erg/s/cm2/A     0       1
    8.525     0.011 nan  mag  TYCHO2.BT         I/239/hip_main             0.036   2.17e-06    1.15e-06  4204.4 2.71829e-12   2.754e-14 erg/s/cm2/A     0       1
    8.309     0.012 nan  mag  TYCHO2.VT         I/239/hip_main             0.036   2.17e-06    1.15e-06 5321.86   1.913e-12 2.11433e-14 erg/s/cm2/A     0       1
     8.02     0.057 nan  mag  COUSINS.I         II/271A/patch2               0.2 -4.983e-05   2.315e-05 7884.05 7.51152e-13 3.94347e-14 erg/s/cm2/A     0       1
    8.287     0.056 nan  mag  JOHNSON.V         II/271A/patch2               0.2 -4.983e-05   2.315e-05 5504.67 1.77773e-12 9.16914e-14 erg/s/cm2/A     0       1
     8.47       nan nan  mag  USNOB1.B1         I/284/out                  0.026   7.17e-06     1.5e-07 4448.06 2.57935e-12 7.73805e-13 erg/s/cm2/A     0       1
     8.19       nan nan  mag  USNOB1.R1         I/284/out                  0.026   7.17e-06     1.5e-07 6939.52 1.02601e-12 3.07803e-13 erg/s/cm2/A     0       1
    8.491     0.012 nan  mag  JOHNSON.B         I/280B/ascc                0.036   2.17e-06    2.15e-06 4448.06 2.53928e-12 2.80651e-14 erg/s/cm2/A     0       1
    8.274     0.013 nan  mag  JOHNSON.V         I/280B/ascc                0.036   2.17e-06    2.15e-06 5504.67 1.79914e-12 2.15419e-14 erg/s/cm2/A     0       1
    7.816     0.023 nan  mag  2MASS.J           II/246/out                 0.118 -2.783e-05   1.715e-05 12412.1 2.28049e-13 4.83093e-15 erg/s/cm2/A     0       1
    7.792     0.021 nan  mag  2MASS.KS          II/246/out                 0.118 -2.783e-05   1.715e-05 21909.2 3.26974e-14 6.32423e-16 erg/s/cm2/A     0       1
    7.825     0.042 nan  mag  2MASS.H           II/246/out                 0.118 -2.783e-05   1.715e-05 16497.1 8.48652e-14 3.28288e-15 erg/s/cm2/A     0       1
    8.272     0.017 nan  mag  GENEVA.V          GCPD                         nan        nan         nan  5482.6 1.88047e-12 2.94435e-14 erg/s/cm2/A     0       1
      1.8     0.004 nan  mag  GENEVA.G-B        GCPD                         nan        nan         nan     nan    0.669837  0.00669837 flux_ratio      1       0
    1.384     0.004 nan  mag  GENEVA.V1-B       GCPD                         nan        nan         nan     nan    0.749504  0.00749504 flux_ratio      1       0
     0.85     0.004 nan  mag  GENEVA.B1-B       GCPD                         nan        nan         nan     nan     1.05773   0.0105773 flux_ratio      1       0
    1.517     0.004 nan  mag  GENEVA.B2-B       GCPD                         nan        nan         nan     nan    0.946289  0.00946289 flux_ratio      1       0
    0.668     0.004 nan  mag  GENEVA.V-B        GCPD                         nan        nan         nan     nan    0.726008  0.00726008 flux_ratio      1       0
    0.599     0.004 nan  mag  GENEVA.U-B        GCPD                         nan        nan         nan     nan     1.13913   0.0113913 flux_ratio      1       0
    7.604 0.0174642 nan  mag  GENEVA.B          GCPD                         nan        nan         nan 4200.85 2.59014e-12 4.16629e-14 erg/s/cm2/A     0       1
    9.404 0.0179165 nan  mag  GENEVA.G          GCPD                         nan        nan         nan 5765.89 1.73497e-12   2.863e-14 erg/s/cm2/A     0       1
    8.988 0.0179165 nan  mag  GENEVA.V1         GCPD                         nan        nan         nan 5395.63 1.94132e-12 3.20351e-14 erg/s/cm2/A     0       1
    8.454 0.0179165 nan  mag  GENEVA.B1         GCPD                         nan        nan         nan 4003.78 2.73968e-12 4.52092e-14 erg/s/cm2/A     0       1
    9.121 0.0179165 nan  mag  GENEVA.B2         GCPD                         nan        nan         nan 4477.56 2.45102e-12  4.0446e-14 erg/s/cm2/A     0       1
    8.203 0.0179165 nan  mag  GENEVA.U          GCPD                         nan        nan         nan 3421.62 2.95052e-12 4.86885e-14 erg/s/cm2/A     0       1
     8.27       nan nan  mag  JOHNSON.V         GCPD                         nan        nan         nan 5504.67 1.80578e-12 1.80578e-13 erg/s/cm2/A     0       1
     0.22       nan nan  mag  JOHNSON.B-V       GCPD                         nan        nan         nan     nan     1.40749    0.140749 flux_ratio      1       0
    -0.66       nan nan  mag  JOHNSON.U-B       GCPD                         nan        nan         nan     nan     1.21491    0.121491 flux_ratio      1       0
     8.49       nan nan  mag  JOHNSON.B         GCPD                         nan        nan         nan 4448.06 2.54162e-12 2.54162e-13 erg/s/cm2/A     0       1
     7.83       nan nan  mag  JOHNSON.U         GCPD                         nan        nan         nan 3641.75 3.08783e-12 3.08783e-13 erg/s/cm2/A     0       1
   -0.035       nan nan  mag  STROMGREN.M1      GCPD                         nan        nan         nan     nan    0.954224   0.0954224 flux_ratio      1       0
    0.031       nan nan  mag  STROMGREN.C1      GCPD                         nan        nan         nan     nan     0.91257    0.091257 flux_ratio      1       0
    0.259       nan nan  mag  STROMGREN.B-Y     GCPD                         nan        nan         nan     nan     1.25605    0.125605 flux_ratio      1       0
   -0.009 0.0478853 nan  mag  2MASS.J-H         II/246/out                 0.118 -2.783e-05   1.715e-05     nan     2.68719    0.118516 flux_ratio      1       0
   -0.033 0.0469574 nan  mag  2MASS.KS-H        II/246/out                 0.118 -2.783e-05   1.715e-05     nan    0.385286   0.0166634 flux_ratio      1       0
    -0.01 0.0412311 nan  mag  JOHNSON.J-H       J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06     nan     2.73514    0.103867 flux_ratio      1       0
    -0.04 0.0412311 nan  mag  JOHNSON.K-H       J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06     nan    0.375718    0.014268 flux_ratio      1       0
    0.209 0.0206155 nan  mag  TYCHO2.BT-VT      I/259/tyc2                 0.042   7.17e-06    1.15e-06     nan     1.43014    0.027155 flux_ratio      1       0

Once a .phot file is written and L{get_photometry} is called again for the same
target, the script will B{not retrieve the photometry from the internet again},
but will use the contents of the file instead. The purpose is minimizing network
traffic and maximizing speed. If you want to refresh the search, simply manually
delete the .phot file.

The content of the .phot file is most easily read using the L{ivs.io.ascii.read2recarry}
function. Be careful, as it contains both absolute fluxes as flux ratios.

>>> data = ascii.read2recarray('HD180642.phot')

Section 2. SED fitting using a grid based approach
==================================================

We make an SED of HD180642 by simply B{exploring a whole grid of Kurucz models}
(constructed via L{fit.generate_grid}, iterated over via L{fit.igrid_search} and
evaluated with L{fit.stat_chi2}). The model with the best parameters is picked
out and we make a full SED with those parameters.

>>> mysed = SED('HD180642')
>>> mysed.get_photometry()

Now we have collected B{all fluxes and colors}, but we do not want to use them
all: first, B{fluxes and colors are not independent}, so you probably want to use
either only absolute fluxes or only colors (plus perhaps one absolute flux per
system to compute the angular diameter) (see L{SED.set_photometry_scheme}). Second,
B{some photometry is better suited for fitting an SED than others}; typically IR
photometry does not add much to the fitting of massive stars, or it can be
contaminated with circumstellar material. Third, B{some photometry is not so
reliable}, i.e. the measurements can have large systematic uncertainties due to
e.g. calibration effects, which are typically not included in the error bars.

Here, we chose to use colors and one absolute flux per system, but exclude IR
photometry (wavelength range above 2.5 micron), and some systems and colors which
we know are not so trustworthy:

>>> mysed.set_photometry_scheme('combo')
>>> mysed.exclude(names=['STROMGREN.HBN-HBW','USNOB1','SDSS','DENIS','COUSINS','ANS','TD1'],wrange=(2.5e4,1e10))

Start the grid based fitting process and show some plots. We use 100000 randomly
distributed points over the grid:

>>> mysed.igrid_search(points=100000)

and make the plot

>>> p = pl.figure()
>>> p = pl.subplot(131);mysed.plot_sed()
>>> p = pl.subplot(132);mysed.plot_grid(limit=None)
>>> p = pl.subplot(133);mysed.plot_grid(x='ebv',y='z',limit=None)

]]include figure]]ivs_sed_builder_example_fitting01.png]

The grid is a bit too coarse for our liking around the minimum, so we zoom in on
the results:

>>> teffrange = mysed.results['igrid_search']['CI']['teffL'],mysed.results['igrid_search']['CI']['teffU']
>>> loggrange = mysed.results['igrid_search']['CI']['loggL'],mysed.results['igrid_search']['CI']['loggU']
>>> ebvrange = mysed.results['igrid_search']['CI']['ebvL'],mysed.results['igrid_search']['CI']['ebvU']
>>> mysed.igrid_search(points=100000,teffrange=teffrange,loggrange=loggrange,ebvrange=ebvrange)

and repeat the plot

>>> p = pl.figure()
>>> p = pl.subplot(131);mysed.plot_sed()
>>> p = pl.subplot(132);mysed.plot_grid(limit=None)
>>> p = pl.subplot(133);mysed.plot_grid(x='ebv',y='z',limit=None)

]]include figure]]ivs_sed_builder_example_fitting02.png]

With the L{SED.plot_grid} function, you can easily visualize the correlation between
effective temperature and E(B-V) for these massive stars:

>>> p = pl.figure()
>>> mysed.plot_grid(x='teff',y='ebv',limit=None)

]]include figure]]ivs_sed_builder_example_fitting03.png]

Subsection 2.1 Saving SED fits
------------------------------

You can save all the data to a multi-extension FITS file via

>>> mysed.save_fits()

This FITS file then contains all B{measurements} (it includes the .phot file),
the B{resulting SED}, and also the B{whole fitted grid}: in the above case, the
extensions of the FITS file contain the following information (we print each
header)::

    EXTNAME = 'DATA    '
    XTENSION= 'BINTABLE'           / binary table extension                         
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                  270 / length of dimension 1                          
    NAXIS2  =                   67 / length of dimension 2                          
    PCOUNT  =                    0 / number of group parameters                     
    GCOUNT  =                    1 / number of groups                               
    TFIELDS =                   18 / number of table fields                         
    TTYPE1  = 'meas    '                                                            
    TFORM1  = 'D       '                                                            
    TTYPE2  = 'e_meas  '                                                            
    TFORM2  = 'D       '                                                            
    TTYPE3  = 'flag    '                                                            
    TFORM3  = '20A     '                                                            
    TTYPE4  = 'unit    '                                                            
    TFORM4  = '30A     '                                                            
    TTYPE5  = 'photband'                                                            
    TFORM5  = '30A     '                                                            
    TTYPE6  = 'source  '                                                            
    TFORM6  = '50A     '                                                            
    TTYPE7  = '_r      '                                                            
    TFORM7  = 'D       '                                                            
    TTYPE8  = '_RAJ2000'                                                            
    TFORM8  = 'D       '                                                            
    TTYPE9  = '_DEJ2000'                                                            
    TFORM9  = 'D       '                                                            
    TTYPE10 = 'cwave   '                                                            
    TFORM10 = 'D       '                                                            
    TTYPE11 = 'cmeas   '                                                            
    TFORM11 = 'D       '                                                            
    TTYPE12 = 'e_cmeas '                                                            
    TFORM12 = 'D       '                                                            
    TTYPE13 = 'cunit   '                                                            
    TFORM13 = '50A     '                                                            
    TTYPE14 = 'color   '                                                            
    TFORM14 = 'L       '                                                            
    TTYPE15 = 'include '                                                            
    TFORM15 = 'L       '                                                            
    TTYPE16 = 'synflux '                                                            
    TFORM16 = 'D       '                                                            
    TTYPE17 = 'mod_eff_wave'                                                        
    TFORM17 = 'D       '                                                            
    TTYPE18 = 'chi2    '                                                            
    TFORM18 = 'D       '                                                            

    EXTNAME = 'MODEL   '      
    XTENSION= 'BINTABLE'           / binary table extension                         
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   24 / length of dimension 1                          
    NAXIS2  =                 1221 / length of dimension 2                          
    PCOUNT  =                    0 / number of group parameters                     
    GCOUNT  =                    1 / number of groups                               
    TFIELDS =                    3 / number of table fields                                                                               
    TTYPE1  = 'wave    '                                                            
    TFORM1  = 'D       '                                                            
    TUNIT1  = 'A       '                                                            
    TTYPE2  = 'flux    '                                                            
    TFORM2  = 'D       '                                                            
    TUNIT2  = 'erg/s/cm2/A'                                                         
    TTYPE3  = 'dered_flux'                                                          
    TFORM3  = 'D       '                                                            
    TUNIT3  = 'erg/s/cm2/A'                                                         
    TEFFL   =    23000.68377498454                                                  
    TEFF    =    23644.49138963689                                                  
    TEFFU   =    29999.33189058337                                                  
    LOGGL   =    3.014445328877565                                                  
    LOGG    =    4.984788506855546                                                  
    LOGGU   =    4.996095055525759                                                  
    EBVL    =   0.4703171900728142                                                  
    EBV     =   0.4933185398871652                                                  
    EBVU    =   0.5645144879211454                                                  
    ZL      =   -2.499929434601112                                                  
    Z       =   0.4332336811329144     
    ZU      =   0.4999776537652627                                                  
    SCALEL  = 2.028305798068863E-20                                                 
    SCALE   = 2.444606991671813E-20                                                 
    SCALEU  = 2.830281842143698E-20                                                 
    E_SCALEL= 8.844041760265867E-22                                                 
    E_SCALE = 1.391012126969951E-21                                                 
    E_SCALEU= 2.366112834925025E-21                                                 
    LABSL   =    250.9613352757437                                                  
    LABS    =     281.771013453664                                                  
    LABSU   =    745.3149766772975                                                  
    CHISQL  =    77.07958742733673                                                  
    CHISQ   =    77.07958742733673                                                  
    CHISQU  =    118.8587169011471                                                      
    CI_RAWL =   0.9999999513255379                                                  
    CI_RAW  =   0.9999999513255379                                                  
    CI_RAWU =   0.9999999999999972                                                  
    CI_RED  =   0.5401112973063139                                                  
    CI_REDL =   0.5401112973063139                                                  
    CI_REDU =   0.9500015229597392                                                  
    
    EXTNAME = 'IGRID_SEARCH' 
    XTENSION= 'BINTABLE'           / binary table extension                         
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   80 / length of dimension 1                          
    NAXIS2  =                99996 / length of dimension 2                          
    PCOUNT  =                    0 / number of group parameters                     
    GCOUNT  =                    1 / number of groups                               
    TFIELDS =                   10 / number of table fields                         
    TTYPE1  = 'teff    '                                                            
    TFORM1  = 'D       '                                                            
    TTYPE2  = 'logg    '                                                            
    TFORM2  = 'D       '                                                            
    TTYPE3  = 'ebv     '                                                            
    TFORM3  = 'D       '                                                            
    TTYPE4  = 'z       '                                                            
    TFORM4  = 'D       '                                                            
    TTYPE5  = 'chisq   '                                                            
    TFORM5  = 'D       '                                                            
    TTYPE6  = 'scale   '                                                            
    TFORM6  = 'D       '                                                            
    TTYPE7  = 'e_scale '                                                            
    TFORM7  = 'D       '                                                            
    TTYPE8  = 'Labs    '                                                            
    TFORM8  = 'D       '                                                            
    TTYPE9  = 'CI_raw  '                                                            
    TFORM9  = 'D       '                                                            
    TTYPE10 = 'CI_red  '                                                            
    TFORM10 = 'D       '                                                            


Subsection 2.2 Loading SED fits
-------------------------------

Once saved, you can load the contents of the FITS file again into an SED object
via

>>> mysed = SED('HD180642')
>>> mysed.load_fits()

and then you can build all the plots again easily. You can of course use the
predefined plotting scripts to start a plot, and then later on change the
properties of the labels, legends etc... for higher quality plots or to better
suit your needs.

"""
import sys
import time
import os
import logging

import pylab as pl
from matplotlib import mlab
import Image
import numpy as np
import scipy.stats
from scipy.interpolate import Rbf
import pyfits

from ivs import config
from ivs.aux import numpy_ext
from ivs.aux.decorators import memoized
from ivs.io import ascii
from ivs.io import fits
from ivs.sed import model
from ivs.sed import filters
from ivs.sed import fit
from ivs.sed import distance
from ivs.sed import extinctionmodels
from ivs.catalogs import crossmatch
from ivs.catalogs import vizier
from ivs.catalogs import mast
from ivs.catalogs import sesame
from ivs.units import conversions
from ivs.units import constants

from ivs.units.uncertainties import unumpy,ufloat
from ivs.units.uncertainties.unumpy import sqrt as usqrt
from ivs.units.uncertainties.unumpy import tan as utan

logger = logging.getLogger("SED.BUILD")
#logger.setLevel(10)

def fix_master(master,e_default=None):
    """
    Clean/extend/fix record array received from C{get_photometry}.
    
    This function does a couple of things:
    
        1. Adds common but uncatalogized colors like 2MASS.J-H if not already
           present. WARNING: these colors can mix values from the same system
           but from different catalogs!
        2. Removes photometry for which no calibration is available
        3. Adds a column 'color' with flag False denoting absolute flux measurement
        and True denoting color.
        4. Adds a column 'include' with flag True meaning the value will be
        included in the fit
        5. Sets some default errors to photometric values, for which we know
        that the catalog values are not trustworthy.
        6. Sets a lower limit to the allowed errors of 1%. Errors below this
        value are untrostworthy because the calibration error is larger than that.
        7. USNOB1 and ANS photometry are set the have a minimum error of 30%
        due to uncertainties in response curves.
        
    To do BUG FIX: only make colors from measurements from the same source
    catalog!
    
    @param master: record array containing photometry. This should have the fields
    'meas','e_meas','unit','photband','source','_r','_RAJ2000','DEJ2000',
    'cmeas','e_cmeas','cwave','cunit'
    @type master: record array
    @param e_default: default error for measurements without errors
    @type e_default: float
    @return: record array extended with fields 'include' and 'color', and with
    rows added (see above description)
    @rtype: numpy recard array
    """
    #-- we recognize uncalibrated stuff as those for which no absolute flux was
    #   obtained, and remove them:
    master = master[-np.isnan(master['cmeas'])]
    
    #-- add common uncatalogized colors:
    columns = list(master.dtype.names)
    for color in ['2MASS.J-H','2MASS.KS-H','TD1.1565-1965','TD1.2365-1965','TD1.2365-2740',
                  'JOHNSON.J-H','JOHNSON.K-H','JOHNSON.B-V','JOHNSON.U-B','GALEX.FUV-NUV',
                  'TYCHO2.BT-VT','ANS.15N-15W','ANS.15W-18','ANS.18-22','ANS.22-25','ANS.25-33']:
        #-- get the filter system and the separate bands for these colors
        system,band = color.split('.')
        band0,band1 = band.split('-')
        band0,band1 = '%s.%s'%(system,band0),'%s.%s'%(system,band1)
        
        #-- if the separate bands are available (and not the color itself),
        #   calculate the colors here
        if band0 in master['photband'] and band1 in master['photband'] and not color in master['photband']:
            #-- where are the bands located?
            index0 = list(master['photband']).index(band0)
            index1 = list(master['photband']).index(band1)
            
            #-- start a new row to add the color
            row = list(master[index1])
            row[columns.index('photband')] = color
            row[columns.index('cwave')] = np.nan
            
            #-- it could be a magnitude difference
            if master['unit'][index0]=='mag':
                row[columns.index('meas')] = master['meas'][index0]-master['meas'][index1]
                row[columns.index('e_meas')] = np.sqrt(master['e_meas'][index0]**2+master['e_meas'][index1]**2)
                #-- error will not always be available...
                try:
                    row[columns.index('cmeas')],row[columns.index('e_cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],row[columns.index('e_meas')],photband=color)
                except AssertionError:
                    row[columns.index('cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],photband=color)
                    row[columns.index('e_cmeas')] = np.nan
            #-- or it could be a flux ratio
            else:
                row[columns.index('meas')] = master['meas'][index0]/master['meas'][index1]
                row[columns.index('e_meas')] = np.sqrt(((master['e_meas'][index0]/master['meas'][index0])**2+(master['e_meas'][index1]/master['meas'][index1])**2)*row[columns.index('meas')]**2)
                row[columns.index('cmeas')],row[columns.index('e_cmeas')] = row[columns.index('meas')],row[columns.index('e_meas')]
            row[columns.index('cunit')] = 'flux_ratio'
            master = numpy_ext.recarr_addrows(master,[tuple(row)])
            
    #-- add an extra column with a flag to distinguish colors from absolute
    #   fluxes, and a column with flags to include/exclude photometry
    #   By default, exclude photometry if the effective wavelength is above
    #   150 mum or if the measurement is nan
    extra_cols = [[filters.is_color(photband) for photband in master['photband']],
                  [(cwave<150e4) or np.isnan(meas) for cwave,meas in zip(master['cwave'],master['cmeas'])]]
    dtypes = [('color',np.bool),('include',np.bool)]
    master = numpy_ext.recarr_addcols(master,extra_cols,dtypes)
    
    #-- set default errors if no errors are available and set really really
    #   small errors to some larger default value
    if e_default is not None:
        no_error = np.isnan(master['e_cmeas'])
        master['e_cmeas'][no_error] = np.abs(e_default*master['cmeas'][no_error])
        small_error = master['e_cmeas']<(0.01*np.abs(master['cmeas']))
        master['e_cmeas'][small_error] = 0.01*np.abs(master['cmeas'][small_error])
    
    #-- the measurements from USNOB1 are not that good because the response
    #   curve is approximated by the JOHNSON filter. Set the errors to min 30%
    #   The same holds for ANS: here, the transmission curves are very uncertain
    for i,photband in enumerate(master['photband']):
        if 'USNOB1' in photband or 'ANS' in photband:
            master['e_cmeas'][i] = max(0.30*master['cmeas'][i],master['e_cmeas'][i])
    
    return master


def decide_phot(master,names=None,wrange=None,ptype='all',include=False):
    """
    Exclude/include photometric passbands containing one of the strings listed in
    photbands.
    
    This function will set the flags in the 'include' field of the extended
    master record array.
    
    Colours also have a wavelength in this function, and it is defined as the
    average wavelength of the included bands (doesn't work for stromgren C1 and M1).
    
    include=False excludes photometry based on the input parameters.
    include=True includes photometry based on the input parameters.
    
    Some examples:
    
        1. Exclude all measurements:
        >>> decide_phot(master,wrange=(-np.inf,+np.inf),ptype='all',include=False)
    
        2. Include all TD1 fluxes and colors:
        >>> decide_phot(master,names=['TD1'],ptype='all',include=True)
    
        3. Include all V band measurements from all systems (but not the colors):
        >>> decide_phot(master,names=['.V'],ptype='abs',include=True)
    
        4. Include all Geneva colors and exclude Geneva magnitudes:
        >>> decide_phot(master,names=['GENEVA'],ptype='col',include=True)
        >>> decide_phot(master,names=['GENEVA'],ptype='abs',include=False)
    
        5. Exclude all infrared measurements beyond 1 micron:
        >>> decide_phot(master,wrange=(1e4,np.inf),ptype='all',include=False)
    
        6. Include all AKARI measurements below 10 micron:
        >>> decide_phot(master,names=['AKARI'],wrange=(-np.inf,1e5),ptype='all',include=True)
    
    @param master: record array containing all photometry
    @type master: numpy record array
    @param names: strings excerpts to match filters
    @type names: list of strings
    @param wrange: wavelength range (most likely angstrom) to include/exclude
    @type wrange: 2-tuple (start wavelength,end wavelength)
    @param ptype: type of photometry to include/exclude: absolute values, colors
    or both
    @type ptype: string, one of 'abs','col','all'
    @param include: flag setting exclusion or inclusion
    @type include: boolean
    @return: master record array with adapted 'include' flags
    @rtype: numpy record array
    """ 
    #-- exclude/include passbands based on their names
    if names is not None:
        for index,photband in enumerate(master['photband']):
            for name in names:
                if name in photband:
                    if ptype=='all' or (ptype=='abs' and -master['color'][index]) or (ptype=='col' and master['color'][index]):
                        master['include'][index] = include
                        break
    #-- exclude/include colors based on their wavelength
    if wrange is not None:
        for index,photband in enumerate(master['photband']):
            system,color = photband.split('.')
            if not '-' in color or ptype=='abs':
                continue
            band0,band1 = color.split('-')
            cwave0 = filters.eff_wave('%s.%s'%(system,band0))
            cwave1 = filters.eff_wave('%s.%s'%(system,band1))
            if (wrange[0]<cwave0<wrange[1]) or (wrange[0]<cwave0<wrange[1]):
                master['include'][index] = include
        #-- exclude/include passbands based on their wavelength
        if not ptype=='col':
            master['include'][(wrange[0]<master['cwave']) & (master['cwave']<wrange[1])] = include    


def photometry2str(master):
    """
    String representation of master record array
    
    @param master: master record array containing photometry
    @type master: numpy record array
    """
    master = master[np.argsort(master['photband'])]
    txt = '%20s %12s %12s %12s %10s %12s %12s %11s %s\n'%('PHOTBAND','MEAS','E_MEAS','UNIT','CWAVE','CMEAS','E_CMEAS','UNIT','SOURCE')
    txt+= '==========================================================================================================================\n'
    for i,j,k,l,m,n,o,p in zip(master['photband'],master['meas'],master['e_meas'],master['unit'],master['cwave'],master['cmeas'],master['e_cmeas'],master['source']):
        txt+='%20s %12g %12g %12s %10.0f %12g %12g erg/s/cm2/A %s\n'%(i,j,k,l,m,n,o,p)
    return txt

@memoized
def get_schaller_grid():
    """
    Download Schaller 1992 evolutionary tracks and return an Rbf interpolation
    function.
    
    @return: Rbf interpolation function
    @rtype: Rbf interpolation function
    """
    #-- translation between table names and masses
    masses = [1,1.25,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,40,60][:-1]
    tables = ['table20','table18','table17','table16','table15','table14',
              'table13','table12','table11','table10','table9','table8',
              'table7','table6','table5','table4','table3'][:-1]
    #-- read in all the tables and compute radii and luminosities.
    all_teffs,all_loggs,all_radii = [],[],[]
    for mass,table in zip(masses,tables):
        data,comms,units = vizier.search('J/A+AS/96/269/%s'%(table))
        all_teffs.append(10**data['logTe'])
        all_radii.append(np.sqrt((10**data['logL']*constants.Lsol_cgs)/(10**data['logTe'])**4/(4*np.pi*constants.sigma_cgs)))
        all_loggs.append(np.log10(constants.GG_cgs*mass*constants.Msol_cgs/(all_radii[-1]**2)))
    all_teffs = np.hstack(all_teffs)
    all_radii = np.hstack(all_radii)/constants.Rsol_cgs
    all_loggs = np.hstack(all_loggs)
    #-- remove low temperature models, the evolutionary tracks are hard to
    #   interpolate there.
    keep = all_teffs>5000
    all_teffs = all_teffs[keep]
    all_radii = all_radii[keep]
    all_loggs = all_loggs[keep]
    #-- make linear interpolation model between all modelpoints
    mygrid = Rbf(np.log10(all_teffs),all_loggs,all_radii,function='linear')
    logger.info('Interpolation of Schaller 1992 evolutionary tracks to compute radii')
    return mygrid


def get_radii(teffs,loggs):
    """
    Retrieve radii from stellar evolutionary tracks from Schaller 1992.
    
    @param teffs: model effective temperatures
    @type teffs: numpy array
    @param loggs: model surface gravities
    @type loggs: numpy array
    @return: model radii (solar units)
    @rtype: numpy array
    """
    mygrid = get_schaller_grid()
    radii = mygrid(np.log10(teffs),loggs)
    return radii


def calculate_distance(plx,gal,teffs,loggs,scales,n=75000):
    """
    Calculate distances and radii of a target given its parallax and location
    in the galaxy.
    
    
    """
    #-- compute distance up to 25 kpc, which is about the maximum distance from
    #   earth to the farthest side of the Milky Way galaxy
    #   rescale to set the maximum to 1
    d = np.logspace(np.log10(0.1),np.log10(25000),100000)
    if plx is not None:
        dprob = distance.distprob(d,gal[1],plx)
        dprob = dprob / dprob.max()
    else:
        dprob = np.ones_like(d)
    #-- compute the radii for the computed models, and convert to parsec
    #radii = np.ones(len(teffs[-n:]))
    radii = get_radii(teffs[-n:],loggs[-n:])
    radii = conversions.convert('Rsol','pc',radii)
    d_models = radii/np.sqrt(scales[-n:])
    #-- we set out of boundary values to zero
    if plx is not None:
        dprob_models = np.interp(d_models,d[-n:],dprob[-n:],left=0,right=0)
    else:
        dprob_models = np.ones_like(d_models)
    #-- reset the radii to solar units for return value
    radii = conversions.convert('pc','Rsol',radii)
    return (d_models,dprob_models,radii),(d,dprob)
    


class SED(object):
    """
    Class that facilitates the use of the ivs.sed module.
    
    This class is meant to be an easy interface to many of the ivs.sed module's
    functionality.
    """
    def __init__(self,ID,photfile=None,plx=None):
        """
        Initialize SED class.
        
        @param plx: parallax (and error) of the object
        @type plx: tuple (plx,e_plx)
        """
        self.ID = ID
        #-- the file containing photometry should have the following name. We
        #   save photometry to a file to avoid having to download photometry
        #   each time from the internet
        if photfile is None:
            self.photfile = '%s.phot'%(ID).replace(' ','_')
        else:
            self.photfile = photfile
        
        #-- keep information on the star from SESAME, but override parallax with
        #   the value from Van Leeuwen's new reduction. Set the galactic
        #   coordinates, which are not available from SESAME.
        try:
            self.info = sesame.search(ID,fix=True)
        except KeyError:
            logger.warning('Star %s not recognised by sesame'%(ID))
            self.info = {}
        if plx is not None:
            if not 'plx' in self.info:
                self.info['plx'] = {}
            self.info['plx']['v'] = plx[0]
            self.info['plx']['e'] = plx[1]
        
        #-- prepare for information on fitting processes
        self.results = {}
    
    #{ Handling photometric data
    def get_photometry(self,radius=None,ra=None,dec=None,include=None,exclude=None):
        """
        Search photometry on the net or from the phot file if it exists.
        
        For bright stars, you can set radius a bit higher...
        
        @param radius: search radius (arcseconds)
        @type radius: float.
        """
        if radius is None:
            if 'mag.V.v' in self.info and self.info['mag.V.v']<6.:
                radius = 60.
            else:
                radius = 10.
        if not os.path.isfile(self.photfile):
            #-- get and fix photometry. Set default errors to 1%, and set
            #   USNOB1 errors to 3%
            if ra is None and dec is None:
                master = crossmatch.get_photometry(ID=self.ID,radius=radius,
                                       include=include,exclude=exclude,
                                       extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            else:
                master = crossmatch.get_photometry(ID=self.ID,ra=ra,dec=dec,radius=radius,
                                       include=include,exclude=exclude,
                                       extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            if 'jradeg' in self.info:
                master['_RAJ2000'] -= self.info['jradeg']
                master['_DEJ2000'] -= self.info['jdedeg']
            
            #-- fix the photometry: set default errors to 2% and print it to the
            #   screen
            self.master = fix_master(master,e_default=0.1)
            print(photometry2str(master))
            
            #-- write to file
            self.save_photometry()
        else:
            self.load_photometry()
    
    def exclude(self,names=None,wrange=None):
        """
        Exclude photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=False,ptype='all')
    
    def exclude_colors(self,names=None,wrange=None):
        """
        Exclude photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=False,ptype='col')
    
    def exclude_abs(self,names=None,wrange=None):
        """
        Exclude photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=False,ptype='abs')
    
    
    def include(self,names=None,wrange=None):
        """
        Include photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=True,ptype='all')
    
    def include_colors(self,names=None,wrange=None):
        """
        Include photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=True,ptype='col')
    
    def include_abs(self,names=None,wrange=None):
        """
        Include photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,include=True,ptype='abs')
    
    def set_photometry_scheme(self,scheme):
        """
        Scheme = 'abs': means excluding all colors, including all absolute values
        scheme = 'color': means including all colors, excluding all absolute values
        scheme = 'combo': means inculding all colors, and one absolute value per
        system (the one with the smallest relative error)
        """
        if 'abs' in scheme.lower():
            self.master['include'][self.master['color']] = False
            self.master['include'][-self.master['color']] = True
            logger.info('Fitting procedure will use only absolute fluxes (%d)'%(sum(self.master['include'])))
        elif 'col' in scheme.lower():
            self.master['include'][self.master['color']] = True
            self.master['include'][-self.master['color']] = False
            logger.info('Fitting procedure will use only colors (%d)'%(sum(self.master['include'])))
        elif 'com' in scheme.lower():
            self.master['include'][self.master['color']] = True
            self.master['include'][-self.master['color']] = False
            systems = np.array([photband.split('.')[0] for photband in self.master['photband']])
            set_systems = sorted(list(set(systems)))
            for isys in set_systems:
                keep = -self.master['color'] & (isys==systems)
                if not sum(keep): continue
                index = np.argmax(self.master['e_cmeas'][keep]/self.master['cmeas'][keep])
                index = np.arange(len(self.master))[keep][index]
                self.master['include'][index] = True
            logger.info('Fitting procedure will use colors + one absolute flux for each system (%d)'%(sum(self.master['include'])))
            
        
    
    def save_photometry(self,photfile=None):
        #-- write to file
        if photfile is not None:
            self.photfile = photfile
        logger.info('Save photometry to file %s'%(self.photfile))
        ascii.write_array(self.master,self.photfile,header=True,auto_width=True,use_float='%g')
    
    def load_photometry(self):
        logger.info('Load photometry from file %s'%(self.photfile))
        self.master = ascii.read2recarray(self.photfile)
    
    def add_photometry_fromarrays(self,meas,e_meas,units,photbands,source):
        """
        Add photometry from a file or arrays
        
        By default unflagged, and at the ra and dec of the star. No color and
        included.
        """
        extra_cols = [np.zeros(len(meas))]
        extra_cols.append(np.zeros(len(meas)))
        extra_cols.append(np.zeros(len(meas)))
        extra_cols.append(np.zeros(len(meas),dtype=str))
        extra_cols.append(np.zeros(len(meas)))
        extra_cols.append(np.zeros(len(meas)))
        extra_cols.append(np.zeros(len(meas)))
        extra_cols.append(np.zeros(len(meas),dtype=bool))
        extra_cols.append(np.zeros(len(meas),dtype=bool))
        for i,(m,e_m,u,p,s) in enumerate(zip(meas,e_meas,units,photbands,source)):
            photsys,photband = p.split('.')
            if photband in ['C1','M1'] or '-' in photband:
                to_unit = 'flux_ratio'
                color = True
            else:
                to_unit = self.master['cunit'][0]
                color = False
            print u,to_unit
            if e_m>0:
                cm,e_cm = conversions.convert(u,to_unit,m,e_m,photband=p)
            else:
                cm,e_cm = conversions.convert(u,to_unit,m,photband=p),np.nan
            eff_wave = filters.eff_wave(p)
            extra_cols[0][i] = cm
            extra_cols[1][i] = e_cm
            extra_cols[2][i] = eff_wave
            extra_cols[3][i] = self.master['cunit'][0]
            extra_cols[7][i] = color
            extra_cols[8][i] = True
            
        extra_cols += [meas,e_meas,units,photbands,source,np.nan*np.zeros(len(meas))]
        #meas     e_meas flag unit photband          source             _r    _RAJ2000   _DEJ2000   cwave       cmeas     e_cmeas cunit       color include
        extra_array = np.rec.fromarrays(extra_cols,names=['cmeas','e_cmeas','cwave','cunit','_r','_RAJ2000','_DEJ2000','color','include','meas','e_meas','unit','photbands','source','flag'])
        print '\n\nblabla'
        print pl.mlab.rec2txt(self.master,precision=6)
        print pl.mlab.rec2txt(extra_array,precision=6)
        print 'blabla\n\n'
        self.master = np.hstack([self.master,extra_array])
        
        print self.master

    #}
    
    #{ Fitting routines
    def igrid_search(self,teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf),
                          ebvrange=(-np.inf,np.inf),zrange=(-np.inf,np.inf),
                          threads='max',iterations=3,increase=1,speed=1,res=1,
                          points=None,compare=True):
        #-- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are include in the fitting process:\n%s'%(photometry2str(self.master[include_grid])))
        
        
        #-- build the grid
        teffs,loggs,ebvs,zs = fit.generate_grid(self.master['photband'][include_grid],teffrange=teffrange,
                                               loggrange=loggrange,ebvrange=ebvrange,zrange=zrange,
                                               res=res,points=points)
        #-- run over the grid and calculate the CHI2
        chisqs,scales,e_scales,lumis = fit.igrid_search(teffs,loggs,ebvs,zs,
                                 self.master['cmeas'][include_grid],
                                 self.master['e_cmeas'][include_grid],
                                 self.master['photband'][include_grid])
        grid_results = np.rec.fromarrays([teffs,loggs,ebvs,zs,chisqs,scales,e_scales,lumis],
                     dtype=[('teff','f8'),('logg','f8'),('ebv','f8'),('z','f8'),
                  ('chisq','f8'),('scale','f8'),('e_scale','f8'),('Labs','f8')])
                
        #-- exclude failures
        grid_results = grid_results[-np.isnan(grid_results['chisq'])]
        
        #-- inverse sort according to chisq: this means the best models are at
        #   the end (mainly for plotting reasons, so that the best models
        #   are on top).
        sa = np.argsort(grid_results['chisq'])[::-1]
        grid_results = grid_results[sa]
        
        #-- do the statistics
        #   degrees of freedom: teff,logg,E(B-V),theta,Z
        N,df = sum(include_grid),5
        k = N-df
        if k<=0:
            logger.warning('Not enough data to compute CHI2: it will not make sense')
            k = 1
        #   rescale if needed and compute confidence intervals
        factor = max(grid_results['chisq'][-1]/k,1)
        logger.warning('CHI2 rescaling factor equals %g'%(factor))
        CI_raw = scipy.stats.distributions.chi2.cdf(grid_results['chisq'],k)
        CI_red = scipy.stats.distributions.chi2.cdf(grid_results['chisq']/factor,k)
        
        #-- add the results to the record array and to the results dictionary
        grid_results = mlab.rec_append_fields(grid_results, 'CI_raw', CI_raw)
        grid_results = mlab.rec_append_fields(grid_results, 'CI_red', CI_red)
        if not 'igrid_search' in self.results:
            self.results['igrid_search'] = {}
        
        self.results['igrid_search']['grid'] = grid_results
        self.results['igrid_search']['factor'] = factor
        self.results['igrid_search']['CI'] = {}
        
        start_CI = np.argmin(np.abs(CI_red-0.95))
        for name in grid_results.dtype.names:
            self.results['igrid_search']['CI'][name+'L'] = grid_results[name][start_CI:].min()
            self.results['igrid_search']['CI'][name] = grid_results[name][-1]
            self.results['igrid_search']['CI'][name+'U'] = grid_results[name][start_CI:].max()
            logger.info('95%% CI %s: %g <= %g <= %g'%(name,self.results['igrid_search']['CI'][name+'L'],
                                                           self.results['igrid_search']['CI'][name],
                                                           self.results['igrid_search']['CI'][name+'U']))
        self.set_best_model()
        if False:
            self.set_distance()
            d,dprob = self.results['igrid_search']['d'] # in parsecs
            dcumprob = np.cumsum(dprob[1:]*np.diff(d))
            dcumprob /= dcumprob.max()
            d_min = d[dcumprob>0.38].min()
            d_best = d[np.argmax(dprob)]
            d_max = d[dcumprob<0.38].max()
            e_d_best = min([abs(d_best-d_min),abs(d_max-d_best)])
            print d_best,e_d_best
            d_best = ufloat((d_best,e_d_best))

            scales = self.results['igrid_search']['grid']['scale']
            R = scales*d_best
            print R
            
            #radii = np.zeros_like()
            #for di,dprobi in zip(d,dprob):
            #    self.results['igrid_search']['grid']['scale']*d_min
            
        
        
        
        
    #}
    
    #{ Interfaces
    
    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004'):
        """
        Get reddenend and unreddened model
        """
        #-- get reddened and unreddened model
        logger.info('Interpolating approximate full SED of best model')
        scale = self.results[mtype]['CI']['scale']
        wave,flux = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                    logg=self.results[mtype]['CI']['logg'],
                                    ebv=self.results[mtype]['CI']['ebv'],
                                    law=law)
        wave_ur,flux_ur = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                          logg=self.results[mtype]['CI']['logg'],
                                          ebv=0,
                                          law=law)
        flux,flux_ur = flux*scale,flux_ur*scale
        
        #-- synthetic flux
        include = self.master['include']
        synflux = np.zeros(len(self.master['photband']))
        keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
        synflux_,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
                                  logg=self.results[mtype]['CI']['logg'],
                                  ebv=self.results[mtype]['CI']['ebv'],
                                  photbands=self.master['photband'][keep])
        synflux[keep] = synflux_
        
        #synflux,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
        #                          logg=self.results[mtype]['CI']['logg'],
        #                          ebv=self.results[mtype]['CI']['ebv'],
        #                          photbands=self.master['photband'])
        synflux[-self.master['color']] *= scale
        chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
        #-- calculate effective wavelengths of the photometric bands via the model
        #   values
        eff_waves = filters.eff_wave(self.master['photband'],model=(wave,flux))
        self.results['model'] = wave,flux,flux_ur
        self.results['synflux'] = eff_waves,synflux,self.master['photband']
        self.results['chi2'] = chi2
    
    def set_distance(self):
        #-- calculate the distance
        cutlogg = (self.results['igrid_search']['grid']['logg']<=4.4) & (self.results['igrid_search']['grid']['CI_red']<=0.95)
        gal = self.info['galpos']
        if 'plx' in self.info:
            plx = self.info['plx']['v'],self.info['plx']['e']
        else:
            plx = None
        (d_models,dprob_models,radii),(d,dprob)\
                   = calculate_distance(plx,self.info['galpos'],self.results['igrid_search']['grid']['teff'][cutlogg],\
                                                   self.results['igrid_search']['grid']['logg'][cutlogg],\
                                                   self.results['igrid_search']['grid']['scale'][cutlogg])
        #-- calculate model extinction
        res = 100
        self.results['igrid_search']['drimmel'] = np.ravel(np.array([extinctionmodels.findext(gal[0], gal[1], model='drimmel', distance=myd) for myd in d[::res]]))
        self.results['igrid_search']['marshall'] = np.ravel(np.array([extinctionmodels.findext(gal[0], gal[1], model='marshall', distance=myd) for myd in d[::res]]))
        self.results['igrid_search']['d_mod'] = (d_models,dprob_models,radii)
        self.results['igrid_search']['d'] = (d,dprob)
    
    #}
    
    #{ Input and output
    def plot_grid(self,x='teff',y='logg',ptype='CI_red',mtype='igrid_search',limit=0.95):
        """
        Plot grid as scatter plot
        
        Possible plot types: 'CI_red','z','ebv'
        """
        if limit is not None:
            region = self.results[mtype]['grid']['CI_red']<limit
        else:
            region = self.results[mtype]['grid']['CI_red']<np.inf
        #-- get the colors and the color scale
        colors = self.results[mtype]['grid'][ptype][region]
        vmin,vmax = colors.min(),colors.max()
        if 'CI' in ptype:
            colors *= 100.
            vmax = 95.
            vmin = colors.min()
        #-- grid scatter plot
        pl.scatter(self.results[mtype]['grid'][x][region],
                   self.results[mtype]['grid'][y][region],
             c=colors,edgecolors='none',cmap=pl.cm.spectral,vmin=vmin,vmax=vmax)
        #-- mark best value
        pl.plot(self.results[mtype]['grid'][x][-1],
                self.results[mtype]['grid'][y][-1],'r+',ms=40,mew=3)
        #-- set the limits to only include the 95 interval
        pl.xlim(self.results[mtype]['grid'][x][region].max(),
                self.results[mtype]['grid'][x][region].min())
        pl.ylim(self.results[mtype]['grid'][y][region].max(),
                self.results[mtype]['grid'][y][region].min())
        cbar = pl.colorbar()
        
        #-- set the x/y labels
        if   x=='teff': pl.xlabel('log (effective temperature [K]) [dex]')
        elif x=='z'   : pl.xlabel('log (Metallicity Z [$Z_\odot$]) [dex]')
        elif x=='logg': pl.xlabel(r'log (surface gravity [cm s$^{-2}$]) [dex]')
        elif x=='ebv' : pl.xlabel('E(B-V) [mag]')
        if   y=='teff': pl.ylabel('log (effective temperature [K]) [dex]')
        elif y=='z'   : pl.ylabel('log (Metallicity Z [$Z_\odot$]) [dex]')
        elif y=='logg': pl.ylabel(r'log (surface gravity [cm s$^{-2}$]) [dex]')
        elif y=='ebv' : pl.ylabel('E(B-V) [mag]')
        
        #-- set the colorbar label
        if 'CI'in ptype:   cbar.set_label('Probability [%]')
        elif ptype=='z':   cbar.set_label('Metallicity Z ($Z_\odot$)')
        elif ptype=='ebv': cbar.set_label('E(B-V) [mag]')
        
        logger.info('Plotted teff-logg diagram of %s'%(ptype))
    
    def plot_data(self):
        """
        Plot only the SED data.
        """
        wave,flux,e_flux = self.master.cwave,self.master.cmeas,self.master.e_cmeas
        iscolor = np.array(self.master.color,bool)
        photbands = self.master.photband
        allsystems = np.array([i.split('.')[0] for i in photbands])
        systems = sorted(set(allsystems))
        
        color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(systems))]
        pl.gca().set_color_cycle(color_cycle)
        pl.gca().set_xscale('log',nonposx='clip')
        pl.gca().set_yscale('log',nonposy='clip')
        
        for system in systems:
            keep = (allsystems==system) & -iscolor
            if keep.sum():
                pl.errorbar(wave[keep],flux[keep],yerr=e_flux[keep],fmt='o',label=system,ms=7)
        pl.legend(prop=dict(size='small'))
        pl.grid()
        pl.title(self.ID)
        pl.ylabel(r'$F_\lambda$ [erg/s/cm$^2/A$]')
        pl.xlabel('wavelength [$\AA$]')
        
    
    
    def plot_sed(self,colors=False,mtype='igrid_search'):
        """
        Plot a fitted SED together with the data.
        """
        
        def plot_sed_getcolors(master,color_dict=None):
            myphotbands = [iphotb.split('.')[1] for iphotb in master['photband'][master['color']]]
            if color_dict is None:
                color_dict = {myphotbands[0]:0}
                for mycol in myphotbands[1:]:
                    if not mycol in color_dict:
                        max_int = max([color_dict[entry] for entry in color_dict])
                        color_dict[mycol] = max_int+1
            x = [color_dict[mycol] for mycol in myphotbands]
            y = master['cmeas']
            e_y = master['e_cmeas']
            return x,y,e_y,color_dict
            
        
        x__,y__,e_y__,color_dict = plot_sed_getcolors(self.master)
        
        systems = np.array([system.split('.')[0] for system in self.master['photband']],str)
        set_systems = sorted(list(set(systems)))
        color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))]
        pl.gca().set_color_cycle(color_cycle)
        
        #-- for plotting reasons, we translate every color to an integer
        
        for system in set_systems:
            keep = (systems==system) & (self.master['color']==colors)
            #-- synthetic:
            pl.gca()._get_lines.count += 1
            if sum(keep):
                pl.gca()._get_lines.count -= 1
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                    y = self.results['synflux'][1][keep]
                else:
                    x = self.results['synflux'][0][keep]
                    y = x*self.results['synflux'][1][keep]
                pl.plot(x,y,'x',ms=10,mew=2,alpha=0.75)
            #-- include:
            keep = (systems==system) & (self.master['color']==colors) & self.master['include']
            if sum(keep):
                pl.gca()._get_lines.count -= 1
                if colors:
                    #-- translate every color to an integer
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results['synflux'][0][keep]
                    y = x*self.master['cmeas'][keep]
                    e_y = x*self.master['e_cmeas'][keep]
                p = pl.errorbar(x,y,yerr=e_y,fmt='o',label=system,
                                capsize=10,ms=7)
            
            #-- exclude:
            label = np.any(keep) and '_nolegend_' or system
            keep = (systems==system) & (self.master['color']==colors) & -self.master['include']
            if sum(keep):
                pl.gca()._get_lines.count -= 1
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results['synflux'][0][keep]
                    y = x*self.master['cmeas'][keep]
                    e_y = x*self.master['e_cmeas'][keep]
                    
                pl.errorbar(x,y,yerr=e_y,fmt='o',label=label,
                            capsize=10,ms=7,mew=2,mfc='1',mec=color_cycle[pl.gca()._get_lines.count%len(color_cycle)])
        
        #-- only set logarithmic scale if absolute fluxes are plotted
        #   and only plot the real model then
        if not colors:
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            pl.gca().set_autoscale_on(False)
        
            #-- the model
            wave,flux,flux_ur = self.results['model']
            pl.plot(wave,wave*flux,'r-')
            pl.plot(wave,wave*flux_ur,'k-')
            pl.ylabel(r'$\lambda F_\lambda$ [erg/s/cm$^2$]')
            pl.xlabel('wavelength [$\AA$]')
        else:
            xlabels = color_dict.keys()
            xticks = [color_dict[key] for key in xlabels]
            pl.xticks(xticks,xlabels,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Color name')
            pl.xlim(min(xticks)-0.5,max(xticks)+0.5)
        
        pl.grid()
        if colors:
            pl.legend(loc='best',prop=dict(size='x-small'))
        else:
            pl.legend(loc='upper right',prop=dict(size='x-small'))
        loc = (0.05,0.05)
        teff = self.results[mtype]['grid']['teff'][-1]
        logg = self.results[mtype]['grid']['logg'][-1]
        ebv = self.results[mtype]['grid']['ebv'][-1]
        scale = self.results[mtype]['grid']['scale'][-1]
        angdiam = conversions.convert('sr','mas',4*np.pi*np.sqrt(scale))
        pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\n$\Theta$=%.3g mas'%(teff,logg,ebv,angdiam),
                    loc,xycoords='axes fraction')
        logger.info('Plotted SED as %s'%(colors and 'colors' or 'absolute fluxes'))
        
    def plot_chi2(self,colors=False,mtype='igrid_search'):
        
        include_grid = self.master['include']
        eff_waves,synflux,photbands = self.results['synflux']
        chi2 = self.results['chi2']
        systems = np.array([system.split('.')[0] for system in self.master['photband'][include_grid]],str)
        set_systems = sorted(list(set(systems)))
        pl.gca().set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0, 1.0, len(set_systems))])
        for system in set_systems:
            pl.gca()._get_lines.count += 1
            keep = systems==system
            if sum(keep) and not colors:
                pl.gca()._get_lines.count -= 1
                try:
                    pl.loglog(eff_waves[include_grid][keep],chi2[include_grid][keep],'o',label=system)
                except:
                    logger.critical('Plotting of CHI2 of absolute values failed')
            elif sum(keep) and colors:
                pl.gca()._get_lines.count -= 1
                pl.semilogy(range(len(eff_waves[include_grid][keep])),chi2[include_grid][keep],'o',label=system)
        pl.legend(loc='upper right',prop=dict(size='x-small'))
        pl.grid()
        k = 5.
        pl.annotate('Exp. val. = %d'%(k),(0.69,0.12),xycoords='axes fraction',color='b')
        pl.annotate('Total $\chi^2$ = %.1f'%(self.results[mtype]['grid']['chisq'][-1]),(0.69,0.065),xycoords='axes fraction',color='r')
        pl.annotate('Error scale = %.2f'%(np.sqrt(self.results[mtype]['factor'])),(0.69,0.03),xycoords='axes fraction',color='k')
        xlims = pl.xlim()
        pl.plot(xlims,[k,k],'b-',lw=2)
        pl.plot(xlims,[self.results[mtype]['grid']['chisq'][-1],self.results[mtype]['grid']['chisq'][-1]],'r-',lw=2)
        pl.xlim(xlims)
        pl.xlabel('wavelength [$\AA$]')
        pl.ylabel(r'Statistic ($\chi^2$)')
        logger.info('Plotted CHI2 of %s'%(colors and 'colors' or 'absolute fluxes'))
        
    def plot_distance(self,mtype='igrid_search'):
        #-- necessary information
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        
        ax_d = pl.gca()
        
        
        gal = self.info['galpos']
        #-- the plot
        dzoom = dprob>1e-4
        pl.plot(d,dprob,'k-')
        pl.grid()
        pl.xlabel('Distance [pc]')
        pl.ylabel('Probability [unnormalized]')
        pl.xlim(d[dzoom].min(),d[dzoom].max())
        xlims = pl.xlim()
        pl.twiny(ax_d)
        pl.xlim(xlims)
        xticks = pl.xticks()
        pl.xticks(xticks[0],['%.2f'%(conversions.convert('pc','Rsol',np.sqrt(self.results['igrid_search']['grid']['scale'][-1])*di)) for di in xticks[0]])
        pl.xlabel('Radius [$R_\odot$]')
        pl.twinx(ax_d)
        res = 100
        d_keep = (xlims[0]<=d[::res]) & (d[::res]<=xlims[1])
        if len(self.results['igrid_search']['drimmel']):
            pl.plot(d[::res][d_keep],self.results['igrid_search']['drimmel'].ravel()[d_keep],'b-',label='Drimmel')
        if len(self.results['igrid_search']['marshall']):
            pl.plot(d[::res][d_keep],self.results['igrid_search']['marshall'].ravel()[d_keep],'b--',label='Marshall')
        ebv = self.results[mtype]['grid']['ebv'][-1]
        pl.plot(xlims,[ebv*3.1,ebv*3.1],'r--',lw=2,label='measured')
        pl.ylabel('Visual extinction $A_v$ [mag]')
        pl.legend(loc='lower right',prop=dict(size='x-small'))
        pl.xlim(xlims)
        logger.info('Plotted distance/reddening')
    
    def plot_grid_model(self,ptype='prob'):
        """
        Grid of models
        """
        if 'spType' in self.info:
            pl.title(self.info['spType'])
        cutlogg = (self.results['igrid_search']['grid']['logg']<=4.4) & (self.results['igrid_search']['grid']['CI_red']<=0.95)
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        gal = self.info['galpos']
        
        n = 75000
        region = self.results['igrid_search']['grid']['CI_red']<0.95
        total_prob = 100-(1-self.results['igrid_search']['grid']['CI_red'][cutlogg][-n:])*d_prob_models*100
        tp_sa = np.argsort(total_prob)[::-1]
        if ptype=='prob':
            pl.scatter(self.results['igrid_search']['grid']['teff'][cutlogg][-n:][tp_sa],self.results['igrid_search']['grid']['logg'][cutlogg][-n:][tp_sa],
                c=total_prob[tp_sa],edgecolors='none',cmap=pl.cm.spectral,
                vmin=total_prob.min(),vmax=total_prob.max())
        elif ptype=='radii':
            pl.scatter(self.results['igrid_search']['grid']['teff'][cutlogg][-n:][tp_sa],self.results['igrid_search']['grid']['logg'][cutlogg][-n:][tp_sa],
                c=radii,edgecolors='none',cmap=pl.cm.spectral,
                vmin=radii.min(),vmax=radii.max())
        pl.xlim(self.results['igrid_search']['grid']['teff'][region].max(),self.results['igrid_search']['grid']['teff'][region].min())
        pl.ylim(self.results['igrid_search']['grid']['logg'][region].max(),self.results['igrid_search']['grid']['logg'][region].min())
        cbar = pl.colorbar()
        pl.xlabel('log (effective temperature [K]) [dex]')
        pl.ylabel(r'log (surface gravity [cm s$^{-2}$]) [dex]')
        
        if ptype=='prob':
            cbar.set_label('Probability (incl. plx) [%]')
        elif ptype=='radii':
            cbar.set_label('Model radii [$R_\odot$]')
        logger.info('Plotted teff-logg diagram of models (%s)'%(ptype))
        
        
    
    def plot_MW_side(self):
        im = Image.open(config.get_datafile('images','NRmilkyway.tif'))
        left,bottom,width,height = 0.0,0.0,1.0,1.0
        startm,endm = 183,-177
        startv,endv = -89,91

        xwidth = startm-endm
        ywidth = 90.
        ratio = ywidth/xwidth

        gal = list(self.info['galpos'])
        if gal[0]>180:
            gal[0] = gal[0] - 360.
        #-- boundaries of ESO image
        pl.imshow(im,extent=[startm,endm,startv,endv],origin='lower')
        pl.plot(gal[0],gal[1],'rx',ms=15,mew=2)
        pl.xlim(startm,endm)
        pl.ylim(startv,endv)
        pl.xticks([])
        pl.yticks([])
    
    def plot_MW_top(self):
        im = Image.open(config.get_datafile('images','topmilkyway.jpg'))
        pl.imshow(im,origin='lower')
        pl.box(on=False)
        pl.xticks([])
        pl.yticks([])
        xlims = pl.xlim()
        ylims = pl.ylim()
        gal = self.info['galpos']
        pl.plot(2800,1720,'ro',ms=10)
        pl.plot([2800,-5000*np.sin(gal[0]/180.*np.pi)+2800],[1720,5000*np.cos(gal[0]/180.*np.pi)+1720],'r-',lw=2)
        
        #-- necessary information
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        
        x = d/10.*1.3
        y = dprob*1000.
        theta = gal[0]/180.*np.pi + np.pi/2.
        x_ = np.cos(theta)*x - np.sin(theta)*y + 2800
        y_ = np.sin(theta)*x + np.cos(theta)*y + 1720
    
        pl.plot(x_,y_,'r-',lw=2)
        index = np.argmax(y)
        pl.plot(np.cos(theta)*x[index] + 2800,np.sin(theta)*x[index] + 1720,'rx',ms=15,mew=2)
    
        pl.xlim(xlims)
        pl.ylim(ylims)
    
    def plot_finderchart(self):
        
        try:
            data,coords,size = mast.get_dss_image(self.ID)
            pl.imshow(data[::-1],extent=[-size[0]/2*60,size[0]/2*60,
                                        -size[1]/2*60,size[1]/2*60],cmap=pl.cm.RdGy_r)#Greys
            xlims,ylims = pl.xlim(),pl.ylim()
            keep_this = -self.master['color'] & (self.master['cmeas']>0)
            toplot = self.master[keep_this]
            systems = np.array([system.split('.')[0] for system in toplot['photband']],str)
            set_systems = sorted(list(set(systems)))
            pl.gca().set_color_cycle([pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))])    
            for system in set_systems:
                keep = systems==system
                if sum(keep):
                    pl.plot(toplot['_RAJ2000'][keep][0]*60,
                            toplot['_DEJ2000'][keep][0]*60,'x',label=system,
                            mew=2.5,ms=15,alpha=0.5)
                else:
                    pl.gca()._get_lines.count += 1
            pl.legend(numpoints=1,prop=dict(size='x-small'))
            pl.xlim(xlims)
            pl.ylim(ylims)
            pl.xlabel(r'Right ascension $\alpha$ [arcmin]')
            pl.ylabel(r'Declination $\delta$ [arcmin]')    
        except:
            pass
        
        if 'pm' in self.info:
            ppm_ra,ppm_de = (self.info['pm']['pmRA'],self.info['pm']['epmRA']),(self.info['pm']['pmDE'],self.info['pm']['epmDE'])
            pl.annotate('',xy=(ppm_ra[0]/50.,ppm_de[0]/50.),
                  xycoords='data',xytext=(0,0),textcoords='data',color='red',
                  arrowprops=dict(facecolor='red', shrink=0.05),
                  horizontalalignment='right', verticalalignment='top')
            pl.annotate('pmRA: %.1f $\pm$ %.1f mas/yr\npmDE: %.1f $\pm$ %.1f mas/yr'%(ppm_ra[0],ppm_ra[1],ppm_de[0],ppm_de[1]),
                        xy=(0.05,0.25),xycoords='axes fraction',color='red')
            (d,dprob) = self.results['igrid_search']['d']
            max_distance = d[np.argmax(dprob)]
            e_max_distance = abs(max_distance - d[np.argmin(np.abs(dprob-0.5*max(dprob)))])
            
            tang_velo = 'Tan. vel. at %.0f+/-%.0f pc: '%(max_distance,e_max_distance)
            
            max_distance = conversions.convert('pc','km',max_distance,e_max_distance)
            ppm_ra = conversions.convert('mas/yr','rad/s',*ppm_ra)
            ppm_de = conversions.convert('mas/yr','rad/s',*ppm_de)
            max_distance = unumpy.uarray([max_distance[0],max_distance[1]])
            x = unumpy.uarray([ppm_ra[0],ppm_ra[1]])
            y = unumpy.uarray([ppm_de[0],ppm_de[1]])
            velocity = max_distance*utan( usqrt(x**2+y**2))
            
            
            pl.annotate(tang_velo + '%s km/s'%(velocity),xy=(0.05,0.2),xycoords='axes fraction',color='red')
            if 'Vel' in self.info and 'v' in self.info['Vel']:
                rad_velo = 'Rad. vel.: %.1f'%(self.info['Vel']['v'])
                if 'e' in self.info['Vel']:
                    rad_velo += '+/-%.1f'%(self.info['Vel']['e'])
                pl.annotate(rad_velo+' km/s',xy=(0.05,0.15),xycoords='axes fraction',color='red')
    
        
        
    def make_plots(self):
        """
        Make all available plots
        """
        pl.figure(figsize=(22,12))
        rows,cols = 3,4
        pl.subplots_adjust(left=0.04, bottom=0.07, right=0.97, top=0.96,
                wspace=0.17, hspace=0.24)
        pl.subplot(rows,cols,1);self.plot_grid(ptype='CI_red')
        pl.subplot(rows,cols,2);self.plot_grid(ptype='ebv')
        pl.subplot(rows,cols,3);self.plot_grid(ptype='z')
        pl.subplot(rows,cols,4);self.plot_distance()
        
        pl.subplot(3,2,3);self.plot_sed(colors=False)
        pl.subplot(3,2,5);self.plot_sed(colors=True)
        
        pl.subplot(rows,cols,7);self.plot_chi2(colors=False)
        pl.subplot(rows,cols,11);self.plot_chi2(colors=True)
        
        pl.subplot(rows,cols,8);self.plot_grid_model(ptype='prob')
        pl.subplot(rows,cols,12);self.plot_grid_model(ptype='radii')
        
        pl.figure(figsize=(12,12))
        pl.axes([0,0.0,1.0,0.5]);self.plot_MW_side()
        pl.axes([0,0.5,0.5,0.5]);self.plot_MW_top()
        pl.axes([0.5,0.5,0.5,0.5]);self.plot_finderchart()
    
    def save_fits(self,filename=None,overwrite=True):
        """
        Save fitting parameters and results to a FITS file.
        """
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.fits'
        if overwrite:
            if os.path.isfile(filename): os.remove(filename)
            
        eff_waves,synflux,photbands = self.results['synflux']
        chi2 = self.results['chi2']
        
        master = mlab.rec_append_fields(self.master, 'synflux',synflux)
        master = mlab.rec_append_fields(master, 'mod_eff_wave',eff_waves)
        master = mlab.rec_append_fields(master, 'chi2',chi2)
        
        results_dict = dict(extname='model')
        keys = sorted(self.results['igrid_search'])
        for key in keys:
            if 'CI' in key:
                for ikey in self.results['igrid_search'][key]:
                    results_dict[ikey] = self.results['igrid_search'][key][ikey]
                    
        fits.write_recarray(master,filename,header_dict=dict(extname='data'))
        fits.write_array(list(self.results['model']),filename,
                                names=('wave','flux','dered_flux'),
                                units=('A','erg/s/cm2/A','erg/s/cm2/A'),
                                header_dict=results_dict)
        
        results_dict['extname'] = 'igrid_search'
        fits.write_recarray(self.results['igrid_search']['grid'],filename,header_dict=results_dict)
    
    def load_fits(self,filename=None):
        """
        Load grid results from a FITS file
        """
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.fits'
        ff = pyfits.open(filename)
        
        #-- grid search results
        self.results['igrid_search'] = {}
        fields = ff[3].columns.names
        master = np.rec.fromarrays([ff[3].data.field(field) for field in fields],names=','.join(fields))
        self.results['igrid_search']['grid'] = master
        
        self.results['model'] = ff[2].data.field('wave'),ff[2].data.field('flux'),ff[2].data.field('dered_flux')
        self.results['chi2'] = ff[1].data.field('chi2')
        self.results['synflux'] = ff[1].data.field('mod_eff_wave'),ff[1].data.field('synflux'),ff[1].data.field('photband')
        
        
        #-- observed photometry
        fields = ff[1].columns.names
        master = np.rec.fromarrays([ff[1].data.field(field) for field in fields],names=','.join(fields))
        self.master = master
        
        ff.close()
    
    #}


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    pl.show()
    sys.exit()
    #-- PCA analysis
    master['include'] = True
    exclude(master,names=['STROMGREN.HBN-HBW','USNOB1'],wrange=(1.5e4,np.inf))
    do_pca = False
    include_pca = master['include']

    if sum(keep)>2:
        do_pca = True
        logger.info("Start of PCA analysis to find fundamental parameters")
        colors,index = np.unique1d(master['photband'][include_pca],return_index=True)
        A,grid = fit.get_PCA_grid(colors,ebvrange=(0,0.5),res=10)
        P,T,(means,stds) = fit.get_PCA(A)
        calib = fit.calibrate_PCA(T,grid,function='linear')
        obsT,e_obsT = master['cmeas'][include_pca][index], master['e_cmeas'][include_pca][index]
        pars = fit.get_PCA_parameters(obsT,calib,P,means,stds,e_obsT=e_obsT,mc=mc)
        teff,logg,ebv = pars[0]
        logger.info("PCA result: teff=%.0f, logg=%.2f, E(B-V)=%.3f"%(teff,logg,ebv))

        #-- find angular diameter
        logger.info('Estimation of angular diameter')
        iflux = model.get_itable(teff=teff,logg=logg,ebv=ebv,photbands=master['photband'][include_grid])
        scale_pca = np.median(master['cmeas'][include_grid]/iflux)
        angdiam = 2*np.arctan(np.sqrt(scale_pca))/np.pi*180*3600*1000
        logger.info('Angular diameter = %.3f mas'%(angdiam))

        #-- get model
        wave_pca,flux_pca = model.get_table(teff=teff,logg=logg,ebv=ebv,law='fitzpatrick1999')
        wave_ur_pca,flux_ur_pca = model.get_table(teff=teff,logg=logg,ebv=0,law='fitzpatrick1999')

    logger.info('...brought to you in %.3fs'%(time.time()-c0))

    pl.figure()
    pl.title(ID)
    toplot = master[-master['color']]
    systems = np.array([system.split('.')[0] for system in toplot['photband']],str)
    set_systems = sorted(list(set(systems)))
    pl.gca().set_color_cycle([pl.cm.spectral(i) for i in np.linspace(0, 1.0, len(set_systems))])
    for system in set_systems:
        keep = systems==system
        pl.errorbar(master['cwave'][keep],master['cmeas'][keep],
                yerr=master['e_cmeas'][keep],fmt='o',label=system)
    pl.gca().set_xscale('log',nonposx='clip')
    pl.gca().set_yscale('log',nonposy='clip')
    pl.gca().set_autoscale_on(False)
    pl.plot(wave_pca,flux_pca*scale_pca,'r-')
    pl.plot(wave_ur_pca,flux_ur_pca*scale_pca,'k-')
    pl.grid()
    pl.legend(loc='lower left',prop=dict(size='x-small'))
    pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\n'%(teff,logg,ebv)+r'$\theta$=%.3f mas'%(angdiam),(0.75,0.80),xycoords='axes fraction')
    pl.xlabel('wavelength [$\mu$m]')
    pl.ylabel(r'$F_\nu$ [Jy]')
    #### -- END SOME FIGURES -- ####


    if mc and do_pca:
        pl.figure()
        pl.subplot(131)
        pl.hexbin(pars[:,0],pars[:,1])
        pl.xlim(pars[:,0].max(),pars[:,0].min())
        pl.ylim(pars[:,1].max(),pars[:,1].min())
        pl.colorbar()
        pl.subplot(132)
        pl.hexbin(pars[:,0],pars[:,2])
        pl.xlim(pars[:,0].max(),pars[:,0].min())
        pl.ylim(pars[:,2].min(),pars[:,2].max())
        pl.colorbar()
        pl.subplot(133)
        pl.hexbin(pars[:,1],pars[:,2])
        pl.xlim(pars[:,1].max(),pars[:,1].min())
        pl.ylim(pars[:,2].min(),pars[:,2].max())
        pl.colorbar()
    pl.show()









    ################################
