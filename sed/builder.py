# -*- coding: utf-8 -*-
"""
SED builder program.

To construct an SED, use the SED class. The functions defined in this module
are mainly convenience functions specifically for that class, but can be used
outside of the SED class if you know what you're doing.

Table of contents:
    
    1. Retrieving and plotting photometry of a target
    2. Where/what is my target?
    3. SED fitting using a grid based approach
        - Saving SED fits
        - Loading SED fits
    4. Accessing the best fitting full SED model
    5. Radii, distances and luminosities
        - Relations between quantities
        - Parallaxes
        - Seismic constraints
        - Reddening constraints
        - Evolutionary constraints

Section 1. Retrieving and plotting photometry of a target
=========================================================

>>> mysed = SED('HD180642')
>>> mysed.get_photometry()
>>> mysed.plot_data()

and call Pylab's C{show} function to show to the screen:

]]include figure]]ivs_sed_builder_example_photometry.png]

Catch IndexErrors and TypeErrors in case no photometry is found.

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
typically slow). If one of the C{gator}, C{vizier} or C{gcpd} is impossibly slow
or the site is down, you can B{include/exclude these sources} via the keywords
C{include} or C{exclude}, which take a list of strings (choose from C{gator},
C{vizier} and/or C{gcpd}). For ViZieR, there is an extra option to change to
another mirror site via

>>> vizier.change_mirror()

The L{vizier.change_mirror} function cycles through all the mirrors continuously,
so sooner or later you will end up with the default one and repeat the cycle.

>>> mysed.get_photometry(exclude=['gator'])

The results will be written to the file B{HD180642.phot}. An example content is::

  #  meas    e_meas flag unit photband          source                        _r   _RAJ2000    _DEJ2000   cwave       cmeas     e_cmeas cunit       color include
  #float64   float64 |S20 |S30 |S30              |S50                     float64    float64     float64 float64     float64     float64 |S50         bool    bool
    7.823      0.02 nan  mag  WISE.W3           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05  123337 4.83862e-17 8.91306e-19 erg/s/cm2/AA     0       1
    7.744     0.179 nan  mag  WISE.W4           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05  222532 4.06562e-18 6.70278e-19 erg/s/cm2/AA     0       1
     7.77     0.023 nan  mag  WISE.W1           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05 33791.9   6.378e-15  1.3511e-16 erg/s/cm2/AA     0       1
    7.803      0.02 nan  mag  WISE.W2           wise_prelim_p3as_psd    0.112931  2.667e-05   2.005e-05   46293 1.82691e-15 3.36529e-17 erg/s/cm2/AA     0       1
    8.505     0.016 nan  mag  TYCHO2.BT         I/259/tyc2                 0.042   7.17e-06    1.15e-06  4204.4 2.76882e-12 4.08029e-14 erg/s/cm2/AA     0       1
    8.296     0.013 nan  mag  TYCHO2.VT         I/259/tyc2                 0.042   7.17e-06    1.15e-06 5321.86 1.93604e-12 2.31811e-14 erg/s/cm2/AA     0       1
     8.27       nan nan  mag  JOHNSON.V         II/168/ubvmeans             0.01    1.7e-07    3.15e-06 5504.67 1.80578e-12 1.80578e-13 erg/s/cm2/AA     0       1
     0.22       nan nan  mag  JOHNSON.B-V       II/168/ubvmeans             0.01    1.7e-07    3.15e-06     nan     1.40749    0.140749 flux_ratio      1       0
    -0.66       nan nan  mag  JOHNSON.U-B       II/168/ubvmeans             0.01    1.7e-07    3.15e-06     nan     1.21491    0.121491 flux_ratio      1       0
     8.49       nan nan  mag  JOHNSON.B         II/168/ubvmeans             0.01    1.7e-07    3.15e-06 4448.06 2.54162e-12 2.54162e-13 erg/s/cm2/AA     0       1
     7.83       nan nan  mag  JOHNSON.U         II/168/ubvmeans             0.01    1.7e-07    3.15e-06 3641.75 3.08783e-12 3.08783e-13 erg/s/cm2/AA     0       1
    2.601       nan nan  mag  STROMGREN.HBN-HBW J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan     1.66181    0.166181 flux_ratio      1       0
   -0.043       nan nan  mag  STROMGREN.M1      J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan    0.961281   0.0961281 flux_ratio      1       0
    8.221       nan nan  mag  STROMGREN.Y       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685 5477.32 1.88222e-12 1.88222e-13 erg/s/cm2/AA     0       1
    0.009       nan nan  mag  STROMGREN.C1      J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan     0.93125    0.093125 flux_ratio      1       0
    0.238       nan nan  mag  STROMGREN.B-Y     J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685     nan     1.28058    0.128058 flux_ratio      1       0
    8.459       nan nan  mag  STROMGREN.B       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685  4671.2 2.41033e-12 2.41033e-13 erg/s/cm2/AA     0       1
    8.654       nan nan  mag  STROMGREN.V       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685 4108.07 2.96712e-12 2.96712e-13 erg/s/cm2/AA     0       1
    8.858       nan nan  mag  STROMGREN.U       J/A+A/528/A148/tables       0.54 -5.983e-05 -0.00013685 3462.92 3.40141e-12 3.40141e-13 erg/s/cm2/AA     0       1
     7.82      0.01 nan  mag  JOHNSON.J         J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06 12487.8 2.36496e-13 2.36496e-15 erg/s/cm2/AA     0       1
     7.79      0.01 nan  mag  JOHNSON.K         J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06 21951.2 3.24868e-14 3.24868e-16 erg/s/cm2/AA     0       1
     7.83      0.04 nan  mag  JOHNSON.H         J/PASP/120/1128/catalog     0.02  1.017e-05    3.15e-06 16464.4 8.64659e-14 3.18552e-15 erg/s/cm2/AA     0       1
   8.3451    0.0065 nan  mag  HIPPARCOS.HP      I/239/hip_main             0.036   2.17e-06    1.15e-06 5275.11 1.91003e-12 1.91003e-14 erg/s/cm2/AA     0       1
    8.525     0.011 nan  mag  TYCHO2.BT         I/239/hip_main             0.036   2.17e-06    1.15e-06  4204.4 2.71829e-12   2.754e-14 erg/s/cm2/AA     0       1
    8.309     0.012 nan  mag  TYCHO2.VT         I/239/hip_main             0.036   2.17e-06    1.15e-06 5321.86   1.913e-12 2.11433e-14 erg/s/cm2/AA     0       1
     8.02     0.057 nan  mag  COUSINS.I         II/271A/patch2               0.2 -4.983e-05   2.315e-05 7884.05 7.51152e-13 3.94347e-14 erg/s/cm2/AA     0       1
    8.287     0.056 nan  mag  JOHNSON.V         II/271A/patch2               0.2 -4.983e-05   2.315e-05 5504.67 1.77773e-12 9.16914e-14 erg/s/cm2/AA     0       1
     8.47       nan nan  mag  USNOB1.B1         I/284/out                  0.026   7.17e-06     1.5e-07 4448.06 2.57935e-12 7.73805e-13 erg/s/cm2/AA     0       1
     8.19       nan nan  mag  USNOB1.R1         I/284/out                  0.026   7.17e-06     1.5e-07 6939.52 1.02601e-12 3.07803e-13 erg/s/cm2/AA     0       1
    8.491     0.012 nan  mag  JOHNSON.B         I/280B/ascc                0.036   2.17e-06    2.15e-06 4448.06 2.53928e-12 2.80651e-14 erg/s/cm2/AA     0       1
    8.274     0.013 nan  mag  JOHNSON.V         I/280B/ascc                0.036   2.17e-06    2.15e-06 5504.67 1.79914e-12 2.15419e-14 erg/s/cm2/AA     0       1
    7.816     0.023 nan  mag  2MASS.J           II/246/out                 0.118 -2.783e-05   1.715e-05 12412.1 2.28049e-13 4.83093e-15 erg/s/cm2/AA     0       1
    7.792     0.021 nan  mag  2MASS.KS          II/246/out                 0.118 -2.783e-05   1.715e-05 21909.2 3.26974e-14 6.32423e-16 erg/s/cm2/AA     0       1
    7.825     0.042 nan  mag  2MASS.H           II/246/out                 0.118 -2.783e-05   1.715e-05 16497.1 8.48652e-14 3.28288e-15 erg/s/cm2/AA     0       1
    8.272     0.017 nan  mag  GENEVA.V          GCPD                         nan        nan         nan  5482.6 1.88047e-12 2.94435e-14 erg/s/cm2/AA     0       1
      1.8     0.004 nan  mag  GENEVA.G-B        GCPD                         nan        nan         nan     nan    0.669837  0.00669837 flux_ratio      1       0
    1.384     0.004 nan  mag  GENEVA.V1-B       GCPD                         nan        nan         nan     nan    0.749504  0.00749504 flux_ratio      1       0
     0.85     0.004 nan  mag  GENEVA.B1-B       GCPD                         nan        nan         nan     nan     1.05773   0.0105773 flux_ratio      1       0
    1.517     0.004 nan  mag  GENEVA.B2-B       GCPD                         nan        nan         nan     nan    0.946289  0.00946289 flux_ratio      1       0
    0.668     0.004 nan  mag  GENEVA.V-B        GCPD                         nan        nan         nan     nan    0.726008  0.00726008 flux_ratio      1       0
    0.599     0.004 nan  mag  GENEVA.U-B        GCPD                         nan        nan         nan     nan     1.13913   0.0113913 flux_ratio      1       0
    7.604 0.0174642 nan  mag  GENEVA.B          GCPD                         nan        nan         nan 4200.85 2.59014e-12 4.16629e-14 erg/s/cm2/AA     0       1
    9.404 0.0179165 nan  mag  GENEVA.G          GCPD                         nan        nan         nan 5765.89 1.73497e-12   2.863e-14 erg/s/cm2/AA     0       1
    8.988 0.0179165 nan  mag  GENEVA.V1         GCPD                         nan        nan         nan 5395.63 1.94132e-12 3.20351e-14 erg/s/cm2/AA     0       1
    8.454 0.0179165 nan  mag  GENEVA.B1         GCPD                         nan        nan         nan 4003.78 2.73968e-12 4.52092e-14 erg/s/cm2/AA     0       1
    9.121 0.0179165 nan  mag  GENEVA.B2         GCPD                         nan        nan         nan 4477.56 2.45102e-12  4.0446e-14 erg/s/cm2/AA     0       1
    8.203 0.0179165 nan  mag  GENEVA.U          GCPD                         nan        nan         nan 3421.62 2.95052e-12 4.86885e-14 erg/s/cm2/AA     0       1
     8.27       nan nan  mag  JOHNSON.V         GCPD                         nan        nan         nan 5504.67 1.80578e-12 1.80578e-13 erg/s/cm2/AA     0       1
     0.22       nan nan  mag  JOHNSON.B-V       GCPD                         nan        nan         nan     nan     1.40749    0.140749 flux_ratio      1       0
    -0.66       nan nan  mag  JOHNSON.U-B       GCPD                         nan        nan         nan     nan     1.21491    0.121491 flux_ratio      1       0
     8.49       nan nan  mag  JOHNSON.B         GCPD                         nan        nan         nan 4448.06 2.54162e-12 2.54162e-13 erg/s/cm2/AA     0       1
     7.83       nan nan  mag  JOHNSON.U         GCPD                         nan        nan         nan 3641.75 3.08783e-12 3.08783e-13 erg/s/cm2/AA     0       1
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
delete the .phot file or set C{force=True} when calling L{get_photometry}. 

The content of the .phot file is most easily read using the L{ivs.io.ascii.read2recarray}
function. Be careful, as it contains both absolute fluxes as flux ratios.

>>> data = ascii.read2recarray('HD180642.phot')

Notice that in the C{.phot} files, also a C{comment} column is added. You can 
find translation of some of the flags here (i.e. upper limit, extended source etc..),
or sometimes just additional remarks on variability etc. Not all catalogs have
this feature implemented, so you are still responsible yourself for checking
the quality of the photometry.

The references to each source are given in the C{bibtex} column. Simply call

>>> mysed.save_bibtex()

to convert those bibcodes to a C{.bib} file.

Using L{SED.plot_MW_side} and L{SED.plot_MW_top}, you can make a picture of where
your star is located with respect to the Milky Way and the Sun. With L{SED.plot_finderchart},
you can check the location of your photometry, and also see if proper motions
etc are available.

Section 2. Where/what is my target?
===================================

To give you some visual information on the target, the following plotting
procedure might be of some help.

To check whether the downloaded photometry is really belonging to the target,
instead of some neighbouring star (don't forget to set C{radius} when looking
for photometry!), you can generate a finderchart with the location of the
downloaded photometry overplotted. On top of that, proper motion data is added
when available, as well as radial velocity data. When a distance is available,
the proper motion velocity will be converted to a true tangential velocity.

>>> p = pl.figure();mysed.plot_finderchart(window_size=1)

]]include figure]]ivs_sed_builder_finderchart.png]

To know the location of your target wrt the Milky Way (assuming your target is
in the milky way), you can call

>>> p = pl.figure();mysed.plot_MW_side()
>>> p = pl.figure();mysed.plot_MW_top()

]]include figure]]ivs_sed_builder_MWside.png]

]]include figure]]ivs_sed_builder_MWtop.png]

Section 3. SED fitting using a grid based approach
==================================================

Subsection 3.1 Single star
--------------------------

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

Currently, four standard schemes are implemented, which you can set via L{SED.set_photometry_scheme}:

    1. C{absolute}: use only absolute values
    2. C{colors}: use only colors (no angular diameter values calculated)
    3. C{combo}: use all colors and one absolute value per photometric system
    4. C{irfm}: (infrared flux method) use colors for wavelengths shorter than
    infrared wavelengths, and absolute values for systems in the infrared. The
    infrared is defined as wavelength longer than 1 micron, but this can be
    customized with the keyword C{infrared=(value,unit)} in
    L{SED.set_photometry_scheme}.

Here, we chose to use colors and one absolute flux per system, but exclude IR
photometry (wavelength range above 2.5 micron), and some systems and colors which
we know are not so trustworthy:

>>> mysed.set_photometry_scheme('combo')
>>> mysed.exclude(names=['STROMGREN.HBN-HBW','USNOB1','SDSS','DENIS','COUSINS','ANS','TD1'],wrange=(2.5e4,1e10))

You can L{include}/L{exclude} photoemtry based on name, wavelength range, source and index,
and only select absolute photometry or colors (L{include_abs},L{include_colors}).
When working in interactive mode, in particular the index is useful. Print the
current set of photometry to the screen with 

>>> print(photometry2str(mysed.master,color=True,index=True))

and you will see in green the included photometry, and in red the excluded photometry.
You will see that each column is preceded by an index, you can use these indices
to select/deselect the photometry.

Speed up the fitting process by copying the model grids to the scratch disk

>>> model.copy2scratch(z='*')

Start the grid based fitting process and show some plots. We use 100000 randomly
distributed points over the grid:

>>> mysed.igrid_search(points=100000)

Delete the model grids from the scratch disk

>>> model.clean_scratch()

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
>>> p = pl.subplot(131);mysed.plot_sed(plot_deredded=True)
>>> p = pl.subplot(132);mysed.plot_grid(limit=None)
>>> p = pl.subplot(133);mysed.plot_grid(x='ebv',y='z',limit=None)

]]include figure]]ivs_sed_builder_example_fitting02.png]

You can automatically make plots of (most plotting functions take C{colors=True/False}
as an argument so you can make e.g. the 'color' SED and 'absolute value' SED):

    1. the grid (see L{SED.plot_grid})
    2. the SED (see L{SED.plot_sed}) 
    3. the fit statistics (see L{SED.plot_chi2})
    
To change the grid, load the L{ivs.sed.model} module and call
L{ivs.sed.model.set_defaults} with appropriate arguments. See that module for
conventions on the grid structure when adding custom grids.

To add arrays manually, i.e. not from the predefined set of internet catalogs,
use the L{SED.add_photometry_fromarrays} function.

B{Warning}: Be careful when interpreting the Chi2 results. In order to always have a
solution, the chi2 is rescaled so that the minimum equals 1, in the case the
probability of the best chi2-model is zero. The Chi2 rescaling factor I{f} mimicks
a rescaling of all errorbars with a factor I{sqrt(f)}, and does not discriminate
between systems (i.e., B{all} errors are blown up). If the errorbars are
underestimated, it could be that the rescaling factor is also wrong, which means
that the true probability region can be larger or smaller!

Subsection 3.2 Binary star 
--------------------------

The SED class can create SEDs for multiple stars as well. There are 2 options 
available, the multiple SED fit which in theory can handle any number of stars,
and the binary SED fit which is for binaries only, and uses the mass of both 
components to restrict the radii when combining two model SEDs.

As an example we take the system PG1104+243, which consists of a subdwarf B star,
and a G2 type mainsequence star. The photometry from the standard catalogues that
are build in this class, is of to low quality, so we use photometry obtained from 
U{the subdwarf database<http://catserver.ing.iac.es/sddb/>}.

>>> mysed = SED('PG1104+243')

>>> meas, e_meas, units, photbands, source = ascii.read2array('pg1104+243_sddb.phot', dtype=str)
>>> meas = np.array(meas, dtype=float)
>>> e_meas = np.array(e_meas, dtype=float)
>>> mysed.add_photometry_fromarrays(meas, e_meas, units, photbands, source)

We use only the absolute fluxes

>>> mysed.set_photometry_scheme('abs')

For the main sequence component we use kurucz models with solar metalicity, and
for the sdB component tmap models. And we copy the model grids to the scratch disk
to speed up the process:

>>> grid1 = dict(grid='kurucz',z=+0.0)
>>> grid2 = dict(grid='tmap')
>>> model.set_defaults_multiple(grid1,grid2)
>>> model.copy2scratch()

The actual fitting. The second fit starts from the 95% probability intervals of
the first fit.

>>> teff_ms = (5000,7000)
>>> teff_sdb = (25000,45000)
>>> logg_ms = (4.00,4.50)
>>> logg_sdb = (5.00,6.50)
>>> mysed.igrid_search(masses=(0.47,0.71) ,teffrange=(teff_ms,teff_fix),loggrange=(logg_ms,logg_sdb), ebvrange=(0.00,0.02), zrange=(0,0), points=2000000, type='binary')
>>> mysed.igrid_search(masses=(0.47,0.71) ,points=2000000, type='binary')

Delete the used models from the scratch disk

>>> model.clean_scratch()

Plot the results

>>> p = pl.figure()
>>> p = pl.subplot(131); mysed.plot_sed(plot_deredded=False)
>>> p = pl.subplot(132); mysed.plot_grid(x='teff', y='logg', limit=0.95)
>>> p = pl.subplot(133); mysed.plot_grid(x='teff-2', y='logg-2', limit=0.95)

]]include figure]]ivs_sed_builder_example_fittingPG1104+243.png]

Subsection 3.3 Saving SED fits
------------------------------

You can save all the data to a multi-extension FITS file via

>>> mysed.save_fits()

This FITS file then contains all B{measurements} (it includes the .phot file),
the B{resulting SED}, the B{confidence intervals of all parameters} and also the
B{whole fitted grid}: in the above case, the extensions of the FITS file contain
the following information (we print part of each header)::

    EXTNAME = 'DATA    '
    XTENSION= 'BINTABLE'           / binary table extension                         
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                  270 / length of dimension 1                          
    NAXIS2  =                   67 / length of dimension 2                          
    TTYPE1  = 'meas    '                                                                                                                   
    TTYPE2  = 'e_meas  '                                                            
    TTYPE3  = 'flag    '                                                            
    TTYPE4  = 'unit    '                                                            
    TTYPE5  = 'photband'                                                            
    TTYPE6  = 'source  '                                                            
    TTYPE7  = '_r      '                                                            
    TTYPE8  = '_RAJ2000'                                                            
    TTYPE9  = '_DEJ2000'                                                            
    TTYPE10 = 'cwave   '                                                            
    TTYPE11 = 'cmeas   '                                                            
    TTYPE12 = 'e_cmeas '                                                            
    TTYPE13 = 'cunit   '                                                            
    TTYPE14 = 'color   '                                                            
    TTYPE15 = 'include '                                                            
    TTYPE16 = 'synflux '                                                            
    TTYPE17 = 'mod_eff_wave'                                                        
    TTYPE18 = 'chi2    '                                                            

    EXTNAME = 'MODEL   '      
    XTENSION= 'BINTABLE'           / binary table extension                         
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   24 / length of dimension 1                          
    NAXIS2  =                 1221 / length of dimension 2                                                                                
    TTYPE1  = 'wave    '                                                            
    TUNIT1  = 'A       '                                                            
    TTYPE2  = 'flux    '                                                            
    TUNIT2  = 'erg/s/cm2/AA'                                                         
    TTYPE3  = 'dered_flux'                                                          
    TUNIT3  = 'erg/s/cm2/AA'                                                         
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
    TTYPE1  = 'teff    '                                                            
    TTYPE2  = 'logg    '                                                            
    TTYPE3  = 'ebv     '                                                            
    TTYPE4  = 'z       '                                                            
    TTYPE5  = 'chisq   '                                                            
    TTYPE6  = 'scale   '                                                            
    TTYPE7  = 'escale '                                                            
    TTYPE8  = 'Labs    '                                                            
    TTYPE9  = 'CI_raw  '                                                            
    TTYPE10 = 'CI_red  '                                                            


Subsection 3.4 Loading SED fits
-------------------------------

Once saved, you can load the contents of the FITS file again into an SED object
via

>>> mysed = SED('HD180642')
>>> mysed.load_fits()

and then you can build all the plots again easily. You can of course use the
predefined plotting scripts to start a plot, and then later on change the
properties of the labels, legends etc... for higher quality plots or to better
suit your needs.

Section 4. Accessing the best fitting full SED model
====================================================

You can access the full SED model that matches the parameters found by the
fitting routine via:

>>> wavelength,flux,deredded_flux = mysed.get_best_model()

Note that this model is retrieved after fitting, and was not in any way used
during the fitting. As a consequence, there could be small differences between
synthetic photometry calculated from this returned model and the synthetic
fluxes stored in C{mysed.results['igrid_search']['synflux'], which is the
synthetic photometry coming from the interpolation of the grid of pre-interpolated
photometry. See the documentation of L{SED.get_model} for more information.

Section 5. Radii, distances and luminosities
============================================

Subsection 4.1. Relations between quantities
--------------------------------------------

Most SED grids don't have the radius as a tunable model parameter. The scale
factor, which is available for all fitted models in the grid when at least one
absolute photometric point is included, is directly propertional to the angular
diameter. The following relations hold::

    >> distance = radius / np.sqrt(scale)
    >> radius = distance * np.sqrt(scale)

Where C{radius} and C{distance} have equal units. The true absolute luminosity
(solar units) is related to the absolute luminosity from the SED (solar units)
via the radius (solar units)::
    
    >> L_abs_true = L_abs * radius**2

Finally, the angular diameter can be computed via::

    >> 2*conversions.convert('sr','mas',scale)

Subsection 4.2. Seismic constraints
-----------------------------------

If the star shows clear solar-like oscillations, you can use the nu_max give an
independent constraint on the surface gravity, given the effective temperature of
the model (you are free to give errors or not, just keep in mind that in the
former case, you will get an array of Uncertainties rather than floats)::
    
    >> teff = 4260.,'K'
    >> nu_max = 38.90,0.86,'muHz'
    >> logg_slo = conversions.derive_logg_slo(teff,nu_max,unit='[cm/s2'])

If you have the large separation (l=0 modes) and the nu_max, you can get an
estimate of the radius::

    >> Deltanu0 = 4.80,0.02,'muHz'
    >> R_slo = conversions.derive_radius_slo(numax,Deltanu0,teff,unit='Rsol')
    
Then from the radius, you can get both the absolute luminosity and distance to
your star.

Subsection 4.3. Parallaxes
--------------------------

From the parallax, you can get an estimate of the distance. This is, however,
dependent on some prior assumptions such as the shape of the galaxy and the
distribution of stars. To estimate the probability density value due to
a measured parallax of a star at particular distance, you can call L{distance.distprob}::
    
    >> gal_latitude = 0.5
    >> plx = 3.14,0.5
    >> d_prob = distance.distprob(d,gal_latitude,plx)

Using this distance, you can get an estimate for the radius of your star, and
thus also the absolute luminosity.

Subsection 4.4: Reddening constraints
-------------------------------------

An estimate of the reddening of a star at a particular distance may be obtained
with the L{extinctionmodels.findext} function. There are currently three
reddening maps available: Drimmel, Marshall and Arenou::

    >> lng,lat = 120.,-0.5
    >> dist = 100. # pc
    >> Rv = 3.1
    >> EBV = extinctionmodels.findext(lng,lat,model='drimmel',distance=dist)/Rv
    >> EBV = extinctionmodels.findext(lng,lat,model='marshall',distance=dist)/Rv
    >> EBV = extinctionmodels.findext(lng,lat,model='arenou',distance=dist)/Rv


"""
import sys
import time
import os
import logging
import itertools
import glob
import json

import pylab as pl
from matplotlib import mlab
import Image
import numpy as np
import scipy.stats
from scipy.interpolate import Rbf
import pyfits

from ivs import config
from ivs.aux import numpy_ext
from ivs.aux import termtools
from ivs.aux.decorators import memoized,clear_memoization
from ivs.io import ascii
from ivs.io import fits
from ivs.sed import model
from ivs.sed import filters
from ivs.sed import fit
from ivs.sed import distance
from ivs.sed import extinctionmodels
from ivs.sed.decorators import standalone_figure
from ivs.spectra import tools
from ivs.catalogs import crossmatch
from ivs.catalogs import vizier
from ivs.catalogs import mast
from ivs.catalogs import sesame
from ivs.catalogs import corot
from ivs.units import conversions
from ivs.units import constants
from ivs.units.uncertainties import unumpy,ufloat
from ivs.units.uncertainties.unumpy import sqrt as usqrt
from ivs.units.uncertainties.unumpy import tan as utan
from ivs.sigproc import evaluate
try:
    from ivs.stellar_evolution import evolutionmodels
except ImportError:
    print("Warning: no evolution models available (probably not important)")

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
    cats = np.array([ph.split('.')[0] for ph in master['source']])
    set_cats = sorted(set(cats))
    #-- add common uncatalogized colors:
    columns = list(master.dtype.names)
    #-- if the separate bands are available (and not the color itself),
    #   calculate the colors here. We need to do this for every separate
    #   system!
    add_rows = []
    for cat in set_cats:
        master_ = master[cats==cat]
    
        for color in ['2MASS.J-H','2MASS.KS-H','TD1.1565-1965','TD1.2365-1965','TD1.2365-2740',
                  'JOHNSON.J-H','JOHNSON.K-H','JOHNSON.B-V','JOHNSON.U-B','JOHNSON.R-V','JOHNSON.I-V',
                  'GALEX.FUV-NUV','TYCHO2.BT-VT','WISE.W1-W2','WISE.W3-W2','WISE.W4-W3',
                  'SDSS.U-G','SDSS.G-R','SDSS.R-I','SDSS.R-Z',
                  'UVEX.U-G','UVEX.G-R','UVEX.R-I','UVEX.HA-I',
                  'DENIS.I-J','DENIS.KS-J',
                  'ANS.15N-15W','ANS.15W-18','ANS.18-22','ANS.22-25','ANS.25-33']:
            #-- get the filter system and the separate bands for these colors
            system,band = color.split('.')
            band0,band1 = band.split('-')
            band0,band1 = '%s.%s'%(system,band0),'%s.%s'%(system,band1)
        
            if band0 in master_['photband'] and band1 in master_['photband'] and not color in master_['photband']:
                #-- where are the bands located?
                index0 = list(master_['photband']).index(band0)
                index1 = list(master_['photband']).index(band1)
                
                #-- start a new row to add the color
                row = list(master_[index1])
                row[columns.index('photband')] = color
                row[columns.index('cwave')] = np.nan
                
                #-- it could be a magnitude difference
                if master_['unit'][index0]=='mag':
                    row[columns.index('meas')] = master_['meas'][index0]-master_['meas'][index1]
                    row[columns.index('e_meas')] = np.sqrt(master_['e_meas'][index0]**2+master_['e_meas'][index1]**2)
                    #-- error will not always be available...
                    try:
                        row[columns.index('cmeas')],row[columns.index('e_cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],row[columns.index('e_meas')],photband=color)
                    except AssertionError:
                        row[columns.index('cmeas')] = conversions.convert('mag_color','flux_ratio',row[columns.index('meas')],photband=color)
                        row[columns.index('e_cmeas')] = np.nan
                #-- or it could be a flux ratio
                else:
                    row[columns.index('meas')] = master_['meas'][index0]/master_['meas'][index1]
                    row[columns.index('e_meas')] = np.sqrt(((master_['e_meas'][index0]/master_['meas'][index0])**2+(master_['e_meas'][index1]/master_['meas'][index1])**2)*row[columns.index('meas')]**2)
                    row[columns.index('cmeas')],row[columns.index('e_cmeas')] = row[columns.index('meas')],row[columns.index('e_meas')]
                row[columns.index('cunit')] = 'flux_ratio'
                add_rows.append(tuple(row))
    master = numpy_ext.recarr_addrows(master,add_rows)
                
    #-- add an extra column with a flag to distinguish colors from absolute
    #   fluxes, and a column with flags to include/exclude photometry
    #   By default, exclude photometry if the effective wavelength is above
    #   150 mum or if the measurement is nan
    if not 'color' in master.dtype.names:
        extra_cols = [[filters.is_color(photband) for photband in master['photband']],
                    [(cwave<150e4) or np.isnan(meas) for cwave,meas in zip(master['cwave'],master['cmeas'])]]
        dtypes = [('color',np.bool),('include',np.bool)]
        master = numpy_ext.recarr_addcols(master,extra_cols,dtypes)
    else:
        iscolor = [filters.is_color(photband) for photband in master['photband']]
        master['color'] = iscolor
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
    
    #-- remove negative fluxes
    master = master[master['cmeas']>0]
    
    return master


def decide_phot(master,names=None,wrange=None,sources=None,indices=None,ptype='all',include=False):
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
    
        1. Exclude all measurements::
        
        >> decide_phot(master,wrange=(-np.inf,+np.inf),ptype='all',include=False)
    
        2. Include all TD1 fluxes and colors::
        
        >> decide_phot(master,names=['TD1'],ptype='all',include=True)
    
        3. Include all V band measurements from all systems (but not the colors)::
        
        >> decide_phot(master,names=['.V'],ptype='abs',include=True)
    
        4. Include all Geneva colors and exclude Geneva magnitudes::
        
        >> decide_phot(master,names=['GENEVA'],ptype='col',include=True)
        >> decide_phot(master,names=['GENEVA'],ptype='abs',include=False)
    
        5. Exclude all infrared measurements beyond 1 micron::
        
        >> decide_phot(master,wrange=(1e4,np.inf),ptype='all',include=False)
    
        6. Include all AKARI measurements below 10 micron::
        
        >> decide_phot(master,names=['AKARI'],wrange=(-np.inf,1e5),ptype='all',include=True)
    
    @param master: record array containing all photometry
    @type master: numpy record array
    @param names: strings excerpts to match filters
    @type names: list of strings
    @param wrange: wavelength range (most likely angstrom) to include/exclude
    @type wrange: 2-tuple (start wavelength,end wavelength)
    @param sources: list of sources
    @type sources: list of strings
    @param indices: list of indices (integers)
    @type indices: list of integers
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
        logger.info('%s photometry based on photband containining one of %s'%((include and 'Include' or "Exclude"),names))
        for index,photband in enumerate(master['photband']):
            for name in names:
                if name in photband:
                    if ptype=='all' or (ptype=='abs' and -master['color'][index]) or (ptype=='col' and master['color'][index]):
                        master['include'][index] = include
                        break
    #-- exclude/include colors based on their wavelength
    if wrange is not None:
        logger.info('%s photometry based on photband wavelength between %s'%((include and 'Include' or "Exclude"),wrange))
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
    #-- exclude/include passbands based on their source
    if sources is not None:
        logger.info('%s photometry based on source catalog in %s'%((include and 'Include' or "Exclude"),sources))
        for index,msource in enumerate(master['source']):
            for source in sources:
                if source==msource:
                    if ptype=='all' or (ptype=='abs' and -master['color'][index]) or (ptype=='col' and master['color'][index]):
                        master['include'][index] = include
                        break
    #-- exclude/include passbands based on their index
    if indices is not None:
        logger.info('%s photometry based on index'%((include and 'Include' or "Exclude")))
        indices = np.array(indices,int)
        if not indices.shape: indices = indices.reshape(1)
        master['include'][indices] = include
    

def photometry2str(master,comment='',sort='photband',color=False,index=False):
    """
    String representation of master record array.
    
    Sorting is disabled when C{index=True}.
    
    @param master: master record array containing photometry
    @type master: numpy record array
    """
    if sort and not index:
        master = master[np.argsort(master[sort])]
    
    templateh = '{:20s} {:>12s} {:>12s} {:12s} {:>10s} {:>12s} {:>12s} {:12s} {:30s}'
    templated = '{:20s} {:12g} {:12g} {:12s} {:10.0f} {:12g} {:12g} {:12s} {:30s}'
    header = ['PHOTBAND','MEAS','E_MEAS','UNIT','CWAVE','CMEAS','E_CMEAS','CUNIT','SOURCE']
    if 'comments' in master.dtype.names:
        templateh += ' {:s}'
        templated += ' {:s}'
        header += ['COMMENTS']
    if index:
        templateh = '{:3s} '+templateh
        templated = '{:3d} '+templated
        header = ['NR']+header
        
    txt = [comment+templateh.format(*header)]
    txt.append(comment+'='*170)
    columns = [master[col.lower()] for col in header if not col=='NR']
    for nr,contents in enumerate(zip(*columns)):
        contents = list(contents)
        if 'comments' in master.dtype.names:
            contents[-1] = contents[-1].replace('_',' ')
        if index:
            contents = [nr] + contents
        line = templated.format(*contents)
        if color:
            mycolor = termtools.green if master['include'][nr] else termtools.red
            line = mycolor(line)
        txt.append(comment + line)
    return "\n".join(txt)

#@memoized
#def get_schaller_grid():
    #"""
    #Download Schaller 1992 evolutionary tracks and return an Rbf interpolation
    #function.
    
    #@return: Rbf interpolation function
    #@rtype: Rbf interpolation function
    #"""
    ##-- translation between table names and masses
    ##masses = [1,1.25,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,40,60][:-1]
    ##tables = ['table20','table18','table17','table16','table15','table14',
    ##          'table13','table12','table11','table10','table9','table8',
    ##          'table7','table6','table5','table4','table3'][:-1]
    ##-- read in all the tables and compute radii and luminosities.
    #data,comms,units = vizier.search('J/A+AS/96/269/models')
    #all_teffs = 10**data['logTe']
    #all_radii = np.sqrt((10**data['logL']*constants.Lsol_cgs)/(10**data['logTe'])**4/(4*np.pi*constants.sigma_cgs))
    #all_loggs = np.log10(constants.GG_cgs*data['Mass']*constants.Msol_cgs/(all_radii**2))
    #all_radii /= constants.Rsol_cgs
    ##-- remove low temperature models, the evolutionary tracks are hard to
    ##   interpolate there.
    #keep = all_teffs>5000
    #all_teffs = all_teffs[keep]
    #all_radii = all_radii[keep]
    #all_loggs = all_loggs[keep]
    ##-- make linear interpolation model between all modelpoints
    #mygrid = Rbf(np.log10(all_teffs),all_loggs,all_radii,function='linear')
    #logger.info('Interpolation of Schaller 1992 evolutionary tracks to compute radii')
    #return mygrid

    
    
    
    


#def get_radii(teffs,loggs):
    #"""
    #Retrieve radii from stellar evolutionary tracks from Schaller 1992.
    
    #@param teffs: model effective temperatures
    #@type teffs: numpy array
    #@param loggs: model surface gravities
    #@type loggs: numpy array
    #@return: model radii (solar units)
    #@rtype: numpy array
    #"""
    #mygrid = get_schaller_grid()
    #radii = mygrid(np.log10(teffs),loggs)
    #return radii


#def calculate_distance(plx,gal,teffs,loggs,scales,n=75000):
    #"""
    #Calculate distances and radii of a target given its parallax and location
    #in the galaxy.
    
    
    #"""
    ##-- compute distance up to 25 kpc, which is about the maximum distance from
    ##   earth to the farthest side of the Milky Way galaxy
    ##   rescale to set the maximum to 1
    #d = np.logspace(np.log10(0.1),np.log10(25000),100000)
    #if plx is not None:
        #dprob = distance.distprob(d,gal[1],plx)
        #dprob = dprob / dprob.max()
    #else:
        #dprob = np.ones_like(d)
    ##-- compute the radii for the computed models, and convert to parsec
    ##radii = np.ones(len(teffs[-n:]))
    #radii = get_radii(teffs[-n:],loggs[-n:])
    #radii = conversions.convert('Rsol','pc',radii)
    #d_models = radii/np.sqrt(scales[-n:])
    ##-- we set out of boundary values to zero
    #if plx is not None:
        #dprob_models = np.interp(d_models,d[-n:],dprob[-n:],left=0,right=0)
    #else:
        #dprob_models = np.ones_like(d_models)
    ##-- reset the radii to solar units for return value
    #radii = conversions.convert('pc','Rsol',radii)
    #return (d_models,dprob_models,radii),(d,dprob)



class SED(object):
    """
    Class that facilitates the use of the ivs.sed module.
    
    This class is meant to be an easy interface to many of the ivs.sed module's
    functionality.
    
    The most important attributes of SED are:
    
        1. C{sed.ID}: star's identification (str)
        2. C{sed.photfile}: name of the file containing all photometry (str)
        3. C{sed.info}: star's information from Simbad (dict)
        4. C{sed.master}: photometry data (record array)
        5. C{sed.results}: results and summary of the fitting process (dict)
    
    After fitting, e.g. via calling L{igrid_search}, you can call L{get_model}
    to retrieve the full SED matching the best fitting parameters (or, rather,
    closely matching them, see the documentation).
        
    """
    def __init__(self,ID=None,photfile=None,plx=None,load_fits=True,label=''):
        """
        Initialize SED class.
        
        Different ways to initialize:
        
        B{1. If no previous data are saved:}
        
        >>> mysed = SED('HD129929')
        
        B{2. If previous data exists:}
        
        >>> #mysed = SED(photfile='HD129929.phot') # will set ID with 'oname' field from SIMBAD
        >>> #mysed = SED(ID='bla', photfile='HD129929.phot') # Sets custom ID
        
        The C{ID} variable is used internally to look up data, so it should be
        something SIMBAD understands and that designates the target.
        
        @param plx: parallax (and error) of the object
        @type plx: tuple (plx,e_plx)
        """
        self.ID = ID
        self.label = label
        self.info = {}
        #-- the file containing photometry should have the following name. We
        #   save photometry to a file to avoid having to download photometry
        #   each time from the internet
        if photfile is None:
            self.photfile = '%s.phot'%(ID).replace(' ','_')
            #-- keep information on the star from SESAME, but override parallax with
            #   the value from Van Leeuwen's new reduction. Set the galactic
            #   coordinates, which are not available from SESAME.
        else:
            self.photfile = photfile
        
        #-- load information from the photometry file if it exists
        if not os.path.isfile(self.photfile):
            try:
                self.info = sesame.search(os.path.basename(ID),fix=True)
            except KeyError:
                logger.warning('Star %s not recognised by SIMBAD'%(os.path.basename(ID)))
                try:
                    self.info = sesame.search(os.path.basename(ID),db='N',fix=True)
                    logger.info('Star %s recognised by NED'%(os.path.basename(ID)))
                except KeyError:
                    logger.warning('Star %s not recognised by NED'%(os.path.basename(ID)))
            #-- final attempt: if it's a CoRoT star, try the catalog
            if 'corot' in ID.lower() and not 'galpos' in self.info:
                corot_id = int("".join([char for char in ID if char.isdigit()]))
                try:
                    self.info['jradeg'],self.info['jdedeg'] = ra,dec = corot.resolve(corot_id)
                    gal = conversions.convert('equatorial','galactic',(ra,dec),epoch='2000')
                    self.info['galpos'] = float(gal[0])/np.pi*180,float(gal[1])/np.pi*180    
                    logger.info('Resolved position via CoRoT EXO catalog')
                except:
                    logger.warning("Star %s not recognised by CoRoT"%(os.path.basename(ID)))
                    
            if plx is not None:
                if not 'plx' in self.info:
                    self.info['plx'] = {}
                self.info['plx']['v'] = plx[0]
                self.info['plx']['e'] = plx[1]
        else:
            self.load_photometry()
            #-- if no ID was given, set the official name as the ID.
            if self.ID is None:
                self.ID = os.path.splitext(os.path.basename(self.photfile))[0]#self.info['oname']
        #--load information from the FITS file if it exists
        self.results = {}
        if load_fits:
            self.load_fits()
            
        #-- prepare for information on fitting processes
        self.CI_limit = 0.95
        
    def __repr__(self):
        """
        Machine readable string representation of an SED object.
        """
        return "SED('{}')".format(self.ID)
    
    def __str__(self):
        """
        Human readable string representation of an SED object.
        """
        txt = []
        #-- object designation
        txt.append("Object identification: {:s}".format(self.ID))
        if hasattr(self,'info') and self.info and 'oname' in self.info:
            txt.append("Official designation: {:s}".format(self.info['oname']))
        #-- additional info
        for key in sorted(self.info.keys()):
            if isinstance(self.info[key],dict):
                txt.append(" {:10s} = ".format(key)+", ".join(["{}: {}".format(i,j) for i,j in self.info[key].iteritems()]))
            else:
                txt.append(" {:10s} = {}".format(key,self.info[key]))

        #txt.append("Included photometry:")
        #if hasattr(self,'master') and self.master is not None:
            #include_grid = self.master['include']
            #txt.append(photometry2str(self.master[include_grid]))
        #txt.append("Excluded photometry:")
        #if hasattr(self,'master') and self.master is not None:
            #include_grid = self.master['include']
            #txt.append(photometry2str(self.master[-include_grid]))
        if hasattr(self,'master'):
            txt.append(photometry2str(self.master,color=True))
        return "\n".join(txt)
        
    #{ Handling photometric data
    def get_photometry(self,radius=None,ra=None,dec=None,
                       include=None,exclude=None,force=False,
                       units='erg/s/cm2/AA'):
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
        if not os.path.isfile(self.photfile) or force:
            #-- get and fix photometry. Set default errors to 1%, and set
            #   USNOB1 errors to 3%
            if ra is None and dec is None:
                master = crossmatch.get_photometry(ID=os.path.basename(self.ID),radius=radius,
                                       include=include,exclude=exclude,to_units=units,
                                       extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            else:
                master = crossmatch.get_photometry(ID=os.path.basename(self.ID),ra=ra,dec=dec,radius=radius,
                                       include=include,exclude=exclude,to_units=units,
                                       extra_fields=['_r','_RAJ2000','_DEJ2000']) # was radius=3.
            if 'jradeg' in self.info:
                master['_RAJ2000'] -= self.info['jradeg']
                master['_DEJ2000'] -= self.info['jdedeg']
            
            #-- fix the photometry: set default errors to 2% and print it to the
            #   screen
            self.master = fix_master(master,e_default=0.1)
            logger.info('\n'+photometry2str(master))
            
            #-- write to file
            self.save_photometry()
    
    def get_spectrophotometry(self,directory=None,force_download=False):
        """
        Retrieve and combine spectra.
        
        B{WARNING:} this function creates FUSE and DIR directories!
        """
        if directory is None and os.path.dirname(self.ID) == '':
            directory = os.getcwd()
        elif os.path.dirname(self.ID) != '':
            directory = os.path.dirname(self.ID)
            
        if not os.path.isdir(directory):
            os.mkdir(directory)
            
        #-- add spectrophotometric filters to the set
        photbands = filters.add_spectrophotometric_filters(R=200,lambda0=950,lambdan=3350)
        if hasattr(self,'master') and self.master is not None and not force_download:
            if any(['BOXCAR' in photband for photband in self.master['photband']]):
                return None
        #-- FUSE spectra
        fuse_direc = os.path.join(directory,'FUSE')
        iue_direc = os.path.join(directory,'IUE')
        if not os.path.isdir(fuse_direc) or force_download:
            out1 = mast.get_FUSE_spectra(ID=os.path.basename(self.ID),directory=fuse_direc,select=['ano'])
            if out1 is None:
                out1 = []
        else:
            out1 = glob.glob(fuse_direc+'/*')
        #-- IUE spectra    
        if not os.path.isdir(iue_direc) or force_download:
            out2 = vizier.get_IUE_spectra(ID=os.path.basename(self.ID),directory=iue_direc,select='lo')
            if out2 is None:
                out2 = []
        else:
            out2 = glob.glob(iue_direc+'/*')
        #-- read them in to combine    
        list_of_spectra = [fits.read_fuse(ff)[:3] for ff in out1]
        list_of_spectra+= [fits.read_iue(ff)[:3] for ff in out2]
        #-- combine
        wave,flux,err,nspec = tools.combine(list_of_spectra)
        #-- add to the master
        N = len(flux)
        units = ['erg/s/cm2/AA']*N
        source = ['specphot']*N
        self.add_photometry_fromarrays(flux,err,units,photbands,source,flags=nspec)
        
        #-- write to file
        self.master = fix_master(self.master,e_default=0.1)
        logger.info('\n'+photometry2str(self.master))
        self.save_photometry()
        
            
            
    
    def exclude(self,names=None,wrange=None,sources=None,indices=None):
        """
        Exclude (any) photometry from fitting process.
        
        If called without arguments, all photometry will be excluded.
        """
        if names is None and wrange is None and sources is None and indices is None:
            wrange = (-np.inf,np.inf)
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=False,ptype='all')
    
    def exclude_colors(self,names=None,wrange=None,sources=None,indices=None):
        """
        Exclude (color) photometry from fitting process.
        """
        if names is None and wrange is None and sources is None and indices is None:
            wrange = (-np.inf,0)
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=False,ptype='col')
    
    def exclude_abs(self,names=None,wrange=None,sources=None,indices=None):
        """
        Exclude (absolute) photometry from fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=False,ptype='abs')
    
    
    def include(self,names=None,wrange=None,sources=None,indices=None):
        """
        Include (any) photometry in fitting process.
        """
        if names is None and wrange is None and sources is None and indices is None:
            wrange = (-np.inf,np.inf)
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=True,ptype='all')
    
    def include_colors(self,names=None,wrange=None,sources=None,indices=None):
        """
        Include (color) photometry in fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=True,ptype='col')
    
    def include_abs(self,names=None,wrange=None,sources=None,indices=None):
        """
        Include (absolute) photometry in fitting process.
        """
        decide_phot(self.master,names=names,wrange=wrange,sources=sources,indices=indices,include=True,ptype='abs')
    
    def set_photometry_scheme(self,scheme,infrared=(1,'micron')):
        """
        Set a default scheme of colors/absolute values to fit the SED.
        
        Possible values:
            
            1. scheme = 'abs': means excluding all colors, including all absolute values
            2. scheme = 'color': means including all colors, excluding all absolute values
            3. scheme = 'combo': means inculding all colors, and one absolute value per
            system (the one with the smallest relative error)
            4. scheme = 'irfm': means mimic infrared flux method: choose absolute values
            in the infrared (define with C{infrared}), and colours in the optical
        
        @param infrared: definition of start of infrared for infrared flux method
        @type infrared: tuple (value <float>, unit <str>)
        """
        #-- only absolute values: real SED fitting
        if 'abs' in scheme.lower():
            self.master['include'][self.master['color']] = False
            self.master['include'][-self.master['color']] = True
            logger.info('Fitting procedure will use only absolute fluxes (%d)'%(sum(self.master['include'])))
        #-- only colors: color fitting
        elif 'col' in scheme.lower():
            self.master['include'][self.master['color']] = True
            self.master['include'][-self.master['color']] = False
            logger.info('Fitting procedure will use only colors (%d)'%(sum(self.master['include'])))
        #-- combination: all colors and one absolute value per system
        elif 'com' in scheme.lower():
            self.master['include'][self.master['color'] & self.master['include']] = True
            systems = np.array([photband.split('.')[0] for photband in self.master['photband']])
            set_systems = sorted(list(set(systems)))
            for isys in set_systems:
                keep = -self.master['color'] & (isys==systems) & self.master['include']
                if not sum(keep): continue
                index = np.argmax(self.master['e_cmeas'][keep]/self.master['cmeas'][keep])
                index = np.arange(len(self.master))[keep][index]
                self.master['include'][keep] = False
                self.master['include'][index] = True
            logger.info('Fitting procedure will use colors + one absolute flux for each system (%d)'%(sum(self.master['include'])))
        #-- infrared flux method scheme
        elif 'irfm' in scheme.lower():
            #-- check which measurements are in the infrared
            for i,meas in enumerate(self.master):
                #-- if measurement is in the infrared and it is a color, remove it (only use absolute values in IR)
                #-- for colors, make sure *both* wavelengths are checked
                if meas['color']:
                    bands = filters.get_color_photband(meas['photband'])
                    eff_waves = filters.eff_wave(bands)
                    is_infrared = any([(conversions.Unit(*infrared)<conversions.Unit(eff_wave,'AA')) for eff_wave in eff_waves])
                else:
                    is_infrared = conversions.Unit(*infrared)<conversions.Unit(meas['cwave'],'AA')
                if is_infrared and meas['color']:
                    self.master['include'][i] = False
                #-- if measurement is in the infrared and it is not color, keep it (only use absolute values in IR)
                elif is_infrared:
                    self.master['include'][i] = True
                #-- if measurement is not in the infrared and it is a color, keep it (only use colors in visible)
                elif not is_infrared and meas['color']:
                    self.master['include'][i] = True
                #-- if measurement is not in the infrared and it is not a color, remove it (only use colors in visible)
                elif not is_infrared:
                    self.master['include'][i] = False
            use_colors = sum(self.master['include'] & self.master['color'])
            use_abs = sum(self.master['include'] & -self.master['color'])
            logger.info('Fitting procedure will use colors (%d) < %.3g%s < absolute values (%d)'%(use_colors,infrared[0],infrared[1],use_abs))
        
    
    def add_photometry_fromarrays(self,meas,e_meas,units,photbands,source,flags=None):
        """
        Add photometry from numpy arrays
        
        By default unflagged, and at the ra and dec of the star. No color and
        included.
        
        @param meas: original measurements (fluxes, magnitudes...)
        @type meas: array
        @param e_meas: error on original measurements in same units as C{meas}
        @type e_meas: array
        @param units: units of original measurements
        @type units: array of strings
        @param photbands: photometric passbands of original measurements
        @type photbands: array of strings
        @param source: source of original measurements
        @type source: array of strings
        """
        if flags is None:
            flags = np.nan*np.ones(len(meas))
        #-- if master record array does not exist, make a new one
        if not hasattr(self,'master') or self.master is None:
            dtypes = [('meas','f8'),('e_meas','f8'),('flag','S20'),('unit','S30'),('photband','S30'),('source','S50'),('_r','f8'),('_RAJ2000','f8'),\
                       ('_DEJ2000','f8'),('cwave','f8'),('cmeas','f8'),('e_cmeas','f8'),('cunit','S50'),('color',bool),('include',bool),\
                       ('bibcode','S20'),('comments','S200')]
            logger.info('No previous measurements available, initialising master record')
            self.master = np.rec.fromarrays(np.array([ [] for i in dtypes]), dtype=dtypes)
            _to_unit = 'erg/s/cm2/AA'
        else:
            _to_unit = self.master['cunit'][0]
        extra_master = np.zeros(len(meas),dtype=self.master.dtype)

        for i,(m,e_m,u,p,s,f) in enumerate(zip(meas,e_meas,units,photbands,source,flags)):
            photsys,photband = p.split('.')
            if filters.is_color(p):
                to_unit = 'flux_ratio'
                color = True
            else:
                to_unit = _to_unit
                color = False
            if e_m>0:
                cm,e_cm = conversions.convert(u,to_unit,m,e_m,photband=p)
            else:
                cm,e_cm = conversions.convert(u,to_unit,m,photband=p),np.nan
            eff_wave = filters.eff_wave(p)
            extra_master['cmeas'][i] = cm
            extra_master['e_cmeas'][i] = e_cm
            extra_master['cwave'][i] = eff_wave
            extra_master['cunit'][i] = to_unit
            extra_master['color'][i] = color
            extra_master['include'][i] = True
            extra_master['meas'][i] = meas[i]
            extra_master['e_meas'][i] = e_meas[i]
            extra_master['unit'][i] = units[i]
            extra_master['photband'][i] = photbands[i]
            extra_master['source'][i] = source[i]
            extra_master['flag'][i] = flags[i]
            if 'bibcode' in extra_master.dtype.names:
                extra_master['bibcode'][i] = '-'
            if 'comments' in extra_master.dtype.names:
                extra_master['comments'][i] = '-'
        
        logger.info('Original measurements:\n%s'%(photometry2str(self.master)))
        logger.info('Appending:\n%s'%(photometry2str(extra_master)))
        self.master = fix_master(np.hstack([self.master,extra_master]))
        #self.master = np.hstack([self.master,extra_array])
        logger.info('Final measurements:\n%s'%(photometry2str(self.master)))

    #}
    #{ Additional information
    
    def is_target(self,name):
        """
        Check if this SED represents the object `name'.
        
        Purpose: solve alias problems. Maybe the ID is 'HD129929', and you are
        checking for "V* V836 Cen", which is the same target.
        
        @param name: object name
        @type name: str
        @return: True if this object instance represent the target "name".
        @rtype: bool
        """
        try:
            info = sesame.search(name)
            oname = info['oname']
        except:
            logger.warning('Did not find {:s} on Simbad}'.format(name))
            return False
        if oname==self.info['oname']:
            return True
    
    def has_photfile(self):
        """
        Check if this SED has a phot file.
        
        @return: True if this object instance has a photfile
        @rtype: bool
        """
        return os.path.isfile(self.photfile)
    
    def get_distance_from_plx(self,plx=None,lutz_kelker=True,unit='pc'):
        """
        Get the probability density function for the distance, given a parallax.
        
        If no parallax is given, the catalogue value will be used. If that one
        is also not available, a ValueError will be raised.
        
        Parallax must be given in mas.
        
        If the parallax is a float without uncertainty, an uncertainty of 0 will
        be assumed.
        
        If C{lutz_kelker=True}, a probability density function will be given for
        the distance, otherwise the distance is computed from simply inverting
        the parallax.
        
        Distance is returned in parsec (pc).
        
        @return: distance
        @rtype: (float,float)
        """
        #-- we need parallax and galactic position
        if plx is None and 'plx' in self.info:
            plx = self.info['plx']['v'],self.info['plx']['e']
        elif plx is None:
            raise ValueError('distance cannot be computed from parallax (no parallax)')
        if not 'galpos' in self.info:
            raise ValueError('distance cannot be computed from parallax (no position)')
        gal = self.info['galpos']
        #-- if the parallax has an uncertainty and the Lutz-Kelker bias needs
        #   to be taken into account, compute probability density function
        if lutz_kelker and hasattr(plx,'__iter__'):
            d = np.logspace(np.log10(0.1),np.log10(25000),100000)
            dprob = distance.distprob(d,gal[1],plx)
            dprob = dprob / dprob.max()
            logger.info('Computed distance to target with Lutz-Kelker bias')
            return d,dprob
        #-- just invert (with or without error)
        else:
            dist = (conversions.Unit(plx,'kpc-1')**-1).convert(unit)
            logger.info('Computed distance to target via parallax inversion: {:s}'.format(dist))
            return dist.get_value()
    
    def get_interstellar_reddening(self,distance=None, Rv=3.1):
        """
        Under construction.
        """
        gal = self.info['galpos']
        if distance is None:
            distance = self.get_distance_from_plx(lutz_kelker=False,unit='pc')[0]
        output = {}
        for model in ['arenou','schlegel','drimmel','marshall']:
            ext = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance)
            if ext is not None:
                output[model] = ext/Rv
        return output
    
    def get_angular_diameter(self):
        """
        Under construction.
        """
        raise NotImplementedError
    
    def compute_distance(self,mtype='igrid_search'):
        """
        Compute distance from radius and angular diameter (scale factor).
        
        The outcome is added to the C{results} attribute::
        
            distance = r/sqrt(scale)
        
        This particularly useful when you added constraints from solar-like
        oscillations (L{add_constraint_slo}).
        
        @return: distance,uncertainty (pc)
        @rtype: (float,float)
        """
        grid = self.results[mtype]['grid']
        radius = grid['radius']
        e_radius = grid['e_radius']
        scale = grid['scale']
        e_scale = grid['escale']
        
        radius = conversions.unumpy.uarray([radius,e_radius])
        scale = conversions.unumpy.uarray([scale,e_scale])
        
        distance = radius/scale**0.5
        distance = conversions.convert('Rsol','pc',distance)
        distance,e_distance = conversions.unumpy.nominal_values(distance),\
                      conversions.unumpy.std_devs(distance)
        self.results[mtype]['grid'] = pl.mlab.rec_append_fields(grid,\
                        ['distance','e_distance'],[distance,e_distance])
        d,e_d = distance.mean(),distance.std() 
        logger.info("Computed distance from R and scale: {0}+/-{1} pc".format(d,e_d))
        return d,e_d
    
    
    
    
    #}
    
    def generate_ranges(self,type='single',start_from='igrid_search',distribution='uniform',**kwargs):
        """
        Generate sensible search range for each parameter.
        """
        limits = {}
        exist_previous = (start_from in self.results and 'CI' in self.results[start_from])
        #-- how many stellar components do we have? It has to be given when no
        #   ranges are given. Otherwise, we can derive it from the ranges.
        components = 1
        if type is None:
            components = max([(isinstance(kwargs[key][0],tuple) and len(kwargs[key]) or 1) for key in kwargs if (kwargs[key] is not None)])
            type = (components==1) and 'single' or 'multiple-{0:d}'.format(components)
        elif type=='binary':
            components = 2
        elif 'multiple' in type:
            components = int(type.split('-')[-1])
        #-- the field names are 'teff', 'logg' etc for the primary component, all
        #   others have a postfix '-n' with n the number of the component
        postfixes = [''] + ['-{0:d}'.format(i) for i in range(2,components+1)]
        #-- run over all parnames, and generate the parameters ranges for each
        #   component:
        #   (1) if a previous parameter search is done and the user has not set
        #       the parameters range, derive a sensible parameter range from the
        #       previous results. 
        #   (2) if no previous search is done, use infinite upper and lower bounds
        #   (3) if ranges are given, use those. Cycle over them, and use the
        #       strategy from (1) if no parameter ranges are given for a specific
        #       component.
        for par_range_name in kwargs:
            #-- when "teffrange" is given, par_range will be the value of the
            #   keyword "teffrange", and parname will be "teff".
            parname = par_range_name.rsplit('range')[0]
            parrange = kwargs.get(par_range_name,None)
            #-- the three cases described above:
            if exist_previous and parrange is None:
                parrange = []
                for postfix in postfixes:
                    lkey = parname+postfix+'_l'
                    ukey = parname+postfix+'_u'
                    #-- if the parameters was not used in the fit, stick to the
                    #   default given value
                    if not lkey in self.results[start_from]['CI']:
                        parrange.append(kwargs[par_range_name])
                    #-- else we can derive a better parameter range
                    else:
                        parrange.append((self.results[start_from]['CI'][lkey],
                                         self.results[start_from]['CI'][ukey]))
            elif parrange is None:
                parrange = [(-np.inf,np.inf) for i in range(components)]
            else:
                if components>1:
                    print parrange,postfixes
                    parrange = [(iparrange is not None) and iparrange or \
                              (exist_previous and (self.results[start_from]['CI'][parname+postfix+'_l'],\
                                                   self.results[start_from]['CI'][parname+postfix+'_u'])\
                                              or (-np.inf,np.inf))\
                               for iparrange,postfix in zip(parrange,postfixes)]
            #-- if there is only one component, make sure the range is a tuple of a tuple
            if not isinstance(parrange[0],tuple):
                parrange = (parrange,)
            #-- now, the ranges are (lower,upper) for the uniform distribution,
            #   and (mu,scale) for the normal distribution
            if distribution=='normal':
                parrange = [((i[1]-i[0])/2.,(i[1]-i[0])/6.) for i in parrange]
            elif distribution!='uniform':
                raise NotImplementedError, 'Any distribution other than "uniform" and "normal" has not been implemented yet!'
            #-- now we *don't* want a tuple if there is only one component:
            if components==1:
                limits[par_range_name] = parrange[0]
            else:
                limits[par_range_name] = parrange
        #-- this returns the kwargs but with filled in limits, and confirms
        #   the type if it was given, or gives the type when it needed to be derived
        logger.info('Parameter ranges calculated for type {0:s}, starting from {1:s} and using distribution {2:s}.'.format(type,start_from,distribution))
        return limits,type
    
    #{ Fitting routines
    
    def clip_grid(self,mtype='igrid_search',CI_limit=None):
        """
        Clip grid on CI limit, to save memory.
        
        @param mtype: type or results to clip
        @type mtype: str
        @param CI_limit: confidence limit to clip on
        @type CI_limit: float (between 0 (clips everything) and 1 (clips nothing))
        """
        if CI_limit is None:
            CI_limit = self.CI_limit
        new_grid = self.results[mtype]['grid']
        new_grid = new_grid[new_grid['ci_red']<=CI_limit]
        self.results[mtype]['grid'] = new_grid
        logger.info("Clipped grid at {:.6f}%".format(CI_limit*100))
    
    def calculate_confidence_intervals(self,mtype='igrid_search',chi2_type='red',CI_limit=None):
        """
        Compute confidence interval of all columns in the results grid.
        
        @param mtype: type of results to compute confidence intervals of
        @type mtype: str
        @param chi2_type: type of chi2 (raw or reduced)
        @type chi2_type: str ('raw' or 'red')
        @param CI_limit: confidence limit to clip on
        @type CI_limit: float (between 0 (clips everything) and 1 (clips nothing))
        """
        #-- get some info
        grid_results = self.results[mtype]['grid']
        if not 'CI' in self.results[mtype]:
            self.results['igrid_search']['CI'] = {}
        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit
        #-- the chi2 grid is ordered: where is the value closest to the CI limit?
        region = self.results[mtype]['grid']['ci_'+chi2_type]<=CI_limit
        if sum(region)==0:
            raise ValueError("No models in the sample have a chi2_{} below the limit {}. Try increasing the number of models, increasing the CI_limit or choose different photometry.".format(chi2_type,CI_limit))
        #-- now compute the confidence intervals
        for name in grid_results.dtype.names:
            self.results['igrid_search']['CI'][name+'_l'] = grid_results[name][region].min()
            self.results['igrid_search']['CI'][name] = grid_results[name][-1]
            self.results['igrid_search']['CI'][name+'_u'] = grid_results[name][region].max()
            logger.info('%i%% CI %s: %g <= %g <= %g'%(CI_limit*100,name,self.results['igrid_search']['CI'][name+'_l'],
                                                           self.results['igrid_search']['CI'][name],
                                                           self.results['igrid_search']['CI'][name+'_u']))
    
    
    
    def igrid_search(self,points=100000,teffrange=None,loggrange=None,ebvrange=None,
                          zrange=None,radiusrange=None,rvrange=None,vradrange=(0,0),
                          compare=True,df=None,CI_limit=None,set_model=True,**kwargs):
        """
        Fit fundamental parameters using a (pre-integrated) grid search.
        
        If called consecutively, the ranges will be set to the CI_limit of previous
        estimations, unless set explicitly.
        
        If called for the first time, the ranges will be +/- np.inf by defaults,
        unless set explicitly.
        """
        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit
        
        #-- set defaults limits
        ranges,rtype = self.generate_ranges(teffrange=teffrange,\
                             loggrange=loggrange,ebvrange=ebvrange,\
                             zrange=zrange,rvrange=rvrange,vradrange=vradrange)
        #-- compute the degrees of freedom
        if df is None:
            df = 1
            df_info = ['theta']
            for range_name in ranges:
                if not np.allclose(ranges[range_name][0],ranges[range_name][1]):
                    df += 1
                    df_info.append(range_name)
            logger.info('Degrees of freedom = {} ({})'.format(df,', '.join(df_info)))
        
                #-- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'%(photometry2str(self.master[include_grid])))
        
        #-- build the grid, run over the grid and calculate the CHI2
        pars = fit.generate_grid_single_pix(self.master['photband'][include_grid],points=points,**ranges) 
        chisqs,scales,escales,lumis = fit.igrid_search_pix(self.master['cmeas'][include_grid],
                             self.master['e_cmeas'][include_grid],
                             self.master['photband'][include_grid],**pars)
        parnames = sorted(pars.keys())
        pardtypes = [(name,'f8') for name in parnames]+[('chisq','f8'),('scale','f8'),('escale','f8'),('labs','f8')]
        pararrays = [pars[name] for name in parnames]+[chisqs,scales,escales,lumis]
        grid_results = np.rec.fromarrays(pararrays,dtype=pardtypes)
                  
        #-- exclude failures
        failures = np.isnan(grid_results['chisq'])
        if sum(failures):
            logger.info('Excluded {0} failed results (nan)'.format(sum(failures)))
            grid_results = grid_results[-failures]
        
        #-- make room for chi2 statistics
        grid_results = mlab.rec_append_fields(grid_results, 'ci_raw', np.zeros(len(grid_results)))
        grid_results = mlab.rec_append_fields(grid_results, 'ci_red', np.zeros(len(grid_results)))
        
        #-- take the previous results into account if they exist:
        if not 'igrid_search' in self.results:
            self.results['igrid_search'] = {}
        elif 'grid' in self.results['igrid_search']:
            logger.info('Appending previous results ({:d}+{:d})'.format(len(self.results['igrid_search']['grid']),len(grid_results)))
            ex_names = grid_results.dtype.names
            ex_grid = np.rec.fromarrays([self.results['igrid_search']['grid'][exname] for exname in ex_names],
                                        names=ex_names)
            grid_results = np.hstack([ex_grid,grid_results])
        
        #-- inverse sort according to chisq: this means the best models are at
        #   the end (mainly for plotting reasons, so that the best models
        #   are on top).
        sa = np.argsort(grid_results['chisq'])[::-1]
        grid_results = grid_results[sa]
        
        #-- do the statistics
        #   degrees of freedom: teff,logg,E(B-V),theta,Z
        N = sum(include_grid)
        k = N-df
        if k<=0:
            logger.warning('Not enough data to compute CHI2: it will not make sense')
            k = 1
        #-- rescale if needed and compute confidence intervals
        factor = max(grid_results['chisq'][-1]/k,1)
        logger.warning('CHI2 rescaling factor equals %g'%(factor))
        CI_raw = scipy.stats.distributions.chi2.cdf(grid_results['chisq'],k)
        CI_red = scipy.stats.distributions.chi2.cdf(grid_results['chisq']/factor,k)
        
        #-- add the results to the record array and to the results dictionary
        grid_results['ci_raw'] = CI_raw
        grid_results['ci_red'] = CI_red
        self.results['igrid_search']['grid'] = grid_results
        self.results['igrid_search']['factor'] = factor
        
        #-- compute the confidence intervals
        self.calculate_confidence_intervals(mtype='igrid_search',chi2_type='red',CI_limit=CI_limit)
        
        #-- remember the best model
        if set_model: self.set_best_model(type=type)
            
        
    def imc(self,teffrange=None,loggrange=None,ebvrange=None,zrange=None,start_from='imc',\
                 distribution='uniform',points=None,fitmethod='fmin',disturb=True):
        
        limits,type = self.generate_ranges(teffrange=teffrange,loggrange=loggrange,\
                                      ebvrange=ebvrange,zrange=zrange,distribution=distribution,\
                                      start_from='imc')
        
        #-- grid search on all include data: extract the best CHI2
        include = self.master['include']
        meas = self.master['cmeas'][include]
        if disturb:
            emeas = self.master['e_cmeas'][include]
        else:
            emeas = self.master['e_cmeas'][include]*1e-6
        photbands = self.master['photband'][include]
        logger.info('The following measurements are included in the fitting process:\n%s'%(photometry2str(self.master[include])))
        
        #-- generate initial guesses
        teffs,loggs,ebvs,zs,radii = fit.generate_grid(self.master['photband'][include],type=type,points=points+25,**limits)         
        NrPoints = len(teffs)>points and points or len(teffs)
        firstoutput = np.zeros((len(teffs)-NrPoints,9))
        output = np.zeros((NrPoints,9))
        
        #-- fit the original data a number of times
        for i,(teff,logg,ebv,z) in enumerate(zip(teffs[NrPoints:],loggs[NrPoints:],ebvs[NrPoints:],zs[NrPoints:])):
            try:
                fittedpars,warnflag = fit.iminimize2(meas,emeas,photbands,teff,logg,ebv,z,fitmethod=fitmethod)
                firstoutput[i,1:] = fittedpars
                firstoutput[i,0] = warnflag
            except IOError:
                firstoutput[i,0] = 3
        
        logger.info("{0}/{1} fits on original data failed (max func call)".format(sum(firstoutput[:,0]==1),firstoutput.shape[0]))
        logger.info("{0}/{1} fits on original failed (max iter)".format(sum(firstoutput[:,0]==2),firstoutput.shape[0]))
        logger.info("{0}/{1} fits on original data failed (outside of grid)".format(sum(firstoutput[:,0]==3),firstoutput.shape[0]))
        
        #-- retrieve the best fitting result and make it the first entry of output
        keep = (firstoutput[:,0]==0) & (firstoutput[:,1]>0)
        best = firstoutput[keep,-2].argmin()
        output[-1,:] = firstoutput[keep][best,:]
        
        # calculate the factor with which to multiply the scale
        #factor = np.sqrt(output[-1,5]/len(meas))
        #print factor
        
        #-- now do the actual Monte Carlo simulation
        for i,(teff,logg,ebv,z) in enumerate(zip(teffs[:NrPoints-1],loggs[:NrPoints-1],ebvs[:NrPoints-1],zs[:NrPoints-1])):
            newmeas = meas + np.random.normal(scale=emeas) #*factor)
            try:
                fittedpars,warnflag = fit.iminimize2(newmeas,emeas,photbands,teff,logg,ebv,z,fitmethod=fitmethod)
                output[i,1:] = fittedpars
                output[i,0] = warnflag
            except IOError:
                output[i,0] = 3
                
        logger.info("{0}/{1} MC simulations failed (max func call)".format(sum(output[:,0]==1),NrPoints))
        logger.info("{0}/{1} MC simulations failed (max iter)".format(sum(output[:,0]==2),NrPoints))
        logger.info("{0}/{1} MC simulations failed (outside of grid)".format(sum(output[:,0]==3),NrPoints))
        
        #-- remove nonsense results
        keep = (output[:,0]==0) & (output[:,1]>0)
        output = output[keep]
        output = np.rec.fromarrays(output[:,1:-1].T,names=['teff','logg','ebv','z','labs','chisq','scale'])
        #-- derive confidence intervals and median values
        #print np.median(output[:,1:],axis=0)
        
        self.results['imc'] = {}
        self.results['imc']['CI'] = {}
        self.results['imc']['grid'] = output
        CI_limit = 0.95
        for name in output.dtype.names:
            sortarr = np.sort(output[name])
            trimarr = scipy.stats.trimboth(sortarr,(1-CI_limit)/2.) # trim 2.5% of top and bottom, to arrive at 95% CI
            self.results['imc']['CI'][name+'_l'] = trimarr.min()
            self.results['imc']['CI'][name] = output[name][-1]#np.median(output[name])
            self.results['imc']['CI'][name+'_u'] = trimarr.max()
            logger.info('%i%% CI %s: %g <= %g <= %g'%(CI_limit*100,name,self.results['imc']['CI'][name+'_l'],
                                                           self.results['imc']['CI'][name],
                                                           self.results['imc']['CI'][name+'_u']))
        self.set_best_model(mtype='imc')

    
    #}
    
    #{ Interfaces
    
    def clear(self):
        """
        Clear the results.
        """
        self.results = {}
    
    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004',type='single'):
        """
        Get reddenend and unreddened model
        """
        logger.info('Interpolating approximate full SED of best model')
        
        #-- synthetic flux
        include = self.master['include']
        synflux = np.zeros(len(self.master['photband']))
        keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
        keep = keep & include
        
        if mtype in ['igrid_search']:
            #-- get the metallicity right
            files = model.get_file(z='*')
            metals = np.array([pyfits.getheader(ff)['Z'] for ff in files])
            metals = metals[np.argmin(np.abs(metals-self.results[mtype]['CI']['z']))]
            scale = self.results[mtype]['CI']['scale']
            #-- get (approximated) reddened and unreddened model
            wave,flux = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                    logg=self.results[mtype]['CI']['logg'],
                                    ebv=self.results[mtype]['CI']['ebv'],
                                    z=metals,
                                    law=law)
            wave_ur,flux_ur = model.get_table(teff=self.results[mtype]['CI']['teff'],
                                        logg=self.results[mtype]['CI']['logg'],
                                        ebv=0,
                                        z=metals,
                                        law=law)
            #-- get synthetic photometry
            synflux_,Labs = model.get_itable(teff=self.results[mtype]['CI']['teff'],
                               logg=self.results[mtype]['CI']['logg'],
                               ebv=self.results[mtype]['CI']['ebv'],
                               z=self.results[mtype]['CI']['z'],
                               photbands=self.master['photband'][keep])
 
            flux,flux_ur = flux*scale,flux_ur*scale
                
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
            self.results[mtype]['model'] = wave,flux,flux_ur
            self.results[mtype]['synflux'] = eff_waves,synflux,self.master['photband']
            self.results[mtype]['chi2'] = chi2
    
    def set_model(self,wave,flux,label='manual'):
        """
        Manually set the best SED model.
        
        This is particularly useful when working with calibrators for which there
        is an observed SED.
        
        The label will be used to store the model in the C{results} attribute.
        
        @param wave: wavelength array (angstrom)
        @type wave: ndarray
        @param flux: flux array (erg/s/cm2/AA)
        @type flux: ndarray
        @param label: key used to store the model and synthetic photometry in C{results}
        @type label:s str
        """
        #-- necessary information
        photbands = self.master['photband']
        is_color = self.master['color']
        include = self.master['include']
        synflux = np.zeros(len(photbands))
        
        #-- compute synthetic photometry
        synflux[(-is_color) & include] = model.synthetic_flux(wave,flux,photbands[(-is_color) & include])
        synflux[is_color & include] = model.synthetic_color(wave,flux,photbands[is_color & include])
        chi2 = (self.master['cmeas']-synflux)**2/self.master['e_cmeas']**2
        eff_waves = filters.eff_wave(photbands)
        
        if not label in self.results:
            self.results[label] = {}
        self.results[label]['model'] = wave,flux,flux
        self.results[label]['synflux'] = eff_waves,synflux,photbands
        self.results[label]['chi2'] = chi2
        logger.debug('Stored model SED in {}'.format(label))
    
    def get_model(self,label='igrid_search'):
        """
        Retrieve the best SED model.
        
        B{Warning}: the grid search interpolation is also done with interpolation
        in metallicity, while this is not the case for the best full SED model.
        Interpolation of the full SED is only done in teff/logg, and the
        metallicity is chosen to be equal to the closest grid point. On the
        other hand, while the reddening value is interpolated in the grid search,
        this is not the case for the best model (the model is simply reddened
        according the found best value). So don't be surprised if you find
        differences between the values of C{self.results[label]['synflux']},
        which are the value sfor the interpolated photometric points of the
        best model, and the synthetic photometry obtained by manual integration
        of the returned full SED model. Those differences should be incredible
        small and mainly due to the metallicity.
        """
        wave,flux,urflux = self.results[label]['model']
        return wave,flux,urflux
    
    def chi2(self,select=None,reduced=False,label='igrid_search'):
        """
        Calculate chi2 of best model.
        
        TEMPORARY!!! API WILL CHANGE!!!
        """
        #-- select photometry
        master = self.master.copy()
        keep = np.ones(len(master),bool)
        if isinstance(select,dict):
            for key in select:
                keep = keep & (master[key]==select[key])
        master = master[keep]
        #-- compute synthetic phtoometry
        photbands = master['photband']
        is_color = master['color']
        synflux = np.zeros(len(photbands))
        
        wave,flux,urflux = self.get_model(label=label)
        synflux[-is_color] = model.synthetic_flux(wave,flux,photbands[-is_color])
        synflux[is_color] = model.synthetic_color(wave,flux,photbands[is_color])
        x2 = (master['cmeas']-synflux)**2/master['e_cmeas']**2
        x2 = x2.mean()
        return x2
        
    #}
    #{ Add constraints
    
    def add_constraint_distance(self,distance=None,mtype='igrid_search',**kwargs):
        """
        Use the distance to compute additional information.
        
        Compute radii, absolute luminosities and masses, and add them to the
        results.
        
        Extra kwargs go to L{get_distance_from_plx}.
        
        B{Warning:} after calling this function, the C{labs} column in the grid
        is actually absolutely calibrated and reflects the true absolute
        luminosity instead of the absolute luminosity assuming 1 solar radius.
        
        @param distance: distance in solar units and error
        @type distance: tuple (float,float)
        @param mtype: type of results to add the information to
        @type mtype: str
        """
        if distance is None:
            kwargs['lutz_kelker'] = False # we can't handle asymmetric error bars
            kwargs['unit'] = 'Rsol'
            distance = self.get_distance_from_plx(**kwargs)
        
        #-- compute the radius, absolute luminosity and mass: note that there is
        #   an uncertainty on the scaling factor *and* on the distance!
        scale = conversions.unumpy.uarray([self.results[mtype]['grid']['scale'],\
                                           self.results[mtype]['grid']['escale']])
        distance = conversions.ufloat(distance)
        radius = distance*np.sqrt(self.results[mtype]['grid']['scale']) # in Rsol
        labs = self.results[mtype]['grid']['labs']*radius**2
        mass = conversions.derive_mass((self.results[mtype]['grid']['logg'],'[cm/s2]'),\
                                      (radius,'Rsol'),unit='Msol')
                                      
        #-- store the results in the grid
        mass,e_mass = conversions.unumpy.nominal_values(mass),\
                      conversions.unumpy.std_devs(mass)
        radius,e_radius = conversions.unumpy.nominal_values(radius),\
                      conversions.unumpy.std_devs(radius)
        labs,e_labs = conversions.unumpy.nominal_values(labs),\
                      conversions.unumpy.std_devs(labs)
        self.results[mtype]['grid']['labs'] = labs # Labs is already there, overwrite
        #-- check if the others are already in there or not:
        labels,data = ['e_labs','radius','e_radius','mass','e_mass'],\
                      [e_labs,radius,e_radius,mass,e_mass]
        if 'e_labs' in self.results[mtype]['grid'].dtype.names:
            for idata,ilabel in zip(data,labels):
                self.results[mtype]['grid'][ilabel] = idata
        else:
            self.results[mtype]['grid'] = pl.mlab.rec_append_fields(self.results[mtype]['grid'],\
                     labels,data)
                     
        #-- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: distance (improved luminosity, added info on radius and mass)')
    
    def add_constraint_slo(self,numax,Deltanu0,mtype='igrid_search',chi2_type='red'):
        """
        Use diagnostics from solar-like oscillations to put additional constraints on the parameters.
        
        If the results are constrained with the distance before L{add_constraint_distance},
        then these results are combined with the SLO constraints.
        
        @param numax: frequency of maximum amplitude
        @type numax: 3-tuple (value,error,unit)
        @param Deltanu0: large separation (l=0)
        @type Deltanu0: 3-tuple (value,error,unit)
        @param chi2_type: type of chi2 (raw or reduced)
        @type chi2_type: str ('raw' or 'red')
        """
        grid = self.results[mtype]['grid']
        #-- we need the teffs, so that we can compute the logg and radius using
        #   the solar-like oscillations information.
        teff = grid['teff']
        logg = grid['logg']
        pvalues = [1-grid['ci_'+chi2_type]]
        names = ['ci_'+chi2_type+'_phot']
        
        #-- we *always* have a logg, that's the way SEDs work:
        logg_slo = conversions.derive_logg_slo((teff,'K'),numax)
        logg_slo,e_logg_slo = conversions.unumpy.nominal_values(logg_slo),\
                      conversions.unumpy.std_devs(logg_slo)
        #   compute the probabilities: scale to standard normal
        logg_prob = scipy.stats.distributions.norm.sf(abs(logg_slo-logg)/e_logg_slo)
        
        pvalues.append(logg_prob)
        names.append('ci_logg_slo')        
        
        #-- but we only sometimes have a radius (viz. when the distance is
        #   known). There is an uncertainty on that radius!
        radius_slo = conversions.derive_radius_slo(numax,Deltanu0,(teff,'K'),unit='Rsol')
        radius_slo,e_radius_slo = conversions.unumpy.nominal_values(radius_slo),\
                                  conversions.unumpy.std_devs(radius_slo)
        if 'radius' in grid.dtype.names:
            radius = grid['radius']
            e_radius = grid['e_radius']
            #total_error = np.sqrt( (e_radius_slo**2+e_radius**2))
            total_error = np.sqrt( (e_radius_slo**2+e_radius**2)/2. + (radius-radius_slo)**2/4.)
            radius_prob = scipy.stats.distributions.norm.sf(abs(radius_slo-radius)/total_error)
            pvalues.append(radius_prob)
            names.append('ci_radius_slo')
            #-- combined standard deviation and mean for two populations with
            #   possibly zero intersection (wiki: standard deviation)
            #total_error = np.sqrt( (e_radius_slo**2+e_radius**2)/2. + (radius-radius_slo)**2/4.)
            total_mean = (radius+radius_slo)/2.
            grid['radius'] = total_mean
            grid['e_radius'] = total_error
            logger.info('Added contraint: combined existing radius estimate with slo estimate')
        #   otherwise, we just add info on the radius
        else:
            labs = grid['labs']*conversions.unumpy.uarray([radius_slo,e_radius_slo])**2
            labs,e_labs = conversions.unumpy.nominal_values(labs),\
                          conversions.unumpy.std_devs(labs)
            grid['labs'] = labs
            mass = conversions.derive_mass((self.results[mtype]['grid']['logg'],'[cm/s2]'),\
                                           (radius_slo,e_radius_slo,'Rsol'),unit='Msol')
            mass,e_mass = conversions.unumpy.nominal_values(mass),\
                          conversions.unumpy.std_devs(mass)
            #-- now we can also derive the distance:
            scale = conversions.unumpy.uarray([self.results[mtype]['grid']['scale'],\
                                           self.results[mtype]['grid']['escale']])
            distance = radius_slo/conversions.sqrt(scale)
            distance,e_distance = conversions.unumpy.nominal_values(distance),\
                                  conversions.unumpy.std_devs(distance)
            grid = pl.mlab.rec_append_fields(grid,['radius','e_radius','e_labs','mass','e_mass','distance','e_distance'],\
                                        [radius_slo,e_radius_slo,e_labs,mass,e_mass,distance,e_distance])
            logger.info('Added constraint: {0:s} via slo, improved luminosity'.format(', '.join(['radius','e_radius','e_labs','mass','e_mass'])))
            
        #-- combine p values using Fisher's method
        combined_pvals = evaluate.fishers_method(pvalues)
        
        #-- add the information to the grid
        self.results[mtype]['grid'] = pl.mlab.rec_append_fields(grid,\
                     names,pvalues)
        
        #-- and replace the original confidence intervals, and re-order
        self.results[mtype]['grid']['ci_'+chi2_type] = 1-combined_pvals
        sa = np.argsort(self.results[mtype]['grid']['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = self.results[mtype]['grid'][sa]
        
        #-- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: {0:s} via slo and replaced ci_{1:s} with combined CI'.format(', '.join(names),chi2_type))
                      
    def add_constraint_reddening(self,distance=None,ebv=None,e_ebv=0.1,Rv=3.1,\
                model=None,mtype='igrid_search',chi2_type='red',upper_limit=False):
        """
        Use reddening maps to put additional constraints on the parameters.
                
        This constraint assumes that, if C{upper_limit=False}, given the
        distance to target, the reddening from the reddening maps is similar to
        the one derived from the SED. All models where this is not the case will
        be deemed improbable.
        
        If C{upper_limit=True}, all models with an E(B-V) value above C{ebv}
        will be considered improbable.
        
        When you don't set C{ebv}, the value will, by default, be derived from
        Drimmel maps if C{upper_limit=False} and from Schlegel if
        C{upper_limit=True}. You can change this behaviour by setting C{model}
        manually.
                
        @param distance: distance and uncertainty in parsec
        @type distance: tuple (float,float)
        @param ebv: E(B-V) reddening in magnitude
        @type ebv: float
        @param e_ebv: error on the reddening in percentage
        @type e_ebv: float
        @param model: model reddening maps
        @type model: str
        @param mtype: type of results to add the information to
        @type mtype: str
        @param upper_limit: consider the E(B-V) value as an upper limit
        @type upper_limit: bool
        """
        #-- for upper limits on E(B-V), we best use Schlegel maps by default,
        #   otherwise we best use Drimmel maps.
        if model is None and upper_limit:
            model = 'schlegel'
        elif model is None:
            model = 'drimmel'
        #-- if I need to figure out the reddening myself, I need the distance!
        if ebv is None and distance is None:
            distance = self.get_distance_from_plx(lutz_kelker=False,unit='pc')
        #-- let me figure out the reddening if you didn't give it:
        if ebv is None:
            gal = self.info['galpos']
            ebv = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance[0])/Rv
            ebv_u = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance[0]-distance[1])/Rv
            ebv_l = extinctionmodels.findext(gal[0], gal[1], model=model, distance=distance[0]+distance[1])/Rv
            e_ebv = max(ebv-ebv_l,ebv_u-ebv)
        else:
            e_ebv = 0.1*ebv
        grid = self.results[mtype]['grid']
        ebvs = grid['ebv']
        #-- probabilities are differently calculated depending on upper limit.
        if not upper_limit:
            ebv_prob = scipy.stats.distributions.norm.cdf(abs(ebv-ebvs)/e_ebv)
            #-- combine p values using Fisher's method
            combined_pvals = evaluate.fishers_method([1-grid['ci_'+chi2_type],ebv_prob])
            grid['ci_'+chi2_type] = 1-combined_pvals
        else:
            grid['ci_'+chi2_type] = np.where(ebvs<=ebv,grid['ci_'+chi2_type],1.)
        
        #-- and replace the original confidence intervals, and re-order
        sa = np.argsort(grid['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = grid[sa]
        
        #-- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: E(B-V)={0}+/-{1}'.format(ebv,e_ebv))
    
    def add_constraint_angular_diameter(self,angdiam):
        raise NotImplementedError
        
    def add_constraint_mass(self,mass,mtype='igrid_search',chi2_type='red'):
        """
        Add constraints on the mass.
        
        C{mass} must be a tuple, if the second element is smaller than the first,
        it is assumed to be (mu,sigma) from a normal distribution. If the second
        element is larger than the first, it is assumed to be (lower, upper)
        from a uniform distribution.
        """
        normal = mass[0]>mass[1]
        grid = self.results[mtype]['grid']
        masses = grid['mass']
        if normal:
            mass_prob = scipy.stats.distributions.norm.cdf(abs(mass[0]-masses)/mass[1])
            #-- combine p values using Fisher's method
            combined_pvals = evaluate.fishers_method([1-grid['ci_'+chi2_type],mass_prob])
            grid['ci_'+chi2_type] = 1-combined_pvals
        else:
            grid['ci_'+chi2_type] = np.where((masses<=mass[0]) | (mass[1]<=masses),1.,grid['ci_'+chi2_type])
        
        #-- and replace the original confidence intervals, and re-order
        sa = np.argsort(grid['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = grid[sa]
        
        #-- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: mass={0}+/-{1}'.format(*mass))
        
    def add_constraint_evolution_models(self,models='siess2000',\
               ylabels=['age','labs','radius'],e_y=None,
               function='linear',mtype='igrid_search',chi2_type='red'):
        """
        Use stellar evolutionary models to put additional constraints on the parameters.
        """
        grid = self.results[mtype]['grid']
        
        #-- make sure y's and ylabels are iterable
        if isinstance(ylabels,str):
            ylabels = [ylabels]        
        #-- cycle over all yvalues and compute the pvalues
        pvalues = [1-grid['ci_'+chi2_type]]
        add_info = []
        #-- create the evolutionary grid and interpolate the stellar evolutinary
        #   grid on the SED-integrated grid points
        output = np.array([evolutionmodels.get_itable(iteff,ilogg,iz) for iteff,ilogg,iz in zip(grid['teff'],grid['logg'],grid['z'])])
        output = np.rec.fromarrays(output.T,names=['age','labs','radius'])
        for label in ylabels:
            y_interpolated = output[label]
            add_info.append(y_interpolated)
            #-- only add the contraint when it is possible to do so.
            if not label in grid.dtype.names:
                logger.info("Cannot put constraint on {} (not a parameter)".format(label))
                continue
            y_computed = grid[label]
            #-- error on the y-value: either it is computed, it is given or it is
            #   assumed it is 10% of the value
            if e_y is None and ('e_'+label) in grid.dtype.names:
                e_y = np.sqrt(grid['e_'+label]**2+(0.1*y_interpolated)**2)
            elif e_y is None:
                e_y = np.sqrt((0.1*y_computed)**2+(0.1*y_interpolated)**2)
            y_prob = scipy.stats.distributions.norm.sf(abs(y_computed-y_interpolated)/e_y)
            pl.figure()
            pl.subplot(221)
            pl.title(label)
            pl.scatter(y_computed,y_interpolated,c=y_prob,edgecolors='none',cmap=pl.cm.spectral)
            pl.plot([pl.xlim()[0],pl.xlim()[1]],[pl.xlim()[0],pl.xlim()[1]],'r-',lw=2)
            pl.xlim(pl.xlim())
            pl.ylim(pl.xlim())
            pl.xlabel('Computed')
            pl.ylabel('Interpolated')
            pl.colorbar()
            pl.subplot(223)
            pl.title('y_interpolated')
            pl.scatter(grid['teff'],grid['logg'],c=y_interpolated,edgecolors='none',cmap=pl.cm.spectral)
            pl.colorbar()
            pl.subplot(224)
            pl.title('y_computed')
            pl.scatter(grid['teff'],grid['logg'],c=y_computed,edgecolors='none',cmap=pl.cm.spectral)
            pl.colorbar()
            pvalues.append(y_prob)
        
        pl.figure()
        pl.subplot(221)
        pl.title('p1')
        sa = np.argsort(pvalues[0])
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=pvalues[0][sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        pl.subplot(222)
        pl.title('p2')
        sa = np.argsort(pvalues[1])
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=pvalues[1][sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        pl.subplot(223)
        pl.title('p3')
        sa = np.argsort(pvalues[2])
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=pvalues[2][sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        
        #-- combine p values using Fisher's method
        combined_pvals = evaluate.fishers_method(pvalues)
        pl.subplot(224)
        pl.title('pcombined')
        sa = np.argsort(combined_pvals)
        pl.scatter(grid['labs'][sa],grid['radius'][sa],c=combined_pvals[sa],edgecolors='none',cmap=pl.cm.spectral)
        pl.colorbar()
        pl.xlabel('labs')
        pl.ylabel('radius')
        pl.show()    
        
        #-- add the information to the grid (assume that if there one label
        #   not in there already, none of them are
        if not 'c_'+ylabels[0] in grid.dtype.names:
            self.results[mtype]['grid'] = pl.mlab.rec_append_fields(grid,\
                        ['c_'+ylabel for ylabel in ylabels],add_info)
        #-- if they are already in there, overwrite!
        else:
            for info,ylabel in zip(add_info,ylabels):
                self.results[mtype]['grid']['c_'+ylabel] = add_info
        
        #-- and replace the original confidence intervals, and re-order
        self.results[mtype]['grid']['ci_'+chi2_type] = 1-combined_pvals
        
        sa = np.argsort(self.results[mtype]['grid']['ci_'+chi2_type])[::-1]
        self.results[mtype]['grid'] = self.results[mtype]['grid'][sa]
        
        #-- update the confidence intervals
        self.calculate_confidence_intervals(mtype=mtype)
        logger.info('Added constraint: {0:s} via stellar models and replaced ci_{1:s} with combined CI'.format(', '.join(ylabels),chi2_type))    
    
    
    #}
    
    #{ Plotting routines
    @standalone_figure
    def plot_grid(self,x='teff',y='logg',ptype='ci_red',mtype='igrid_search',limit=0.95,d=None,**kwargs):
        """
        Plot grid as scatter plot
        
        PrameterC{ptype} sets the colors of the scattered points (e.g., 'ci_red','z','ebv').
        
        Example usage:
        
        First set the SED:
        
        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()
        
        Then make the plots:
        
        >>> p = pl.figure()
        >>> p = pl.subplot(221);mysed.plot_grid(limit=None)
        >>> p = pl.subplot(222);mysed.plot_grid(x='ebv',y='z',limit=None)
        >>> p = pl.subplot(223);mysed.plot_grid(x='teff',y='ebv',limit=None)
        >>> p = pl.subplot(224);mysed.plot_grid(x='logg',y='z',limit=None)
        
        ]]include figure]]ivs_sed_builder_plot_grid_01.png]
        
        """
        #-- if no distance is set, derive the most likely distance from the plx:
        #if d is None:
            #try:
                #d = self.get_distance_from_plx()
            #except ValueError:
                #d = None
                #logger.info('Distance to {0} unknown'.format(self.ID))
            #if isinstance(d,tuple):
                #d = d[0][np.argmax(d[1])]
        ##-- it's possible that we still don't have a distance
        #if d is not None:
            #logger.info('Assumed distance to {0} = {1:.3e} pc'.format(self.ID,d))
            #radius  = d*np.sqrt(self.results[mtype]['grid']['scale'])
            #radius  = conversions.convert('pc','Rsol',radius) # in Rsol
            #labs = np.log10(self.results[mtype]['grid']['labs']*radius**2) # in [Lsol]
            #mass = conversions.derive_mass((self.results[mtype]['grid']['logg'].copy(),'[cm/s2]'),\
                                           #(radius,'Rsol'),unit='Msol')
        #-- compute angular diameter
        theta = 2*conversions.convert('sr','mas',self.results[mtype]['grid']['scale'])
        
        if limit is not None:
            region = self.results[mtype]['grid']['ci_red']<=limit
        else:
            region = self.results[mtype]['grid']['ci_red']<np.inf
        #-- get the colors and the color scale
        if d is not None and ptype=='labs':
            colors = locals()[ptype][region]
        elif ptype in self.results[mtype]['grid'].dtype.names:
            colors = self.results[mtype]['grid'][ptype][region]
        elif isinstance(ptype,str):
            colors = locals()[ptype][region]
        else:
            colors = ptype[region]
            ptype = 'custom variable'
        
        if 'ci_' in ptype.lower():
            colors *= 100.
            vmin = colors.min()
            vmax = 95.
        else:
            vmin = kwargs.pop('vmin',colors.min())
            vmax = kwargs.pop('vmax',colors.max())
        
        #-- define abbrevation for plotting
        if x in self.results[mtype]['grid'].dtype.names:
            X = self.results[mtype]['grid'][x]
        else:
            X = locals()[x]
            
        if y in self.results[mtype]['grid'].dtype.names:
            Y = self.results[mtype]['grid'][y]
        else:
            Y = locals()[y]
        
        #-- for setting the x/y/color labels
        label_dict = dict(teff='Effective temperature [K]',\
                          z='log (Metallicity Z [$Z_\odot$]) [dex]',\
                          logg=r'log (surface gravity [cm s$^{-2}$]) [dex]',\
                          ebv='E(B-V) [mag]',\
                          ci_raw='Raw probability [%]',\
                          ci_red='Reduced probability [%]',\
                          #labs=r'log (Absolute Luminosity [$L_\odot$]) [dex]',\
                          labs=r'Absolute Luminosity [$L_\odot$]',\
                          radius=r'Radius [$R_\odot$]',\
                          mass=r'Mass [$M_\odot$]',
                          mc=r'MC [Nr. points in hexagonal bin]',
                          rv=r'Extinction parameter $R_v$')
        
        #-- make the plot
        if mtype == 'imc':
            pl.hexbin(X,Y,mincnt=1,cmap=pl.cm.spectral)  #bins='log'
            ptype = 'mc'
                        
            #-- set the limits
            pl.xlim(X.max(),X.min())
            pl.ylim(Y.max(),Y.min())
            cbar = pl.colorbar()
        else:
            #if limit is not None:
                #region = self.results[mtype]['grid']['ci_red']<limit
            #else:
                #region = self.results[mtype]['grid']['ci_red']<np.inf
            ##-- get the colors and the color scale
            #if d is not None and ptype=='labs':
                #colors = locals()[ptype][region]
            #elif ptype in self.results[mtype]['grid'].dtype.names:
                #colors = self.results[mtype]['grid'][ptype][region]
            #else:
                #colors = locals()[ptype][region]
            
            #if 'ci' in ptype:
                #colors *= 100.
                #vmin = colors.min()
                #vmax = 95.
            #else:
                #vmin = kwargs.pop('vmin',colors.min())
                #vmax = kwargs.pop('vmax',colors.max())
            
            #-- grid scatter plot
            pl.scatter(X[region],Y[region],
                 c=colors,edgecolors='none',cmap=pl.cm.spectral,vmin=vmin,vmax=vmax)
            #-- set the limits to only include the 95 interval
            pl.xlim(X[region].max(),X[region].min())
            pl.ylim(Y[region].max(),Y[region].min())
            cbar = pl.colorbar()
            
        #-- mark best value
        pl.plot(X[-1],Y[-1],'r+',ms=40,mew=3)
        
        pl.xlabel(label_dict[x.rstrip('-2')])
        pl.ylabel(label_dict[y.rstrip('-2')])
        if ptype in label_dict:
            cbar.set_label(label_dict[ptype])
        else:
            cbar.set_label(ptype)
        
        logger.info('Plotted %s-%s diagram of %s'%(x,y,ptype))


    @standalone_figure
    def plot_data(self,colors=False, plot_unselected=True,
                  unit_wavelength='angstrom',unit_flux=None,**kwargs):
        """
        Plot only the SED data.
        
        Extra kwargs are passed to plotting functions.
        
        The decorator provides an keyword C{savefig}. When set to C{True}, an
        image name is generated automatically (very long!), and the figure is
        closed. When set to a string, the image is saved with that name, and
        the figure is closed. When C{savefig} is not given, the image will
        stay open so that the user can still access all plot elements and do
        further enhancements.
        
        @param colors: if False, plot absolute values, otherwise plot colors
        (flux ratios)
        @type colors: boolean
        @param plot_unselected: if True, all photometry is plotted, otherwise
        only those that are selected
        """
        if not plot_unselected:
            master = self.master[self.master['include']]
        else:
            master = self.master
        if unit_flux is None:
            unit_flux = master['cunit'][0]
        wave,flux,e_flux = master['cwave'],master['cmeas'],master['e_cmeas']
        sources = master['source']
        iscolor = np.array(master['color'],bool)
        photbands = master['photband']
        indices = np.arange(len(master))
        
        allsystems = np.array([i.split('.')[0] for i in photbands])
        systems = sorted(set(allsystems))
        color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(systems))]
        
        if not colors:
            color_cycle = itertools.cycle(color_cycle)
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            wave = conversions.convert('angstrom',unit_wavelength,wave)
            flux,e_flux = conversions.convert(master['cunit'][0],unit_flux,flux,e_flux,wave=(wave,unit_wavelength))
            mf = []
            # plot photometric systems in different colors
            for system in systems:
                keep = (allsystems==system) & -iscolor
                if keep.sum():
                    # plot each photometric points separately, so that we could
                    # use it interactively. Label them all with a unique ID
                    # and make them pickable.
                    color = color_cycle.next()
                    #for i in range(sum(keep)):
                        #label = system if i==0 else '_nolegend_'
                        #pltlin,caplins,barlincs = pl.errorbar(wave[keep][i],flux[keep][i],yerr=e_flux[keep][i],fmt='o',label=label,ms=7,picker=5,color=color,**kwargs)
                        #pltlin.sed_index = indices[keep][i]
                        #caplins[0].sed_index = indices[keep][i]
                        #caplins[1].sed_index = indices[keep][i]
                        #barlincs[0].sed_index = indices[keep][i]
                    pl.errorbar(wave[keep],flux[keep],yerr=e_flux[keep],fmt='o',label=system,ms=7,color=color,**kwargs)
                    mf.append(flux[keep])
            if keep.sum():
                pl.ylabel(conversions.unit2texlabel(unit_flux,full=True))
            pl.xlabel('Wavelength [{0}]'.format(conversions.unit2texlabel(unit_wavelength)))
            #-- scale y-axis (sometimes necessary for data with huge errorbars)
            mf = np.log10(np.hstack(mf))
            lmin,lmax = np.nanmin(mf),np.nanmax(mf)
            lrange = np.abs(lmin-lmax)
            pl.ylim(10**(lmin-0.1*lrange),10**(lmax+0.1*lrange))
        else:
            pl.gca().set_color_cycle(color_cycle)
            names = []
            start_index = 1
            for system in systems:
                keep = (allsystems==system) & iscolor
                if keep.sum():
                    pl.errorbar(range(start_index,start_index+keep.sum()),flux[keep],yerr=e_flux[keep],fmt='o',label=system,ms=7,**kwargs)
                    names += [ph.split('.')[1] for ph in photbands[keep]]
                start_index += keep.sum()
            pl.xticks(range(1,len(names)+1),names,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Index')
        leg = pl.legend(prop=dict(size='small'),loc='best',fancybox=True)
        leg.get_frame().set_alpha(0.5)
        pl.grid()
        pl.title(self.ID)
            
    
    @standalone_figure
    def plot_sed(self,colors=False,mtype='igrid_search',plot_deredded=False,
            plot_unselected=True,wave_units='AA',flux_units='erg/s/cm2',**kwargs):
        """
        Plot a fitted SED together with the data.
        
        Example usage:
        
        First set the SED:
        
        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()
        
        Then make the plots:
        
        >>> p = pl.figure()
        >>> p = pl.subplot(121)
        >>> mysed.plot_sed(colors=False)
        >>> p = pl.subplot(122)
        >>> mysed.plot_sed(colors=True)
        
        ]]include figure]]ivs_sed_builder_plot_sed_01.png]
        """
        
        def plot_sed_getcolors(master,color_dict=None):
            myphotbands = [iphotb.split('.')[1] for iphotb in master['photband'][master['color']]]
            if not myphotbands:  #-- If there are no colours none can be returned (added by joris 30-01-2012)
                return [],[],[],None
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
        
        #-- get the color cycle
        systems = np.array([system.split('.')[0] for system in self.master['photband']],str)
        set_systems = sorted(list(set(systems)))
        color_cycle = itertools.cycle([pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(set_systems))])
        
        #-- for plotting reasons, we translate every color to an integer
        for system in set_systems:
            color = color_cycle.next()
            keep = (systems==system) & (self.master['color']==colors)
            if not plot_unselected:
                keep = keep & self.master['include']
            #-- synthetic:
            if sum(keep) and mtype in self.results and 'synflux' in self.results[mtype]:
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                    y = self.results[mtype]['synflux'][1][keep]
                else:
                    x = self.results[mtype]['synflux'][0][keep]
                    y = self.results[mtype]['synflux'][1][keep]
                    #-- convert to correct units
                    y = conversions.convert('erg/s/cm2/AA',flux_units,y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)
                pl.plot(x,y,'x',ms=10,mew=2,alpha=0.75,color=color,**kwargs)
            #-- include:
            keep = (systems==system) & (self.master['color']==colors) & self.master['include']
            if sum(keep):
                if colors:
                    #-- translate every color to an integer
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    if mtype in self.results and 'synflux' in self.results[mtype]:
                        x = self.results[mtype]['synflux'][0][keep]
                    else:
                        x = self.master['cwave'][keep]
                    y = self.master['cmeas'][keep]
                    e_y = self.master['e_cmeas'][keep]
                    y,e_y = conversions.convert('erg/s/cm2/AA',flux_units,y,e_y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)    
                    
                    
                p = pl.errorbar(x,y,yerr=e_y,fmt='o',label=system,
                                capsize=10,ms=7,color=color,**kwargs)
            
            #-- exclude:
            label = np.any(keep) and '_nolegend_' or system
            keep = (systems==system) & (self.master['color']==colors) & -self.master['include']
            if sum(keep) and plot_unselected:
                if colors:
                    x,y,e_y,color_dict = plot_sed_getcolors(self.master[keep],color_dict)
                else:
                    x = self.results[mtype]['synflux'][0][keep]
                    if np.any(np.isnan(x)):
                        x = self.master['cwave'][keep]
                    y = self.master['cmeas'][keep]
                    e_y = self.master['e_cmeas'][keep]
                    y,e_y = conversions.convert('erg/s/cm2/AA',flux_units,y,e_y,wave=(x,'AA'))
                    x = conversions.convert('AA',wave_units,x)    
                pl.errorbar(x,y,yerr=e_y,fmt='o',label=label,
                            capsize=10,ms=7,mew=2,color=color,mfc='1',mec=color,**kwargs)
        
        #-- only set logarithmic scale if absolute fluxes are plotted
        #   and only plot the real model then
        if not colors:
            pl.gca().set_xscale('log',nonposx='clip')
            pl.gca().set_yscale('log',nonposy='clip')
            pl.gca().set_autoscale_on(False)
        
            #-- the model
            if mtype in self.results and 'model' in self.results[mtype]:
                wave,flux,flux_ur = self.results[mtype]['model']
                
                flux = conversions.convert('erg/s/cm2/AA',flux_units,flux,wave=(wave,'AA'))
                flux_ur = conversions.convert('erg/s/cm2/AA',flux_units,flux_ur,wave=(wave,'AA'))
                wave = conversions.convert('AA',wave_units,wave)    
                
                pl.plot(wave,flux,'r-',**kwargs)
                if plot_deredded:
                    pl.plot(wave,flux_ur,'k-',**kwargs)
            pl.ylabel(conversions.unit2texlabel(flux_units,full=True))
            pl.xlabel('wavelength [{0}]'.format(conversions.unit2texlabel(wave_units)))
        else:
            xlabels = color_dict.keys()
            xticks = [color_dict[key] for key in xlabels]
            pl.xticks(xticks,xlabels,rotation=90)
            pl.ylabel(r'Flux ratio')
            pl.xlabel('Color name')
            pl.xlim(min(xticks)-0.5,max(xticks)+0.5)
        
        pl.grid()
        if colors:
            leg = pl.legend(loc='best',prop=dict(size='x-small'))
        else:
            leg = pl.legend(loc='upper right',prop=dict(size='x-small'))
            leg.get_frame().set_alpha(0.5)
        loc = (0.45,0.05)
        if mtype in self.results and 'grid' in self.results[mtype]:
            teff = self.results[mtype]['grid']['teff'][-1]
            logg = self.results[mtype]['grid']['logg'][-1]
            ebv = self.results[mtype]['grid']['ebv'][-1]
            scale = self.results[mtype]['grid']['scale'][-1]
            angdiam = 2*conversions.convert('sr','mas',scale)
            try:
                teff2 = self.results[mtype]['CI']['teff-2']
                logg2 = self.results[mtype]['CI']['logg-2']
                radii = self.results[mtype]['CI']['rad-2']/self.results[mtype]['CI']['rad']
                pl.annotate('Teff=%i   %i K\nlogg=%.2f   %.2f cgs\nE(B-V)=%.3f mag\nr2/r1=%.2f\n$\Theta$=%.3g mas'%(teff,teff2,logg,logg2,ebv,radii,angdiam),
                        loc,xycoords='axes fraction')
            except:
                pl.annotate('Teff=%d K\nlogg=%.2f cgs\nE(B-V)=%.3f mag\n$\Theta$=%.3g mas'%(teff,logg,ebv,angdiam),
                        loc,xycoords='axes fraction')
        logger.info('Plotted SED as %s'%(colors and 'colors' or 'absolute fluxes'))
        
    @standalone_figure    
    def plot_chi2(self,colors=False,mtype='igrid_search'):
        """
        Plot chi2 statistic for every datapoint included in the fit.
        
        To plot the statistic from the included absolute values, set
        C{colors=False}.
        
        To plot the statistic from the included colors, set C{colors=True}.
        
        Example usage:
        
        First set the SED:
        
        >>> mysed = SED('HD180642')
        >>> mysed.load_fits()
        
        Then make the plots:
        
        >>> p = pl.figure()
        >>> p = pl.subplot(121)
        >>> mysed.plot_chi2(colors=False)
        >>> p = pl.subplot(122)
        >>> mysed.plot_chi2(colors=True)
        
        ]]include figure]]ivs_sed_builder_plot_chi2_01.png]
        
        @param colors: flag to distinguish between colors and absolute values
        @type colors: boolean
        """
        
        include_grid = self.master['include']
        systems = np.array([system.split('.')[0] for system in self.master['photband'][include_grid]],str)
        set_systems = sorted(list(set(systems)))
        color_cycle = itertools.cycle([pl.cm.spectral(i) for i in np.linspace(0, 1.0, len(set_systems))])
        if mtype in self.results:
            eff_waves,synflux,photbands = self.results[mtype]['synflux']
            chi2 = self.results[mtype]['chi2']
            for system in set_systems:
                color = color_cycle.next()
                keep = systems==system
                if sum(keep) and not colors:
                    try:
                        pl.loglog(eff_waves[include_grid][keep],chi2[include_grid][keep],'o',label=system,color=color)
                    except:
                        logger.critical('Plotting of CHI2 of absolute values failed')
                elif sum(keep) and colors:
                    pl.semilogy(range(len(eff_waves[include_grid][keep])),chi2[include_grid][keep],'o',label=system,color=color)
            pl.legend(loc='upper right',prop=dict(size='x-small'))
            pl.grid()
            pl.annotate('Total $\chi^2$ = %.1f'%(self.results[mtype]['grid']['chisq'][-1]),(0.59,0.120),xycoords='axes fraction',color='r')
            pl.annotate('Total Reduced $\chi^2$ = %0.2f'%(sum(chi2)),(0.59,0.075),xycoords='axes fraction',color='r')
            if 'factor' in self.results[mtype]:
                pl.annotate('Error scale = %.2f'%(np.sqrt(self.results[mtype]['factor'])),(0.59,0.030),xycoords='axes fraction',color='k')
            xlims = pl.xlim()
            pl.plot(xlims,[self.results[mtype]['grid']['chisq'][-1],self.results[mtype]['grid']['chisq'][-1]],'r-',lw=2)
            pl.xlim(xlims)
            pl.xlabel('wavelength [$\AA$]')
            pl.ylabel(r'Reduced ($\chi^2$)')
            logger.info('Plotted CHI2 of %s'%(colors and 'colors' or 'absolute fluxes'))
        else:
            logger.info('%s not in results, no plot could be made.'%(mtype))
        
        
    @standalone_figure    
    def plot_distance(self,mtype='igrid_search'):
        try:
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
        except KeyError:
            logger.info('No distance/reddening plotted due to KeyError.')
    
    @standalone_figure
    def plot_grid_model(self,ptype='prob'):
        """
        Grid of models
        """
        if 'spType' in self.info:
            pl.title(self.info['spType'])
        cutlogg = (self.results['igrid_search']['grid']['logg']<=4.4) & (self.results['igrid_search']['grid']['ci_red']<=0.95)
        (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
        (d,dprob) = self.results['igrid_search']['d']
        gal = self.info['galpos']
        
        n = 75000
        region = self.results['igrid_search']['grid']['ci_red']<0.95
        total_prob = 100-(1-self.results['igrid_search']['grid']['ci_red'][cutlogg][-n:])*d_prob_models*100
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
        
        
    @standalone_figure
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
    
    @standalone_figure
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
        if 'igrid_search' in self.results and 'd_mod' in self.results['igrid_search']:
            (d_models,d_prob_models,radii) = self.results['igrid_search']['d_mod']
            (d,dprob) = self.results['igrid_search']['d']
        else:
            d = np.linspace(0,1000,2)
            dprob = np.zeros(len(d))
        
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
    
    @standalone_figure
    def plot_finderchart(self,cmap_photometry=pl.cm.spectral,window_size=5.):
        """
        Size is x and y width in arcminutes
        """
        try:
            dec = 'jdedeg' in self.info and self.info['jdedeg'] or None
            ra = 'jradeg' in self.info and self.info['jradeg'] or None
            data,coords,size = mast.get_dss_image(self.ID,ra=ra,dec=dec)
            pl.imshow(data[::-1],extent=[-size[0]/2*60,size[0]/2*60,
                                        -size[1]/2*60,size[1]/2*60],cmap=pl.cm.RdGy_r)#Greys
            pl.xlim(-window_size/2.,+window_size/2.)
            pl.ylim(-window_size/2.,+window_size/2.)
            xlims,ylims = pl.xlim(),pl.ylim()
            keep_this = -self.master['color'] & (self.master['cmeas']>0)
            toplot = self.master[keep_this]
            systems = np.array([system.split('.')[0] for system in toplot['photband']],str)
            set_systems = sorted(list(set(systems)))
            color_cycle = itertools.cycle([cmap_photometry(j) for j in np.linspace(0, 1.0, len(set_systems))])    
            for system in set_systems:
                color = color_cycle.next()
                keep = systems==system
                if sum(keep):
                    pl.plot(toplot['_RAJ2000'][keep][0]*60,
                            toplot['_DEJ2000'][keep][0]*60,'x',label=system,
                            mew=2.5,ms=15,alpha=0.5,color=color)
            leg = pl.legend(numpoints=1,prop=dict(size='x-small'),loc='best',fancybox=True)
            leg.get_frame().set_alpha(0.75)
            pl.xlim(xlims)
            pl.ylim(ylims)
            pl.xlabel(r'Right ascension $\alpha$ [arcmin]')
            pl.ylabel(r'Declination $\delta$ [arcmin]')    
        except:
            logger.warning('No image found of %s'%(self.ID))
            pass
        
        if 'pm' in self.info:
            logger.info("Found proper motion info")
            ppm_ra,ppm_de = (self.info['pm']['pmRA'],self.info['pm']['epmRA']),(self.info['pm']['pmDE'],self.info['pm']['epmDE'])
            pl.annotate('',xy=(ppm_ra[0]/50.,ppm_de[0]/50.),
                  xycoords='data',xytext=(0,0),textcoords='data',color='red',
                  arrowprops=dict(facecolor='red', shrink=0.05),
                  horizontalalignment='right', verticalalignment='top')
            pl.annotate('pmRA: %.1f $\pm$ %.1f mas/yr\npmDE: %.1f $\pm$ %.1f mas/yr'%(ppm_ra[0],ppm_ra[1],ppm_de[0],ppm_de[1]),
                        xy=(0.05,0.25),xycoords='axes fraction',color='red')
            if 'igrid_search' in self.results and 'd' in self.results['igrid_search']:
                (d,dprob) = self.results['igrid_search']['d']
                max_distance = d[np.argmax(dprob)]
                e_max_distance = abs(max_distance - d[np.argmin(np.abs(dprob-0.5*max(dprob)))])
            elif 'plx' in self.info and 'v' in self.info['plx'] and 'e' in self.info['plx']:
                plx,eplx = self.info['plx']['v'],self.info['plx']['e']
                dist = 1000./ufloat((plx,eplx))
                max_distance,e_max_distance = dist.nominal_value,dist.std_dev()
            else:
                max_distance = 1000.
                e_max_distance = 100.
            
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
        pl.subplot(rows,cols,1);self.plot_grid(ptype='ci_red')
        pl.subplot(rows,cols,2);self.plot_grid(ptype='ebv')
        pl.subplot(rows,cols,3);self.plot_grid(ptype='z')
        pl.subplot(rows,cols,4);self.plot_distance()
        
        pl.subplot(3,2,3);self.plot_sed(colors=False)
        pl.subplot(3,2,5);self.plot_sed(colors=True)
        
        pl.subplot(rows,cols,7);self.plot_chi2(colors=False)
        pl.subplot(rows,cols,11);self.plot_chi2(colors=True)
        
        #pl.subplot(rows,cols,8);self.plot_grid_model(ptype='prob')
        #pl.subplot(rows,cols,12);self.plot_grid_model(ptype='radii')
        
        pl.figure(figsize=(12,12))
        pl.axes([0,0.0,1.0,0.5]);self.plot_MW_side()
        pl.axes([0,0.5,0.5,0.5]);self.plot_MW_top()
        pl.axes([0.5,0.5,0.5,0.5]);self.plot_finderchart()
    
    #}
    
    #{Input and output
    
    def save_photometry(self,photfile=None):
        """
        Save master photometry to a file.
        
        @param photfile: name of the photfile. Defaults to C{starname.phot}.
        @type photfile: str
        """
        #-- write to file
        if photfile is not None:
            self.photfile = photfile
        logger.info('Save photometry to file %s'%(self.photfile))
        #-- add some comments
        if self.ID:
            if not 'bibcode' in self.master.dtype.names:
                self.master = crossmatch.add_bibcodes(self.master)
            if not 'comments' in self.master.dtype.names:
                self.master = vizier.quality_check(self.master,self.ID)
        ascii.write_array(self.master,self.photfile,header=True,auto_width=True,use_float='%g',comments=['#'+json.dumps(self.info)])
    
    def load_photometry(self,photfile=None):
        """
        Load the contents of the photometry file to the master record.
        
        @param photfile: name of the photfile. Defaults to the value of C{self.photfile}.
        @type photfile: str
        """
        if photfile is not None:
            self.photfile = photfile
        logger.info('Load photometry from file %s'%(self.photfile))
        self.master,comments = ascii.read2recarray(self.photfile,return_comments=True)
        self.info = json.loads(comments[-3])
    
    
    
    def save_fits(self,filename=None,overwrite=True):
        """
        Save content of SED object to a FITS file. 
        
        The .fits file will contain the following extensions if they are present in the object:
           1) data (all previously collected photometric data = content of .phot file)
           2) model_igrid_search (full best model SED table from grid_search)
           3) igrid_search (CI from grid search = actual grid samples)
           4) synflux_igrid_search (integrated synthetic fluxes for best model from grid_search)
           5) model_imc (full best model SED table from monte carlo)
           6) imc (CI from monte carlo = actual monte carlo samples)
           7) synflux_imc (integrated synthetic fluxes for best model from monte carlo)
        
        @param filename: name of SED FITS file
        @type filename: string
        @param overwrite: overwrite old FITS file if true
        @type overwrite: boolean
        
        Example usage:
        
        >>> #mysed.save_fits()
        >>> #mysed.save_fits(filename='myname.fits')
        """
        if filename is None:
            filename = str(os.path.splitext(self.photfile)[0]+'.fits')
        if overwrite:
            if os.path.isfile(filename):
                os.remove(filename)
                logger.info('Old FITS file removed')
                
        #-- write primary header
        #prim_header = {}
        #for key in self.info:
            #if not (isinstance(self.info[key],float) or isinstance(self.info[key],str)):
                #continue
            #prim_header[key] = self.info[key]
        #fits.write_recarray(np.array([[0]]),filename,header_dict=prim_header,ext=0)
        
        #-- write master data
        master = self.master.copy()
        fits.write_recarray(master,filename,header_dict=dict(extname='data'))
        
        #-- write the rest
        for mtype in self.results:#['igrid_search','imc']:
            if mtype in self.results:
                eff_waves,synflux,photbands = self.results[mtype]['synflux']
                chi2 = self.results[mtype]['chi2']
                
                results_modeldict = dict(extname='model_'+mtype)
                results_griddict = dict(extname=mtype)
                keys = sorted(self.results[mtype])
                for key in keys:
                    if 'CI' in key:
                        for ikey in self.results[mtype][key]:
                            if '_l' not in ikey and '_u' not in ikey and ikey != 'chisq':
                                results_modeldict[ikey] = self.results[mtype][key][ikey]
                            results_griddict[ikey] = self.results[mtype][key][ikey]    
                    if key=='factor':
                        results_griddict[key] = self.results[mtype][key]
            
                fits.write_array(list(self.results[mtype]['model']),filename,
                                names=('wave','flux','dered_flux'),
                                units=('AA','erg/s/cm2/AA','erg/s/cm2/AA'),
                                header_dict=results_modeldict)
                if 'grid' in self.results[mtype]:
                    fits.write_recarray(self.results[mtype]['grid'],filename,header_dict=results_griddict)
                
                results = np.rec.fromarrays([synflux,eff_waves,chi2],dtype=[('synflux','f8'),('mod_eff_wave','f8'),('chi2','f8')])
                
                fits.write_recarray(results,filename,header_dict=dict(extname='synflux_'+mtype))
        
        logger.info('Results saved to FITS file: %s'%(filename))
        
    
    def load_fits(self,filename=None):
        """
        Load a previously made SED FITS file. Only works for SEDs saved with 
        the save_fits function after 14.06.2012.
        
        The .fits file can contain the following extensions:
           1) data (all previously collected photometric data = content of .phot file)
           2) model_igrid_search (full best model SED table from grid_search)
           3) igrid_search (CI from grid search = actual grid samples)
           4) synflux_igrid_search (integrated synthetic fluxes for best model from grid_search)
           5) model_imc (full best model SED table from monte carlo)
           6) imc (CI from monte carlo = actual monte carlo samples)
           7) synflux_imc (integrated synthetic fluxes for best model from monte carlo)
        
        @param filename: name of SED FITS file
        @type filename: string
        @rtype: bool
        @return: true if Fits file could be loaded
        """
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.fits'
        if not os.path.isfile(filename):
            logger.warning('No previous results saved to FITS file {:s}'.format(filename))
            return False
        ff = pyfits.open(filename)
        
        #-- observed photometry
        fields = ff['data'].columns.names
        master = np.rec.fromarrays([ff['data'].data.field(field) for field in fields],names=','.join(fields))
        #-- remove the whitespace in columns with strings added by fromarrays
        for i,name in enumerate(master.dtype.names):
            if master.dtype[i].str.count('S'):
                for j,element in enumerate(master[name]):
                    master[name][j] = element.strip()
        self.master = master
        
        #-- add dictionary that will contain the results
        if not hasattr(self,'results'):
            self.results = {}
        
        #-- grid search and MC results
        mtypes = [ext.header['extname'] for ext in ff[1:]]
        mtypes = list(set(mtypes) - set(['DATA']))
        for mtype in mtypes:
            mtype = mtype.lower().lstrip('synflux_').lstrip('model_')
            if not mtype in self.results:
                    self.results[mtype] = {}
            self.results[mtype]['model'] = np.array(ff['model_'+mtype].data.field('wave'),dtype='float64'),np.array(ff['model_'+mtype].data.field('flux'),dtype='float64'),np.array(ff['model_'+mtype].data.field('dered_flux'),dtype='float64')
            self.results[mtype]['chi2'] = np.array(ff['synflux_'+mtype].data.field('chi2'),dtype='float64')
            self.results[mtype]['synflux'] = np.array(ff['synflux_'+mtype].data.field('mod_eff_wave'),dtype='float64'),np.array(ff['synflux_'+mtype].data.field('synflux'),dtype='float64'),self.master['photband']
            
                
        for mtype in ['igrid_search','imc']:
            try:
                fields = ff[mtype].columns.names
                master = np.rec.fromarrays([ff[mtype].data.field(field) for field in fields],names=','.join(fields))
                if not mtype in self.results:
                    self.results[mtype] = {}
                self.results[mtype]['grid'] = master
                if 'factor' in ff[mtype].header:
                    self.results[mtype]['factor'] = np.array([ff[mtype].header['factor']])[0]
                
                headerkeys = ff[mtype].header.ascardlist().keys()
                for key in headerkeys[::-1]:
                    for badkey in ['xtension','bitpix','naxis','pcount','gcount','tfields','ttype','tform','tunit','factor','extname']:
                        if key.lower().count(badkey):
                            headerkeys.remove(key)
                            continue
                self.results[mtype]['CI'] = {}
                for key in headerkeys:
                    #-- we want to have the same types as the original: numpy.float64 --> np.array([..])[0]
                    self.results[mtype]['CI'][key.lower()] = np.array([ff[mtype].header[key]])[0]
            except KeyError:
                continue
        
        #self.results['igrid_search'] = {}
        #fields = ff['igrid_search'].columns.names
        #master = np.rec.fromarrays([ff['igrid_search'].data.field(field) for field in fields],names=','.join(fields))
        #self.results['igrid_search']['grid'] = master
        #self.results['igrid_search']['factor'] = ff['igrid_search'].header['factor']
        
        #self.results['model'] = ff[2].data.field('wave'),ff[2].data.field('flux'),ff[2].data.field('dered_flux')
        #self.results['chi2'] = ff[4].data.field('chi2')
        #self.results['synflux'] = ff[4].data.field('mod_eff_wave'),ff[4].data.field('synflux'),ff[1].data.field('photband')
                
        ff.close()
        
        logger.info('Loaded previous results from FITS')
        return filename
    
    def save_bibtex(self):
        """
        Convert the bibcodes in a phot file to a bibtex file.
        
        The first line in the bibtex file contains a \citet command citing
        all photometry.
        """
        filename = os.path.splitext(self.photfile)[0]+'.bib'
        crossmatch.make_bibtex(self.master,filename=filename)
    
    def save_summary(self,filename=None,CI_limit=None,method='igrid_search',chi2type='ci_red'):
        """
        Save a summary of the results to an ASCII file.
        """
        #-- open the summary file to write the results
        if filename is None:
            filename = os.path.splitext(self.photfile)[0]+'.sum'
        
        if CI_limit is None:
            CI_limit = self.CI_limit
        
        #-- gather the results:
        grid_results = self.results[method]['grid']
        start_CI = np.argmin(np.abs(grid_results[chi2type]-self.CI_limit))
        factor = self.results[method]['factor']
        names = ['factor','chi2_type','ci_limit']
        results = [factor,chi2type,CI_limit*100]
        for name in grid_results.dtype.names:
            lv,cv,uv = grid_results[name][start_CI:].min(),\
                       grid_results[name][-1],\
                       grid_results[name][start_CI:].max()
            names += [name+'_l',name,name+'_u']
            results += [lv,cv,uv]
        #-- write the used photometry to a file
        include_grid = self.master['include']
        photbands = ":".join(self.master[include_grid]['photband'])
        references = ",".join(self.master[include_grid]['bibcode'])
        used_photometry = photometry2str(self.master[include_grid],comment='#')
        used_atmosphere = '#'+model.defaults2str()+'\n'
        used_photbands = '#'+photbands+'\n'
        used_references = '#'+references
        comments = used_photometry+used_atmosphere+used_photbands+used_references
        
        contents = np.array([results]).T
        contents = np.rec.fromarrays(contents,names=names)
        ascii.write_array(contents,filename,auto_width=True,header=True,
                          comments=comments.split('\n'),mode='a',use_float='%g')
        
        logger.info('Saved summary to {0}'.format(filename))
        
        
    
    #}

class BinarySED(SED):
    
    def igrid_search(self,points=100000,teffrange=None,loggrange=None,ebvrange=None,
                          zrange=None,radiusrange=None,rvrange=None,vradrange=(0,0),
                          masses=None,compare=True,df=5,CI_limit=None,
                          primary_hottest=False, gr_diff=None,set_model=True,**kwargs):
        """
        Fit fundamental parameters using a (pre-integrated) grid search.
        
        If called consecutively, the ranges will be set to the CI_limit of previous
        estimations, unless set explicitly.
        
        If called for the first time, the ranges will be +/- np.inf by defaults,
        unless set explicitly.
        """
        if CI_limit is None or CI_limit > 1.0:
            CI_limit = self.CI_limit
        #-- set defaults limits
        ranges,rtype = self.generate_ranges(teffrange=teffrange,\
                             loggrange=loggrange,ebvrange=ebvrange,\
                             zrange=zrange,rvrange=rvrange,vradrange=vradrange)
        
        if masses is None:
            # if masses are not given it doesn`t make much sense to use a binary grid...
            logger.warning('Using igridsearch with type binary, but no masses are provided for the components! Using masses=(1,1)')
            masses = (1,1)
        
        #-- grid search on all include data: extract the best CHI2
        include_grid = self.master['include']
        logger.info('The following measurements are included in the fitting process:\n%s'%(photometry2str(self.master[include_grid])))
        
        #-- build the grid, run over the grid and calculate the CHI2
        teffs,loggs,ebvs,zs,radii = fit.generate_grid(self.master['photband'][include_grid],teffrange=teffrange,
                   loggrange=loggrange,ebvrange=ebvrange, zrange=zrange, radiusrange=radiusrange, masses=masses,
                   points=points,res=res, type=type,primary_hottest=primary_hottest, gr_diff=gr_diff) 
        chisqs,scales,escales,lumis = fit.igrid_search(self.master['cmeas'][include_grid],
                            self.master['e_cmeas'][include_grid],
                            self.master['photband'][include_grid],
                            teffs,loggs,ebvs,zs,radii,threads=threads,model_func=model.get_itable_multiple)
        grid_results = np.rec.fromarrays([teffs[:,0],loggs[:,0],ebvs[:,0],zs[:,0],radii[:,0],
                                         teffs[:,1],loggs[:,1],ebvs[:,1],zs[:,1],radii[:,1],
                                         chisqs,scales,escales,lumis],
                                         dtype=[('teff','f8'),('logg','f8'),('ebv','f8'),('z','f8'),('rad','f8'),
                                         ('teff-2','f8'),('logg-2','f8'),('ebv-2','f8'),('z-2','f8'),('rad-2','f8') ,
                                         ('chisq','f8'),('scale','f8'),('escale','f8'),('labs','f8')])
        
        #-- exclude failures
        failures = np.isnan(grid_results['chisq'])
        if sum(failures):
            logger.info('Excluded {0} failed results (nan)'.format(sum(failures)))
            grid_results = grid_results[-failures]
        
        #-- make room for chi2 statistics
        grid_results = mlab.rec_append_fields(grid_results, 'ci_raw', np.zeros(len(grid_results)))
        grid_results = mlab.rec_append_fields(grid_results, 'ci_red', np.zeros(len(grid_results)))
        
        #-- take the previous results into account if they exist:
        if not 'igrid_search' in self.results:
            self.results['igrid_search'] = {}
        elif 'grid' in self.results['igrid_search']:
            logger.info('Appending previous results ({:d}+{:d})'.format(len(self.results['igrid_search']['grid']),len(grid_results)))
            ex_names = grid_results.dtype.names
            ex_grid = np.rec.fromarrays([self.results['igrid_search']['grid'][exname] for exname in ex_names],
                                        names=ex_names)
            grid_results = np.hstack([ex_grid,grid_results])
        
        #-- inverse sort according to chisq: this means the best models are at
        #   the end (mainly for plotting reasons, so that the best models
        #   are on top).
        sa = np.argsort(grid_results['chisq'])[::-1]
        grid_results = grid_results[sa]
        
        #-- do the statistics
        #   degrees of freedom: teff,logg,E(B-V),theta,Z
        N = sum(include_grid)
        k = N-df
        if k<=0:
            logger.warning('Not enough data to compute CHI2: it will not make sense')
            k = 1
        #-- rescale if needed and compute confidence intervals
        factor = max(grid_results['chisq'][-1]/k,1)
        logger.warning('CHI2 rescaling factor equals %g'%(factor))
        CI_raw = scipy.stats.distributions.chi2.cdf(grid_results['chisq'],k)
        CI_red = scipy.stats.distributions.chi2.cdf(grid_results['chisq']/factor,k)
        
        #-- add the results to the record array and to the results dictionary
        grid_results['ci_raw'] = CI_raw
        grid_results['ci_red'] = CI_red
        self.results['igrid_search']['grid'] = grid_results
        self.results['igrid_search']['factor'] = factor
        
        #-- compute the confidence intervals
        self.calculate_confidence_intervals(mtype='igrid_search',chi2_type='red',CI_limit=CI_limit)
        
        #-- remember the best model
        if set_model: self.set_best_model(type=type)
    
    
    def set_best_model(self,mtype='igrid_search',law='fitzpatrick2004',type='single'):
        """
        Get reddenend and unreddened model
        """
        logger.info('Interpolating approximate full SED of best model')
        
        #-- synthetic flux
        include = self.master['include']
        synflux = np.zeros(len(self.master['photband']))
        keep = (self.master['cwave']<1.6e6) | np.isnan(self.master['cwave'])
        keep = keep & include
        
        if mtype in ['igrid_search']:
            scale = self.results[mtype]['CI']['scale']
            
            #-- get (approximated) reddened and unreddened model
            wave,flux = model.get_table_multiple(teff=(self.results[mtype]['CI']['teff'],self.results[mtype]['CI']['teff-2']),
                                    logg=(self.results[mtype]['CI']['logg'],self.results[mtype]['CI']['logg-2']),
                                    ebv=(self.results[mtype]['CI']['ebv'],self.results[mtype]['CI']['ebv-2']),
                                    radius=(self.results[mtype]['CI']['rad'],self.results[mtype]['CI']['rad-2']),
                                    law=law)
            wave_ur,flux_ur = model.get_table_multiple(teff=(self.results[mtype]['CI']['teff'],self.results[mtype]['CI']['teff-2']),
                                      logg=(self.results[mtype]['CI']['logg'],self.results[mtype]['CI']['logg-2']),
                                      ebv=(0,0),
                                      radius=(self.results[mtype]['CI']['rad'],self.results[mtype]['CI']['rad-2']),
                                      law=law)
            #-- get synthetic photometry
            synflux_,Labs = model.get_itable_multiple(teff=(self.results[mtype]['CI']['teff'],self.results[mtype]['CI']['teff-2']),
                              logg=(self.results[mtype]['CI']['logg'],self.results[mtype]['CI']['logg-2']),
                              ebv=(self.results[mtype]['CI']['ebv'],self.results[mtype]['CI']['ebv-2']),
                              z=(self.results[mtype]['CI']['z'],self.results[mtype]['CI']['z-2']),
                              radius=(self.results[mtype]['CI']['rad'],self.results[mtype]['CI']['rad-2']),
                              photbands=self.master['photband'][keep])
            flux,flux_ur = flux*scale,flux_ur*scale
                
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
            self.results[mtype]['model'] = wave,flux,flux_ur
            self.results[mtype]['synflux'] = eff_waves,synflux,self.master['photband']
            self.results[mtype]['chi2'] = chi2
    
class Calibrator(SED):
    """
    Convenience class for a photometric standard star or calibrator.
    """
    def __init__(self,ID=None,photfile=None,plx=None,load_fits=True,label='',library='calspec'):
        names,fits_files,phot_files = model.read_calibrator_info(library='ngsl')
        index = names.index(ID)
        #--retrieve fitsfile information
        fits_file = pyfits.open(fits_files[index])
        wave = fits_file[1].data.field('wavelength')
        flux = fits_file[1].data.field('flux')
        fits_file.close()
        #--photfile:
        photfile = phot_files[index]
        super(Calibrator,self).__init__(photfile=photfile,plx=plx,load_fits=load_fits,label=label)
        self.set_model(wave,flux)
        

class SampleSEDs(object):
    """
    Class representing a list of SEDs.
    """
    def __init__(self,targets,**kwargs):
        """
        Initialize a sample.
        
        This can be done either with a list of IDs or with a list of SED
        instances. The SED must exist! That is, for each ID, there must be a
        phot or FITS file.
        """
        #-- we don't want to load the FITS files by default because they
        #   can eat up memory
        kwargs.setdefault('load_fits',False)
        self.targets = []
        self.seds = []
        #-- create all SEDs
        for target in targets:
            #-- perhaps we gave an ID: in that case, create the SED instance
            if not isinstance(target,SED):
                mysed = SED(target,**kwargs)
                if not mysed.has_photfile():
                    raise ValueError("No phot file found for {}".format(target))
            #-- perhaps already an SED object: then do nothing
            else:
                mysed = target
            self.seds.append(mysed)
            self.targets.append(mysed.ID)
        logger.info("Loaded {} SEDs".format(len(self.seds)))
    
    def __iter__(self):
        """
        Allow iteration over the SED instances.
        """
        for sed in self.seds:
            yield sed
    
    def __len__(self):
        """
        The length of a SampleSEDs instance is the number SED instances.
        """
        return len(self.seds)
    
    def __getitem__(self,key):
        """
        Implements various ways to get individual seds.
        
        Allows integer indexing, slicing, indexing with integer and boolean arrays.
        """
        #-- via slicing
        if isinstance(key,slice):
            return SampleSEDs([self[ii] for ii in range(*key.indices(len(self)))])
        #-- via an integer
        elif isinstance(key,int):
            return self.seds[key]
        else:
            #-- try to make the input an array
            try:
                key = np.array(key)
            except:
                raise TypeError("Cannot use instance of type {} for indexing".format(type(key)))
            #-- integer array slicing
            if key.dtype==np.dtype(int):
                return SampleSEDs([self[ii] for ii in key])
            #-- boolean array slicing
            elif key.dtype==np.dtype(bool):
                return SampleSEDs([self[ii] for ii in range(len(key)) if key[ii]])
            #-- that's all I can come up with
            else:
                raise TypeError("Cannot use arrays of type {} for indexing".format(key.dtype))
    
    def summarize(self):
        #-- collect the names of all the different sources
        sources = {}
        for sed in self.seds:
            #-- collect the different sources
            these_sources = list(set(sed.master['source']))
            for source in these_sources:
                if not source in sources:
                    sources[source] = []
                #-- now for each source, collect the names of the passbands
                keep = sed.master['source']==source
                sources[source] += list(set(sed.master[keep]['photband']))
                sources[source] = sorted(list(set(sources[source])))
        #-- next step: for each source, create a record array where the columns
        #   are the different photbands. Fill in the values for the photbands
        #   for each target when possible.
        summary = []
        source_names = sorted(sources.keys())
        for source in source_names:
            #-- create the record array
            data = np.zeros((len(self.targets),2*len(sources[source])))
            #   but remember to have errors with the photbands
            names = []
            for photband in sources[source]:
                names += [photband,'e_'+photband]
            data = np.rec.fromarrays(data.T,names=names)
            #-- fill in the values
            for i,sed in enumerate(self.seds):
                for photband in sources[source]:
                    keep = (sed.master['source']==source) & (sed.master['photband']==photband)
                    #-- fill in nans for value and error when not present
                    if not sum(keep):
                        data[photband][i] = np.nan
                        data['e_'+photband][i] = np.nan
                    #-- otherwise give the first value you've found
                    else:
                        data[photband][i] = sed.master[keep]['cmeas'][0]
                        data['e_'+photband][i] = sed.master[keep]['e_cmeas'][0]
                        #if not sum(keep): logger.warning('multiple defined photband ({}) and source ({})'.format(photband,source))
            summary.append(data)
        dtypes = [(source,summary[i].dtype) for i,source in enumerate(source_names)]
        output = np.zeros(len(self.targets),dtype=np.dtype(dtypes))
        for i,name in enumerate(output.dtype.names):
            output[name] = summary[i]
        return output
    
    def get_data(self,source,photband,label=None):
        """
        Get all data on a particular passband from a particular source.
        
        If label is not None, synthetic flux from a model will be added to the
        columns.
        """
        records = []
        if label is not None:
            synflux_label = [('synflux','f8')]
        else:
            synflux_label = []
        for targetname,sed in zip(self.targets,self.seds):
            #-- for the source, collect the names of the passbands
            keep = (sed.master['source']==source) & (sed.master['photband']==photband)
            #-- is there a model SED for which we want to retrieve synthetic photometry?
            if label is not None:
                if sum(keep)==0:
                    synflux = [0]
                else:
                    synflux = [sed.results[label]['synflux'][1][keep][0]]
            else:
                synflux = []
            #-- if no data on the matter is present, put zero values everywhere
            if sum(keep)==0:
                records.append([targetname]+list(np.zeros(len(sed.master.dtype.names)))+synflux)
            else:
                records.append([targetname]+list(sed.master[keep][0])+synflux)
        #-- make the thing into a record array
        dtypes = np.dtype([('targetname','|S25')] + sed.master.dtype.descr + synflux_label)
        output = np.rec.fromrecords(records,names=dtypes.names)
        output = np.array(output,dtype=dtypes)
        return output
    
    def get_confidence_interval(self,parameter='teff',mtype='igrid_search'):
        values = np.zeros((len(self),3))
        for i,sed in enumerate(self):
            sed.load_fits()
            values[i,0] = sed.results[mtype]['CI'][parameter+'_l']
            values[i,1] = sed.results[mtype]['CI'][parameter]
            values[i,2] = sed.results[mtype]['CI'][parameter+'_u']
            sed.clear()
        return values
        



if __name__ == "__main__":
    import sys
    import doctest
    import pprint
    from ivs.aux import loggers
    
    if not sys.argv[1:]:
        doctest.testmod()
        pl.show()
    else:
        name = " ".join([string for string in sys.argv[1:] if not '=' in string])
        units = [string.split('=')[1] for string in sys.argv[1:] if 'units=' in string]
        if not units:
            units = 'erg/s/cm2/AA'
        else:
            units = units[0]
        logger = loggers.get_basic_logger("")
        mysed = SED(name)
        pprint.PrettyPrinter(indent=4).pprint(mysed.info)
        mysed.get_photometry(units=units)
        mysed.plot_data()
        pl.show()
        answer = raw_input('Keep photometry file %s (y/N)'%(mysed.photfile))
        if not 'y' in answer.lower():
            os.unlink(mysed.photfile)
            logger.info('Removed %s'%(mysed.photfile))
    raise SystemExit
    
    #-- clean up
    if os.path.isfile('HD180642.fits'):
        os.remove('HD180642.fits')
    raise SystemExit
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
