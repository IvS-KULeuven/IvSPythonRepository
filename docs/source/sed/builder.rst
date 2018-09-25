To construct an SED, use the SED class. The functions defined in this module
are mainly convenience functions specifically for that class, but can be used
outside of the SED class if you know what you're doing.

1. Retrieving and plotting photometry of a target
=================================================

>>> mysed = SED('HD180642')
>>> mysed.get_photometry()
>>> mysed.plot_data()

and call Pylab's ``show`` function to show to the screen:

]]include figure]]ivs_sed_builder_example_photometry.png]

Catch IndexErrors and TypeErrors in case no photometry is found.

You can give a **search radius** to ``get_photometry`` via the keyword ``radius``.
The default value is 10 arcseconds for stars dimmer than 6th magnitude, and 60
arcseconds for brighter stars. The best value of course depends on the density
of the field.

>>> mysed.get_photometry(radius=5.)

If your star's name is not recognised by any catalog, you can give coordinates
to look for photometry. In that case, the ID of the star given in the ``SED``
command will not be used to search photometry (only to save the phot file):

>>> mysed.get_photometry(ra=289.31167983,dec=1.05941685)

Note that ``ra`` and ``dec`` are given in **degrees**.

You best **switch on the logger** (see L{ivs.aux.loggers.get_basic_logger}) to see the progress:
sometimes, access to catalogs can take a long time (the GATOR sources are
typically slow). If one of the ``gator``, ``vizier`` or ``gcpd`` is impossibly slow
or the site is down, you can **include/exclude these sources** via the keywords
``include`` or ``exclude``, which take a list of strings (choose from ``gator``,
``vizier`` and/or ``gcpd``). For ViZieR, there is an extra option to change to
another mirror site via

>>> vizier.change_mirror()

The L{vizier.change_mirror} function cycles through all the mirrors continuously,
so sooner or later you will end up with the default one and repeat the cycle.

>>> mysed.get_photometry(exclude=['gator'])

The results will be written to the file **HD180642.phot**. An example content is:

.. literalinclude:: HD180642.phot

Once a .phot file is written and L{get_photometry} is called again for the same
target, the script will **not retrieve the photometry from the internet again**,
but will use the contents of the file instead. The purpose is minimizing network
traffic and maximizing speed. If you want to refresh the search, simply manually
delete the .phot file or set ``force=True`` when calling L{get_photometry}.

The content of the .phot file is most easily read using the L{ivs.inout.ascii.read2recarray}
function. Be careful, as it contains both absolute fluxes as flux ratios.

>>> data = ascii.read2recarray('HD180642.phot')

Notice that in the ``.phot`` files, also a ``comment`` column is added. You can
find translation of some of the flags here (i.e. upper limit, extended source etc..),
or sometimes just additional remarks on variability etc. Not all catalogs have
this feature implemented, so you are still responsible yourself for checking
the quality of the photometry.

The references to each source are given in the ``bibtex`` column. Simply call

>>> mysed.save_bibtex()

to convert those bibcodes to a ``.bib`` file.

Using L{SED.plot_MW_side} and L{SED.plot_MW_top}, you can make a picture of where
your star is located with respect to the Milky Way and the Sun. With L{SED.plot_finderchart},
you can check the location of your photometry, and also see if proper motions
etc are available.

2. Where/what is my target?
===========================

To give you some visual information on the target, the following plotting
procedure might be of some help.

To check whether the downloaded photometry is really belonging to the target,
instead of some neighbouring star (don't forget to set ``radius`` when looking
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

3. SED fitting using a grid based approach
==========================================

3.1 Single star
---------------

We make an SED of HD180642 by simply **exploring a whole grid of Kurucz models**
(constructed via L{fit.generate_grid}, iterated over via L{fit.igrid_search} and
evaluated with L{fit.stat_chi2}). The model with the best parameters is picked
out and we make a full SED with those parameters.

>>> mysed = SED('HD180642')
>>> mysed.get_photometry()

Now we have collected **all fluxes and colors**, but we do not want to use them
all: first, **fluxes and colors are not independent**, so you probably want to use
either only absolute fluxes or only colors (plus perhaps one absolute flux per
system to compute the angular diameter) (see L{SED.set_photometry_scheme}). Second,
**some photometry is better suited for fitting an SED than others**; typically IR
photometry does not add much to the fitting of massive stars, or it can be
contaminated with circumstellar material. Third, **some photometry is not so
reliable**, i.e. the measurements can have large systematic uncertainties due to
e.g. calibration effects, which are typically not included in the error bars.

Currently, four standard schemes are implemented, which you can set via L{SED.set_photometry_scheme}:

1. ``absolute``: use only absolute values
2. ``colors``: use only colors (no angular diameter values calculated)
3. ``combo``: use all colors and one absolute value per photometric system
4. ``irfm``: (infrared flux method) use colors for wavelengths shorter than
   infrared wavelengths, and absolute values for systems in the infrared. The
   infrared is defined as wavelength longer than 1 micron, but this can be
   customized with the keyword ``infrared=(value,unit)`` in
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

You can automatically make plots of (most plotting functions take ``colors=True/False``
as an argument so you can make e.g. the 'color' SED and 'absolute value' SED):

1. the grid (see L{SED.plot_grid})
2. the SED (see L{SED.plot_sed})
3. the fit statistics (see L{SED.plot_chi2})

To change the grid, load the L{ivs.sed.model} module and call
L{ivs.sed.model.set_defaults} with appropriate arguments. See that module for
conventions on the grid structure when adding custom grids.

To add arrays manually, i.e. not from the predefined set of internet catalogs,
use the L{SED.add_photometry_fromarrays} function.

**Warning**: Be careful when interpreting the Chi2 results. In order to always have a
solution, the chi2 is rescaled so that the minimum equals 1, in the case the
probability of the best chi2-model is zero. The Chi2 rescaling factor I{f} mimicks
a rescaling of all errorbars with a factor I{sqrt(f)}, and does not discriminate
between systems (i.e., **all** errors are blown up). If the errorbars are
underestimated, it could be that the rescaling factor is also wrong, which means
that the true probability region can be larger or smaller!

3.2 Binary star - ### W.I.P ###
-------------------------------

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

3.3 Saving SED fits
-------------------

You can save all the data to a multi-extension FITS file via

>>> mysed.save_fits()

This FITS file then contains all **measurements** (it includes the .phot file),
the **resulting SED**, the **confidence intervals of all parameters** and also the
**whole fitted grid**: in the above case, the extensions of the FITS file contain
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


3.4 Loading SED fits  ### BROKEN ###
------------------------------------

Unfortunately this is not yet working properly!

Once saved, you can load the contents of the FITS file again into an SED object
via

>>> mysed = SED('HD180642')
>>> mysed.load_fits()

and then you can build all the plots again easily. You can of course use the
predefined plotting scripts to start a plot, and then later on change the
properties of the labels, legends etc... for higher quality plots or to better
suit your needs.

4. Accessing the best fitting full SED model
============================================

You can access the full SED model that matches the parameters found by the
fitting routine via:

>>> wavelength,flux,deredded_flux = mysed.get_best_model()

Note that this model is retrieved after fitting, and was not in any way used
during the fitting. As a consequence, there could be small differences between
synthetic photometry calculated from this returned model and the synthetic
fluxes stored in ``mysed.results['igrid_search']['synflux']``, which is the
synthetic photometry coming from the interpolation of the grid of pre-interpolated
photometry. See the documentation of L{SED.get_model} for more information.

5. Radii, distances and luminosities
====================================

5.1. Relations between quantities
---------------------------------

Most SED grids don't have the radius as a tunable model parameter. The scale
factor, which is available for all fitted models in the grid when at least one
absolute photometric point is included, is directly propertional to the angular
diameter. The following relations hold::

>> distance = radius / np.sqrt(scale)
>> radius = distance * np.sqrt(scale)

Where ``radius`` and ``distance`` have equal units. The true absolute luminosity
(solar units) is related to the absolute luminosity from the SED (solar units)
via the radius (solar units)::

>> L_abs_true = L_abs * radius**2

Finally, the angular diameter can be computed via::

>> 2*conversions.convert('sr','mas',scale)

5.2. Seismic constraints
------------------------

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

5.3. Parallaxes
---------------

From the parallax, you can get an estimate of the distance. This is, however,
dependent on some prior assumptions such as the shape of the galaxy and the
distribution of stars. To estimate the probability density value due to
a measured parallax of a star at particular distance, you can call L{distance.distprob}::

>> gal_latitude = 0.5
>> plx = 3.14,0.5
>> d_prob = distance.distprob(d,gal_latitude,plx)

Using this distance, you can get an estimate for the radius of your star, and
thus also the absolute luminosity.

5.4: Reddening constraints
--------------------------

An estimate of the reddening of a star at a particular distance may be obtained
with the L{extinctionmodels.findext} function. There are currently three
reddening maps available: Drimmel, Marshall and Arenou::

>> lng,lat = 120.,-0.5
>> dist = 100. # pc
>> Rv = 3.1
>> EBV = extinctionmodels.findext(lng,lat,model='drimmel',distance=dist)/Rv
>> EBV = extinctionmodels.findext(lng,lat,model='marshall',distance=dist)/Rv
>> EBV = extinctionmodels.findext(lng,lat,model='arenou',distance=dist)/Rv

Reference/API
=============

.. automodule:: sed.builder
    :members:
