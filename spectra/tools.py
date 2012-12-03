"""
Manipulate or extract information from spectra.

Subsection 1. Example usage
===========================

In this example, we:

    1. retrieve a synthetic spectrum
    2. apply a doppler shift of 20 km/s
    3. rotationally broaden the synthetic spectrum according to a vsini of 66 km/s,
    a solar-like limb darkening, and an instrumental profile with FWHM=0.25 angstrom
    4. add noise to the synthetic spectrum
    5. calculate the vsini with the Fourier Transform Method and check the influence
    of the limb darkening coefficients

Retrieve a synthetic spectrum with effective temperature and surface gravity
resembling a main sequence star of spectral type A0. For this purpose, we do
not need the whole spectrum but just cut out a piece

>>> from ivs.spectra import model
>>> wave = np.linspace(4570,4574,300)
>>> wave,flux,cont = model.get_table(teff=10000,logg=4.0,wave=wave)
>>> clam = wave[np.argmin(flux)]

Apply a doppler shift of 20 km/s (the star travels away from us)

>>> wave_ = doppler_shift(wave,20.)
>>> clam_shift = doppler_shift(clam,20.)

Rotationally broaden the spectrum and convolve with an instrumental profile

>>> wave_,flux_ = rotational_broadening(wave_,flux,66.,stepr=-1,fwhm=0.25,stepi=-1,epsilon=0.6)

Add some noise

>>> fluxn_ = flux_ + np.random.normal(size=len(flux_),scale=0.01)

Calculate the vsini with the Fourier transform method. To compare the results,
first compute the FT of the synthetic broadened spectrum without noise:

>>> pergram,minima,vsinis,error = vsini(wave_,flux_,epsilon=0.6,clam=clam_shift,window=(4571,4573.5),df=1e-4)

Make a plot of what we already have:

>>> p = pl.figure()
>>> p = pl.subplot(121)
>>> p = pl.plot(wave,flux,'k-',label='Original template')
>>> p = pl.plot(wave_,fluxn_,'b-',label='Spectrum + noise')
>>> p = pl.plot(wave_,flux_,'r-',lw=2,label='Broadened')
>>> p = pl.legend(loc='lower right',prop=dict(size='small'))

Set the color cycle of the Fourier Transform plot to spectral

>>> p = pl.subplot(122)
>>> color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, 10)]
>>> p = pl.gca().set_color_cycle(color_cycle)
>>> p = pl.plot(pergram[0],pergram[1],'r-',lw=2,label='Correct')
>>> p = pl.gca().set_yscale('log')

Now compute the vsini of the noisy spectrum, assuming different limb darkening
parameters

>>> for epsilon in np.linspace(0.0,1.0,10):
...   pergram,minima,vsinis,error = vsini(wave_,fluxn_,epsilon=epsilon,clam=clam_shift,window=(4571,4573.5),df=1e-4)
...   p = pl.plot(pergram[0],pergram[1],label='$\epsilon$=%.2f: vsini = %.1f km/s'%(epsilon,vsinis[0]))

Set the xticks to vsini values for clarity:

>>> xticks = np.array([25,35,50.,66,80,100,500][::-1])
>>> p = pl.xticks(1/xticks,['%.0f'%(i) for i in xticks])
>>> p = pl.xlim(0,0.04);p = pl.grid()
>>> p = pl.ylim(5e-4,1)
>>> p = pl.legend(loc='lower left',prop=dict(size='small'))

]include figure]]ivs_spectra_tools_vsini01.png]

]include figure]]ivs_spectra_tools_vsini02.png]

"""
from ivs.spectra import pyrotin4
import numpy as np
import logging
from numpy import pi,sin,cos,sqrt
import scipy.stats
from ivs.timeseries import pergrams
from ivs.units import conversions
from ivs.units import constants

from ivs import config
from ivs.io import fits
from scipy.interpolate import interp1d

logger = logging.getLogger("SPEC.TOOLS")

def doppler_shift(wave,vrad,vrad_units='km/s',flux=None):
    """
    Shift a spectrum with towards the red or blue side with some radial velocity.
    
    You can give units with the extra keywords C{vrad_units} (units of
    wavelengths are not important). The shifted wavelengths will be in the same
    units as the input wave array.
    
    If units are not supplied, the radial velocity is assumed to be in km/s.
    
    If you want to apply a barycentric (or orbital) correction, you'd probabl
    want to reverse the sign of the radial velocity!
    
    When the keyword C{flux} is set, the spectrum will be interpolated onto
    the original wavelength grid (so the original wavelength grid will not
    change). When the keyword C{flux} is not set, the wavelength array will be
    changed (but the fluxes not, obviously).
    
    When C{flux} is set, fluxes will be returned.
    
    When C{flux} is not set, wavelengths will be returned.
    
    Example usage: shift a spectrum to the blue ('left') with 20 km/s.
    
    >>> wave = np.linspace(3000,8000,1000)
    >>> wave_shift1 = doppler_shift(wave,20.)
    >>> wave_shift2 = doppler_shift(wave,20000.,vrad_units='m/s')
    >>> print(wave_shift1[0],wave_shift1[-1])
    (3000.200138457119, 8000.5337025523177)
    >>> print(wave_shift2[0],wave_shift2[-1])
    (3000.200138457119, 8000.5337025523177)
    
    @param wave: wavelength array
    @type wave: ndarray
    @param vrad: radial velocity (negative shifts towards blue, positive towards red)
    @type vrad: float (units: km/s) or tuple (float,'units')
    @param vrad_units: units of radial velocity (default: km/s)
    @type vrad_units: str (interpretable for C{units.conversions.convert})
    @return: shifted wavelength array/shifted flux
    @rtype: ndarray
    """ 
    cc = constants.cc
    cc = conversions.convert('m/s',vrad_units,cc)
    wave_out = wave * (1+vrad/cc)
    if flux is not None:
        flux = np.interp(wave,wave_out,flux)
        return flux
    else:
        return wave_out

def vsini(wave,flux,epsilon=0.6,clam=None,window=None,**kwargs):
    """
    Deterimine vsini of an observed spectrum via the Fourier transform method.
    
    According to Simon-Diaz (2006) and Carroll (1933):
    
    vsini = 0.660 * c/ (lambda * f1)
    
    But more general (see Reiners 2001, Dravins 1990)
    
    vsini = q1 * c/ (lambda*f1)
    
    where f1 is the first minimum of the Fourier transform.
    
    The error is estimated as the Rayleigh limit of the Fourier Transform
    
    Example usage and tests: Generate some data. We need a central wavelength (A),
    the speed of light in angstrom/s, limb darkening coefficients and test
    vsinis:
    
    >>> clam = 4480.
    >>> c = conversions.convert('m/s','A/s',constants.cc)
    >>> epsilons = np.linspace(0.,1.0,10)
    >>> vsinis = np.linspace(50,300,10)
    
    We analytically compute the shape of the Fourier transform in the following
    domain (and need the C{scipy.special.j1} for this)
    
    >>> x = np.linspace(0,30,1000)[1:]
    >>> from scipy.special import j1
    
    Keep track of the calculated and predicted q1 values:
    
    >>> q1s = np.zeros((len(epsilons),len(vsinis)))
    >>> q1s_pred = np.zeros((len(epsilons),len(vsinis)))
    
    Start a figure and set the color cycle
    
    >>> p= pl.figure()
    >>> p=pl.subplot(131)
    >>> color_cycle = [pl.cm.spectral(j) for j in np.linspace(0, 1.0, len(epsilons))]
    >>> p = pl.gca().set_color_cycle(color_cycle)
    >>> p=pl.subplot(133);p=pl.title('Broadening kernel')
    >>> p = pl.gca().set_color_cycle(color_cycle)
    
    Now run over all epsilons and vsinis and determine the q1 constant:
    
    >>> for j,epsilon in enumerate(epsilons):
    ...    for i,vsini in enumerate(vsinis):
    ...       vsini = conversions.convert('km/s','A/s',vsini)
    ...       delta = clam*vsini/c
    ...       lambdas = np.linspace(-5.,+5.,20000)
    ...       #-- analytical rotational broadening kernel
    ...       y = 1-(lambdas/delta)**2
    ...       G = (2*(1-epsilon)*np.sqrt(y)+pi*epsilon/2.*y)/(pi*delta*(1-epsilon/3))
    ...       lambdas = lambdas[-np.isnan(G)]
    ...       G = G[-np.isnan(G)]
    ...       G /= max(G)
    ...       #-- analytical FT of rotational broadening kernel
    ...       g = 2. / (x*(1-epsilon/3.)) * ( (1-epsilon)*j1(x) +  epsilon* (sin(x)/x**2 - cos(x)/x))
    ...       #-- numerical FT of rotational broadening kernel
    ...       sigma,g_ = pergrams.deeming(lambdas-lambdas[0],G,fn=2e0,df=1e-3,norm='power')
    ...       myx = 2*pi*delta*sigma
    ...       #-- get the minima
    ...       rise = np.diff(g_[1:])>=0
    ...       fall = np.diff(g_[:-1])<=0
    ...       minima = sigma[1:-1][rise & fall]
    ...       minvals = g_[1:-1][rise & fall]
    ...       q1s[j,i] =  vsini / (c/clam/minima[0])
    ...       q1s_pred[j,i] = 0.610 + 0.062*epsilon + 0.027*epsilon**2 + 0.012*epsilon**3 + 0.004*epsilon**4
    ...    p= pl.subplot(131)
    ...    p= pl.plot(vsinis,q1s[j],'o',label='$\epsilon$=%.2f'%(epsilon));pl.gca()._get_lines.count -= 1
    ...    p= pl.plot(vsinis,q1s_pred[j],'-')
    ...    p= pl.subplot(133)
    ...    p= pl.plot(lambdas,G,'-')
    
    And plot the results:
    
    >>> p= pl.subplot(131)
    >>> p= pl.xlabel('v sini [km/s]');p = pl.ylabel('q1')
    >>> p= pl.legend(prop=dict(size='small'))
    
    
    >>> p= pl.subplot(132);p=pl.title('Fourier transform')
    >>> p= pl.plot(x,g**2,'k-',label='Analytical FT')
    >>> p= pl.plot(myx,g_/max(g_),'r-',label='Computed FT')
    >>> p= pl.plot(minima*2*pi*delta,minvals/max(g_),'bo',label='Minima')
    >>> p= pl.legend(prop=dict(size='small'))
    >>> p= pl.gca().set_yscale('log')
    
    ]include figure]]ivs_spectra_tools_vsini_kernel.png]
    
    Extra keyword arguments are passed to L{pergrams.deeming}
    
    @param wave: wavelength array in Angstrom
    @type wave: ndarray
    @param flux: normalised flux of the profile
    @type flux: ndarray
    @rtype: (array,array),(array,array),array,float
    @return: periodogram (freqs in s/km), extrema (weird units), vsini values (km/s), error (km/s)
    """
    cc = conversions.convert('m/s','AA/s',constants.cc)
    #-- clip the wavelength and flux region if needed:
    if window is not None:
        keep = (window[0]<=wave) & (wave<=window[1])
        wave,flux = wave[keep],flux[keep]
    #-- what is the central wavelength? If not given, it's the middle of the
    #   wavelength range
    if clam is None: clam = ((wave[0]+wave[-1])/2)
    #-- take care of the limb darkening
    q1 = 0.610 + 0.062*epsilon + 0.027*epsilon**2 + 0.012*epsilon**3 + 0.004*epsilon**4
    #-- do the Fourier transform and normalise
    #flux = flux / (np.median(np.sort(flux)[-5:]))
    freqs,ampls = pergrams.deeming(wave,(1-flux),**kwargs)
    ampls = ampls/max(ampls)
    error = 1./wave.ptp()
    #-- get all the peaks
    rise = np.diff(ampls[1:])>=0
    fall = np.diff(ampls[:-1])<=0
    minima = freqs[1:-1][rise & fall]
    minvals = ampls[1:-1][rise & fall]
    #-- compute the vsini and convert to km/s
    freqs = freqs*clam/q1/cc
    freqs = conversions.convert('s/AA','s/km',freqs,wave=(clam,'AA'))
    vsini_values = cc/clam*q1/minima
    vsini_values = conversions.convert('AA/s','km/s',vsini_values)#,wave=(clam,'AA'))
    #-- determine the error as the rayleigh limit
    error = error*clam/q1/cc
    error = conversions.convert('s/AA','s/km',error,wave=(clam,'AA'))
    return (freqs,ampls),(minima,minvals),vsini_values,error




def rotational_broadening(wave_spec,flux_spec,vrot,fwhm=0.25,epsilon=0.6,
                         chard=None,stepr=0,stepi=0,alam0=None,alam1=None,
                         irel=0,cont=None):
    """
    Apply rotational broadening to a spectrum assuming a linear limb darkening
    law.
    
    This function is based on the ROTIN program. See Fortran file for
    explanations of parameters.
    
    Limb darkening law is linear, default value is epsilon=0.6
    
    Possibility to normalize as well by giving continuum in 'cont' parameter.
    
    Section 1. Parameters for rotational convolution 
    ================================================

    C{VROT}: v sin i (in km/s):
    
        -  if VROT=0 - rotational convolution is 
                 a) either not calculated,
                 b) or, if simultaneously FWHM is rather large
                 (vrot/c*lambda < FWHM/20.),
                 vrot is set to  FWHM/20*c/lambda;
        -  if VROT >0 but the previous condition b) applies, the
        value of VROT is changed as  in the previous case
        -  if VROT<0 - the value of abs(VROT) is used regardless of
        how small compared to FWHM it is
     
    C{CHARD}: characteristic scale of the variations of unconvolved stellar
    spectrum (basically, characteristic distance between two neighbouring
    wavelength points) - in A:
     
        - if =0 - program sets up default (0.01 A)
        
    C{STEPR}: wavelength step for evaluation rotational convolution;
     
        - if =0, the program sets up default (the wavelength
        interval corresponding to the rotational velocity
        devided by 3.)                           
        - if <0, convolved spectrum calculated on the original
        (detailed) SYNSPEC wavelength mesh


    Section 2. parameters for instrumental convolution
    ==================================================

    C{FWHM}: WARNING: this is not the full width at half maximum for Gaussian
    instrumental profile, but the sigma (FWHM = 2.3548 sigma).
    
    C{STEPI}: wavelength step for evaluating instrumental convolution
          - if =0, the program sets up default (FWHM/10.)
          - if <0, convolved spectrum calculated with the previous
          wavelength mesh:
          either the original (SYNSPEC) one if vrot=0,
          or the one used in rotational convolution (vrot > 0)


    Section 3. wavelength interval and normalization of spectra
    ===========================================================

    C{ALAM0}: initial wavelength
    C{ALAM1}: final wavelength
    C{IREL}: for =1 relative spectrum, =0 absolute spectrum
    
    @return: wavelength,flux
    @rtype: array, array
    """
    #-- set arguments
    if alam0 is None: alam0 = wave_spec[0]
    if alam1 is None: alam1 = wave_spec[-1]
    if cont is None: cont = (np.ones(1),np.ones(1))
    contw,contf = cont
    if chard is None:
        chard = np.diff(wave_spec).mean()
    
    #-- apply broadening
    w3,f3,ind = pyrotin4.pyrotin(wave_spec,flux_spec,contw,contf,
                  vrot,chard,stepr,fwhm,stepi,alam0,alam1,irel,epsilon)
    logger.info('ROTIN rot.broad. with vrot=%.3f (epsilon=%.2f)'%(vrot,epsilon))
    
    return w3[:ind],f3[:ind]

def combine(list_of_spectra,R=200.,lambda0=(950.,'AA'),lambdan=(3350.,'AA')):
    """
    Combine and weight-average spectra on a common wavelength grid.
    
    C{list_of_spectra} should be a list of lists/arrays. Each element in the
    main list should be (wavelength,flux,error).
    
    If you have FUSE fits files, use L{ivs.fits.read_fuse}.
    If you have IUE FITS files, use L{ivs.fits.read_iue}.
    
    After Peter Woitke.
    
    @param R: resolution
    @type R: float
    @param lambda0: start wavelength, unit
    @type lambda0: tuple (float,str)
    @param lambdan: end wavelength, unit
    @type lambdan: tuple (float,str)
    @return: binned spectrum (wavelengths,flux, error)
    @rtype: array, array, array
    """
    l0 = conversions.convert(lambda0[1],'AA',lambda0[0])
    ln = conversions.convert(lambdan[1],'AA',lambdan[0])
    #-- STEP 1: define wavelength bins
    Delta = np.log10(1.+1./R)
    x = np.arange(np.log10(l0),np.log10(ln)+Delta,Delta)
    x = 10**x
    lamc_j = 0.5*(np.roll(x,1)+x)

    #-- STEP 2: rebinning of data onto newly defined wavelength bins
    Ns = len(list_of_spectra)
    Nw = len(lamc_j)-1
    binned_fluxes = np.zeros((Ns,Nw))
    binned_errors = np.inf*np.ones((Ns,Nw))

    for snr,(wave,flux,err) in enumerate(list_of_spectra):
        wave0 = np.roll(wave,1)
        wave1 = np.roll(wave,-1)
        lam_i0_dc = 0.5*(wave0+wave)
        lam_i1_dc = 0.5*(wave1+wave)
        dlam_i = lam_i1_dc-lam_i0_dc
        
        for j in range(Nw):
            A = np.min(np.vstack([lamc_j[j+1]*np.ones(len(wave)),lam_i1_dc]),axis=0)
            B = np.max(np.vstack([lamc_j[j]*np.ones(len(wave)),lam_i0_dc]),axis=0)
            overlaps = scipy.stats.threshold(A-B,threshmin=0)
            norm = np.sum(overlaps)
            binned_fluxes[snr,j] = np.sum(flux*overlaps)/norm
            binned_errors[snr,j] = np.sqrt(np.sum((err*overlaps)**2))/norm
    
    #-- STEP 3: all available spectra sets are co-added, using the inverse
    #   square of the bin uncertainty as weight
    binned_fluxes[np.isnan(binned_fluxes)] = 0
    binned_errors[np.isnan(binned_errors)] = 1e300
    weights = 1./binned_errors**2
    
    totalflux = np.sum(weights*binned_fluxes,axis=0)/np.sum(weights,axis=0)
    totalerr = np.sqrt(np.sum((weights*binned_errors)**2,axis=0))/np.sum(weights,axis=0)
    totalspec = np.sum(binned_fluxes>0,axis=0)
    
    #-- that's it!
    return x[:-1],totalflux,totalerr,totalspec

def get_response(instrument='hermes'):
    """
    Returns the response curve of the given instrument. Up till now only a HERMES response cruve is
    available. This response curve is based on 25 spectra of the single sdB star Feige 66, and has
    a wavelenght range of 3800 to 8000 A.
    
    @param instrument: the instrument of which you want the response curve
    @type instrument: string
    """
        
    basedir = 'spectables/responses/'
    
    if instrument == 'hermes':
        basename = 'response_Feige66.fits'
        
    response = config.get_datafile(basedir,basename)
    
    wave, flux = fits.read_spectrum(response)
        
    return wave, flux

def remove_response(wave, flux, instrument='hermes'):
    """
    Divides the provided spectrum by the response curve of the given instrument. Up till now only
    a HERMES response curve is available.
    
    @param wave: wavelenght array
    @type wave: numpy array
    @param flux: flux array
    @type flux: numpy array
    @param instrument: the instrument of which you want the response curve
    @type instrument: string
    
    @return: the new spectrum (wave, flux)
    @rtype: (array, array)
    """
    
    #-- get the response curve
    wr, fr = get_response(instrument=instrument)
    fr_ = interp1d(wr, fr, kind='linear')
    
    #-- select usefull part of the spectrum
    flux = flux[(wave >= wr[0]) & (wave <= wr[-1])]
    wave = wave[(wave >= wr[0]) & (wave <= wr[-1])]
    
    flux = flux / fr_(wave)
    
    return wave, flux
    

if __name__=="__main__":
    import doctest
    import pylab as pl
    doctest.testmod()
    pl.show()
