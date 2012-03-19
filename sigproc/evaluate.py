"""
Evaluate model functions.

Most functions accept a domain of interest and a set of parameters, and return
a model function (e.g., a sine, gaussian etc...). Parameters are usually a
record array containing named fields for the parameters.

Some nonlinear optimizers expect and return parameters as a 1D array of the
independent parameters. Therefore, all functions also accept these 1D numpy
arrays, and they are converted by L{check_input} to record arrays. The latter
function actually accepts both record arrays and 1D arrays, and acts as a
translator between the two parameter representations.
"""
import logging
import functools

import numpy as np
from numpy import pi,cos,sin,sqrt,tan,arctan
from scipy.interpolate import splev
from scipy.special import erf,jn
from ivs.timeseries import keplerorbit

logger = logging.getLogger('SIGPROC.EVAL')

def check_input(fctn):
    """
    Decorator to check input of evaluating model, and transforming input from
    record array to 1D numpy array for optimization, and vice versa.
    """
    @functools.wraps(fctn)
    def check(*args,**kwargs):
        args = list(args)
        parameters = np.asarray(args[1])
        if not parameters.dtype.names:
            parameters = globals()[fctn.__name__+'_preppars'](parameters)
        args[1] = parameters
        return fctn(*args,**kwargs)
    return check

#{ Timeseries

def phasediagram(time, signal, nu0, D=0, t0=None, return_indices=False,
                 chronological=False):
    """
    Construct a phasediagram, using frequency nu0.
    
    Possibility to include a linear frequency shift D.
    
    If needed, the indices can be returned. To 'unphase', just do:
    
    Example usage:
    
    >>> original = np.random.normal(size=10)
    >>> indices = np.argsort(original)
    >>> inv_indices = np.argsort(indices)
    >>> all(original == np.take(np.take(original,indices),inv_indices))
    True
    
    Optionally, the zero time point can be given, it will be subtracted from all
    time points
    
    Joris De Ridder, Pieter Degroote
    
    @param time: time points [0..Ntime-1]
    @type time: ndarray
    @param signal: observed data points [0..Ntime-1]
    @type signal: ndarray
    @param nu0: frequency to compute the phasediagram with (inverse unit of time)
    @type nu0: float
    @param t0: zero time point, if None: t0 = time[0]
    @type t0: float
    @param D: frequency shift
    @type D: float
    @param chronological: set to True if you want to have a list of every phase
    separately
    @type chronological: boolean
    @param return_indices: flag to return indices
    @type return_indices: boolean
    @return: phase points (sorted), corresponding observations
    @rtype: ndarray,ndarray(, ndarray)
    """
    if (t0 is None): t0 = 0.
    phase = np.fmod(nu0 * (time-t0) + D/2. * (time-t0)**2,1.0)
    phase = np.where(phase<0,phase+1,phase)
    
    #-- keep track of the phase number, and prepare lists to add the points
    if chronological:
        chainnr = np.floor((time-t0)*nu0)
        myphase = [[]]
        mysignl = [[]]
        currentphase = chainnr[0]
        #-- divide the signal into separate phase bins
        for i in xrange(len(signal)):
            if chainnr[i]==currentphase:
                myphase[-1].append(phase[i])
                mysignl[-1].append(signal[i])
            else:
                myphase[-1] = np.array(myphase[-1])
                mysignl[-1] = np.array(mysignl[-1])
                currentphase = chainnr[i]
                myphase.append([phase[i]])
                mysignl.append([signal[i]])
        myphase[-1] = np.array(myphase[-1])
        mysignl[-1] = np.array(mysignl[-1])
    #-- if we 
    else:
        indices = np.argsort(phase)
    
    #-- possibly only return phase and signal
    if not return_indices and not chronological:
        return phase[indices], signal[indices]
    #-- maybe return phase,signal and indices
    elif not chronological:
        return phase[indices], signal[indices], indices
    #-- or return seperate phase bins
    else:
        return myphase,mysignl



@check_input
def sine(times, parameters):
    """
    Creates a harmonic function based on given parameters.
    
    Parameter fields: C{const, ampl, freq, phase}.
    
    This harmonic function is of the form
    
    f(p1,...,pm;x1,...,x_n) = S{sum} c_i + a_i sin(2 S{pi} (f_i+ S{phi}_i))
    
    where p_i are parameters and x_i are observation times.
    
    Parameters can be given as a record array containing columns 'ampl', 'freq'
    and 'phase', and optionally also 'const' (the latter array is summed up so
    to reduce it to one number).
    
    C{pars = np.rec.fromarrays([consts,ampls,freqs,phases],names=('const','ampl','freq','phase'))}
    
    Parameters can also be given as normal 1D numpy array. In that case, it will
    be transformed into a record array by the C{sine_preppars} function. The
    array you give must be 1D, optionally contain a constant as the first
    parameter, followed by 3 three numbers per frequency:
        
    C{pars = [const, ampl1, freq1, phase1, ampl2, freq2, phase2...]}
    
    Example usage: We construct a synthetic signal with two frequencies.
    
    >>> times = np.linspace(0,100,1000)
    >>> const = 4.
    >>> ampls = [0.01,0.004]
    >>> freqs = [0.1,0.234]
    >>> phases = [0.123,0.234]
    
    The signal can be made via a record array
    
    >>> parameters = np.rec.fromarrays([[const,0],ampls,freqs,phases],names = ('const','ampl','freq','phase'))
    >>> mysine = sine(times,parameters)
    
    or a normal 1D numpy array
    
    >>> parameters = np.hstack([const,np.ravel(np.column_stack([ampls,freqs,phases]))])
    >>> mysine = sine(times,parameters)
    
    Of course, the latter is easier if you have your parameters in a list:
    
    >>> freq1 = [0.01,0.1,0.123]
    >>> freq2 = [0.004,0.234,0.234]
    >>> parameters = np.hstack([freq1,freq2])
    >>> mysine = sine(times,parameters)
    
    @param times: observation times
    @type times: numpy array
    @param parameters: record array containing amplitudes ('ampl'), frequencies ('freq')
    phases ('phase') and optionally constants ('const'), or 1D numpy array (see above).
    B{Warning:} record array should have shape (n,)!
    @type parameters: record array or 1D array
    @return: sine signal (same shape as C{times})
    @rtype: array
    """
    if 'const' in parameters.dtype.names:
        total_fit = parameters['const'].sum()
    else:
        total_fit = 0.
    for i in xrange(len(parameters)):
        total_fit += parameters['ampl'][i]*sin(2*pi*(parameters['freq'][i]*times+parameters['phase'][i]))
    return total_fit






def sine_freqshift(times,parameters,t0=None):
    """
    Creates a sine function with a linear frequency shift.
    
    Parameter fields: C{const, ampl, freq, phase, D}.
    
    Similar to C{sine}, but with extra column 'D', which is the linear frequency
    shift parameter.
    
    Only accepts record arrays as parameters (for the moment).
    
    @param times: observation times
    @type times: numpy array
    @param parameters: record array containing amplitudes ('ampl'), frequencies ('freq'),
    phases ('phase'), frequency shifts ('D') and optionally constants ('const')
    @type parameters: record array
    @return: sine signal with frequency shift (same shape as C{times})
    @rtype: array
    """
    #-- take epoch into account
    if t0 is None: t0 = times[0]
    
    #-- fit the signal
    signal = np.zeros(len(times))
    if 'const' in parameters.dtype.names:
        signal += np.sum(parameters['const'])
        
    for par in parameters:
        frequency = (par['freq'] + par['D']/2.*(times-t0)) * (times-t0)
        signal += par['ampl'] * sin(2*pi*(frequency + par['phase']))
    
    return signal
    
    
def sine_orbit(times,parameters,t0=None,nmax=10):
    """
    Creates a sine function with a sinusoidal frequency shift.
    
    Parameter fields: C{const, ampl, freq, phase, asini, omega}.
    
    Similar to C{sine}, but with extra column 'asini' and 'forb', which are the
    orbital parameters. forb in cycles/day or something similar, asini in au.
    
    For eccentric orbits, add longitude of periastron 'omega' (radians) and
    'ecc' (eccentricity).
    
    Only accepts record arrays as parameters (for the moment).
    
    @param times: observation times
    @type times: numpy array
    @param parameters: record array containing amplitudes ('ampl'), frequencies ('freq'),
    phases ('phase'), frequency shifts ('D') and optionally constants ('const')
    @type parameters: record array
    @return: sine signal with frequency shift (same shape as C{times})
    @rtype: array
    """
    #-- take epoch into account
    if t0 is None: t0 = times[0]
    
    def ane(n,e): return 2.*np.sqrt(1-e**2)/e/n*jn(n,n*e)    
    def bne(n,e): return 1./n*(jn(n-1,n*e)-jn(n+1,n*e))
    
    #-- fit the signal
    signal = np.zeros(len(times))
    names = parameters.dtype.names
    if 'const' in names:
        signal += np.sum(parameters['const'])
        
    cc = 173.144632674 # speed of light in AU/d
    for par in parameters:
        alpha = par['freq']*par['asini']/cc
        if not 'ecc' in names:# or par['ecc']==0:
            frequency = par['freq']*(times-t0) + \
              alpha*(sin(2*pi*par['forb']*times) - sin(2*pi*par['forb']*t0))
        else:
            e,omega = par['ecc'],par['omega']
            ns = np.arange(1,nmax+1,1)
            ans,bns = np.array([[ane(n,e),bne(n,e)] for n in ns]).T
            ksins = sqrt(ans**2*cos(omega)**2+bns**2*sin(omega)**2)
            thns = arctan(bns/ans*tan(omega))
            tau = -np.sum(bns*sin(omega))
            frequency = par['freq']*(times-t0) + \
               alpha*(np.sum(np.array([ksins[i]*sin(2*pi*ns[i]*par['forb']*(times-t0)+thns[i]) for i in range(nmax)]),axis=0)+tau)
        signal += par['ampl'] * sin(2*pi*(frequency + par['phase']))
    
    return signal




def periodic_spline(times,parameters,t0=None):
    """
    Evaluate a periodic spline function
    
    Parameter fields: C{freq, knots, coeffs, degree}, as returned from the
    corresponding fitting function L{fit.periodic_spline}.
    
    @param times: observation times
    @type times: numpy array
    @param parameters: [frequency,"cj" spline coefficients]
    @return: sine signal with frequency shift (same shape as C{times})
    @rtype: array
    """
    #-- get time offset
    if t0 is None:
        t0 = times[0]
        
    mytrend = 0.
    #-- reconstruct entire time series, instead of phased version
    for i in range(len(parameters)):
        frequency = parameters[i]['freq']
        #-- reconstruct coefficients array
        coefficients = (parameters[i]['knots'],parameters[i]['coeffs'],parameters[i]['degree'])
        phases = np.fmod((times-t0) * frequency,1.0)
        mytrend += splev(phases,coefficients)
    
    return mytrend
    


@check_input
def kepler(times,parameters,itermax=8):
    """
    Construct a radial velocity curve due to Kepler orbit(s).
    
    Parameter fields: C{gamma, P, T0, e, omega, K}
    
    @param times: observation times
    @type times: numpy array
    @param parameters: record array containing periods ('P' in days),
    times of periastron ('T0' in days), eccentricities ('e'),
    longitudes of periastron ('omega' in radians), amplitude ('K' in km/s) and
    systemic velocity ('gamma' in km/s)
    phases ('phase'), frequency shifts ('D') and optionally constants ('const')
    @type parameters: record array
    @param itermax: maximum number of iterations for solving Kepler equation
    @type itermax: integer
    @return: radial velocity curve (same shape as C{times})
    @rtype: array
    """
    if 'gamma' in parameters.dtype.names:
        RVfit = parameters['gamma'].sum()
    else:
        RVfit = 0
    for pars in parameters:
        p = [pars['P'],pars['T0'],pars['e'],pars['omega'],pars['K'],0]
        RVfit += keplerorbit.radial_velocity(p,times=times,itermax=itermax)
    return RVfit


@check_input
def box(times,parameters,t0=None):
    """
    Evaluate a box transit model.
    
    Parameter fields: C{cont, freq, ingress, egress, depth}
    
    Parameters [[frequency,depth, fractional start, fraction end, continuum]]
    @rtype: array
    """
    if t0 is None: t0 = 0.
    
    #-- continuum
    parameters = np.asarray(parameters)
    model = np.ones(len(times))*sum(parameters['cont'])
    
    for parameter in parameters:
        #-- set up the model, define which point is where in a
        #   phasediagram
        phase = np.fmod((times - t0) * parameter['freq'], 1.0)
        phase = np.where(phase<0,phase+1,phase)
        transit_place = (parameter['ingress']<=phase) & (phase<=parameter['egress'])
        model = np.where(transit_place,model-parameter['depth'],model)
    return model



def spots(times,spots,t0=None,**kwargs):
    """
    Spotted model from Lanza (2003).
    Extract from Carrington Map.
    
    Included:
        - presence of faculae through Q-factor (set constant to 5 or see 
        Chapman et al. 1997 or Solank & Unruh for a scaling law, since this
        factor decreases strongly with decreasing rotation period.
        - differential rotation through B-value (B/period = 0.2 for the sun)
        - limb-darkening through coefficients ap, bp and cp
        - non-equidistant time series
        - individual starspot evolution and lifetime following a Gaussian law,
        (see Hall & Henry (1993) for a relation between stellar radius and
        sunspot lifetime), through an optional 3rd (time of maximum) and
        4th parameter (decay time)
    
    Not included:
        - distinction between umbra's and penumbra
        - sunspot structure
    
    Parameters of global star:
        - i_sun: inclination angle, i=pi/2 is edge-on, i=0 is pole on.
        - l0: 'epoch angle', if L0=0, spot with coordinate (0,0)
        starts in center of disk, if L0=0, spot is in center of
        disk after half a rotation period
        - period: equatorial period
        - b: differential rotation b = Omega_pole - Omega_eq where Omega is in radians per time unit (angular velocity)
        - ap,bp,cp: limb darkening
    
    Parameters of spots:
        - lambda: 'greenwich' coordinate (longitude) (pi means push to back of star if l0=0)
        - theta: latitude (latitude=0 means at equator, latitude=pi/2 means on pole (B{not
        colatitude})
        - A: size
        - maxt0 (optional): time of maximum size
        - decay (optional): exponential decay speed
        
    
    Note: alpha = (Omega_eq - Omega_p) / Omega_eq
          alpha = -b /2*pi * period
    
    Use for fitting!
    """
    #-- set zeropoint
    if t0 is None: t0 = 0.
    
    #-- simulated light curve template
    signal = np.zeros(len(times))
            
    #-- run through all timepoints
    this_signal = 0
    for spot in spots:
        #-- compute limb-darkening law
        C = 1. / (0.25 * spots[0]['ap'] + 2./3.*spots[0]['bp'] + 1./2.*spots[0]['cp'])
        #-- position and size of spot
        lambda_i = spot['lambda']
        A_i      = spot['A']
        if 'decay' in spot.dtype.names:
            #-- include evolving spot
            A_i  *= np.exp(-((spot['maxt0'] - times)/spot['decay'])**2)
        #-- find the velocity at this position (possible dependent on theta)
        Omega_theta = 2*pi/spots[0]['P_eq'] + spots[0]['b'] * sin(spot['theta'])**2
        #-- Correct the position of lambda_i according to the position
        lambda_i += Omega_theta*(times-t0)
        #-- Calculate mu_i = cos(psi_i) with psi_i the angle between the normal
        #   to the ith active region and the line of sight
        mu_i = cos(spots[0]['i_sun']) * sin(spot['theta']) + sin(spots[0]['i_sun'])*cos(spot['theta'])*cos(lambda_i - spots[0]['l0'])
        #-- Calculate total irradiance for this spot and add it to the signal
        irradiance_spot = C*spots[0]['s0'] * A_i * mu_i * \
                      (spots[0]['ap'] + spots[0]['bp']*mu_i + spots[0]['cp']*mu_i**2) *\
                      ((spots[0]['cs']-1) + spots[0]['q']*(spots[0]['cf'] + spots[0]['cf_']*mu_i - 1))
        #-- only take visible spots into account
        this_signal = where(mu_i<0,this_signal,this_signal+irradiance_spot)
    this_signal += spots[0]['s0'] + spots[0]['sr']
    return this_signal



#}


#{ General purpose

@check_input
def gauss(x, parameters):
    """
    Evaluate a gaussian fit.
    
    Parameter fields: C{'ampl','mu','sigma'} and optionally C{'const'}.
    
    Evaluating one Gaussian:
    
    >>> x = np.linspace(-20,20,10000)
    >>> pars = [0.5,0.,3.]
    >>> y = gauss(x,pars)
    >>> p = pl.plot(x,y,'k-')
    
    Evaluating multiple Gaussian, with and without constants:
    
    >>> pars = [0.5,0.,3.,0.1,-10,1.,.1]
    >>> y = gauss(x,pars)
    >>> p = pl.plot(x,y,'r-')
    
    @parameter x: domain
    @type x: ndarray
    @parameter parameters: record array.
    @type x: numpy record array
    @return: fitted signal
    @rtype: array
    """
    if 'const' in parameters.dtype.names:
        C = parameters['const'].sum()
    else:
        C = 0.
    
    y = C
    for pars in parameters:
        y += pars['ampl'] * np.exp( -(x-pars['mu'])**2 / (2.*pars['sigma']**2))
    
    return y

@check_input
def lorentz(x,parameters):
    """
    Evaluate  Lorentz profile
    
    P(f) = A / ( (x-mu)**2 + gamma**2)
    
    Parameters fields: C{ampl,mu,gamma} and optionally C{const}.
    
    Evaluating one Lorentzian:
    
    >>> x = np.linspace(-20,20,10000)
    >>> pars = [0.5,0.,3.]
    >>> y = lorentz(x,pars)
    >>> p = pl.figure()
    >>> p = pl.plot(x,y,'k-')
    
    Evaluating multiple Lorentzians, with and without constants:
    
    >>> pars = [0.5,0.,3.,0.1,-10,1.,.1]
    >>> y = lorentz(x,pars)
    >>> p = pl.plot(x,y,'r-')
    """
    if 'const' in parameters.dtype.names:
        y = parameters['const'].sum()
    else:
        y = 0.
    for par in parameters:
        y += par['ampl'] / ((x-par['mu'])**2 + par['gamma']**2)
    return y
    


@check_input
def voigt(x,parameters):
    """
    Evaluate a Voigt profile.
    
    Field parameters: C{ampl, mu, gamma, sigma} and optionally C{const}
    
    Function::
    
        z = (x + gamma*i) / (sigma*sqrt(2))
        V = A * Real[cerf(z)] / (sigma*sqrt(2*pi))
    
    >>> x = np.linspace(-20,20,10000)
    >>> pars1 = [0.5,0.,2.,2.]
    >>> pars2 = [0.5,0.,1.,2.]
    >>> pars3 = [0.5,0.,2.,1.]
    >>> y1 = voigt(x,pars1)
    >>> y2 = voigt(x,pars2)
    >>> y3 = voigt(x,pars3)
    >>> p = pl.figure()
    >>> p = pl.plot(x,y1,'b-')
    >>> p = pl.plot(x,y2,'g-')
    >>> p = pl.plot(x,y3,'r-')
    
    Multiple voigt profiles:
    
    >>> pars4 = [0.5,0.,2.,1.,0.1,-3,0.5,0.5,0.01]
    >>> y4 = voigt(x,pars4)
    >>> p = pl.plot(x,y4,'c-')
    
    @rtype: array
    @return: Voigt profile
    """
    if 'const' in parameters.dtype.names:
        y = parameters['const'].sum()
    else:
        y = 0.
    for par in parameters:
        x = x-par['mu']
        z = (x+1j*par['gamma'])/(par['sigma']*sqrt(2))
        y += par['ampl']*_complex_error_function(z).real/(par['sigma']*sqrt(2*pi))
    return y

@check_input
def power_law(x,parameters):
    """
    Evaluate a power law.
    
    Parameter fields: C{A, B, C}, optionally C{const}
    
    Function:
    
    P(f) = A / (1+ Bf)**C + const
    
    >>> x = np.linspace(0,10,1000)
    >>> pars1 = [0.5,1.,1.]
    >>> pars2 = [0.5,1.,2.]
    >>> pars3 = [0.5,2.,1.]
    >>> y1 = power_law(x,pars1)
    >>> y2 = power_law(x,pars2)
    >>> y3 = power_law(x,pars3)
    >>> p = pl.figure()
    >>> p = pl.plot(x,y1,'b-')
    >>> p = pl.plot(x,y2,'g-')
    >>> p = pl.plot(x,y3,'r-')
    
    Multiple power laws:
    >>> pars4 = pars1+pars2+pars3
    >>> y4 = power_law(x,pars4)
    >>> p = pl.plot(x,y4,'c-')
    
    @parameter x: domain
    @type x: ndarray
    @parameter parameters: record array.
    @type x: numpy record array
    @return: fitted signal
    @rtype: array
    """
    if 'const' in parameters.dtype.names:
        y = par['const']
    else:
        y = 0.
    for par in parameters:
        y += par['A'] / (1+ (par['B']*x)**par['C'])
    return y

#}



#{ Convert record arrays to flat arrays and back


def sine_preppars(pars):
    """
    Prepare sine function parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        # with constant
        if 'const' in pars.dtype.names:
            converted_pars = np.zeros(3*len(pars)+1)
            converted_pars[0] = pars['const'].sum()
            converted_pars[1::3] = pars['ampl']
            converted_pars[2::3] = pars['freq']
            converted_pars[3::3] = pars['phase']
        # without constant
        else:
            converted_pars = np.zeros(3*len(pars))
            converted_pars[0::3] = pars['ampl']
            converted_pars[1::3] = pars['freq']
            converted_pars[2::3] = pars['phase']
    #-- from flat to record array
    else:
        # with constant
        if len(pars)%3==1:
            converted_pars = np.zeros((4,(len(pars)-1)/3))
            converted_pars[0,0] = pars[0]
            converted_pars[1] = pars[1::3]
            converted_pars[2] = pars[2::3]
            converted_pars[3] = pars[3::3]
            names = ['const','ampl','freq','phase']
        # without constant
        else:
            converted_pars = np.zeros((3,len(pars)/3))
            converted_pars[0] = pars[0::3]
            converted_pars[1] = pars[1::3]
            converted_pars[2] = pars[2::3]
            names = ['ampl','freq','phase']
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars


def kepler_preppars(pars):
    """
    Prepare Kepler orbit parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        # with systemic velocity
        if 'gamma' in pars.dtype.names:
            converted_pars = np.zeros(5*len(pars)+1)
            converted_pars[0] = pars['gamma'].sum()
            converted_pars[1::5] = pars['P']
            converted_pars[2::5] = pars['T0']
            converted_pars[3::5] = pars['e']
            converted_pars[4::5] = pars['omega']
            converted_pars[5::5] = pars['K']
        else:
            converted_pars = np.zeros(5*len(pars))
            converted_pars[0::5] = pars['P']
            converted_pars[1::5] = pars['T0']
            converted_pars[2::5] = pars['e']
            converted_pars[3::5] = pars['omega']
            converted_pars[4::5] = pars['K']
    #-- from flat to record array
    else:
        # with constant
        if len(pars)%5==1:
            converted_pars = np.zeros((6,(len(pars)-1)/5))
            converted_pars[0,0] = pars[0]
            converted_pars[1] = pars[1::5]
            converted_pars[2] = pars[2::5]
            converted_pars[3] = pars[3::5]
            converted_pars[4] = pars[4::5]
            converted_pars[5] = pars[5::5]
            names = ['gamma','P','T0','e','omega','K']
        # without constant
        else:
            converted_pars = np.zeros((5,len(pars)/5))
            converted_pars[0] = pars[0::5]
            converted_pars[1] = pars[1::5]
            converted_pars[2] = pars[2::5]
            converted_pars[3] = pars[3::5]
            converted_pars[4] = pars[4::5]
            names = ['P','T0','e','omega','K']
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars
    
def box_preppars(pars):
    """
    Prepare sine function parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        converted_pars = np.ones(4*len(pars)+1)
        converted_pars[0] += sum(pars['cont'])
        converted_pars[1::4] = pars['freq']
        converted_pars[2::4] = pars['depth']
        converted_pars[3::4] = pars['ingress']
        converted_pars[4::4] = pars['egress']
    #-- from flat to record array
    else:
        converted_pars = np.ones((5,(len(pars)-1)/4))
        converted_pars[0,0] = pars[0]
        converted_pars[1] = pars[1::4]
        converted_pars[2] = pars[2::4]
        converted_pars[3] = pars[3::4]
        converted_pars[4] = pars[4::4]
        names = ['cont','freq','depth','ingress','egress']
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars    
    

def gauss_preppars(pars):
    """
    Prepare gauss function parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        if 'const' in pars.dtype.names:
            converted_pars = np.zeros(3*len(pars)+1)
            converted_pars[-1] = pars['const'].sum()
        else:
            converted_pars = np.zeros(3*len(pars))
        converted_pars[0::3] = pars['ampl']
        converted_pars[1::3] = pars['mu']
        converted_pars[2::3] = pars['sigma']
    #-- from flat to record array
    else:
        if len(pars)%3==0:
            converted_pars = np.zeros((3,len(pars)/3))
            names = ['ampl','mu','sigma']
        else:
            converted_pars = np.zeros((4,(len(pars)-1)/3))
            converted_pars[3,0] = pars[-1]
            names = ['ampl','mu','sigma','const']
        converted_pars[0] = pars[0:-1:3]
        converted_pars[1] = pars[1::3]
        converted_pars[2] = pars[2::3]
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars

def lorentz_preppars(pars):
    """
    Prepare Lorentz function parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        if 'const' in pars.dtype.names:
            converted_pars = np.zeros(3*len(pars)+1)
            converted_pars[-1] = pars['const'].sum()
        else:
            converted_pars = np.zeros(3*len(pars))
        converted_pars[0::3] = pars['ampl']
        converted_pars[1::3] = pars['mu']
        converted_pars[2::3] = pars['gamma']
    #-- from flat to record array
    else:
        if len(pars)%3==0:
            converted_pars = np.zeros((3,len(pars)/3))
            names = ['ampl','mu','gamma']
        else:
            converted_pars = np.zeros((4,(len(pars)-1)/3))
            converted_pars[3,0] = pars[-1]
            names = ['ampl','mu','gamma','const']
        converted_pars[0] = pars[0:-1:3]
        converted_pars[1] = pars[1::3]
        converted_pars[2] = pars[2::3]
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars

def voigt_preppars(pars):
    """
    Prepare voigt function parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        if 'const' in pars.dtype.names:
            converted_pars = np.zeros(4*len(pars)+1)
            converted_pars[-1] = pars['const'].sum()
        else:
            converted_pars = np.zeros(4*len(pars))
        converted_pars[0::4] = pars['ampl']
        converted_pars[1::4] = pars['mu']
        converted_pars[2::4] = pars['sigma']
        converted_pars[3::4] = pars['gamma']
    #-- from flat to record array
    else:
        if len(pars)%4==0:
            converted_pars = np.zeros((4,len(pars)/4))
            names = ['ampl','mu','sigma','gamma']
        else:
            converted_pars = np.zeros((5,(len(pars)-1)/4))
            converted_pars[4,0] = pars[-1]
            names = ['ampl','mu','sigma','gamma','const']
        converted_pars[0] = pars[0:-1:4]
        converted_pars[1] = pars[1::4]
        converted_pars[2] = pars[2::4]
        converted_pars[3] = pars[3::4]
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars


def power_law_preppars(pars):
    """
    Prepare gauss function parameters in correct form for evaluating/fitting.
    
    If you input a record array, this function will output a 1D numpy array
    containing only the independent parameters for use in nonlinear fitting
    algorithms.
    
    If you input a 1D numpy array, it will output a record array.
    
    @param pars: input parameters in record array or normal array form
    @type pars: record array or normal numpy array
    @return: input parameters in normal array form or record array
    @rtype: normal numpy array or record array
    """
    #-- from record array to flat
    if pars.dtype.names:
        if 'const' in pars.dtype.names:
            converted_pars = np.zeros(3*len(pars)+1)
            converted_pars[-1] = pars['const'].sum()
        else:
            converted_pars = np.zeros(3*len(pars))
        converted_pars[0::3] = pars['A']
        converted_pars[1::3] = pars['B']
        converted_pars[2::3] = pars['C']
    #-- from flat to record array
    else:
        if len(pars)%3==0:
            converted_pars = np.zeros((3,len(pars)/3))
            names = ['A','B','C']
        else:
            converted_pars = np.zeros((4,(len(pars)-1)/3))
            converted_pars[3,0] = pars[-1]
            names = ['A','B','C','const']
        converted_pars[0] = pars[0:-1:3]
        converted_pars[1] = pars[1::3]
        converted_pars[2] = pars[2::3]
        converted_pars = np.rec.fromarrays(converted_pars,names=names)
    return converted_pars
#}

#{ Helper functions

def _complex_error_function(x):
    """
    Complex error function
    """
    cef_value = np.exp(-x**2)*(1-erf(-1j*x))
    if sum(np.isnan(cef_value))>0:
        logger.warning("Complex Error function: NAN encountered, results are biased")
        noisnans = np.compress(1-np.isnan(cef_value),cef_value)
        try:
            last_value = noisnans[-1]
        except:
            last_value = 0
            logger.warning("Complex Error function: all values are NAN, results are wrong")
        cef_value = np.where(np.isnan(cef_value),last_value,cef_value)

    return cef_value

#}


if __name__=="__main__":
    import doctest
    import pylab as pl
    import sys
    doctest.testmod()
    pl.show()

