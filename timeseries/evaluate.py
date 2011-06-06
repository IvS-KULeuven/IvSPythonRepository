"""
Evaluate model functions of timeseries.
"""
import logging

import numpy as np
from numpy import pi,cos,sin
from scipy.interpolate import splev

from ivs.binary import keplerorbit

logger = logging.getLogger('TS.EVAL')
            
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
    if (t0 is None): t0 = time[0]
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




def sine(times, parameters):
    """
    Creates a harmonic function based on given parameters.
    
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
    
    >>> parameters = np.rec.fromarrays([const,ampls,freqs,phases],names = ('const','ampl','freq','phase'))
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
    phases ('phase') and optionally constants ('const'), or 1D numpy array (see above)
    @type parameters: record array or 1D array
    @return: sine signal (same shape as C{times})
    @rtype: array
    """
    if not parameters.dtype.names:
        parameters = sine_preppars(parameters)
    
    if 'const' in parameters.dtype.names:
        total_fit = parameters['const'].sum()
    for i in xrange(len(parameters)):
        total_fit += parameters['ampl'][i]*sin(2*pi*(parameters['freq'][i]*times+parameters['phase'][i]))
    return total_fit






def sine_freqshift(times,parameters,t0=None):
    """
    Creates a sine function with a linear frequency shift.
    
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
    if t0 is None:
        t0 = times[0]
    
    #-- fit the signal
    signal = np.zeros(len(times))
    if 'const' in parameters.dtype.names:
        signal += np.sum(parameters['const'])
        
    for par in parameters:
        frequency = (par['freq'] + par['D']/2.*(times-t0)) * (times-t0)
        signal += par['ampl'] * sin(2*pi*(frequency + par['phase']))
    
    return signal
    
    





def periodic_spline(times,parameters,t0=None):
    """
    Evaluate a periodic spline function
    
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
    



def kepler(times,parameters,itermax=8):
    """
    Construct a radial velocity curve due to Kepler orbit(s).
    
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
    if not parameters.dtype.names:
        parameters = kepler_preppars(parameters)
    
    if 'gamma' in parameters.dtype.names:
        RVfit = parameters['gamma'].sum()
    else:
        RVfit = 0
    for pars in parameters:
        p = [pars['P'],pars['T0'],pars['e'],pars['omega'],pars['K'],0]
        RVfit += keplerorbit.radial_velocity(p,times=times,itermax=itermax)
    return RVfit




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
#}