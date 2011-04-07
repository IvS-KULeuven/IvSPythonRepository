from numpy import *
from numpy.fft import *
from numpy.linalg import *
from numpy.random import normal
  



def DFTpower(time,signal,startfreq,stopfreq,stepfreq):

    """
    # PURPOSE: Discrete Fourier Transform. Computes the modulus square 
    #          of the fourier transform. Unit: square of unit of signal.
    #
    # INPUT: . time[0..Ntime-1]: time points
    #        . signal[0..Ntime-1]: signal
    #        . startfreq: freq = arange(startfreq,stopfreq,stepfreq)
    #        . stopfreq:  freq = arange(startfreq,stopfreq,stepfreq)
    #        . stepfreq:  freq = arange(startfreq,stopfreq,stepfreq)
    # 
    # OUTPUT: . power[0..Nfreq-1]
    #
    # REMARKS: . time, signal are assumed to be numpy arrays.
    #          . time points need *not* to be equidistant
    #          . normalisation is such that a signal A*sin(2*pi*nu_0*t)
    #            gives power A^2 at nu=nu_0
    """
  
    Ntime = len(time)
    Nfreq = int(ceil((stopfreq-startfreq)/stepfreq))
  
    A = exp(1j*2.*pi*startfreq*time) * signal
    B = exp(1j*2.*pi*stepfreq*time)
    ft = zeros(Nfreq, complex) 
    ft[0] = A.sum()

    for k in range(1,Nfreq):
        A *= B
        ft[k] = sum(A)
    
    return (ft.real**2 + ft.imag**2) * 4.0 / Ntime**2
  
  
  
  
def SFTpower(time, signal, freq):

    """
    # PURPOSE: Slow Fourier Transform. Computes the modulus square 
    #          of the fourier transform. Unit: square of unit of signal
    #
    # INPUT: . time[0..Ntime-1]: time points
    #        . signal[0..Ntime-1]: signal
    #        . freq[0..Nfreq-1]: frequency points: unit: inverse of time unit
    # 
    # OUTPUT: . power[0..Nfreq-1]
    #
    # REMARKS: . time, signal and freq are assumed to be numarrays.
    #          . time points need *not* to be equidistant
    #          . normalisation is such that a signal A*sin(2*pi*nu_0*t)
    #            gives power A^2 at nu=nu_0
    """
  
    Nfreq = len(freq)
    Ntime = len(time)
    result = zeros(Nfreq)

    for i in range(Nfreq):
        arg = 2.0*pi*freq[i]*time
        result[i] = sum(signal * cos(arg))**2 + sum(signal * sin(arg))**2

    result = result * 4.0 / Ntime**2

    return(result)



  
  
  
  
def FFTpower(signal, timestep):

    """
    # PURPOSE: computes the power of an equidistant time series 'signal',
    #          using the FFT algorithm
    #
    # INPUT: . signal[0..Ntime-1]: signal
    #        . timestep: time step of the equidistante time series
    #
    # OUTPUT: . freq[0..Nfreq-1]: frequencies from 0 up to the Nyquist freq. 
    #         . power[0..Nfreq-1]: the power...
    #
    # REMARKS: . Ntime need not be a power of 2
    #            (internally, zero padding is automatic)
    #          . normalisation is such that a signal A*sin(2*pi*nu_0*t)
    #            gives power A^2 at nu=nu_0  (IF nu_0 is in the 'freq' array)
    """
  
    # Compute the FFT of a real-valued signal. If N is the number 
    # of points of the original signal, 'Nfreq' is (N/2+1).
  
    fourier = refft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)
  
    # Compute the power
  
    power = abs(fourier)**2 * 4.0 / Ntime**2
  
    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep
  
    freq = arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep
  
    # That's it!
  
    return (freq, power)





  
  
  
def FFTpowerdensity(signal,timestep):
  
    """
    # PURPOSE: computes the power density of an equidistant time series 'signal',
    #          using the FFT algorithm
    #
    # INPUT: . signal[0..Ntime-1]: signal
    #        . timestep: time step of the equidistante time series
    #
    # OUTPUT: . freq[0..Nfreq-1]: frequencies from 0 up to the Nyquist freq.
    #         . powerdensity[0..Nfreq-1]: the power density...
    #
    # REMARKS: . Ntime need not be a power of 2
    #            (internally, zero padding is automatic)
    """
  
    # Compute the FFT of a real-valued signal. If N is the number 
    # of points of the original signal, 'Nfreq' is (N/2+1).
  
    fourier = rfft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)
  
    # Compute the power density
  
    powerdensity = abs(fourier)**2 / Ntime * timestep
  
    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep
  
    freq = arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep
  
    # That's it!
  
    return (freq, powerdensity)

  


  
  



def weightedpower(time,signal,weight,freq):

    """
    # PURPOSE: compute the weighted power spectrum
    #
    # INPUT: . time[0..Ntime-1]:   time points
    #        . signal[0..Ntime-1]: observations
    #        . weight[0..Ntime-1]:  1/sigma_i^2 of observation
    #        . freq[0..Nfreq]: frequencies to compute the power
    #
    # OUTPUT: . weighted power [0..Nfreq]
    #
    """

    result = zeros(len(freq))

    for i in range(len(freq)):
        if (freq[i] != 0.0):
            sine   = sin(2.0*pi*freq[i]*time)
            cosine = cos(2.0*pi*freq[i]*time)
            a11= sum(weight*sine*sine)
            a12 = sum(weight*sine*cosine)
            a21 = a12
            a22 = sum(weight*cosine*cosine)
            b1 = sum(weight*signal*sine)
            b2 = sum(weight*signal*cosine)
            denominator = a11*a22-a12*a21
            A = (b1*a22-b2*a12)/denominator
            B = (b2*a11-b1*a21)/denominator
            result[i] = A*A+B*B
        else:
            result[i] = sum(signal)/len(signal)

    return(result)









def Zwavelet(time,signal,freq,position,sigma=10.0):

    """
    # PURPOSE: It computes "Weighted Wavelet Z-Transform" as defined
    #          by G. Foster, 1996, Astronomical Journal, 112, 1709 
    #          The result is to be plotted as a colorimage, with 
    #          the x-axis the frequency, and y-axis the position.
    #
    # INPUT: . time[0..Ntime-1]: time points
    #        . signal[0..Ntime-1]: observed data points
    #        . freq[0..Nfreq-1]: frequencies (omega/2pi in the paper)
    #        . position[0..Npos-1]: tau in the paper
    #
    # OUTPUT: . Z[0..Npos-1, 0..Nfreq-1]: the Z-transform.
    #
    # REMARK: . Z can be visualized with pylab with:
    #           >>> imshow(Z,extent=[min(freq),max(freq),max(position),min(position)])
    #           Mind the max, min order for position. Alternatively:
    #           >>> imshow(Z[::-1],extent=[min(freq),max(freq),min(position),max(position)])
    #           Also try log(Z), and different sigmas.
    """
  
    y = zeros(3)
    S = zeros([3,3])
    Z = zeros([len(position),len(freq)])
  
    for i in range(len(position)):
        
        print i,"/",len(position)   
        tau = position[i]
        arg = 2.0*pi*(time-tau)
 
        for j in range(len(freq)):
        
            nu = freq[j]
      
            # Compute statistical weights akin the Morlet wavelet
      
            weight = exp(-(time-tau)**2 * (nu / 2./sigma)**2)
            W = sum(weight)
      
            # Compute the base functions. A 3rd base function is the constant 1.
      
            cosine = cos(arg*nu)
            sine = sin(arg*nu)
          
            # Compute the innerproduct of the base functions
            # phi_0 = 1 (constant), phi_1 = cosine, phi_2 = sine
      
            S[0,0] = 1.0
            S[0,1] = S[1,0] = sum(weight * 1.0 * cosine) / W
            S[0,2] = S[2,0] = sum(weight * 1.0 * sine) / W
            S[1,1] = sum(weight * cosine * cosine) / W
            S[1,2] = S[2,1] = sum(weight * cosine * sine) / W
            S[2,2] = sum(weight * sine * sine) / W
            invS = inverse(S)
      
            # Determine the best-fit coefficients y_k of the base functions
      
            for k in range(3):
                y[k] =   sum(weight * 1.0 * signal) / W * invS[k,0]     \
                       + sum(weight * cosine * signal) / W * invS[k,1]  \
                       + sum(weight * sine * signal) / W * invS[k,2]      
      
            # Compute the best-fit model
      
            model = y[0] + y[1] * cosine + y[2] * sine
      
            # Compute the weighted variation of the signal and the model functions
      
            Vsignal = sum(weight * signal**2) / W -  (sum(weight * signal) / W)**2
            Vmodel = sum(weight * model**2) / W -  (sum(weight * model) / W)**2
      
            # Calculate the weighted Wavelet Z-Transform
      
            Neff = W**2 / sum(weight**2)
            Z[i,j] = (Neff - 3) * Vmodel / 2. / (Vsignal - Vmodel)
      
    # That's it!
  
    return Z
  
  
  
  
           
           
           
           
           
           
def windowfunction(time, freq):

    """
    # PURPOSE: computes the power spectrum of the window function
    #          of a time series 'time' at the frequencies 'freq'
    #
    # INPUT: . time[0..Ntime-1]: time points, not necessarily equidistant
    #        . freq[0..Nfreq-1]: frequency points. Units: inverse unit of 'time'
    #
    # OUTPUT: . winkernel[0..Nfreq-1]: window function |W(freq)|^2
    """
  
    Ntime = len(time)
    Nfreq = len(freq)
    winkernel = empty(Nfreq)

    for i in range(Nfreq):
        winkernel[i] = sum(cos(2.0*pi*freq[i]*time))**2     \
                     + sum(sin(2.0*pi*freq[i]*time))**2

    # Normalise such that winkernel(nu = 0.0) = 1.0 

    return winkernel/Ntime**2







def pdm(time, signal, freq, Nbin=10, Ncover=5):

    """
    # PURPOSE: compute the theta-statistics to do a Phase Dispersion Minimisation.
    #          (cf. Stellingwerf R.F., 1978, ApJ, 224, 953)
    #
    # INPUT: . time[0..Ntime-1]: array of time points
    #        . signal[0..Ntime-1]: array of observations
    #        . freq[0..Nfreq-1]: array of frequencies to compute theta for
    #        . Nbin: the number of phase bins (with length 1 / Nbin)
    #        . Ncover: the number of covers (i.e. bin shifts)
    #
    # OUTPUT: . theta[0..Nfreq-1]: array containing the theta-statistic for each
    #                              frequency
    """
  
    Ntime = len(time)
    Nfreq = len(freq)
  
    binsize = 1.0 / Nbin
    covershift = 1.0 / (Nbin * Ncover)
  
    theta = zeros(Nfreq)
  
    for i in range(Nfreq):
  
        # Compute the phases in [0,1[ for all time points
    
        phase = fmod((time - time[0]) * freq[i], 1.0)
    
        # Reset the number of (shifted) bins without datapoints
    
        Nempty = 0
    
        # Loop over all Nbin * Ncover (shifted) bins
    
        for k in range(Nbin):
            for n in range(Ncover):
        
            # Determine the left and right boundary of one such bin
            # Note that due to the modulo, right may be < left. So instead
            # of 0-----+++++------1, the bin might be 0++-----------+++1 .
        
            left = fmod(k * binsize + n * covershift, 1.0) 
            right = fmod((k+1) * binsize + n * covershift, 1.0) 

            # Select all data points in that bin
        
            if (left < right):
                bindata = compress((left <= phase) & (phase < right), signal)
            else:
                bindata = compress(~((right <= phase) & (phase < left)), signal)

            # Compute the contribution of that bin to the theta-statistics  
          
            if (len(bindata) != 0):
                theta[i] += (len(bindata) - 1) * bindata.var()
            else:
                Nempty += 1
  
        # Normalize the theta-statistics        

        theta[i] /= Ncover * Ntime - (Ncover * Nbin - Nempty)     
    
    # Normalize the theta-statistics again
  
    theta /= signal.var()  
    
    # That's it!
 
    return theta
    
  
  



  

def phasediagram(time, signal, nu0, t0 = None):

    """
    # PURPOSE: compute the phases of a phase diagram 
    # 
    # INPUT: . time[0..Ntime-1]: time points
    #        . signal[0..Ntime-1]: observations
    #        . nu0: frequency (inverse unit of time)
    #        . t0: zero time point, if None: t0 = time[0]
    #
    # OUTPUT: . phase[0..Ntime-1]: phases, sorted from small to large
    #         . xxxxx[0..Ntime-1]: corresponding observations
    """

    if (t0 is None): t0 = time[0]
    phase = fmod((time - t0) * nu0, 1.0)
    indices = argsort(phase)
    
    return phase[indices], signal[indices]
  
            
  
  
  
  
  
  
  
  
def autocor(freq, spectrum, maxlag, threshold = 0., interval = ()):

    """
    # PURPOSE: computes the autocorrelation of a signal
    #
    # INPUT: . freq[0..Nfreq-1]: frequencies, equidistant
    #        . spectrum[0..Nfreq-1]: power or amplitude spectrum
    #        . maxlag: the maximum lag value for the crosscorrelation
    #        . threshold: peaks in the spectrum with a top below this threshold
    #                     are not taken into account.
    #        . interval: = (lower,upper), lower and upper boundaries of the
    #                      spectrum interval to be used. If nothing given, the
    #                      entire spectrum will be used.
    #
    # OUTPUT: . lag: from 0 to maxlag with the same step as the spectrum
    #         . cor: the autocorrelation
    #
    # REMARK: . Peaks with a top above the threshold are taken fully into account,
    #           i.e. the peak profile all the way down to 0.
    """

    Nfreq = len(freq)
    if interval is not ():
        lower,upper = interval
    else:
        lower,upper = freq[0],freq[-1]

    low = empty(Nfreq, float)
    low.fill(lower)
    up = empty(Nfreq, float)
    up.fill(upper)

    #-- Extract the interesting piece of the spectrum

    subfreq = compress((freq >= low) & (freq <= up), freq)
    subspec = compress((freq >= low) & (freq <= up), spectrum)

    #-- Recompute the number of freqs, 'cause the array changed in size

    Nfreq = len(subfreq)
    print "Selected", Nfreq, "frequency points in specified range"

    #-- Extract all peaks larger than the treshold

    subspec = thresholdpeaks(subspec, threshold)

    #-- Compute the correlation

    Nlag = int(maxlag / (freq[1] - freq[0]))
    lag = arange(Nlag)*(freq[1] - freq[0])
    correl = zeros(Nlag, float)

    for i in range(Nlag):
        for j in range(Npeak):
            for k in range(left[j],right[j]+1):
                correl[i] = correl[i] + subspec[k] * subspec[(k+i)%Nfreq]

    # That's it!

    return lag, correl









def peaksAboveThreshold(spec, threshold):

    """
    # PURPOSE: Select peaks higher than a given threshold. Don't do this just
    #          by zeroing everything below this threshold because this would
    #          also cut part of the wings of the peaks. Rather follow the wings
    #          of those peaks with a top higher than threshold, until they reach
    #          a minimum. Keep those peaks intact and zero everything else.
    #
    # INPUT: . spec[0..Nspec-1]: amplitude or power spectrum
    #        . threshold: minimum height of the peak if it's not to be zeroed.
    #
    # OUTPUT: . subspec[0..Nspec-1]: 
    #
    """
  
    Nspec = len(spec)
  
    #-- Extract all peaks larger than the treshold

    peakindex = [ ]
    for i in range(1,Nspec-1):
        if ((spec[i]>spec[i-1]) & (spec[i]>spec[i+1]) & (spec[i] > threshold)) :
            peakindex.append(i)

    Npeak = len(peakindex)
    print "Found", Npeak, "peaks above treshold"

    #-- Find the left and right hand border of each peak

    left = zeros(Npeak)
    right = zeros(Npeak)

    for i in range(Npeak):

        left[i] = peakindex[i]
        while ((spec[left[i]-1] < spec[left[i]]) & (left[i] > 0)):
            left[i] -= 1

        right[i] = peakindex[i]
        while ((spec[right[i]+1] < spec[right[i]]) & (right[i] < Nspec-2)):
            right[i] += 1


    #-- Construct a mask that contains ones in peaks and zeros everywhere else

    mask = zeros(Nspec)
    for i in range(Npeak):
        mask[left[i]:right[i]] = 1


    #-- Zero everything outside the peaks

    return where(mask == 1, spec, 0.0)  
  
  







def peakTops(freq, power):

    # Return only the peaks up of a power spectrum

    Nfreq = len(freq)

    subfreq = []
    subpower = []

    for i in range(1, Nfreq-2):
        if ((power[i] > power[i-1]) & (power[i] > power[i+1])):
            subfreq.append(freq[i])
            subpower.append(power[i])

    return array(subfreq),array(subpower)













  
  
  
  
  
  
  
  
  
  
  
#------------------------------------------------------------------------------------
def harmonicfit(time, signal, sigma, freq, constant = None, error = None, t0 = None):
#------------------------------------------------------------------------------------

    """
    # PURPOSE: Given a time series and frequencies, fit a function of the form
    #            C + \sum_j A_j sin(2pi nu_j (t-t0) + phi_j)
    #          where the presence of the constant term is an option. The error
    #          bars on the fitted parameters can also be requested.
     #
    # INPUT: . time[0..Ntime-1]: the time points
    #        . signal[0..Ntime-1]: the observations
    #        . sigma[0..Ntime-1]: the standard error of the observations
    #        . freq[0..Nfreq-1]: the frequencies of the harmonics
    #        . constant: flag, if not None, also fit a constant
    #        . error: flag, if not None, also compute the errorbars
    #        . t0: time zero point. If none given, t0 = time[0]
    #
    # OUTPUT:  Depending on the flags:
    #          either: amplitude, phase
    #          or    : amplitude, phase, constant
    #          or    : amplitude, phase, err_ampl, err_phase
    #          or    : amplitude, phase, constant, err_ampl, err_phase, err_const
    #          with the phase in radians
    #
    # REMARK: . All input arrays should be numpy arrays
    """
  
    # Determine the number of fitparameters

    Ndata = len(time)
    Nfreq = len(freq)
    if constant is None:
        Nparam = 2*Nfreq
    else:
        Nparam = 2*Nfreq+1

    # Subtract the zero point from the time points. If none is given, use the
    # first time point.
  
    if t0 is None:
        t = time - time[0]
    else:
        t = time - t0
    
    # The fit function used is of the form
    #    C + \sum_j a_j sin(2pi\nu_j t_i) + b_j cos(2pi\nu_j t_i)
    # which is linear in its fit parameters. These parameters p can therefore be
    # solved by minimizing ||b - A p||, where b are the observations, and A is the
    # basisfunction matrix, i.e. A[i,j] is the j-th base function evaluated in 
    # the i-th timepoint. The first Nfreq columns are the amplitudes of the sine,
    # the second Nfreq columns are the amplitudes of the cosine, and if requested,
    # the last column belongs to the constant

    A = empty((Ndata,Nparam), float)
    for j in range(Nfreq):
        A[:,j]       = sin(2.0*pi*freq[j]*t) / sigma
        A[:,Nfreq+j] = cos(2.0*pi*freq[j]*t) / sigma
    if constant is not None:
        A[:,2*Nfreq] = ones(Ndata, float) / sigma
    b = signal / sigma
  
    # Solve using SVD decomposition

    fitparam, chisq, rank, s = lstsq(A,b)
  
     # Compute the amplitudes and phases: A_j sin(2pi*\nu_j t_i + phi_j)
  
    amplitude = zeros(Nfreq, float)
    phase = zeros(Nfreq, float)
    for j in range(Nfreq):
        amplitude[j] = sqrt(fitparam[j]**2 + fitparam[Nfreq+j]**2)
        phase[j] = arctan2(fitparam[Nfreq+j], fitparam[j])
    if constant is not None: const = fitparam[2*Nfreq]  

    # If no error bars are needed, we are finished here
  
    if error is None:
        if constant is None:
            return amplitude, phase
        else:
            return amplitude, phase, const
  
    # If error bars are needed, we compute them here.
    # First the derivative matrix: the 1st Nfreq columns the derivative w.r.t.
    # the amplitude, the 2nd Nfreq columns w.r.t. the phase, and the if relevant
    # the last column w.r.t. the constant. From this the covariance matrix.
  
    F = zeros((Ndata, Nparam), float)   
    for i in range(Ndata):
        F[i,0:Nfreq] = sin(2.0*pi*freq*t[i] + phase)
        F[i,Nfreq:2*Nfreq] = amplitude * cos(2.0*pi*freq*t[i] + phase)
        if constant is not None:
            F[i,2*Nfreq] = 1.0
      
    covariance = inverse(dot(F.T, F))
    covariance *= chisq / (Ndata - Nparam)
  
    error_ampl = sqrt(covariance.diagonal()[0:Nfreq])
    error_phase = sqrt(covariance.diagonal()[Nfreq:2*Nfreq])
    if constant is not None:
        error_const = sqrt(covariance[2*Nfreq,2*Nfreq])
    
    # And we return the estimates + their error bars
    # And that's it!
  
    if constant is None:
        return amplitude, phase, error_ampl, error_phase
    else:
        return amplitude, phase, const, error_ampl, error_phase, error_const  
  






  
    
    
    
    
    
    
    
    
  
  
  






  
  
  
  
  
  
  
  
  





#--------------------------------------------------------------------
def combresponse(freq, spec, separation, nu0, weight, threshold=0.0):
#--------------------------------------------------------------------

  """
  # PURPOSE: compute comb response of a spectrum S around a central peak nu_0:
  #           C(nu_0,Delta) = sum_i=0^N (S(nu_0+i*Delta) * S(nu_0-i*Delta))^w_i
  #          where Delta is the separation between two adjacent comb peaks,
  #          and w_i are the weights of the peaks.
  #
  # INPUT: . freq[0..Nfreq-1]: equidistant frequency array of power spectrum
  #        . spec[0..Nfreq-1]: power spectrum S()
  #        . separation[0..Nsep-1]: freq separations Delta, one per comb
  #        . nu0: freq of peak that *must* exactly be in the comb
  #        . weight[0..Nweight-1]: e.g. weight[1] determines the weight for the
  #                                +/- 1 peaks. Put ones, if no weight wanted.
  #        . threshold: only part of spectrum >= threshold is taken into account
  #
  # OUTPUT: response[0..Nsep-1]: response of each of the combs
  """
  
  # Figure out the frequency resolution
  
  step = freq[1] - freq[0]

  # Compute the comb responses. sep  is the frequency separation
  
  response = zeros(len(separation), float)
  
  for i in range(len(separation)):
          
    # Construct two arrays with the freq positions of the peaks in 
    # the comb: {nu0, nu0+sep, nu0+2sep,...} and {nu0, nu0-sep, nu0-2sep,...}           
    # Then determine the indices of these frequencies in the freq array
    
    poscomb = array([nu0 + k * separation[i] for k in range(len(weight))])
    negcomb = array([nu0 - k * separation[i] for k in range(len(weight))])
    
    posindex = around((poscomb-freq[0])/step).astype(Int)
    negindex = around((negcomb-freq[0])/step).astype(Int)
    
    # Threshold the power spectrum
    
    spec = where(spec >= threshold, spec, 0.0)

    # Compute the weighted product of the peaks
    
    for j in range(len(posindex)):
      response[i] += (spec[posindex[j]] * spec[negindex[j]])**weight[j]
        
    response[i] /= len(poscomb)
    
  # That's it!

  return response
  
