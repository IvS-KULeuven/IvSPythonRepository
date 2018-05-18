"""
Compute the moments of a line profile.
"""
import numpy as np
import scipy.integrate
import pylab as pl
from ivs import config

def moments(velo,flux,SNR=200.,max_mom=3):
    """
    Compute the moments from a line profile.
    """
    mymoms,e_mymoms = np.zeros(max_mom+1),np.zeros(max_mom+1)

    #-- compute zeroth moment (equivalent width)
    m0 = scipy.integrate.simps( (1-flux),x=velo)

    #-- to calculate the uncertainties, we need the error in each velocity
    #   bin, and the error on the zeroth moment
    sigma_i = 1./SNR / np.sqrt(flux)
    Delta_v0 = np.sqrt(scipy.integrate.simps( sigma_i,x=velo)**2)
    e_m0 = np.sqrt(2)*(Delta_v0/m0)

    mymoms[0] = m0
    e_mymoms[0] = e_m0

    #-- calculate other moments
    for n in range(1,max_mom+1):
        mymoms[n] = scipy.integrate.simps((1-flux)*velo**n,x=velo)/ m0
        #-- and the uncertainty on it
        Delta_vn = np.sqrt(scipy.integrate.simps(sigma_i*velo**n)**2)
        e_mymoms[n] = np.sqrt(  (Delta_vn/m0)**2 + (Delta_v0/m0 * mymoms[n])**2 )

    return mymoms,e_mymoms


def moments_fromfiles(filelist,read_func,max_mom=3,**kwargs):
    """
    Compute the moments from a list of files containing line profiles.

    The C{read_func} should return 2 arrays and a float, representing velocities,
    normalised fluxes and the SNR of the line profile.

    m0: equivalent width
    m1: radial velocity
    m2: variance
    m3: skewness

    @param filelist: list of filenames
    @type filelist: list of strings
    @param read_func: function which reads in a file and returns velocities (array),
    the line profile (array) and the SNR of the whole profile (float)
    @type read_func: Python function
    @param max_mom: maximum moment to compute
    @type max_mom: integer
    @return: a list containing the moments and a list containing the errors
    @rtype: [max_mom x array],[max_mom x array]

    Example usage:

    >>> filelist = np.loadtxt(  config.get_datafile('Spectra_testfiles','profiles.list')   ,str)

    >>> output,errors = moments_fromfiles(filelist,get_profile_from_file)
    """

    output = [np.zeros(len(filelist)) for i in range(max_mom+1)]
    errors = [np.zeros(len(filelist)) for i in range(max_mom+1)]
    for i,f in enumerate(filelist):
        #-- read in the profile
        velo,flux,SNR = read_func(f,**kwargs)
        mymoms,e_mymoms = moments(velo,flux,SNR)
        for n in range(len(mymoms)):
            output[n][i] = mymoms[n]
            errors[n][i] = e_mymoms[n]
        pl.plot(velo,flux+i*0.01)
    return output,errors

def profiles_fromfiles(filelist,read_func,max_mom=3,**kwargs):
    """
    Compute an average profile from a file list of profiles.

    The C{read_func} should return 2 arrays and a float, representing velocities,
    normalised fluxes and the SNR of the line profile.

    m0: equivalent width
    m1: radial velocity
    m2: variance
    m3: skewness

    @param filelist: list of filenames
    @type filelist: list of strings
    @param read_func: function which reads in a file and returns velocities (array),
    the line profile (array) and the SNR of the whole profile (float)
    @type read_func: Python function
    @param max_mom: maximum moment to compute
    @type max_mom: integer
    @return: a list containing the moments and a list containing the errors
    @rtype: velo,av_prof,fluxes

    Example usage:

    >>> filelist = np.loadtxt(  config.get_datafile('Spectra_testfiles','profiles.list')   ,str)

    >>> templatevelo,av_prof,fluxes = profiles_fromfiles(filelist,get_profile_from_file)
    """

    output = [np.zeros(len(filelist)) for i in range(max_mom+1)]
    errors = [np.zeros(len(filelist)) for i in range(max_mom+1)]
    for i,f in enumerate(filelist):
        #-- read in the profile
        velo,flux,SNR = read_func(f,**kwargs)
        if i==0:
            template_velo = velo
            fluxes = np.zeros((len(filelist),len(flux)))
            fluxes[0] = flux
        else:
            fluxes[i] = np.interp(template_velo,velo,flux)

    av_prof = fluxes.mean(axis=0)
    fluxes -= av_prof

    return template_velo,av_prof,fluxes



if __name__=="__main__":
    import doctest
    """
    The doctest examples are un-normalised profiles in wavelength space.
    Correct science usage should be normalised profiles in velocity space.
    """
    def get_profile_from_file(filename):
        """
        Used for doctest purposes: import line profile from a given test file
        """
        filedata = np.loadtxt(config.get_datafile('Spectra_testfiles',filename),delimiter = '\t').T
        #Calculate SNR of profiles
        SNR = np.nanmean(filedata[1])/np.nanstd(filedata[1])
        #returns the profile velocity, flux, and SNR
        return filedata[0],filedata[1],SNR

    doctest.testmod()
    pl.show()
