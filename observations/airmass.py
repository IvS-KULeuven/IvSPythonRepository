"""
Calculate airmass according to various models.
"""

import numpy as np

def _deg2rad(degrees):
   return degrees*np.pi/180.

def _rad2deg(radians):
   return radians*180./np.pi

def _sec(radians):
    return 1./np.cos(radians)

def get_modelnames():
  return(np.array(['YoungAndIrvine1967', 'Hardie1962', 'Rozenberg1966', 'KastenAndYoung1989', 'Young1994', 'Pickering2002']))

def airmass(zz, model='Pickering2002'):
  """
  Calculate the airmass from an interpolative model
  
  Input:
    @param z: array with zenith angles (90-altitude in degrees) 
    @type z: ndarray
    @param model: name of the model
    @type model: string
  
  Output:
    @return: airmass
    @rtype: ndarray
    
    Many formulas have been developed to fit tabular values of air mass:
  + The first one by Young and Irvine (1967) included a simple corrective term:
      
      X = sec(z)*(1.-0.0012*(sec(z)**2. - 1.))
    
    where zt is the true zenith angle. This gives usable results up to about
    80 degrees, but the accuracy degrades rapidly at greater zenith angles. The
    calculated air mass reaches a maximum of 11.13 at 86.6 degrees, becomes zero
    at 88 degrees, and approaches negative infinity at the horizon. The plot of this
    formula on the accompanying graph includes a correction for atmospheric
    refraction so that the calculated air mass is for apparent rather than true
    zenith angle.
    
  + Hardie (1962) introduced a polynomial in:
    
      X = sec(z) - 0.0018167*(sec(z) - 1.) - 0.002875*(sec(z) - 1.)**2. - 0.0008083*(sec(z)-1.)**3.
    
    which gives usable results for zenith angles of up to perhaps 85 degrees.
    As with the previous formula, the calculated air mass reaches a maximum,
    and then approaches negative infinity at the horizon.
    
  + Rozenberg (1966) suggested:
    
      X = 1./(cos(z) + 0.025*exp(-11.*cos(z)))
    
    which gives reasonable results for high zenith angles, with a horizon
    airmass of 40.
    
  + Kasten and Young (1989) developed:
    
      X = 1./(cos(z) + 0.50572*(96.07995 - zd)**(-1.6364))
    
    which gives reasonable results for zenith angles of up to 90 degrees, with
    an airmass of approximately 38 at the horizon. Here the zd term is z in
    degrees.
    
  + Young (1994) developed:
    
      X = (1.002432*cos(z)**2. + 0.148386*cos(z) + 0.0096467)/(cos(z)**3. + 0.149864*cos(z)**2. + 0.0102963*cos(z) + 0.000303978)
      
    in terms of the true zenith angle zt, for which he claimed a maximum error
    (at the horizon) of 0.0037 air mass.
    
  + Finally, Pickering (2002) developed
    
      X = 1./sin(h+244/(165+47*h**1.1))

  where h is apparent altitude (90.-z) in degrees. Pickering claimed his
  equation to have a tenth the error of Schaefer (1998) near the horizon.
  
  Respectively, the keyword should be:
  'YoungAndIrvine1967', 'Hardie1962', 'Rozenberg1966', 'KastenAndYoung1989',
  'Young1994', 'Pickering2002'
  
  Example usage:
    get the different airmasses for all models between 70 and 90 degrees
    >>> zz     = np.arange(70.000, 89.999, .001)
    >>> models = get_modelnames()
    >>> airmasses = np.zeros((len(zz), len(models)))
    >>> for i, model in enumerate(models): airmasses[:,i] = airmass(zz, model)
    >>> print(airmasses[-10,:])
    [ -1.69573409e+08  -1.14232330e+08   3.97784409e+01   3.77562570e+01
       3.16224873e+01   3.85395375e+01]
  """
  # define the standard quantities
  zrad = _deg2rad(zz)
  if model.lower() == 'youngandirvine1967':
    out = _sec(zrad)*(1.-0.0012*(_sec(zrad)**2. - 1.))
  elif model.lower() == 'hardie1962':
    out = _sec(zrad) - 0.0018167*(_sec(zrad) - 1.) - 0.002875*(_sec(zrad) - 1.)**2. - 0.0008083*(_sec(zrad) - 1.)**3.
  elif model.lower() == 'rozenberg1966':
    out = 1./(np.cos(zrad) + 0.025*np.exp(-11.*np.cos(zrad)))
  elif model.lower() == 'kastenandyoung1989':
    out = 1./(np.cos(zrad) + 0.50572*(96.07995 - zz)**(-1.6364))
  elif model.lower() == 'young1994':
    out = (1.002432*np.cos(zrad)**2. + 0.148386*np.cos(zrad) + 0.0096467)/(np.cos(zrad)**3. + 0.149864*np.cos(zrad)**2. + 0.0102963*np.cos(zrad) + 0.000303978)
  elif model.lower() == 'pickering2002':
    hh   = 90.-zz
    out = 1./np.sin(_deg2rad(hh+244/(165+47*hh**1.1)))
  return (out)