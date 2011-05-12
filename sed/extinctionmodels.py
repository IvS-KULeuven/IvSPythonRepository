# -*- coding: utf-8 -*-
# Extinction models of Arenou, Drimmel and Marschall
# +++ Uniformly rewritten by K. Smolders

from ivs.catalogs import vizier
from ivs.sed.reddening import get_law
from ivs.misc.decorators import memoized
from ivs.misc import loggers
from ivs import config

import numpy  as np
from numpy import (abs, arange, array, ceil, cos, dot, floor, int, logical_and,
                   logical_or, max, min, ones, pi, sin, sqrt, where, zeros)
import scipy  as sc
import pyfits as pf
import logging

logger = logging.getLogger("SED.EXT")
logger.addHandler(loggers.NullHandler)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#{ Top level wrapper
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def findext(lng, lat, model='drimmel', distance=None, **kwargs):
  """
  Get the "model" extinction at a certain longitude and latitude.
  
  Find the predicted V-band extinction (Av) based on three
  dimensional models of the galactic interstellar extinction.
  The user can choose between different models by setting the
  model keyword:
  
  1) "arenou": model from Arenou et al. (1992).
  2) "schlegel": model from Schlegel et al. (1998)
  3) "drimmel": model from Drimmel et al. (2003)
  4) "marshall": model from Marshall et al. (2006)
    
  example useage:
  
    1. Find the total galactic extinction for a star at galactic lattitude 10.2
       and longitude 59.0 along the line of sight, as given by the model of
       Arenou et al. (1992)
       
        >>> lng = 10.2
        >>> lat = 59.0
        >>> av = findext(lng, lat, model='arenou')
        >>> print("Av at lng = %.2f, lat = %.2f is %.2f magnitude" %(lng, lat, av))
        Av at lng = 10.20, lat = 59.00 is 0.05 magnitude
    
    2. Find the extinction for a star at galactic lattitude 107.05 and
       longitude -34.93 and a distance of 144.65 parsecs, as given by the
       model of Arenou et al. (1992)
       
        >>> lng = 107.05
        >>> lat = -34.93
        >>> dd  = 144.65
        >>> av = findext(lng, lat, distance = dd, model='arenou')
        >>> print("Av at lng = %.2f, lat = %.2f and distance = %.2f parsecs is %.2f magnitude" %(lng, lat, dd, av))
        Av at lng = 107.05, lat = -34.93 and distance = 144.65 parsecs is 0.15 magnitude
        
    3. Find the Marschall extinction for a star at galactic longitude 10.2 and
       latitude 9.0. If the distance is not given, we return the
       complete extinction along the line of sight (i.e. put the star somewhere
       out of the galaxy).
       
        >>> lng = 10.2
        >>> lat = 9.0
        >>> av = findext(lng, lat, model='marshall')
        >>> print("Av at lng = %.2f, lat = %.2f is %.2f magnitude" %(lng, lat, av))
        Av at lng = 10.20, lat = 9.00 is 10.67 magnitude
        
    4. Find the Marshall extinction for a star at galactic lattitude 271.05 and
       longitude -4.93 and a distance of 144.65 parsecs, but convert to Av using
       Rv = 2.5 instead of Rv = 3.1
       
        >>> lng = 271.05
        >>> lat = -4.93
        >>> dd  = 144.65
        >>> av = findext(lng, lat, distance = dd, model='marshall', Rv=2.5)
        >>> print("Av at lng = %.2f, lat = %.2f and distance = %.2f parsecs is %.2f magnitude" %(lng, lat, dd, av))
        Av at lng = 271.05, lat = -4.93 and distance = 144.65 parsecs is 13.95 magnitude
        
    5. Find the extinction for a star at galactic longitude 10.2 and
       latitude 9.0, using the schlegel model, using Rv=2.9 instead of
       Rv=3.1
       
        >>> lng = 58.2
        >>> lat = 24.0
        >>> distance = 848.
        >>> av = findext(lng, lat, distance=distance)
        >>> print("Av at lng = %.2f, lat = %.2f is %.2f magnitude" %(lng, lat, av))
        Av at lng = 58.20, lat = 24.00 is 0.12 magnitude
        
  @param lng: Galactic Longitude (in degrees)
  @type lng: float
  @param lat: Galactic Lattitude (in degrees)
  @type lat: float
  @keyword model: the name of the extinction model: ("arenou", "schlegel", "drimmel" or "marshall"; if none given, the program uses "drimmel")
  @type model: str
  @keyword distance: Distance to the source (in parsecs), if the distance is not given, the total galactic extinction along the line of sight is calculated
  @type distance: float
  @return: The extinction in Johnson V-band
  @rtype: float
  
  REMARKS:
    a) Schlegel actually returns E(B-V), this value is then converted to Av (the desired value for Rv can be set as a keyword; standard sets Rv=3.1)
    b) Schlegel is very dubious for latitudes between -5 and 5 degrees
    c) Marschall actually returns Ak, this value is then converted to Av (the reddening law and Rv can be set as keyword; standard sets Rv=3.1, redlaw='cardelli1989')
    d) Marschall is only available for certain longitudes and latitudes:
       0 < lng < 100 or 260 < lng < 360 and -10 < lat < 10
  """
  
  if model.lower() == 'drimmel':
    av = findext_drimmel(lng, lat, distance=distance, **kwargs)
  elif model.lower() == 'marshall' or model.lower() == 'marschall':
    av = findext_marshall(lng, lat, distance=distance, **kwargs)
  elif model.lower() == 'arenou':
    av = findext_arenou(lng, lat, distance=distance, **kwargs)
  elif model.lower() == 'schlegel':
    av = findext_schlegel(lng, lat, distance=distance, **kwargs)
  return(av)

#}  
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#{ Arenou 3D extinction model
# (based on Arenou et al, "A tridimensional model of the 
# galactic interstellar extinction" published in Astronomy and
# Astrophysics (ISSN 0004-6361), vol. 258, no. 1, p. 104-111.)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def findext_arenou(ll, bb, distance=None):
  """
  Find the predicted V-band extinction (Av) according to the
  3D model for galactic extinction of Arenou et al, "Modelling
  the Galactic interstellar extinction distribution in three
  dimensions", Arenou et al, "A tridimensional model of the 
  galactic interstellar extinction" published in Astronomy and
  Astrophysics (ISSN 0004-6361), vol. 258, no. 1, p. 104-111, 1992
  
  Example usage:
    
    1. Find the extinction for a star at galactic lattitude 10.2 and
       longitude 59.0. If the distance is not given, we use the
       maximal distance r0 for that line of sight, as given in the
       Appendix of Arenou et al. (1992)
       
        >>> lng = 10.2
        >>> lat = 59.0
        >>> av = findext_arenou(lng, lat)
        >>> print("Av at lng = %.2f, lat = %.2f is %.2f magnitude" %(lng, lat, av))
        Av at lng = 10.20, lat = 59.00 is 0.05 magnitude
    
    2. Find the extinction for a star at galactic lattitude 107.05 and
       longitude -34.93 and a distance of 144.65 parsecs
    
        >>> lng = 107.05
        >>> lat = -34.93
        >>> dd  = 144.65
        >>> av = findext_arenou(lng, lat, distance = dd)
        >>> print("Av at lng = %.2f, lat = %.2f and distance = %.2f parsecs is %.2f magnitude" %(lng, lat, dd, av))
        Av at lng = 107.05, lat = -34.93 and distance = 144.65 parsecs is 0.15 magnitude
        
  @param ll: Galactic Longitude (in degrees)
  @type ll: float
  @param bb: Galactic Lattitude (in degrees)
  @type bb: float
  @keyword distance: Distance to the source (in parsecs)
  @type distance: float
  @return: The extinction in Johnson V-band
  @rtype: float
  """
  # make sure that the values for b and l are within the correct range
  if bb < -90. or bb > 90:
    logger.error("galactic lattitude outside [-90,90] degrees")
  elif ll < 0. or ll > 360:
    logger.error("galactic longitude outside [0,360] degrees")
  elif distance < 0 and distance is not None:
    logger.error("distance is negative")
    
  # find the Arenou paramaters in the Appendix of Arenou et al. (1992)
  alpha, beta, gamma, rr0, saa = _getarenouparams(ll, bb)
  logger.info("Arenou params: alpha = %.2f, beta = %.2f, gamma = %.2f, r0 = %.2f and saa = %.2f" %(alpha, beta, gamma, rr0, saa))
  
  # compute the visual extinction from the Arenou paramaters using Equation 5
  # and 5bis
  if distance is None:
    av = alpha*rr0 + beta*rr0**2.
  else:
    distance = distance/1e3 # to kparsec
    if distance <= rr0:
      av = alpha*distance + beta*distance**2.
    else:
      av = alpha*rr0 + beta*rr0**2. + (distance-rr0)*gamma
  return av

def _getarenouparams(ll,bb):
  """
  Input: galactic coordinates
  Output: Arenou 1992 alpha, beta, gamma, rr0, saa
  """
  if -90 < bb < -60:
    gamma = 0
    if   0   <= ll < 29:
      alpha = 2.22538 ; beta = -6.00212 ; rr0 = 0.052 ; saa = 13.
    elif 29  <= ll < 57:
      alpha = 3.35436 ; beta = -14.74567 ; rr0 = 0.035 ; saa = 40
    elif 57  <= ll < 85:
      alpha = 2.77637 ; beta = -9.62706 ; rr0 = 0.042 ; saa = 15
    elif 85  <= ll < 110:
      alpha = 4.44190 ; beta = -19.92097 ; rr0 = 0.025 ; saa = 36
    elif 110 <= ll < 150:
      alpha = 4.46685 ; beta = -26.07305 ; rr0 = 0.026 ; saa = 28
    elif 150 <= ll < 180:
      alpha = 7.63699 ; beta = -46.10856  ; rr0 = 0.014 ; saa = 38
    elif 180 <= ll < 210:
      alpha = 2.43412 ; beta = -8.69913 ; rr0 = 0.050 ; saa = 36
    elif 210 <= ll < 240:
      alpha = 3.34481 ; beta = -13.93228 ; rr0 = 0.035 ; saa = 38
    elif 240 <= ll < 270:
      alpha = 1.40733 ; beta = -13.43418  ; rr0 = 0.091 ; saa = 30
    elif 270 <= ll < 300:
      alpha = 1.64466 ; beta = -3.97380 ; rr0 = 0.074 ; saa = 28
    elif 300 <= ll < 330:
      alpha = 2.12696 ; beta = -6.05682 ; rr0 = 0.056 ; saa = 14
    elif 330 <= ll < 360:
      alpha = 2.34636 ; beta = -8.17784 ; rr0 = 0.052 ; saa = 16
  
  if -60 < bb < -45:
    gamma = 0
    if   0  <= ll < 30:  
      alpha = 2.77060 ; beta = -9.52310 ; rr0 = 0.145 ; saa = 16
    elif 30 <= ll < 60:
      alpha = 1.96533 ; beta = -9.52310 ; rr0 = 0.174 ; saa = 06
    elif 60 <= ll < 110:
      alpha = 1.93622 ; beta = -13.31757 ; rr0 = 0.073 ; saa = 26
    elif 110 <= ll < 180:
      alpha = 1.05414 ; beta = -2.236540 ; rr0 = 0.223 ; saa = 74
    elif 180 <= ll < 210:
      alpha = 1.39990 ; beta = -1.35325 ; rr0 = 0.252 ; saa = 10
    elif 210 <= ll < 240:
      alpha = 2,73481 ; beta = -11.70266 ; rr0 = 0.117 ; saa = 8
    elif 240 <= ll < 270:
      alpha = 2.99784 ; beta = -11.64272 ; rr0 = 0.129 ; saa = 3
    elif 270 <= ll < 300:
      alpha = 3.23802 ; beta = -11.63810 ; rr0 = 0.139 ; saa = 7
    elif 300 <= ll < 330:
      alpha = 1.72679 ; beta = -6.05085 ; rr0 = 0.143 ; saa = 7
    elif 330 <= ll < 360:
      alpha = 1.88890 ; beta = -5.51861 ; rr0 = 0.171 ; saa = 14
  
  if -45 < bb < -30:
    gamma = 0
    if 0  <= ll < 30:
      alpha = 1.98973 ; beta = -4.86206  ; rr0 = 0.205 ; saa = 6 
    elif 30  <= ll < 60:
      alpha = 1.49901 ; beta = -3.75837  ; rr0 = 0.199 ; saa = 16 
    elif 60  <= ll < 90:
      alpha = 0.90091 ; beta = -1.30459  ; rr0 = 0.329 ; saa = 73 
    elif 90  <= ll < 120:
      alpha = 1.94200 ; beta = -6.26833  ; rr0 = 0.155 ; saa = 18 
    elif 120 <= ll < 160:
      alpha = -0.37804 ; beta = 10.77372  ; rr0 = 0.210 ; saa = 100
    elif 160 <= ll < 200:
      alpha = -0.15710 ; beta = 5.03190   ; rr0 = 0.294 ; saa = 24
    elif 200 <= ll < 235:
      alpha = 3.20162 ; beta = -10.59297 ; rr0 = 0.151 ; saa = 9
    elif 235 <= ll < 265:
      alpha = 1.95079 ; beta = 4.73280   ; rr0 = 0.206 ; saa = 21
    elif 265 <= ll < 300:
      alpha = 1.91200 ; beta = 4.97640   ; rr0 = 0.192 ; saa = 17
    elif 300 <= ll < 330:
      alpha = 2.50487 ; beta = -9.63106  ; rr0 = 0.145 ; saa = 28
    elif 330 <= ll < 360:
      alpha = 2.44394 ; beta = -9.17612  ; rr0 = 0.133 ; saa = 7
  
  if -30 < bb < -15:
    gamma = 0
    if 0 <= ll < 20:
      alpha = 2.82440 ; beta = -4.78217 ; rr0 = 0.295 ; saa = 32
    elif 20  <= ll < 40:
      alpha = 3.84362 ; beta = -8.04690 ; rr0 = 0.239 ; saa = 46 
    elif 40  <= ll < 80:
      alpha = 0.60365 ; beta = 0.07893  ; rr0 = 0.523 ; saa = 22 
    elif 80  <= ll < 100: 
      alpha = 0.58307 ; beta = -0.21053 ; rr0 = 0.523 ; saa = 53 
    elif 100 <= ll < 120: 
      alpha = 2.03861 ; beta = -4.40843 ; rr0 = 0.231 ; saa = 60 
    elif 120 <= ll < 140:
      alpha = 1.14271 ; beta = -1.35635 ; rr0 = 0.421 ; saa = 34 
    elif 140 <= ll < 160:
      alpha = 0.79908 ; beta = 1.48074  ; rr0 = 0.513 ; saa = 61 
    elif 160 <= ll < 180:
      alpha = 0.94260 ; beta = 8.16346  ; rr0 = 0.441 ; saa = 42 
    elif 180 <= ll < 200:
      alpha = 1.66398 ; beta = 0.26775  ; rr0 = 0.523 ; saa = 42 
    elif 200 <= ll < 220: 
      alpha = 1.08760 ; beta = -1.02443 ; rr0 = 0.523 ; saa = 45 
    elif 220 <= ll < 240:
      alpha = 1.20087 ; beta = -2.45407 ; rr0 = 0.245 ; saa = 6 
    elif 240 <= ll < 260:
      alpha = 1.13147 ; beta = -1.87916 ; rr0 = 0.301 ; saa = 16 
    elif 260 <= ll < 280:
      alpha = 0.97804 ; beta = -2.92838 ; rr0 = 0.338 ; saa = 21 
    elif 290 <= ll < 300:
      alpha = 1.40086 ; beta = -1.12403 ; rr0 = 0.523 ; saa = 19 
    elif 300 <= ll < 320:
      alpha = 2.06355 ; beta = -3.68278 ; rr0 = 0.280 ; saa = 42 
    elif 320 <= ll < 340: 
      alpha = 1.59260 ; beta = -2.18754 ; rr0 = 0.364 ; saa = 15 
    elif 310 <= ll < 360:
      alpha = 1.45589 ; beta = -1.90598 ; rr0 = 0.382 ; saa = 21
  
  if -15 < bb < -5:
    gamma = 0
    if 0   <= ll < 10:
      alpha = 2.56438 ; beta = -2.31586 ; rr0 = 0.554 ; saa = 37 
    elif 10  <= ll < 20:
      alpha = 3.24095 ; beta = -2.78217 ; rr0 = 0.582 ; saa = 38 
    elif 20  <= ll < 30:
      alpha = 2.95627 ; beta = -2.57422 ; rr0 = 0.574 ; saa = 41 ; gamma = 0.08336
    elif 30  <= ll < 40:
      alpha = 1.85158 ; beta = -0.67716 ; rr0 = 1.152 ; saa = 4 
    elif 40  <= ll < 50:
      alpha = 1.60770 ; beta = 0.35279  ; rr0 = 0.661 ; saa = 24 
    elif 50  <= ll < 60:
      alpha = 0.69920 ; beta = -0.09146 ; rr0 = 0.952 ; saa = 2  ; gamma = 0.12839
    elif 60  <= ll < 80:
      alpha = 1.36189 ; beta = -1.05290 ; rr0 = 0.647 ; saa = 45 ; gamma = 0.16258
    elif 80  <= ll < 90:
      alpha = 0.33179 ; beta = 0.37629  ; rr0 = 1.152 ; saa = 62 
    elif 90  <= ll < 100:
      alpha = 1.70303 ; beta = -0.75246 ; rr0 = 1.132 ; saa = 31 
    elif 100 <= ll < 110:
      alpha = 1.97414 ; beta = -1.59784 ; rr0 = 0.618 ; saa = 35 ; gamma = 0.12847
    elif 110 <= ll < 120:
      alpha = 1.07407 ; beta = -0.40066 ; rr0 = 1.152 ; saa = 14 ; gamma = 0.17698
    elif 120 <= ll < 130:
      alpha = 1.69495 ; beta = -1.00071 ; rr0 = 0.847 ; saa = 28 ; gamma = 0.08567
    elif 130 <= ll < 140:
      alpha = 1.51449 ; beta = -0.08441 ; rr0 = 0.897 ; saa = 12 
    elif 140 <= ll < 150:
      alpha = 1.87894 ; beta = -0.73314 ; rr0 = 1.152 ; saa = 23 
    elif 150 <= ll < 160:
      alpha = 1.43670 ; beta = 0.67706  ; rr0 = 0.778 ; saa = 3  ; gamma = 0.42624
    elif 160 <= ll < 180:
      alpha = 6.84802 ; beta = -5.06864 ; rr0 = 0.676 ; saa = 50 
    elif 180 <= ll < 190:
      alpha = 4.16321 ; beta = -5.80016 ; rr0 = 0.359 ; saa = 51 ; gamma = 0.60387
    elif 190 <= ll < 200:
      alpha = 0.78135 ; beta = -0.27826 ; rr0 = 1.152 ; saa = 4 
    elif 200 <= ll < 210:
      alpha = 0.85535 ; beta = 0.20848  ; rr0 = 1.152 ; saa = 17 
    elif 210 <= ll < 220:
      alpha = 0.52521 ; beta = 0.65726  ; rr0 = 1.152 ; saa = 7 
    elif 220 <= ll < 230:
      alpha = 0.88376 ; beta = -0.44519 ; rr0 = 0.993 ; saa = 40 ; gamma = 0.06013
    elif 230 <= ll < 240:
      alpha = 0.42228 ; beta = 0.26304  ; rr0 = 0.803 ; saa = 26 
    elif 240 <= ll < 250:
      alpha = 0.71318 ; beta = -0.67229 ; rr0 = 0.530 ; saa = 55 ; gamma = 0.20994
    elif 250 <= ll < 260:
      alpha = 0.99606 ; beta = -0.70103 ; rr0 = 0.710 ; saa = 48 ; gamma = 0.01323
    elif 260 <= ll < 270:
      alpha = 0.91519 ; beta = -0.39690 ; rr0 = 1.152 ; saa = 48 ; gamma = 0.01961
    elif 270 <= ll < 280:
      alpha = 0.85791 ; beta = -0.29115 ; rr0 = 1.152 ; saa = 19 
    elif 280 <= ll < 290:
      alpha = 1.44226 ; beta = -1.09775 ; rr0 = 0.657 ; saa = 39 
    elif 290 <= ll < 300:
      alpha = 2.55486 ; beta = -1.68293 ; rr0 = 0.759 ; saa = 31 
    elif 300 <= ll < 310:
      alpha = 3.18047 ; beta = -2.69796 ; rr0 = 0.589 ; saa = 40 
    elif 210 <= ll < 320:
      alpha = 2.11235 ; beta = -1.77506 ; rr0 = 0.595 ; saa = 29 
    elif 320 <= ll < 330:
      alpha = 1.75884 ; beta = -1.38574 ; rr0 = 0.635 ; saa = 25 
    elif 330 <= ll < 340:
      alpha = 1.97257 ; beta = -1.55545 ; rr0 = 0.634 ; saa = 34 ; gamma = 0.00043
    elif 340 <= ll < 350:
      alpha = 1.41497 ; beta = -1.05722 ; rr0 = 0.669 ; saa = 46 ; gamma = 0.03264
    elif 350 <= ll < 360:
      alpha = 1.17795 ; beta = -0.95012 ; rr0 = 0.620 ; saa = 46 ; gamma = 0.03339
  
  if -5 < bb < 5:
    gamma = 0
    if 0   <= ll < 10:
      alpha = 2.62556 ; beta = -1.11097 ; rr0 = 1.182 ; saa = 57 ; gamma = 0.00340
    elif 10  <= ll < 20:
      alpha = 3.14461 ; beta = -1.01140 ; rr0 = 1.555 ; saa = 42
    elif 20  <= ll < 30:
      alpha = 4.26624 ; beta = -1.61242 ; rr0 = 1.323 ; saa = 34  
    elif 30  <= ll < 40:
      alpha = 2.54447 ; beta = -0.12771 ; rr0 = 1.300 ; saa = 30  
    elif 40  <= ll < 50:
      alpha = 2.27030 ; beta = -0.68720 ; rr0 = 1.652 ; saa = 52 ; gamma = 0.04928 
    elif 50  <= ll < 60:
      alpha = 1.34359 ; beta = -0.05416 ; rr0 = 2.000 ; saa = 32  
    elif 60  <= ll < 70:
      alpha = 1.76327 ; beta = -0.26407 ; rr0 = 2.000 ; saa = 37  
    elif 70  <= ll < 80:
      alpha = 2.20666 ; beta = -0.41651 ; rr0 = 2.000 ; saa = 36  
    elif 80  <= ll < 90:
      alpha = 1.50130 ; beta = -0.09855 ; rr0 = 1.475 ; saa = 45  
    elif 90  <= ll < 100:
      alpha = 2.43965 ; beta = -0.82128 ; rr0 = 1.485 ; saa = 36 ; gamma = 0.01959 
    elif 100 <= ll < 110:
      alpha = 3.35775 ; beta = -1.16400 ; rr0 = 0.841 ; saa = 35 ; gamma = 0.00298 
    elif 110 <= ll < 120:
      alpha = 2.60621 ; beta = -0.68687 ; rr0 = 1.897 ; saa = 36  
    elif 120 <= ll < 130:
      alpha = 2.90112 ; beta = -0.97988 ; rr0 = 1.480 ; saa = 32  
    elif 130 <= ll < 140:
      alpha = 2.55377 ; beta = -0.71214 ; rr0 = 1.793 ; saa = 38  
    elif 140 <= ll < 150:
      alpha = 3.12598 ; beta = -1.21437 ; rr0 = 1.287 ; saa = 23 ; gamma = 0.15298 
    elif 150 <= ll < 160:
      alpha = 3.66930 ; beta = -2.29731 ; rr0 = 0.799 ; saa = 40 ; gamma = 0.33473 
    elif 160 <= ll < 170:
      alpha = 2.15465 ; beta = -0.70690 ; rr0 = 1.524 ; saa = 37 ; gamma = 0.14017 
    elif 170 <= ll < 180:
      alpha = 1.82465 ; beta = -0.60223 ; rr0 = 1.515 ; saa = 29 ; gamma = 0.20730 
    elif 180 <= ll < 190:
      alpha = 1.76269 ; beta = -0.35945 ; rr0 = 2.000 ; saa = 28 ; gamma = 0.08052 
    elif 190 <= ll < 200:
      alpha = 1.06085 ; beta = -0.14211 ; rr0 = 2.000 ; saa = 48  
    elif 200 <= ll < 210:
      alpha = 1.21333 ; beta = -0.23225 ; rr0 = 2.000 ; saa = 57  
    elif 210 <= ll < 220:
      alpha = 0.58326 ; beta = -0.06097 ; rr0 = 2.000 ; saa = 30 ; gamma = 0.36962 
    elif 220 <= ll < 230:
      alpha = 0.74200 ; beta = -0.19293 ; rr0 = 1.923 ; saa = 48 ; gamma = 0.07459  
    elif 230 <= ll < 240:
      alpha = 0.67520 ; beta = -0.21041 ; rr0 = 1.604 ; saa = 49 ; gamma = 0.16602 
    elif 240 <= ll < 250:
      alpha = 0.62609 ; beta = -0.25312 ; rr0 = 1.237 ; saa = 73 ; gamma = 0.14437 
    elif 250 <= ll < 260:
      alpha = 0.61415 ; beta = -0.13788 ; rr0 = 2.000 ; saa = 42 ; gamma = 0.26859 
    elif 260 <= ll < 270:
      alpha = 0.58108 ; beta = 0.01195  ; rr0 = 2.000 ; saa = 40 ; gamma = 0.07661 
    elif 270 <= ll < 280:
      alpha = 0.68352 ; beta = -0.10743 ; rr0 = 2.000 ; saa = 50 ; gamma = 0.00849  
    elif 280 <= ll < 290:
      alpha = 0.61747 ; beta = 0.02675  ; rr0 = 2,000 ; saa = 49  
    elif 290 <= ll < 300:
      alpha = 0.06827 ; beta = -0.26290 ; rr0 = 2.000 ; saa = 44  
    elif 300 <= ll < 310:
      alpha = 1.53631 ; beta = -0.36833 ; rr0 = 2.000 ; saa = 37 ; gamma = 0.02960  
    elif 210 <= ll < 320:
      alpha = 1.94300 ; beta = -0.71445 ; rr0 = 1.360 ; saa = 36 ; gamma = 0.15643  
    elif 320 <= ll < 330:
      alpha = 1.22185 ; beta = -0.40185 ; rr0 = 1.520 ; saa = 48 ; gamma = 0.07354  
    elif 330 <= ll < 340:
      alpha = 1.79429 ; beta = -0.48657 ; rr0 = 1.844 ; saa = 43 
    elif 340 <= ll < 350:
      alpha = 2.29545 ; beta = -0.84096 ; rr0 = 1.365 ; saa = 32  
    elif 350 <= ll < 360:
      alpha = 2.07408 ; beta = -0.64745 ; rr0 = 1.602 ; saa = 36 ; gamma = 0.12750
  
  
  if 5 < bb < 15:
    gamma = 0
    if 0   <= ll < 10:
      alpha = 2.94213 ; beta = -2.09158 ; rr0 = 0.703 ; saa = 41 ; gamma = 0.05490 
    elif 10  <= ll < 30:
      alpha = 3.04627 ; beta = 7.71159  ; rr0 = 0.355 ; saa = 45 
    elif 30  <= ll < 40:
      alpha = 3.78033 ; beta = -3.91956 ; rr0 = 0.482 ; saa = 42
    elif 40  <= ll < 50:
      alpha = 2.18119 ; beta = -2.4050  ; rr0 = 0.453 ; saa = 27 
    elif 50  <= ll < 60:
      alpha = 1.45372 ; beta = -0.49522 ; rr0 = 1.152 ; saa = 31 
    elif 60  <= ll < 70:
      alpha = 1.05051 ; beta = -1.01704 ; rr0 = 0.516 ; saa = 2
    elif 70  <= ll < 80:
      alpha = 0.48416 ; beta = -0.27182 ; rr0 = 0.891 ; saa = 94 ; gamma = 0.08639 
    elif 80  <= ll < 90:
      alpha = 0.61963 ; beta = 0.41697  ; rr0 = 1.152 ; saa = 35 ; gamma = 0.47171 
    elif 90  <= ll < 100:
      alpha = 4.40348 ; beta = -2.95011 ; rr0 = 0.745 ; saa = 52  
    elif 100 <= ll < 120:
      alpha = 2.50938 ; beta = -0.56541 ; rr0 = 1.152 ; saa = 27  
    elif 120 <= ll < 130:
      alpha = 0.44180 ; beta = 1.58923  ; rr0 = 0.949 ; saa = 4  
    elif 130 <= ll < 140:
      alpha = 3.96081 ; beta = -3.37374 ; rr0 = 0.587 ; saa = 40 ; gamma = 0.34109 
    elif 140 <= ll < 160:
      alpha = 2.53335 ; beta = -0.40541 ; rr0 = 1.152 ; saa = 38  
    elif 160 <= ll < 170:
      alpha = 2.03760 ; beta = -0.66136 ; rr0 = 1.152 ; saa = 23 
    elif 170 <= ll < 200:
      alpha = 1.06946 ; beta = -0.87395 ; rr0 = 0.612 ; saa = 29 ; gamma = 0.29230
    elif 200 <= ll < 210:
      alpha = 0.86348 ; beta = -0.65870 ; rr0 = 0.655 ; saa = 79 ; gamma = 0.09089
    elif 210 <= ll < 230:
      alpha = 0.30117 ; beta = -0.16136 ; rr0 = 0.933 ; saa = 17 ; gamma = 0.07495 
    elif 230 <= ll < 240:
      alpha = 0.75171 ; beta = -0.57143 ; rr0 = 0.658 ; saa = 12 ; gamma = 0.00534  
    elif 240 <= ll < 250:
      alpha = 1.97427 ; beta = -2.02654 ; rr0 = 0.487 ; saa = 67  
    elif 250 <= ll < 260:
      alpha = 1.25208 ; beta = -1.47763 ; rr0 = 0.424 ; saa = 19 ; gamma = 0.09089  
    elif 260 <= ll < 270:
      alpha = 0.89448 ; beta = -0.43870 ; rr0 = 1.019 ; saa = 5  
    elif 270 <= ll < 280:
      alpha = 0.81141 ; beta = -0.51001 ; rr0 = 0.795 ; saa = 27 ; gamma = 0.03505 
    elif 280 <= ll < 290:
      alpha = 0.83781 ; beta = -0.44138 ; rr0 = 0.949 ; saa = 50 ; gamma = 0.02820 
    elif 290 <= ll < 300:
      alpha = 1.10600 ; beta = -0.86263 ; rr0 = 0.641 ; saa = 28 ; gamma = 0.03402 
    elif 300 <= ll < 310:
      alpha = 1.37040 ; beta = -1.02779 ; rr0 = 0.667 ; saa = 28 ; gamma = 0.05608 
    elif 310 <= ll < 320:
      alpha = 1.77590 ; beta = -1.26951 ; rr0 = 0.699 ; saa = 37 ; gamma = 0.06972 
    elif 320 <= ll < 330:
      alpha = 1.20865 ; beta = -0.70679 ; rr0 = 0.855 ; saa = 35 ; gamma = 0.02902 
    elif 330 <= ll < 340:
      alpha = 2.28830 ; beta = -1.71890 ; rr0 = 0.666 ; saa = 42 ; gamma = 0.22887 
    elif 340 <= ll < 350:
      alpha = 3.26278 ; beta = -0.94181 ; rr0 = 1.152 ; saa = 38 
    elif 350 <= ll < 360:
      alpha = 2.58100 ; beta = -1.69237 ; rr0 = 0.763 ; saa = 53
      
  if 15 < bb < 30:
    gamma = 0
    if 0   <= ll < 20:
      alpha = 6.23279  ; beta = -10.30384 ; rr0 = 0.302 ; saa = 42  
    elif 20  <= ll < 40:
      alpha = -4.47693 ; beta = -7.28366  ; rr0 = 0.307 ; saa = 29
    elif 40  <= ll < 60 :
      alpha =  1.22938 ; beta = -1.19030  ; rr0 = 0.516 ; saa = 5  
    elif 60  <= ll < 80 :
      alpha =  0.84291 ; beta = -1.59338  ; rr0 = 0.265 ; saa = 4  
    elif 80  <= ll < 100 :
      alpha =  0.23996 ; beta = 0.06304   ; rr0 = 0.523 ; saa = 32 
    elif 100 <= ll < 140 :
      alpha =  0.40062 ; beta = -1.75628  ; rr0 = 0.114 ; saa = 16 
    elif 140 <= ll < 180 :
      alpha =  0.56898 ; beta = -0.53331  ; rr0 = 0.523 ; saa = 41 
    elif 180 <= ll < 200 :
      alpha = -0.95721 ; beta = 11.6917   ; rr0 = 0.240 ; saa = 2 
    elif 200 <= ll < 220 :
      alpha = -0.19051 ; beta = 1.45670   ; rr0 = 0.376 ; saa = 1
    elif 220 <= ll < 240 :
      alpha =  2.31305 ; beta = -7.82531  ; rr0 = 0.148 ; saa = 95 
    elif 240 <= ll < 260:
      alpha =  1.39169 ; beta = -1.72984  ; rr0 = 0.402 ; saa = 6 
    elif 260 <= ll < 260:
      alpha =  1.59418 ; beta = -1.28296  ; rr0 = 0.523 ; saa = 36 
    elif 280 <= ll < 300 :
      alpha =  1.57082 ; beta = -197295   ; rr0 = 0.398 ; saa = 10  
    elif 300 <= ll < 320 :
      alpha =  1.95998 ; beta = -3.26159  ; rr0 = 0.300 ; saa = 11 
    elif 320 <= ll < 340:
      alpha =  2.59567 ; beta = -4.84133  ; rr0 = 0.168 ; saa = 37  
    elif 340 <= ll < 360:
      alpha =  5.30273 ; beta = -7.43033  ; rr0 = 0.357 ; saa = 37
  
  if 30 < bb < 45:
    gamma = 0
    if 0   <= ll < 20:
      alpha =  2.93960 ; beta  = -6.48019  ; rr0 = 0.227 ; saa = 77  
    elif 20  <= ll < 50:
      alpha =  1.65864 ; beta  = -9.99317  ; rr0 = 0.083 ; saa = 99  
    elif 50  <= ll < 80:
      alpha =  1.71831 ; beta  = -7.25286  ; rr0 = 0.118 ; saa = 28  
    elif 80  <= ll < 110:
      alpha =  1.33617 ; beta  = -10.39799 ; rr0 = 0.064 ; saa = 99  
    elif 110 <= ll < 160:
      alpha = -0.31330 ; beta  = 1.35622   ; rr0 = 0.329 ; saa = 24
    elif 160 <= ll < 190:
      alpha =  1.51984 ; beta  = -8.69502  ; rr0 = 0.087 ; saa = 99 
    elif 190 <= ll < 220:
      alpha = -0.50758 ; beta  = 4.73320   ; rr0 = 0.250 ; saa = 78 
    elif 220 <= ll < 250:
      alpha =  1.25864 ; beta  = -12.59627 ; rr0 = 0.050 ; saa = 70  
    elif 250 <= ll < 280:
      alpha =  1.54243 ; beta  = -3.75065  ; rr0 = 0.205 ; saa = 10 
    elif 280 <= ll < 320:
      alpha =  2.72258 ; beta  = -7.47806  ; rr0 = 0.182 ; saa = 5 
    elif 320 <= ll < 340:
      alpha =  2.81435 ; beta  = -5.52139  ; rr0 = 0.255 ; saa = 10 
    elif 340 <= ll < 360:
      alpha =  2.23818 ; beta  = 0.81772   ; rr0 = 0.329 ; saa = 19 
  
  if 45 < bb < 60:
    gamma = 0
    if 0   <= ll < 60:
      alpha = 1.38587 ; beta  = -9.06536  ; rr0 = 0.076 ; saa = 3 
    elif 60  <= ll < 90:
      alpha = 2.28570 ; beta  = -9.88812  ; rr0 = 0.116 ; saa = 3 
    elif 90  <= ll < 110:
      alpha = 1.36385 ; beta  = -8.10127  ; rr0 = 0.084 ; saa = 4 
    elif 110 <= ll < 170:
      alpha = 0.05943 ; beta  = -1.08126  ; rr0 = 0.027 ; saa = 50  
    elif 170 <= ll < 200:
      alpha = 1.40171 ; beta  = -3.21783  ; rr0 = 0.218 ; saa = 99  
    elif 200 <= ll < 230:
      alpha = 0.14718 ; beta  = 3.92670   ; rr0 = 0.252 ; saa = 14  
    elif 230 <= ll < 290:
      alpha = 0.57124 ; beta  = -4.30242  ; rr0 = 0.066 ; saa = 10 
    elif 290 <= ll < 330:
      alpha = 3.69891 ; beta  = -19.62204 ; rr0 = 0.094 ; saa = 5 
    elif 330 <= ll < 360:
      alpha = 1.19563 ; beta  = -0.45043  ; rr0 = 0.252 ; saa = 9
  
  if 60 < bb < 90:
    gamma = 0
    if 0   <= ll < 30:
      alpha = 0.69443 ;  beta = -0.27600  ; rr0 = 0.153 ; saa = 99
    elif 30  <= ll < 60:
      alpha = 1.11811 ;  beta = 0.71179   ; rr0 = 0.085 ; saa = 73  
    elif 60  <= ll < 90:
      alpha = 1.10427 ;  beta = -2.37654  ; rr0 = 0.123 ; saa = 99 
    elif 90  <= ll < 120:
      alpha = -0.42211 ; beta = 5.24037   ; rr0 = 0.184 ; saa = 12 
    elif 120 <= ll < 150:
      alpha = 0.87576 ;  beta = -4.38033  ; rr0 = 0.100 ; saa = 35  
    elif 150 <= ll < 180:
      alpha = 1.27477 ;  beta = -4.98307  ; rr0 = 0.128 ; saa = 72 
    elif 180 <= ll < 210:
      alpha = 1.19512 ;  beta = -6.58464  ; rr0 = 0.091 ; saa = 49 
    elif 210 <= ll < 240:
      alpha = 0.97581 ;  beta = -4.89869  ; rr0 = 0.100 ; saa = 95 
    elif 240 <= ll < 270:
      alpha = 0.54379 ;  beta = -0.84403  ; rr0 = 0.207 ; saa = 35 
    elif 270 <= ll < 300:
      alpha = 0.85054 ;  beta = 13.01249  ; rr0 = 0.126 ; saa = 39
    elif 300 <= ll < 330:
      alpha = 0.74347 ;  beta = 1.39825   ; rr0 = 0.207 ; saa = 10
    elif 330 <= ll < 360:
      alpha = 0.77310 ;  beta = -4.45005  ; rr0 = 0.087 ; saa = 16 

  return alpha, beta, gamma, rr0, saa

#}
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#{ Marshall 3D extinction model
# (Marshall et al. (2006) published in Astronomy and Astrophysics,
#  Volume 453, Issue 2, July II 2006, pp.635-651
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@memoized
def get_marshall_data():
  # Read in the Marshall data
  data_ma, units_ma, comments_ma = vizier.search("J/A+A/453/635")
  return data_ma, units_ma, comments_ma

def findext_marshall(ll, bb, distance=None, redlaw='cardelli1989', Rv=3.1,  **kwargs):
  """
  Find the V-band extinction according to the reddening model of
  Marshall et al. (2006) published in Astronomy and Astrophysics,
  Volume 453, Issue 2, July II 2006, pp.635-651
  
  Example usage:
    
    1. Find the extinction for a star at galactic longitude 10.2 and
       latitude 9.0. If the distance is not given, we return the
       complete extinction along the line of sight (i.e. put the star somewhere
       out of the galaxy).
       
        >>> lng = 10.2
        >>> lat = 9.0
        >>> ak = findext_marshall(lng, lat)
        >>> print("Ak at lng = %.2f, lat = %.2f is %.2f magnitude" %(lng, lat, ak))
        Ak at lng = 10.20, lat = 9.00 is 0.15 magnitude
    
    2. Find the extinction for a star at galactic lattitude 271.05 and
       longitude -4.93 and a distance of 144.65 parsecs
    
        >>> lng = 271.05
        >>> lat = -4.93
        >>> dd  = 144.65
        >>> ak = findext_marshall(lng, lat, distance = dd)
        >>> print("Ak at lng = %.2f, lat = %.2f and distance = %.2f parsecs is %.2f magnitude" %(lng, lat, dd, ak))
        Ak at lng = 271.05, lat = -4.93 and distance = 144.65 parsecs is 0.02 magnitude
        
    3. The extinction for galactic longitude 10.2 and latitude 59.0 can not
         be calculated using the Marshall models, because they are only valid at
         certain lng and lat (0 < lng < 100 or 260 < lng < 360, -10 < lat < 10)
       
        >>> lng = 10.2
        >>> lat = 59.0
        >>> ak = findext_marshall(lng, lat)
        >>> print(av)
        None
        

  @param ll: Galactic Longitude (in degrees) should be between 0 and 100 or 260 and 360 degrees
  @type ll: float
  @param bb: Galactic Lattitude (in degrees) should be between -10 and 10 degrees
  @type bb: float
  @keyword distance: Distance to the source (in parsecs)
  @type distance: float
  @keyword redlaw: the used reddening law (standard: 'cardelli1989')
  @type redlaw: str
  @keyword Rv: Av/E(B-V) (standard: 3.1)
  @type Rv: float
  @return: The extinction in K-band
  @rtype: float
  """
  
  # get Marshall data
  data_ma, units_ma, comments_ma = get_marshall_data()
  
  # Check validity of the coordinates 
  if (ll > 100. and ll < 260.) or ll < 0 or ll > 360:
    Av = None
    logger.error("Galactic longitude invalid")
    return Av
  elif bb > 10. or bb < -10.:
    Av = None
    logger.error("Galactic lattitude invalid")
    return Av
  
  # Find the galactic lattitude and longitude of the model, closest to your star
  dist = np.sqrt((data_ma.GLAT - bb)**2. + (data_ma.GLON - ll)**2.)
  kma  = np.where(dist == dist.min())
  
  if min(dist) > .5:
    logger.error("Could not find a good model value")
    return Av
  
  # find the correct index for the distance
  nb    = data_ma.nb[kma]
  rr = [] ; e_rr = [] ; ext = [] ; e_ext = []
  for i in np.arange(nb)+1.:
    rr.append(data_ma["r%i"%i][kma])
    e_rr.append(data_ma["e_r%i"%i][kma])
    ext.append(data_ma["ext%i"%i][kma])
    e_ext.append(data_ma["e_ext%i"%i][kma])
  rr = np.array(rr)[:,0] ; e_rr = np.array(e_rr)[:,0] ; ext = np.array(ext)[:,0] ; e_ext = np.array(e_ext)[:,0]  
  
  # Interpolate linearly in distance. If beyond furthest bin, keep that value.
  if distance is None:
    logger.info("No distance given")
    dd = max(rr)
    ak = ext[-1]
  else:
    dd = distance/1e3
    
  if dd < min(rr):
    logger.info("Distance below distance to first bin")
    ak = (dd/rr[0])*ext[0]
  elif dd > max(rr):
    logger.info("Distance more than distance to last bin")
    ak = ext[-1]
    dd = max(rr)
  else:
    logger.info("Interpolating linearly")
    ak = sc.interp(dd, rr, ext)
  
  redwave, redflux = get_law(redlaw,Rv=Rv,wave_units='micron',norm='Av', wave=array([0.54,2.22]))
  
  return ak*redflux[0]/redflux[1]

#}
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#{ Drimmel 3D extinction model presented in Drimmel et al. in
# Astronomy and Astrophysics, v.409, p.205-215 (2003) and
# Proceedings of the Gaia Symposium "The Three-Dimensional
# Universe with Gaia"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# DISCLAIMER: This software is based on the IDL scripts written at the
# Cosmology Data Analysis Center in support of the Cosmic Background Explorer
# (COBE) Project under NASA contract number NAS5-30750.
# This software may be used, copied, modified or redistributed so long as it
# is not sold and this disclaimer is distributed along with the software. If
# you modify the software please indicate your modifications in a prominent
# place in the source code. All routines are provided "as is" without any
# express or implied warranties whatsoever. All routines are distributed
# without guarantee of support. If errors are found in this code it is
# requested that you contact us by sending email to the address below to
# report the errors but we make no claims regarding timely fixes. This
# software has been used for analysis of COBE data but has not been validated
# and has not been used to create validated data sets of any type.
# Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.

avdisk      = pf.getdata(config.get_datafile('drimmel',"avdisk.fits"      ))
avloc       = pf.getdata(config.get_datafile('drimmel',"avloc.fits"       ))
avspir      = pf.getdata(config.get_datafile('drimmel',"avspir.fits"      ))
avdloc      = pf.getdata(config.get_datafile('drimmel',"avdloc.fits"      ))
avori       = pf.getdata(config.get_datafile('drimmel',"avori.fits"       ))
coordinates = pf.getdata(config.get_datafile('drimmel',"coordinates.fits" ))
avgrid      = pf.getdata(config.get_datafile('drimmel',"avgrid.fits"      ))
avori2      = pf.getdata(config.get_datafile('drimmel',"avori2.fits"      ))
rf_allsky   = pf.getdata(config.get_datafile('drimmel',"rf_allsky.fits"   ))
glat        = rf_allsky.glat
glng        = rf_allsky.glng
ncomp       = rf_allsky.ncomp
pnum        = rf_allsky.pnum
rfac        = rf_allsky.rfac

g2e         = array([[-0.054882486, -0.993821033, -0.096476249], [0.494116468, -0.110993846,  0.862281440], [-0.867661702, -0.000346354,  0.497154957]])

def findext_drimmel(lng, lat, distance=None, rescaling=True,**kwargs):
  """
  procedure to retrieve the absorption in V from three-dimensional 
  grids, based on the Galactic dust distribution of Drimmel & Spergel.
  
  Example usage:
    
    1. Find the extinction for a star at galactic longitude 10.2 and
       latitude 9.0. If the distance is not given, we return the
       complete extinction along the line of sight (i.e. put the star somewhere
       out of the galaxy).
       
        >>> lng = 10.2
        >>> lat = 9.0
        >>> av = findext_drimmel(lng, lat)
        >>> print("Av at lng = %.2f, lat = %.2f is %.2f magnitude" %(lng, lat, av))
        Av at lng = 10.20, lat = 9.00 is 1.24 magnitude
    
    2. Find the extinction for a star at galactic lattitude 107.05 and
       longitude -34.93 and a distance of 144.65 parsecs
    
        >>> lng = 271.05
        >>> lat = -4.93
        >>> dd  = 144.65
        >>> ak = findext_marshall(lng, lat, distance = dd)
        >>> print("Ak at lng = %.2f, lat = %.2f and distance = %.2f parsecs is %.2f magnitude" %(lng, lat, dd, ak))
        Ak at lng = 271.05, lat = -4.93 and distance = 144.65 parsecs is 0.02 magnitude


  @param lng: Galactic Longitude (in degrees)
  @type lng: float
  @param lat: Galactic Lattitude (in degrees)
  @type lat: float
  @keyword distance: Distance to the source (in parsecs)
  @type distance: float
  @keyword rescaling: Rescaling needed or not?
  @type rescaling: boolean
  @return: extinction in V band with/without rescaling
  @rtype: float
  """
  # Constants
  deg2rad = pi/180. # convert degrees to rads
  nsky    = 393216  # number of COBE pixels
  
  # Sun's coordinates (get from dprms)
  xsun    = -8.0
  zsun    = 0.015
  
  # if distance is not given, make it large and put it in kiloparsec
  if distance is None:
    d = 1e10
  else:
    d = distance/1e3
  
  # build skymaps of rescaling parameters for each component
  dfac       = ones(nsky)
  sfac       = ones(nsky)
  lfac       = ones(nsky)
  indx       = where(ncomp == 1)
  dfac[indx] = rfac[indx]
  indx       = where(ncomp == 2)
  sfac[indx] = rfac[indx]
  indx       = where(ncomp == 3)
  lfac[indx] = rfac[indx]
  
  # define abs
  num     = array(d).size
  out     = zeros(num)
  avloc   = zeros(num)
  abspir  = zeros(num)
  absdisk = zeros(num)
  
  l = lng*deg2rad # [radians]
  b = lat*deg2rad # [radians]
  
  # Now for UIDL code:
  # -find the index of the corresponding COBE pixel
  # dimensions of sixpack = 768 x 512 (= 393216)
  vectoarr       = arange(768*512)
  vectoarr.shape = (512,768) # necessary because I reduced the arrays to vectors.
  vectoarr       = vectoarr.transpose()
  res            = 9
  pxindex        = _ll2pix(lng, lat, res)
  xout, yout     = _pix2xy(pxindex, res, sixpack=True)
  tblindex       = vectoarr[xout, yout] # calculate the maximum distance in the grid
  dmax           = ones(num)*100.
  if b != 0.:
    dmax = .49999/abs(sin(b)) - zsun/sin(b)
  if cos(l) != 0.:
    dmax = min([dmax,(14.9999/abs(cos(l)) - xsun/cos(l))])
  if sin(l) != 0.:
    dmax = min([dmax, 14.9999/abs(sin(l))])
  
  # replace distance with dmax when greater
  r    = array([d])
  indx = where(array([d]) > dmax)[0]
  if len(indx) > 0:
    r[indx] = dmax
  
  # heliocentric cartesian coordinates
  x = array([r*cos(b)*cos(l)])
  y = array([r*cos(b)*sin(l)])
  z = array([r*sin(b) + zsun])

  # for stars in Solar neighborhood
  i  = where(logical_and(abs(x) < 1.,abs(y) < 2.))[0]
  ni = len(i)
  j  = where(logical_or(abs(x) >= 1., abs(y) >= 2.))[0]
  nj = len(j)
  
  if ni > 0:
    # define the local grid
    dx = 0.02
    dy = 0.02
    dz = 0.02
    nx = 101
    ny = 201
    nz = 51

    # grid indices
    xi = x[i]/dx + float(nx - 1)/2.
    yj = y[i]/dy + float(ny - 1)/2.
    zk = z[i]/dz + float(nz - 1)/2.
    
    # interpolate
    avloc[i] = _trint(avori2, xi, yj, zk, missing=0.)
  
  # for stars in Solar neighborhood
  k = where(logical_and(abs(x) < 0.75, abs(y) < 0.75))[0]
  nk = len(k)
  m = where(logical_or(abs(x) >= 0.75, abs(y) >= 0.75))[0]
  nm = len(m)
  
  if nk > 0:
  
    # define the local grid
    dx = 0.05
    dy = 0.05
    dz = 0.02
    nx = 31
    ny = 31
    nz = 51

    # grid indices
    xi = x[k]/dx + float(nx - 1)/2.
    yj = y[k]/dy + float(ny - 1)/2.
    zk = z[k]/dz + float(nz - 1)/2.
    
    # trilinear_interpolate
    absdisk[k] = _trint(avdloc,xi,yj,zk,missing=0.)
  
  # galacto-centric cartesian
  x = x + xsun

  #larger orion arm grid 
  if nj > 0:
    # calculate the allowed maximum distance for larger orion grid
    dmax = 100.
    if b != 0.:
      dmax = .49999/abs(sin(b)) - zsun/sin(b)
    if cos(l) > 0.:
      dmax = min([dmax, 2.374999/abs(cos(l))])
    if cos(l) < 0.:
      dmax = min([dmax, 1.374999/abs(cos(l))])
    if sin(l) != 0.:
      dmax = min([dmax, 3.749999/abs(sin(l))])
      
    # replace distance with dmax when greater
    r1    = array([d])
    indx  = where(array([d]) >= dmax)[0]
    n     = len(indx)
    if n > 0:
      r1[indx] = dmax
    
    # galactocentric centric cartesian coordinates
    x1 = array([r1*cos(b)*cos(l) + xsun])
    y1 = array([r1*cos(b)*sin(l)])
    z1 = array([r1*sin(b) + zsun])
    
    # define the grid
    dx = 0.05
    dy = 0.05
    dz = 0.02
    nx = 76
    ny = 151
    nz = 51
    
    # grid indices
    xi = x1[j]/dx + 2.5*float(nx - 1)
    yj = y1[j]/dy + float(ny - 1)/2.
    zk = z1[j]/dz + float(nz - 1)/2.
    
    # trilinear_interpolate
    avloc[j] = _trint(avori,xi,yj,zk,missing=0.)
  
  # define the grid
  dx = 0.2
  dy = 0.2
  dz = 0.02
  nx = 151
  ny = 151
  nz = 51

  # grid indices
  xi = x/dx + float(nx - 1)/2.
  yj = y/dy + float(ny - 1)/2.
  zk = z/dz + float(nz - 1)/2.

  # trilinear_interpolate
  abspir = _trint(avspir,xi,yj,zk,missing=0.)
  if nm > 0:
    absdisk[m] = _trint(avdisk,xi[m],yj[m],zk[m],missing=0.)

  # apply rescaling factors or not
  if rescaling:
    out = (dfac[tblindex]*absdisk + sfac[tblindex]*abspir + lfac[tblindex]*avloc).flatten()
  else:
    out = (absdisk + abspir + avloc).flatten()
  
  return(out)

def _pix2xy(pixel, resolution, sixpack=0):
  """
  _pix2xy creates a raster image (sky cube or face) given a pixelindex and a
  resolution  The data array can be either a vector or two-dimensional array.
  In the latter case, the data for each raster image can be stored in either
  the columns or rows.  The procedure also returns the x and y raster
  coordinates of each pixel. Only right oriented, ecliptic coordinate images
  are built.
  SMOLDERS SEAL OF APPROVAL
  """
  # reform pixel to get rid of all length-1 dimensions
  pixel      = array(pixel)
  resolution = int(resolution)
  newshape   = []
  for s in pixel.shape:
    if s > 1:
      newshape.append()
  pixel.shape = newshape
  pix_sz = (pixel)
  
  if max(pixel) > 6*4**(resolution-1):
    raise('Maximum pixel number too large for resolution')
  
  # set up flag values for RASTR
  data = [-1]
  face = -1
  bad_pixval = 0.0
  # call rasterization routine
  xout, yout = _rastr(pixel,resolution)
  return(xout, yout)

def _rastr(pixel,resolution):
  """
  SMOLDERS SEAL OF APPROVAL
  """
  pixel     = array(pixel)
  n         = pixel.size
  i0        = 3
  j0        = 2
  offx      = [0,0,1,2,2,1]
  offy      = [1,0,0,0,1,1]
  fij       = _pix2fij(pixel,resolution)
  cube_side = 2**(resolution-1)
  lenc      = i0*cube_side
  if n > 1:
    x_out = offx[fij[0,:]] * cube_side + fij[1,:]
    x_out = lenc - (x_out+1)
    y_out = offy[fij[0,:]] * cube_side + fij[2,:]
  else:
    x_out = offx[fij[0]] * cube_side + fij[1]
    x_out = lenc - (x_out+1)
    y_out = offy[fij[0]] * cube_side + fij[2]
  return(x_out, y_out)

def _pix2fij(pixel,resolution):
  """
  This function takes an n-element pixel array and generates an n by 3 element
  array containing the corresponding face, column, and row number (the latter
  two within the face).
  SMOLDERS SEAL OF APPROVAL
  """
  # get number of pixels
  pixel        = array(pixel)
  n            = pixel.size
  output       = array(zeros((3,n)), int)
  res1         = resolution - 1
  num_pix_face = 4**res1
  #
  face = pixel/num_pix_face
  fpix = pixel-num_pix_face*face
  output[0,:] = face
  pow_2       = 2**arange(16)
  ii          = array(zeros(n), int)
  jj          = array(zeros(n), int)
  #
  if n > 1:
    for i in arange(n):
      for bit in arange(res1):
        ii[i] = ii[i] | (pow_2[bit]*(1 & fpix))
        fpix  = fpix >> 1
        jj[i] = jj[i] | (pow_2[bit]*(1 & fpix))
        fpix  = fpix >> 1
  else:
    for bit in arange(res1):
        ii    = ii | (pow_2[bit]*(1 & fpix))
        fpix  = fpix >> 1
        jj    = jj | (pow_2[bit]*(1 & fpix))
        fpix  = fpix >> 1
  output[1,:] = ii
  output[2,:] = jj
  return output

def _incube(alpha, beta):
  gstar =  1.37484847732
  g     = -0.13161671474
  m     =  0.004869491981
  w1    = -0.159596235474
  c00   =  0.141189631152
  c10   =  0.0809701286525
  c01   = -0.281528535557
  c11   =  0.15384112876
  c20   = -0.178251207466
  c02   =  0.106959469314
  d0    =  0.0759196200467
  d1    = -0.0217762490699
  r0    =  0.577350269
  aa    = alpha**2.
  bb    = beta**2
  a4    = aa**2
  b4    = bb**2
  onmaa = 1.-aa
  onmbb = 1.-bb
  x = alpha*(gstar + aa*(1.-gstar) + onmaa*(bb*(g+(m-g)*aa + onmbb*(c00+c10*aa+c01*bb+c11*aa*bb+c20*a4+c02*b4)) + aa*(w1-onmaa*(d0+d1*aa))))
  y = beta *(gstar + bb*(1.-gstar) + onmbb*(aa*(g+(m-g)*bb + onmaa*(c00+c10*bb+c01*aa+c11*bb*aa+c20*b4+c02*a4)) + bb*(w1-onmbb*(d0+d1*bb))))
  return x, y

def _ll2uv(lng, lat):
  """
  Convert longitude, latitude to unit vectors
  SMOLDERS SEAL OF APPROVAL
  
  @param input : galactic longitude
  @type  input : float 
  @param infmt : galactic lattitude
  @type  infmt : float
  @return      : unitvector
  @rtype       : ndarray
  """
  vector = zeros(3)
  d2r    = pi/180
  lng = lng * d2r
  lat = lat * d2r
  vector[0] = cos(lat) * cos(lng)
  vector[1] = cos(lat) * sin(lng)
  vector[2] = sin(lat)
  return vector

def _galvec2eclvec(in_uvec):
  """
  Unitvector for galactic coordinates to unitvector for ecliptic coordinates
  SMOLDERS SEAL OF APPROVAL
  """
  out_uvec = dot(in_uvec, g2e)
  return out_uvec

def _axisxy(vector):
  """
  converts unitvector into nface number (0-5) and X,Y in range 0-1
  SMOLDERS SEAL OF APPROVAL
  """
  vector = array(vector)
  vs   = vector.shape
  if len(vs) > 1:
    vec0 = array(vector[:,0])
    vec1 = array(vector[:,1])
    vec2 = array(vector[:,2])
  else:
    vec0 = array([vector[0]])
    vec1 = array([vector[1]])
    vec2 = array([vector[2]])
  abs_yx = abs(vec1/vec0)
  abs_zx = abs(vec2/vec0)
  abs_zy = abs(vec2/vec1)
  #
  nface = zeros(len(vec0))
  for i in arange(len(vec0)):
    nface[i] = (0 * (abs_zx[i] >= 1 and abs_zy[i] >= 1 and vec2[i] >= 0) +
                5 * (abs_zx[i] >= 1 and abs_zy[i] >= 1 and vec2[i] <  0) +
                1 * (abs_zx[i] <  1 and abs_yx[i] <  1 and vec0[i] >= 0) +
                3 * (abs_zx[i] <  1 and abs_yx[i] <  1 and vec0[i] <  0) +
                2 * (abs_zy[i] <  1 and abs_yx[i] >= 1 and vec1[i] >= 0) +
                4 * (abs_zy[i] <  1 and abs_yx[i] >= 1 and vec1[i] <  0))
  #
  nface_0 = (nface == 0)*1.
  nface_1 = (nface == 1)*1.
  nface_2 = (nface == 2)*1.
  nface_3 = (nface == 3)*1.
  nface_4 = (nface == 4)*1.
  nface_5 = (nface == 5)*1.
  #
  eta = (vec2/vec0)*(nface_1 - nface_3) + (vec2/vec1)*(nface_2 - nface_4) - (vec0/vec2)*(nface_0 + nface_5)
  xi  = (vec1/vec0)*(nface_1 + nface_3) - (vec0/vec1)*(nface_2 + nface_4) + (vec1/vec2) * (nface_0 - nface_5)
  x, y = _incube(xi,eta)
  x = (x+1.)/2.
  y = (y+1.)/2.
  return x,y,array(nface, dtype=int)

def _fij2pix(fij,res):
  """
  This function takes an n by 3 element vector containing the face, column, and
  row number (the latter two within the face) of a pixel and converts it into an
  n-element pixel array for a given resolution.
  """
  #get the number of pixels
  fij     = array(fij)
  n       = fij.size/3
  # generate output pixel and intermediate pixel arrays
  pixel   = array(zeros(n), dtype=int)
  pixel_1 = array(zeros(n), dtype=int)
  # get input column and row numbers
  if n > 1:
    ff = array(fij[0,:], dtype=int)
    ii = array(fij[1,:], dtype=int)
    jj = array(fij[2,:], dtype=int)
  else:
    ff = int(fij[0])
    ii = int(fij[1])
    jj = int(fij[2])
  # calculate the number of pixels in a face
  num_pix_face = 4**(res-1)
  pow_2        = 2**arange(16)
  # if col bit set then set corresponding even bit in pixel_l
  # if row bit set then set corresponding odd bit in pixel_l
  if n > 1:
    for i in arange(n):
      for bit in arange(res-1):
        pixel_1[i] = pixel_1[i] | ((pow_2[bit] & ii[i]) << bit)
        pixel_1[i] = pixel_1[i] | ((pow_2[bit] & jj[i]) << (bit+1))
  else:
    for bit in arange(res-1):
      pixel_1 = pixel_1 | ((pow_2[bit] & ii) << bit)
      pixel_1 = pixel_1 | ((pow_2[bit] & jj) << (bit+1))
  # add face number offset
  pixel = ff*num_pix_face + pixel_1
  return pixel

def _uv2pix(vector, resolution):
  """
  Routine returns pixel number given unit vector pointing to center of pixel
  resolution of the cube.
  SMOLDERS SEAL OF APPROVAL
  """
  two          = 2
  vector       = array(vector)
  n_vec        = int(vector.size/3)
  res1         = resolution - 1
  num_pix_side = int(two**res1)
  x, y, face   = _axisxy(vector)
  ia           = array(x*num_pix_side, dtype=int)
  ib           = array([2.**res1 - 1]*n_vec, dtype=int)
  ja           = array(y*num_pix_side, dtype=int)
  jb           = array([2.**res1 - 1]*n_vec, dtype=int)
  i            = zeros(n_vec)
  j            = zeros(n_vec)
  for k in arange(n_vec):
    i[k] = min([ia[k],ib[k]])
    j[k] = min([ja[k],jb[k]])
  pixel        = _fij2pix(array([face,i,j]),resolution)
  return pixel

def _ll2pix(lng, lat, res=9):
  """
  _ll2pix is a python function to convert galactic coordinates to DIRBE pixels
  (coorconv([lng,lat], infmt='L', outfmt='P', inco='G', outco='R9'))
  
  Input
  @param input     : galactic longitude
  @type  input     : float 
  @param infmt     : galactic lattitude
  @type  infmt     : float
  @return out_coor : output coordinate array
  @rtype  our_coor : ndarray
  
  (Note: The default coordinate system for pixels is ecliptic.)
    """
  # COMMON sky_conv_com, e2g,e2q,g2e,g2q,q2e,q2g,load_flag
  uvec = _ll2uv(lng, lat)
  uvec = _galvec2eclvec(uvec)
  output = _uv2pix(uvec,res)
  return output

def _trint(mm,x,y,z,missing=None):
  """
  Pixel-based trilinear interpolation, where mm is a 3d data cube.
  """
  xmi = array(floor(x),int)
  xd  = abs(xmi-x)
  xma = array(ceil(x),int)
  ymi = array(floor(y), int)
  yd  = abs(ymi-y)
  yma = array(ceil(y),int)
  zmi = array(floor(z), int)
  zd  = abs(zmi-z)
  zma = array(ceil(z), int)
  i1  = mm[zmi,ymi,xmi]*(1.-zd) + mm[zma,ymi,xmi]*(zd)
  i2  = mm[zmi,yma,xmi]*(1.-zd) + mm[zma,yma,xmi]*(zd)
  j1  = mm[zmi,ymi,xma]*(1.-zd) + mm[zma,ymi,xma]*(zd)
  j2  = mm[zmi,yma,xma]*(1.-zd) + mm[zma,yma,xma]*(zd)
  w1  = i1*(1-yd)  + i2*yd
  w2  = j1*(1-yd)  + j2*yd
  iv  = w1*(1.-xd) + w2*xd
  return iv

#}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#{ Schlegel 3D extinction model presented in Schlegel et al.
# The Astrophysical Journal, 500:525È553, 1998 June 20,
# "MAPS OF DUST INFRARED EMISSION FOR USE IN ESTIMATION
# OF REDDENING AND COSMIC MICROWAVE BACKGROUND RADIATION
# FOREGROUNDS"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@memoized
def get_schlegel_data_south():
  # Read in the Schlegel data of the southern hemisphere
  dustname = config.get_datafile('schlegel',"SFD_dust_4096_sgp.fits")
  maskname = config.get_datafile('schlegel',"SFD_mask_4096_sgp.fits")
  data     = pf.getdata(dustname)
  mask     = pf.getdata(maskname)
  return data, mask

def get_schlegel_data_north():
  # Read in the Schlegel data of the northern hemisphere
  dustname = config.get_datafile('schlegel',"SFD_dust_4096_ngp.fits")
  maskname = config.get_datafile('schlegel',"SFD_mask_4096_ngp.fits")
  data     = pf.getdata(dustname)
  mask     = pf.getdata(maskname)
  return data, mask

def _lb2xy_schlegel(ll, bb):
  """
  Converts coordinates in lattitude and longitude to coordinates
  in x and y pixel coordinates.
  
  The x and y coordinates you find with these formulas are not
  the coordinates you can read in ds9, but 1 smaller. Hence
  (x+1, y+1) are the coordinates you find in DS9.
  
  Input
  @param ll     : galactic longitude
  @type  ll     : float 
  @param bb     : galactic lattitude
  @type  bb     : float
  @return out_coor : output coordinate array
  @rtype: ndarray
  """
  deg2rad = pi/180. # convert degrees to rads

  if bb <= 0  : hs = -1.
  elif bb >  0: hs = +1.
  
  yy =  2048 * sqrt(1. - hs * sin(bb*deg2rad)) * cos(ll*deg2rad) + 2047.5
  xx = -2048 * hs * sqrt(1 - hs * sin(bb*deg2rad)) * sin(ll*deg2rad) + 2047.5
  
  return xx, yy

def findext_schlegel(ll, bb, distance = None, Rv=3.1,**kwargs):
  """
  Get the "Schlegel" extinction at certain longitude and latitude
  
  This function returns the E(B-V) maps of Schlegel et al. (1998),
  depending on wether the distance is given or not, the E(B-V) value
  is corrected. If distance is set we use the distance-corrected
  values:

    E(B-V) = E(B-V)_maps * (1 - exp(-10 * r * sin(|b|)))
    where E(B-V) is the value to be used and E(B-V)_maps the value
    as found with the Schlegel dust maps
  
  Then we convert the E(B-V) to Av. Standard we use Av = E(B-V)*Rv with Rv=3.1, but the value of Rv can be given as a keyword.

  ! WARNING: the schlegel maps are not usefull when |b| < 5 degrees !
  """
  deg2rad = pi/180. # convert degrees to rads
  dd = distance
  if distance is not None:
    dd      = distance/1.e3 # convert to kpc
  
  
  # first get the right pixel coordinates
  xx, yy = _lb2xy_schlegel(ll,bb)
  
  # read in the right map:
  if bb <= 0:
    data, mask = get_schlegel_data_south()
  elif bb >  0:
    data, mask = get_schlegel_data_north()

  if abs(bb) < 10.:
    logger.warning("Schlegel is not good for lattitudes > 10 degrees")
    
  # the xy-coordinates are:
  xl = floor(xx)
  yl = floor(yy)
  xh = xl + 1.
  yh = yl + 1.
  
  # the weights are just the distances to the points
  w1 = (xl-xx)**2 + (yl-yy)**2
  w2 = (xl-xx)**2 + (yh-yy)**2
  w3 = (xh-xx)**2 + (yl-yy)**2
  w4 = (xh-xx)**2 + (yh-yy)**2
  
  # the values of these points are:
  v1 = data[xl, yl]
  v2 = data[xl, yh]
  v3 = data[xh, yl]
  v4 = data[xh, yh]  
  f1 = mask[xl, yl]
  f2 = mask[xl, yh]
  f3 = mask[xh, yl]
  f4 = mask[xh, yh]
  
  ebv = (w1*v1 + w2*v2 + w3*v3 + w4*v4) / (w1 + w2 + w3 + w4)
  
  # Check flags at the right pixels
  logger.info("flag of pixel 1 is: %i" %f1)
  logger.info("flag of pixel 2 is: %i" %f2)
  logger.info("flag of pixel 3 is: %i" %f3)
  logger.info("flag of pixel 4 is: %i" %f4)
  
  if dd is not None:
    ebv = ebv * (1. - np.exp(-10. * dd * sin(abs(bb*deg2rad))))
  
  # if Rv is given, by definition we find Av = Ebv*Rv
  av = ebv*Rv
  
  return av
#}