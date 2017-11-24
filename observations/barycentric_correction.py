# -*- coding: utf-8 -*-
"""
Compute the influence of earth's orbit around the sun, and absolute stellar velocities
"""
#from numpy import *
import numpy as np
from ivs.units.conversions import convert

def baryvel(dje, deq=0):
  """
  Calculates heliocentric and barycentric velocity components of Earth.

  This function takes into account the Earth-Moon motion, and is useful for
  radial velocity work to an accuracy of  ~1 m/s.
  
  dvel_hel, dvel_bary = baryvel(dje, deq)
   
  Inputs:
  @param dje: Julian ephemeris date.
  @type dje: array
  @param deq: epoch of mean equinox of dvelh and dvelb. If deq = 0 then deq is assumed to be equal to dje
  @type deq: float
  @return: heliocentric and barycentric velocity components in km/s (The 3-vectors dvelh and dvelb are given in a right-handed coordinate system with the +X axis toward the Vernal Equinox, and +Z axis toward the celestial pole.)
  @rtype: array (2X3)          
  
  Functions called:
    premat() -- computes precession matrix
    
  NOTES:
    Algorithm taken from FORTRAN program of Stumpff (1980, A&A Suppl, 41,1)
    Stumpf claimed an accuracy of 42 cm/s for the velocity. A comparison with
    the JPL FORTRAN planetary ephemeris program PLEPH found agreement to within
    about 65 cm/s between 1986 and 1994
   
  REVISION HISTORY:
    Jeff Valenti,  U.C. Berkeley    Translated BARVEL.FOR to IDL.
    W. Landsman, Cleaned up program sent by Chris McCarthy (SfSU) June 1994
    Converted to IDL V5.0   W. Landsman   September 1997
    Added /JPL keyword  W. Landsman   July 2001
    Documentation update W. Landsman Dec 2005
    Converted to Python S. Koposov 2009-2010
    Adapted for IvS library
  """
  
  #Define constants
  dc2pi = 2 * np.pi
  cc2pi = 2 * np.pi
  dc1 = 1.0e0
  dcto = 2415020.0e0
  dcjul = 36525.0e0                     #days in Julian year
  dcbes = 0.313e0
  dctrop = 365.24219572e0               #days in tropical year (...572 insig)
  dc1900 = 1900.0e0
  au = 1.4959787e8
  
  #Constants dcfel(i,k) of fast changing elements.
  dcfel = np.array([1.7400353e00, 6.2833195099091e02, 5.2796e-6, 6.2565836e00, 6.2830194572674e02, -2.6180e-6, 4.7199666e00, 8.3997091449254e03, -1.9780e-5, 1.9636505e-1, 8.4334662911720e03, -5.6044e-5, 4.1547339e00, 5.2993466764997e01, 5.8845e-6, 4.6524223e00, 2.1354275911213e01, 5.6797e-6, 4.2620486e00, 7.5025342197656e00, 5.5317e-6, 1.4740694e00, 3.8377331909193e00, 5.6093e-6])
  dcfel = np.reshape(dcfel, (8, 3))
  
  #constants dceps and ccsel(i,k) of slowly changing elements.
  dceps = np.array([4.093198e-1, -2.271110e-4, -2.860401e-8])
  ccsel = np.array([1.675104e-2, -4.179579e-5, -1.260516e-7, 2.220221e-1, 2.809917e-2, 1.852532e-5, 1.589963e00, 3.418075e-2, 1.430200e-5, 2.994089e00, 2.590824e-2, 4.155840e-6, 8.155457e-1, 2.486352e-2, 6.836840e-6, 1.735614e00, 1.763719e-2, 6.370440e-6, 1.968564e00, 1.524020e-2, -2.517152e-6, 1.282417e00, 8.703393e-3, 2.289292e-5, 2.280820e00, 1.918010e-2, 4.484520e-6, 4.833473e-2, 1.641773e-4, -4.654200e-7, 5.589232e-2, -3.455092e-4, -7.388560e-7, 4.634443e-2, -2.658234e-5, 7.757000e-8, 8.997041e-3, 6.329728e-6, -1.939256e-9, 2.284178e-2, -9.941590e-5, 6.787400e-8, 4.350267e-2, -6.839749e-5, -2.714956e-7, 1.348204e-2, 1.091504e-5, 6.903760e-7, 3.106570e-2, -1.665665e-4, -1.590188e-7])
  ccsel = np.reshape(ccsel, (17, 3))
  
  #Constants of the arguments of the short-period perturbations.
  dcargs = np.array([5.0974222e0, -7.8604195454652e2, 3.9584962e0, -5.7533848094674e2, 1.6338070e0, -1.1506769618935e3, 2.5487111e0, -3.9302097727326e2, 4.9255514e0, -5.8849265665348e2, 1.3363463e0, -5.5076098609303e2, 1.6072053e0, -5.2237501616674e2, 1.3629480e0, -1.1790629318198e3, 5.5657014e0, -1.0977134971135e3, 5.0708205e0, -1.5774000881978e2, 3.9318944e0, 5.2963464780000e1, 4.8989497e0, 3.9809289073258e1, 1.3097446e0, 7.7540959633708e1, 3.5147141e0, 7.9618578146517e1, 3.5413158e0, -5.4868336758022e2])
  dcargs = np.reshape(dcargs, (15, 2))
  
  #Amplitudes ccamps(n,k) of the short-period perturbations.
  ccamps = np.array([-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5, -2.490817e-7, -3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5, -1.823138e-7, 6.593466e-7, 1.322572e-5, 9.258695e-6, -4.674248e-7, -3.646275e-7, 1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7, 9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7, 7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7, -2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6, -1.655307e-7, -3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6, -3.736225e-7, 3.442177e-7, 2.671323e-6, 1.832858e-6, -2.394688e-7, -3.478444e-7, 8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8, -1.488378e-6, -1.251789e-5, 5.226868e-7, -2.049301e-7, 0.e0, -8.043059e-6, -2.991300e-6, 1.473654e-7, -3.154542e-7, 0.e0, 3.699128e-6, -3.316126e-6, 2.901257e-7, 3.407826e-7, 0.e0, 2.550120e-6, -1.241123e-6, 9.901116e-8, 2.210482e-7, 0.e0, -6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.e0])
  ccamps = np.reshape(ccamps, (15, 5))
  
  #Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
  ccsec3 = -7.757020e-8
  ccsec = np.array([1.289600e-6, 5.550147e-1, 2.076942e00, 3.102810e-5, 4.035027e00, 3.525565e-1, 9.124190e-6, 9.990265e-1, 2.622706e00, 9.793240e-7, 5.508259e00, 1.559103e01])
  ccsec = np.reshape(ccsec, (4, 3))
  
  #Sidereal rates.
  dcsld = 1.990987e-7 #sidereal rate in longitude
  ccsgd = 1.990969e-7 #sidereal rate in mean anomaly
  
  #Constants used in the calculation of the lunar contribution.
  cckm = 3.122140e-5
  ccmld = 2.661699e-6
  ccfdi = 2.399485e-7
  
  #Constants dcargm(i,k) of the arguments of the perturbations of the motion
  # of the moon.
  dcargm = np.array([5.1679830e0, 8.3286911095275e3, 5.4913150e0, -7.2140632838100e3, 5.9598530e0, 1.5542754389685e4])
  dcargm = np.reshape(dcargm, (3, 2))
  
  #Amplitudes ccampm(n,k) of the perturbations of the moon.
  ccampm = np.array([1.097594e-1, 2.896773e-7, 5.450474e-2, 1.438491e-7, -2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8, 1.148966e-2, 5.658888e-8, 8.249439e-3, 4.063015e-8])
  ccampm = np.reshape(ccampm, (3, 4))
  
  #ccpamv(k)=a*m*dl,dt (planets), dc1mme=1-mass(earth+moon)
  ccpamv = np.array([8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12])
  dc1mme = 0.99999696e0
  
  #Time arguments.
  dt = (dje - dcto) / dcjul
  tvec = np.array([1e0, dt, dt * dt])
  
  #Values of all elements for the instant dje.
  temp = (np.transpose(np.dot(np.transpose(tvec), np.transpose(dcfel)))) % dc2pi
  dml = temp[0]
  forbel = temp[1:8]
  g = forbel[0]      #old fortran equivalence
  
  deps = (tvec * dceps).sum() % dc2pi
  sorbel = (np.transpose(np.dot(np.transpose(tvec), np.transpose(ccsel)))) % dc2pi
  e = sorbel[0]     #old fortran equivalence
  
  #Secular perturbations in longitude.
  dummy = np.cos(2.0)
  sn = np.sin((np.transpose(np.dot(np.transpose(tvec[0:2]), np.transpose(ccsec[:,1:3])))) % cc2pi)
  
  #Periodic perturbations of the emb (earth-moon barycenter).
  pertl = (ccsec[:,0] * sn).sum() + dt * ccsec3 * sn[2]
  pertld = 0.0
  pertr = 0.0
  pertrd = 0.0
  for k in range(0, 15):
    a = (dcargs[k,0] + dt * dcargs[k,1]) % dc2pi
    cosa = np.cos(a)
    sina = np.sin(a)
    pertl = pertl + ccamps[k,0] * cosa + ccamps[k,1] * sina
    pertr = pertr + ccamps[k,2] * cosa + ccamps[k,3] * sina
    if k < 11:   
      pertld = pertld + (ccamps[k,1] * cosa - ccamps[k,0] * sina) * ccamps[k,4]
      pertrd = pertrd + (ccamps[k,3] * cosa - ccamps[k,2] * sina) * ccamps[k,4]
      
  #Elliptic part of the motion of the emb.
  phi = (e * e / 4e0) * (((8e0 / e) - e) * np.sin(g) + 5 * np.sin(2 * g) + (13 / 3e0) * e * np.sin(3 * g))
  f = g + phi
  sinf = np.sin(f)
  cosf = np.cos(f)
  dpsi = (dc1 - e * e) / (dc1 + e * cosf)
  phid = 2 * e * ccsgd * ((1 + 1.5 * e * e) * cosf + e * (1.25 - 0.5 * sinf * sinf))
  psid = ccsgd * e * sinf / np.sqrt(dc1 - e * e)
  
  #Perturbed heliocentric motion of the emb.
  d1pdro = dc1 + pertr
  drd = d1pdro * (psid + dpsi * pertrd)
  drld = d1pdro * dpsi * (dcsld + phid + pertld)
  dtl = (dml + phi + pertl) % dc2pi
  dsinls = np.sin(dtl)
  dcosls = np.cos(dtl)
  dxhd = drd * dcosls - drld * dsinls
  dyhd = drd * dsinls + drld * dcosls
  
  #Influence of eccentricity, evection and variation on the geocentric
  # motion of the moon.
  pertl = 0.0
  pertld = 0.0
  pertp = 0.0
  pertpd = 0.0
  for k in range(0, 3):
    a = (dcargm[k,0] + dt * dcargm[k,1]) % dc2pi
    sina = np.sin(a)
    cosa = np.cos(a)
    pertl = pertl + ccampm[k,0] * sina
    pertld = pertld + ccampm[k,1] * cosa
    pertp = pertp + ccampm[k,2] * cosa
    pertpd = pertpd - ccampm[k,3] * sina
    
  #Heliocentric motion of the earth.
  tl = forbel[1] + pertl
  sinlm = np.sin(tl)
  coslm = np.cos(tl)
  sigma = cckm / (1.0 + pertp)
  a = sigma * (ccmld + pertld)
  b = sigma * pertpd
  dxhd = dxhd + a * sinlm + b * coslm
  dyhd = dyhd - a * coslm + b * sinlm
  dzhd = -sigma * ccfdi * np.cos(forbel[2])
  
  #Barycentric motion of the earth.
  dxbd = dxhd * dc1mme
  dybd = dyhd * dc1mme
  dzbd = dzhd * dc1mme
  for k in range(0, 4):
    plon = forbel[k + 3]
    pomg = sorbel[k + 1]
    pecc = sorbel[k + 9]
    tl = (plon + 2.0 * pecc * np.sin(plon - pomg)) % cc2pi
    dxbd = dxbd + ccpamv[k] * (np.sin(tl) + pecc * np.sin(pomg))
    dybd = dybd - ccpamv[k] * (np.cos(tl) + pecc * np.cos(pomg))
    dzbd = dzbd - ccpamv[k] * sorbel[k + 13] * np.cos(plon - sorbel[k + 5])
    
  #Transition to mean equator of date.
  dcosep = np.cos(deps)
  dsinep = np.sin(deps)
  dyahd = dcosep * dyhd - dsinep * dzhd
  dzahd = dsinep * dyhd + dcosep * dzhd
  dyabd = dcosep * dybd - dsinep * dzbd
  dzabd = dsinep * dybd + dcosep * dzbd
  
  #Epoch of mean equinox (deq) of zero implies that we should use
  # Julian ephemeris date (dje) as epoch of mean equinox.
  if deq == 0:   
    dvelh = au * (np.array([dxhd, dyahd, dzahd]))
    dvelb = au * (np.array([dxbd, dyabd, dzabd]))
    return (dvelh,dvelb)
  
  #General precession from epoch dje to deq.
  deqdat = (dje - dcto - dcbes) / dctrop + dc1900
  prema = premat(deqdat, deq, fk4=True)
  dvelh = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxhd, dyahd, dzahd])))))
  dvelb = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxbd, dyabd, dzabd])))))
  
  return (dvelh, dvelb)

def bprecess(ra0, dec0, mu_radec=None, parallax=None, rad_vel=None, epoch=None):
  """
  Precess positions from J2000.0 (FK5) to B1950.0 (FK4)
  
  Calculates the mean place of a star at B1950.0 on the FK4 system from the mean
  place at J2000.0 on the FK5 system.
  
  Input:
    @param ra0: J2000 right ascension in degrees
    @type ra0: N array
    @param dec0: J2000 declination in degrees
    @type dec0: N array
    @keyword mu_radec: array containing the proper motion in seconds of arc per tropical *century* in right ascension and declination.
    @type mu_radec: 2xN array
    @keyword parallax: array giving stellar parallax (seconds of arc)
    @type parallax: N array
    @keyword rad_vel: array giving radial velocity in km/s
    @type rad_vel: N array
    @keyword epoch: scalar giving epoch of original observations, default 2000. This keyword value is only used if mu_radec is not given
    @type epoch:
    @return: The corresponding B1950 right ascension and declination in *degrees*.
    @rtype: N array, N array
    
  NOTES:
    The algorithm is taken from the Explanatory Supplement to the Astronomical Almanac 1992, page 186. Also see Aoki et al (1983), A&A, 128,263
    
    BPRECESS distinguishes between the following two cases:
    (1) The proper motion is known and non-zero
    (2) the proper motion is unknown or known to be exactly zero (i.e.
        extragalactic radio sources).   In this case, the reverse of the
        algorithm in Appendix 2 of Aoki et al. (1983) is used to ensure that the
        output proper motion is exactly zero. Better precision can be achieved
        in this case by inputting the epoch of the original observations.
        
  REVISION HISTORY:
    Written (W. Landsman, october 1992)
    Vectorized (W. Landsman, february, 1994)
    Treat case where proper motion not known or exactly zero  (november 1994)
    Handling of arrays larger than 32767   (Lars L. Christensen, march 1995)
    Converted to IDL V5.0   (W. Landsman, september 1997)
    Fixed bug where A term not initialized for vector input (W. Landsman, february 2000)
    Converted to python (Sergey Koposov, july 2010)
  """
  scal = True
  if isinstance(ra0, np.ndarray):
    ra   = ra0
    dec  = dec0
    n    = ra.size
    scal = False
  else:
    n   = 1
    ra  = np.array([ra0])
    dec = np.array([dec0])
    
  if rad_vel is None:   
    rad_vel = np.zeros(n)
  else:
    if not isinstance(rad_vel, np.ndarray):
      rad_vel = np.array([rad_vel],dtype=float)
    if rad_vel.size != n:   
      raise Exception('ERROR - RAD_VEL keyword vector must be of the same length as RA and DEC')
      
  if (mu_radec is not None):   
    if (np.array(mu_radec).size != 2 * n):   
      raise Exception('ERROR - MU_RADEC keyword (proper motion) be dimensioned (2,' + strtrim(n, 2) + ')')
    mu_radec = mu_radec * 1.
   
  if parallax is None:   
    parallax = np.zeros(n)
  else:   
    if not isinstance(parallax, np.ndarray):
      parallax = np.array([parallax],dtype=float)
      
  if epoch is None:   
    epoch = 2000.0e0
    
  radeg = 180.e0 / np.pi
  sec_to_radian = lambda x : np.deg2rad(x/3600.)
  
  m = np.array([np.array([+0.9999256795e0, -0.0111814828e0, -0.0048590040e0, -0.000551e0, -0.238560e0, +0.435730e0]),
  np.array([+0.0111814828e0, +0.9999374849e0, -0.0000271557e0, +0.238509e0, -0.002667e0, -0.008541e0]),
  np.array([+0.0048590039e0, -0.0000271771e0, +0.9999881946e0, -0.435614e0, +0.012254e0, +0.002117e0]),
  np.array([-0.00000242389840e0, +0.00000002710544e0, +0.00000001177742e0, +0.99990432e0, -0.01118145e0, -0.00485852e0]),
  np.array([-0.00000002710544e0, -0.00000242392702e0, +0.00000000006585e0, +0.01118145e0, +0.99991613e0, -0.00002716e0]),
  np.array([-0.00000001177742e0, +0.00000000006585e0, -0.00000242404995e0, +0.00485852e0, -0.00002717e0, +0.99996684e0])])
  
  a_dot = 1e-3 * np.array([1.244e0, -1.579e0, -0.660e0])           #in arc seconds per century
  
  ra_rad = np.deg2rad(ra)
  dec_rad = np.deg2rad(dec)
  cosra = np.cos(ra_rad)
  sinra = np.sin(ra_rad)
  cosdec = np.cos(dec_rad)
  sindec = np.sin(dec_rad)
  dec_1950 = dec * 0.
  ra_1950 = ra * 0.
  
  for i in range(n):
    # Following statement moved inside loop in Feb 2000.
    a = 1e-6 * np.array([-1.62557e0, -0.31919e0, -0.13843e0])        #in radians
    r0 = np.array([cosra[i] * cosdec[i], sinra[i] * cosdec[i], sindec[i]])
    if (mu_radec is not None):   
      mu_a = mu_radec[i,0]
      mu_d = mu_radec[i,1]
      r0_dot = np.array([-mu_a * sinra[i] * cosdec[i] - mu_d * cosra[i] * sindec[i], mu_a * cosra[i] * cosdec[i] - mu_d * sinra[i] * sindec[i], mu_d * cosdec[i]]) + 21.095e0 * rad_vel[i] * parallax[i] * r0
    else:   
      r0_dot = np.array([0.0e0, 0.0e0, 0.0e0])
      
    r_0 = np.concatenate((r0, r0_dot))
    r_1 = np.transpose(np.dot(np.transpose(m), np.transpose(r_0)))
    
    # Include the effects of the E-terms of aberration to form r and r_dot.
    r1     = r_1[0:3]
    r1_dot = r_1[3:6]
    if mu_radec is None:
      r1 = r1 + sec_to_radian ( r1_dot * (epoch - 1950.0e0) / 100. )
      a  = a + sec_to_radian ( a_dot * (epoch - 1950.0e0) / 100. )

    x1     = r_1[0]
    y1     = r_1[1]
    z1     = r_1[2]
    rmag   = np.sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2)
    s1     = r1 / rmag
    s1_dot = r1_dot / rmag
    s      = s1
    for j in np.arange(0, 3):
      r = s1 + a - ((s * a).sum()) * s
      s = r / rmag

    x    = r[0]
    y    = r[1]
    z    = r[2]
    r2   = x ** 2 + y ** 2 + z ** 2
    rmag = np.sqrt(r2)
    
    if mu_radec is not None:   
      r_dot         = s1_dot + a_dot - ((s * a_dot).sum()) * s
      x_dot         = r_dot[0]
      y_dot         = r_dot[1]
      z_dot         = r_dot[2]
      mu_radec[i,0] = (x * y_dot - y * x_dot) / (x ** 2 + y ** 2)
      mu_radec[i,1] = (z_dot * (x ** 2 + y ** 2) - z * (x * x_dot + y * y_dot)) / (r2 * np.sqrt(x ** 2 + y ** 2))
    
    dec_1950[i] = np.arcsin(z / rmag)
    ra_1950[i]  = np.arctan2(y, x)
    
    if parallax[i] > 0.:   
      rad_vel[i]  = (x * x_dot + y * y_dot + z * z_dot) / (21.095 * parallax[i] * rmag)
      parallax[i] = parallax[i] / rmag
  
  neg = (ra_1950 < 0)
  if neg.any() > 0:   
    ra_1950[neg] = ra_1950[neg] + 2.e0 * np.pi
  
  ra_1950  = np.rad2deg(ra_1950)
  dec_1950 = np.rad2deg(dec_1950)
  
  # Make output scalar if input was scalar
  if scal:
     return ra_1950[0],dec_1950[0]
  else:
     return ra_1950, dec_1950



## COLE START HERE
def convolve(image, psf, ft_psf=None, ft_image=None, no_ft=None, correlate=None, auto_correlation=None):
  """
  Convolution of an image with a Point Spread Function (PSF)
  The default is to compute the convolution using a product of Fourier transforms (for speed).
  
  Inputs:
    @param image: 2-D array (matrix) to be convolved with psf
    @type image: array
    @param psf: the Point Spread Function, (size < or = to size of image).
    @type psf array
    @keyword ft_psf: passes out/in the Fourier transform of the PSF, (so that it can be re-used the next time function is called).
    @type ft_psc: bool
    @param ft_image: passes out/in the Fourier transform of image.
    @type ft_image: bool
    @param correlate: uses the conjugate of the Fourier transform of PSF to compute the cross-correlation of image and PSF (equivalent to IDL function convol() with NO rotation of PSF)
    @type correlate: bool
    @param auto_corr: computes the auto-correlation function of image using FFT.
    @type auto_corr: bool
  
  METHOD:
    When using FFT, PSF is centered & expanded to size of image.
  
  HISTORY:
    Written (Frank Varosi, 1992)
    Appropriate precision type for result depending on input image (Markus Hundertmark, february 2006)
    Fix the bug causing the recomputation of FFT(psf) and/or FFT(image) (Sergey Koposov, december 2006)
  """
  from numpy.fft import fft2, ifft2
  n_params = 2
  psf_ft   = ft_psf
  imft     = ft_image
  noft     = no_ft
  auto     = auto_correlation
  sp       = np.array(np.shape(psf_ft)) 
  sif      = np.array(np.shape(imft))
  sim      = np.array(np.shape(image))
  sc       = sim / 2
  npix     = np.array(image, copy=0).size
  
  if image.ndim!=2 or noft!=None:   
    if (auto is not None):   
      message("auto-correlation only for images with FFT", inf=True)
      return image
    else:   
      if (correlate is not None):   
        return convol(image, psf)
      else:
        return convol(image, rotate(psf, 2))
        
  if imft==None or (imft.ndim!=2) or imft.shape!=im.shape: #add the type check
    imft = ifft2(image)
  if (auto is not None):   
    return np.roll(np.roll(npix * real(fft2(imft * conjugate(imft))), sc[0], 0),sc[1],1)
    
  if (ft_psf==None or ft_psf.ndim!=2 or ft_psf.shape!=image.shape or ft_psf.dtype!=image.dtype):
    sp     = np.array(shape(psf))
    loc    = np.maximum((sc - sp / 2), 0)         #center PSF in new array,
    s      = np.maximum((sp / 2 - sc), 0)         #handle all cases: smaller or bigger
    l      = np.minimum((s + sim - 1), (sp - 1))
    psf_ft = np.conjugate(image) * 0 #initialise with correct size+type according to logic of conj and set values to 0 (type of ft_psf is conserved)
    psf_ft[loc[1]:loc[1]+l[1]-s[1]+1,loc[0]:loc[0]+l[0]-s[0]+1] = psf[s[1]:(l[1])+1,s[0]:(l[0])+1]
    psf_ft = ifft2(psf_ft)
    
  if (correlate is not None):   
    conv = npix * np.real(fft2(imft * np.conjugate(psf_ft)))
  else:   
    conv = npix * np.real(fft2(imft * psf_ft))
    
  sc = sc + (sim % 2)   #shift correction for odd size images.
  return np.roll(np.roll(conv, sc[0],0), sc[1],1)

def cv_coord(a,b,c,fr=None,to=None,degr=False):
  #import numpy
  if degr:
    degrad = np.deg2rad
    raddeg = np.rad2deg
  else:
    degrad = lambda x: x
    raddeg = lambda x: x
  if fr=='sph':
    cosa = np.cos(degrad(a))
    sina = np.sin(degrad(a))
    cosb = np.cos(degrad(b))
    sinb = np.sin(degrad(b))
    x=c*cosa*cosb
    y=c*sina*cosb
    z=c*sinb
  elif fr=='rect':
    x=a
    y=b
    z=c
  elif fr is None:
    raise Exception('You must specify the input coordinate system')
  else:
    raise Exception('Unknown input coordinate system')
  if to=='rect':
    return (x,y,z)
  elif to=='sph':
    ra = raddeg(np.arctan2(y,x))
    dec = raddeg(np.arctan2(z,np.sqrt(x**2+y**2)))
    rad = np.sqrt(x**2+y**2+z**2)
    return (ra,dec,rad)
  elif to is None:
    raise Exception('You must specify the output coordinate system')
  else:
    raise Exception('Unknown output coordinate system')

def daycnv(xjd):
  """
  Converts Julian dates to Gregorian calendar dates
  
  @param xjd: Julian date, positive double precision scalar or vector
  @type: array
  
 OUTPUTS:
   YR: Year (Integer)
   MN: Month (Integer)
   DAY: Day (Integer)
   HR: Hours and fractional hours (Real). If XJD is a vector, then YR,MN,DAY and HR will be vectors of the same length.
   
 REVISION HISTORY:
   Converted to IDL from Yeoman's Comet Ephemeris Generator, B. Pfarr, STX, 6/16/88
   Converted to IDL V5.0   W. Landsman   September 1997
  """
  # Adjustment needed because Julian day starts at noon, calendar day at midnight
  
  jd = np.array(xjd).astype(int)                         #Truncate to integral day
  frac = np.array(xjd).astype(float) - jd + 0.5          #Fractional part of calendar day
  after_noon = (frac >= 1.0)
  
  if after_noon.any(): #Is it really the next calendar day?
    if frac.ndim>0: # proper array
      frac[after_noon] = frac[after_noon] - 1.0
      jd[after_noon] = jd[after_noon] + 1
    else:  # scalar
      frac = frac - 1.0
      jd = jd + 1
  hr  = frac * 24.0
  l   = jd + 68569
  n   = 4 * l / 146097
  l   = l - (146097 * n + 3) / 4
  yr  = 4000 * (l + 1) / 1461001
  l   = l - 1461 * yr / 4 + 31        #1461 = 365.25 * 4
  mn  = 80 * l / 2447
  day = l - 2447 * mn / 80
  l   = mn / 11
  mn  = mn + 2 - 12 * l
  yr  = 100. * (n - 49) + yr + l
  return (yr, mn, day, hr)

def euler(ai, bi, select=1, fk4=False):
  """
  Transform between Galactic, celestial, and ecliptic coordinates.
  
  Input:
    @param ai: Input longitude in degrees
    @type ai: float
    @param bi: Input latitude in degrees
    @type bi: float
    @keyword select: Integer (1-6) specifying type of coordinate transformation.
    @type select: int 1-6
    @keyword fk4: set the fk4, otherwise use J2000
    @type fk4: bool
    
  SELECT   From          To        |   SELECT      From            To
   1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
   2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
   3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
  """
  
  twopi = 2.0e0 * np.pi
  fourpi = 4.0e0 * np.pi
  
  if fk4:   
    equinox = '(B1950)'
    psi     = np.array ([0.57595865315e0, 4.9261918136e0, 0.00000000000e0, 0.0000000000e0, 0.11129056012e0, 4.7005372834e0])
    stheta  = np.array ([0.88781538514e0, -0.88781538514e0, 0.39788119938e0, -0.39788119938e0, 0.86766174755e0, -0.86766174755e0])
    ctheta  = np.array([0.46019978478e0, 0.46019978478e0, 0.91743694670e0, 0.91743694670e0, 0.49715499774e0, 0.49715499774e0])
    phi     = np.array([4.9261918136e0, 0.57595865315e0, 0.0000000000e0, 0.00000000000e0, 4.7005372834e0, 0.11129056012e0])
  else:   
    equinox = '(J2000)'
    psi     = np.array([0.57477043300e0, 4.9368292465e0, 0.00000000000e0, 0.0000000000e0, 0.11142137093e0, 4.71279419371e0])
    stheta  = np.array([0.88998808748e0, -0.88998808748e0, 0.39777715593e0, -0.39777715593e0, 0.86766622025e0, -0.86766622025e0])
    ctheta  = np.array([0.45598377618e0, 0.45598377618e0, 0.91748206207e0, 0.91748206207e0, 0.49714719172e0, 0.49714719172e0])
    phi     = np.array([4.9368292465e0, 0.57477043300e0, 0.0000000000e0, 0.00000000000e0, 4.71279419371e0, 0.11142137093e0])
    
  i    = select - 1
  a    = np.deg2rad(ai) - phi[i]
  b    = np.deg2rad(bi)
  sb   = np.sin(b)
  cb   = np.cos(b)
  cbsa = cb * np.sin(a)
  b    = -stheta[i] * cbsa + ctheta[i] * sb
  bo   = np.rad2deg(arcsin(minimum(b, 1.0)))
  del b
  a    = np.arctan2(ctheta[i] * cbsa + stheta[i] * sb, cb * np.cos(a))
  del cb, cbsa, sb
  ao   = np.rad2deg(((a + psi[i] + fourpi) % twopi) )
  return (ao, bo)

def gal_uvw(distance=None, lsr=None, ra=None, dec=None, pmra=None, pmdec=None, vrad=None, plx=None):
  """
  Calculate the Galactic space velocity (U,V,W) of star
  
  Calculates the Galactic space velocity U, V, W of star given its
  (1) coordinates, (2) proper motion, (3) distance (or parallax), and (4) radial velocity.

  OUTPUT PARAMETERS:
    U - Velocity (km/s) positive toward the Galactic *anti*center
    V - Velocity (km/s) positive in the direction of Galactic rotation
    W - Velocity (km/s) positive toward the North Galactic Pole
  INPUT KEYWORDS:
    User must supply a position, proper motion,radial velocity and distance
    (or parallax).    Either scalars or vectors can be supplied.
   (1) Position:
    RA - Right Ascension in *Degrees*
    Dec - Declination in *Degrees*
   (2) Proper Motion
    PMRA = Proper motion in RA in arc units (typically milli-arcseconds/yr)
    PMDEC = Proper motion in Declination (typically mas/yr)
   (3) Radial Velocity
    VRAD = radial velocity in km/s
   (4) Distance or Parallax
    DISTANCE - distance in parsecs
               or
    PLX - parallax with same distance units as proper motion measurements
          typically milliarcseconds (mas)
  
  REVISION HISTORY:
    Written, W. Landsman                       December   2000
    fix the bug occuring if the input arrays are longer than 32767
      and update the Sun velocity           Sergey Koposov June 2008
      vectorization of the loop -- performance on large arrays
      is now 10 times higher                Sergey Koposov December 2008
  """

  n_params = 3
  
  if ra is None or dec is None:   
     raise Exception('ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)')
  if plx is None and distance is None:
     raise Exception('ERROR - Either a parallax or distance must be specified')
  if distance is not None:
     if np.any(distance == 0):
        raise Exception('ERROR - All distances must be > 0')
     plx = 1e3 / distance
  if plx is not None and any(plx==0):   
     raise Exception('ERROR - Parallaxes must be > 0')
  
  cosd = np.cos(np.deg2rad(dec))
  sind = np.sin(np.deg2rad(dec))
  cosa = np.cos(np.deg2rad(ra))
  sina = np.sin(np.deg2rad(ra))
  k    = 4.74047
  a_g  = np.array([[0.0548755604, +0.4941094279, -0.8676661490], [0.8734370902, -0.4448296300, -0.1980763734], [0.4838350155, 0.7469822445, +0.4559837762]])  
  vec1 = vrad
  vec2 = k * pmra / plx
  vec3 = k * pmdec / plx
  u    = (a_g[0,0] * cosa * cosd + a_g[1,0] * sina * cosd + a_g[2,0] * sind) * vec1 + (-a_g[0,0] * sina + a_g[1,0] * cosa) * vec2 + (-a_g[0,0] * cosa * sind - a_g[1,0] * sina * sind + a_g[2,0] * cosd) * vec3
  v    = (a_g[0,1] * cosa * cosd + a_g[1,1] * sina * cosd + a_g[2,1] * sind) * vec1 + (-a_g[0,1] * sina + a_g[1,1] * cosa) * vec2 + (-a_g[0,1] * cosa * sind - a_g[1,1] * sina * sind + a_g[2,1] * cosd) * vec3
  w    = (a_g[0,2] * cosa * cosd + a_g[1,2] * sina * cosd + a_g[2,2] * sind) * vec1 + (-a_g[0,2] * sina + a_g[1,2] * cosa) * vec2 + (-a_g[0,2] * cosa * sind - a_g[1,2] * sina * sind + a_g[2,2] * cosd) * vec3
  lsr_vel = np.array([-10.00, 5.25, 7.17])
  if (lsr is not None):   
    u = u + lsr_vel[0]
    v = v + lsr_vel[1]
    w = w + lsr_vel[2]
  return (u,v,w)

def helcorr(ra2000, dec2000, jd, obs_long=None, obs_lat=None, obs_alt=None, debug=False):
  """
  calculates heliocentric Julian date, baricentric and heliocentric radial
  velocity corrections from:

  INPUT:
  @param ra2000: Right ascension of object for epoch 2000.0 (hours)
  @type ra2000: float
  @param dec2000: declination of object for epoch 2000.0 (degrees)
  @type dec2000: float
  @param jd: julian date for the middle of exposure
  @type jd: float
  @keyword obs_long: longitude of observatory (degrees, western direction is positive)
  @type obs_long: float
  @keyword obs_lat: latitude of observatory (degrees)
  @type obs_lat: float
  @keyword obs_alt: altitude of observatory (meters)
  @type obs_alt: float
  
  Algorithms used are taken from the IRAF task noao.astutils.rvcorrect
  and some procedures of the IDL Astrolib are used as well.
  Accuracy is about 0.5 seconds in time and about 1 m/s in velocity.
  
  History:
  written by Peter Mittermayer, Nov 8,2003
  2005-January-13   Kudryavtsev   Made more accurate calculation of the sideral time.
                                  Conformity with MIDAS compute/barycorr is checked.
  2005-June-20      Kochukhov Included precession of RA2000 and DEC2000 to current epoch
  """
  _radeg = 180.0 / np.pi
  
  # If observatory not given, set La Palma
  if obs_long is None:
    obs_long = 17.878333
  if obs_lat is None:
    obs_lat = 28.762222 
  if obs_alt is None:
    obs_alt = 2333.
  
  #covert JD to Gregorian calendar date
  xjd = jd*1.0
  jd  = jd-2400000.0
  year, month, dayraw = convert("JD", "CD", xjd)
  day = dayraw-dayraw%24.
  ut  = 24.*(dayraw%24)
  
  #current epoch
  epoch = year + month / 12. + day / 365.
  
  #precess ra2000 and dec2000 to current epoch
  ra,dec=precess(ra2000*15., dec2000, 2000.0, epoch)
  
  #calculate heliocentric julian date
  hjd = np.array(helio_jd(jd, ra, dec)).astype(float)
  
  #DIURNAL VELOCITY (see IRAF task noao.astutil.rvcorrect)
  dlat = -(11. * 60. + 32.743) * np.sin(2 * obs_lat / _radeg) + 1.1633 * np.sin(4 * obs_lat / _radeg) - 0.0026 * np.sin(6 * obs_lat / _radeg)
  lat = obs_lat + dlat / 3600
  
  #calculate distance of observer from earth center
  r = 6378160.0 * (0.998327073 + 0.001676438 * np.cos(2 * lat / _radeg) - 0.00000351 * np.cos(4 * lat / _radeg) + 0.000000008 * np.cos(6 * lat / _radeg)) + obs_alt
  
  #calculate rotational velocity (perpendicular to the radius vector) in km/s
  #23.934469591229 is the siderial day in hours for 1986
  v = 2. * np.pi * (r / 1000.) / (23.934469591229 * 3600.)
  
  #calculating local mean siderial time (see astronomical almanach)
  tu = (jd - 51545.0) / 36525
  gmst = 6.697374558 + ut + (236.555367908 * (jd - 51545.0) + 0.093104 * tu ** 2 - 6.2e-6 * tu ** 3) / 3600
  lmst = (gmst - obs_long / 15) % 24
  
  #projection of rotational velocity along the line of sight
  vdiurnal = v * np.cos(lat / _radeg) * np.cos(dec / _radeg) * np.sin((ra - lmst * 15) / _radeg)
  
  #BARICENTRIC and HELIOCENTRIC VELOCITIES
  vh,vb=baryvel(xjd, 0)
  
  #project to line of sight
  vbar = vb[0] * np.cos(dec / _radeg) * np.cos(ra / _radeg) + vb[1] * np.cos(dec / _radeg) * np.sin(ra / _radeg) + vb[2] * np.sin(dec / _radeg)
  vhel = vh[0] * np.cos(dec / _radeg) * np.cos(ra / _radeg) + vh[1] * np.cos(dec / _radeg) * np.sin(ra / _radeg) + vh[2] * np.sin(dec / _radeg)
  corr = (vdiurnal + vbar) #using baricentric velocity for correction 
  
  return (corr, hjd)

def helio_jd(date, ra, dec, b1950=False, time_diff=False):
  """
  Convert geocentric (reduced) Julian date to heliocentric Julian date
  
  This procedure correct for the extra light travel time between the Earth and the Sun.
  An online calculator for this quantity is available at http://www.physics.sfasu.edu/astro/javascript/hjd.html
  
  Input:
    @param date: Reduced Julian date (= JD - 2400000)
    @type date: float
    @param ra: J2000 right ascension in DEGREES (unless B1950 keyword is set)
    @type ra: float
    @param dec: J2000 declination in DEGREES (unless B1950 keyword is set)
    @type dec: float
    @keyword b1950: if set, then input coordinates are assumed to be in equinox B1950 coordinates
    @type b1950: bool
    @keyword time_diff: if set, then HELIO_JD() returns the time difference (heliocentric JD - geocentric JD ) in seconds
    @type time_diff: bool
    @return: heliocentric reduced Julian date.  If time_diff is True, then helio_jd() instead returns the time difference in seconds between the geocentric and heliocentric Julian date.
  """
  #Because XYZ uses default B1950 coordinates, we'll convert everything to B1950
  if not b1950:
    ra1, dec1 = bprecess(ra, dec)
  else:   
    ra1 = ra
    dec1 = dec
    
  delta_t     = (np.array(date).astype(float) - 33282.42345905e0) / 36525.0e0
  epsilon_sec = np.poly1d([44.836e0, -46.8495, -0.00429, 0.00181][::-1])(delta_t)
  epsilon     = np.deg2rad(23.433333e0 + epsilon_sec / 3600.0e0)
  ra1         = np.deg2rad(ra1)
  dec1        = np.deg2rad(dec1)
  x, y, z, tmp, tmp, tmp = xyz(date)
  
  #Find extra distance light must travel in AU, multiply by 1.49598e13 cm/AU,
  #and divide by the speed of light, and multiply by 86400 second/year
  time = -499.00522e0 * (np.cos(dec1) * np.cos(ra1) * x + (np.tan(epsilon) * np.sin(dec1) + np.cos(dec1) * np.sin(ra1)) * y)
  if time_diff:
    return time
  else:
    return np.array(date).astype(float) + time / 86400.0e0
   
def precess(ra0, dec0, equinox1, equinox2, doprint=False, fk4=False, radian=False):
  """
  Precess coordinates from EQUINOX1 to EQUINOX2.

  For interactive display, one can use the procedure ASTRO which calls PRECESS
  or use the /PRINT keyword. The default (RA,DEC) system is FK5 based on epoch
  J2000.0 but FK4 based on B1950.0 is available via the /FK4 keyword.
  
  Input:
    @param ra0: Input right ascension in DEGREES
    @type ra0: float
    @param dec0: Input declination in degrees
    @type dec0: float
    @param equinox1: first equinox
    @type equinox1: float
    @param equinox2: second equinox
    @type equinox2: float
    @keyword fk4: If this keyword is set and non-zero, the FK4 (B1950.0) system will be used otherwise FK5 (J2000.0) will be used instead.
    @type fk4: float
    @keyword radian: If this keyword is set, input is in radian instead of degrees
    @type radian
    
  Restrictions:
    Accuracy of precession decreases for declination values near 90 degrees.
    PRECESS should not be used more than 2.5 centuries from 2000 on the FK5 system (1950.0 on the FK4 system).
  
  Procedure:
    Algorithm from Computational Spherical Astronomy by Taff (1983), p24. (FK4).
    FK5 constants from "Astronomical Almanac Explanatory Supplement 1992, page 104 Table 3.211.1.
  
  Revision history
    Written, Wayne Landsman, STI Corporation  August 1986
    Correct negative output RA values   February 1989
    Provided FK5 (J2000.0)  I. Freedman   January 1994
    Precession Matrix computation now in PREMAT   W. Landsman June 1994
    Work for arrays, not just vectors (W. Landsman, september 2003)
    Convert to Python (S. Koposov, july 2010)
    Converted for use at IvS (K. Smolders)
  """
  scal = True
  if isinstance(ra0, np.ndarray):
    ra   = ra0.copy()  
    dec  = dec0.copy()
    scal = False
  else:
    ra  = np.array([ra0])
    dec = np.array([dec0])
  npts = ra.size 
  
  if not radian:   
    ra_rad  = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)
  else:   
    ra_rad  = ra
    dec_rad = dec
  
  a      = np.cos(dec_rad)
  x      = np.zeros((npts, 3))
  x[:,0] = a * np.cos(ra_rad)
  x[:,1] = a * np.sin(ra_rad)
  x[:,2] = np.sin(dec_rad)
  
  # Use PREMAT function to get precession matrix from Equinox1 to Equinox2
  r       = premat(equinox1, equinox2, fk4=fk4)
  x2      = np.transpose(np.dot(np.transpose(r), np.transpose(x)))
  ra_rad  = np.zeros(npts) + np.arctan2(x2[:,1], x2[:,0])
  dec_rad = np.zeros(npts) + np.arcsin(x2[:,2])
  
  if not radian:   
    ra  = np.rad2deg(ra_rad)
    ra  = ra + (ra < 0.) * 360.e0
    dec = np.rad2deg(dec_rad)
  else:   
    ra  = ra_rad
    dec = dec_rad
    ra  = ra + (ra < 0.) * 2.0e0 * np.pi
  
  if doprint:   
    print 'Equinox (%.2f): %f,%f' % (equinox2, ra, dec)
  if scal:
    ra, dec = ra[0], dec[0]
  return ra, dec


def precess_xyz(x, y, z, equinox1, equinox2):
  """
  Precess equatorial geocentric rectangular coordinates.
  
  Input:
  @param x: heliocentric rectangular x-coordinates
  @type x: float
  @param y: heliocentric rectangular y-coordinates
  @type y: float
  @param z: heliocentric rectangular z-coordinates
  @type z: float
  @param equinox1: first equinox
  @type equinox1: float
  @param equinox2: second equinox
  @type equinox2: float

  NOTES:
    The equatorial geocentric rectangular coords are converted to RA and Dec,
    precessed in the normal way, then changed back to x, y and z using unit
    vectors.
  """
  #take input coords and convert to ra and dec (in radians)   
  ra   = np.arctan2(y, x)
  _del = np.sqrt(x * x + y * y + z * z)  #magnitude of distance to Sun
  dec  = np.arcsin(z / _del)
  
  #   precess the ra and dec
  ra,dec = precess(ra, dec, equinox1, equinox2, radian=True)
  
  #convert back to x, y, z
  xunit = np.cos(ra) * np.cos(dec)
  yunit = np.sin(ra) * np.cos(dec)
  zunit = np.sin(dec)  
  x = xunit * _del
  y = yunit * _del
  z = zunit * _del
  return x,y,z

def premat(equinox1, equinox2, fk4=False):
  """
  Return the precession matrix needed to go from EQUINOX1 to EQUINOX2.
  
  This matrix is used by the procedures PRECESS and BARYVEL to precess astronomical coordinates
  @param equinox1: Original equinox of coordinates, numeric scalar.
  @type equinox1: float
  @param equinox2: Equinox of precessed coordinates.
  @type equinox2: float
  @return: double precision 3 x 3 precession matrix, used to precess equatorial rectangular coordinates
  @rtype: 3 by 3 array
  @keyword fk4: If this keyword is set, the FK4 (B1950.0) system precession angles are used to compute the precession matrix. The default is to use FK5 (J2000.0) precession angles
  
  Revision history:
    Written, Wayne Landsman, HSTX Corporation, June 1994
    Converted to IDL V5.0   W. Landsman   September 1997
  """
  deg_to_rad = np.pi / 180.0e0
  sec_to_rad = deg_to_rad / 3600.e0
  t          = 0.001e0 * (equinox2 - equinox1)
  
  if not fk4:   
    st = 0.001e0 * (equinox1 - 2000.e0)
    #  Compute 3 rotation angles
    a = sec_to_rad * t * (23062.181e0 + st * (139.656e0 + 0.0139e0 * st) + t * (30.188e0 - 0.344e0 * st + 17.998e0 * t))
    b = sec_to_rad * t * t * (79.280e0 + 0.410e0 * st + 0.205e0 * t) + a
    c = sec_to_rad * t * (20043.109e0 - st * (85.33e0 + 0.217e0 * st) + t * (-42.665e0 - 0.217e0 * st - 41.833e0 * t))
  else:   
    st = 0.001e0 * (equinox1 - 1900.e0)
    #  Compute 3 rotation angles
    a = sec_to_rad * t * (23042.53e0 + st * (139.75e0 + 0.06e0 * st) + t * (30.23e0 - 0.27e0 * st + 18.0e0 * t))
    b = sec_to_rad * t * t * (79.27e0 + 0.66e0 * st + 0.32e0 * t) + a
    c = sec_to_rad * t * (20046.85e0 - st * (85.33e0 + 0.37e0 * st) + t * (-42.67e0 - 0.37e0 * st - 41.8e0 * t))
    
  sina   = np.sin(a)
  sinb   = np.sin(b)
  sinc   = np.sin(c)
  cosa   = np.cos(a)
  cosb   = np.cos(b)
  cosc   = np.cos(c)
  r      = np.zeros((3, 3))
  r[0,:] = np.array([cosa * cosb * cosc - sina * sinb, sina * cosb + cosa * sinb * cosc, cosa * sinc])
  r[1,:] = np.array([-cosa * sinb - sina * cosb * cosc, cosa * cosb - sina * sinb * cosc, -sina * sinc])
  r[2,:] = np.array([-cosb * sinc, -sinb * sinc, cosc])
  
  return r


def sphdist (ra1, dec1, ra2, dec2):
  """
  Measures the spherical distance in degrees. The input has to be in degrees too
  """
  dec1_r = np.deg2rad(dec1)
  dec2_r = no.deg2rad(dec2)
  return 2 * np.rad2deg(np.arcsin(np.sqrt((np.sin((dec1_r - dec2_r) / 2))**2+np.cos(dec1_r) * np.cos(dec2_r) *(np.sin((np.deg2rad(ra1 - ra2)) / 2))**2)))

def xyz(date, equinox=None):
  """
  Calculate geocentric X,Y, and Z  and velocity coordinates of the Sun
  
  Calculates geocentric X,Y, and Z vectors and velocity coordinates
  (dx, dy and dz) of the Sun.   (The positive X axis is directed towards the
  equinox, the y-axis, towards the point on the equator at right ascension 6h,
  and the z axis toward the north pole of the equator). Typical position
  accuracy is <1e-4 AU (15000 km).
  
  Input:
  @param date: reduced julian date (=JD - 2400000)
  @type date: float
  @keyword equinox: equinox of output. Default is 1950.
  @type equinox: float
  
  Output:
    x,y,z: scalars or vectors giving heliocentric rectangular coordinates
           (in A.U) for each date supplied. Note that sqrt(x^2 + y^2 + z^2)
           gives the Earth-Sun distance for the given date.
    xvel, yvel, zvel: velocity vectors corresponding to X, Y and Z.
  """
  picon = np.pi / 180.0e0
  t     = (date - 15020.0e0) / 36525.0e0         #Relative Julian century from 1900
  
  # NOTE: longitude arguments below are given in *equinox* of date.
  #   Precess these to equinox 1950 to give everything an even footing.
  #   Compute argument of precession from equinox of date back to 1950
  pp = (1.396041e0 + 0.000308e0 * (t + 0.5e0)) * (t - 0.499998e0)
  
  # Compute mean solar longitude, precessed back to 1950
  el = 279.696678e0 + 36000.76892e0 * t + 0.000303e0 * t * t - pp
  
  # Compute Mean longitude of the Moon
  c = 270.434164e0 + 480960.e0 * t + 307.883142e0 * t - 0.001133e0 * t * t - pp
  
  # Compute longitude of Moon's ascending node
  n = 259.183275e0 - 1800.e0 * t - 134.142008e0 * t + 0.002078e0 * t * t - pp
  
  # Compute mean solar anomaly
  g = 358.475833e0 + 35999.04975e0 * t - 0.00015e0 * t * t
  
  # Compute the mean jupiter anomaly
  j = 225.444651e0 + 2880.0e0 * t + 154.906654e0 * t * t
  
  # Compute mean anomaly of Venus
  v = 212.603219e0 + 58320.e0 * t + 197.803875e0 * t + 0.001286e0 * t * t
  
  # Compute mean anomaly of Mars
  m = 319.529425e0 + 19080.e0 * t + 59.8585e0 * t + 0.000181e0 * t * t
  
  # Convert degrees to radians for trig functions
  el = el * picon
  g = g * picon
  j = j * picon
  c = c * picon
  v = v * picon
  n = n * picon
  m = m * picon
  
  # Calculate X,Y,Z using trigonometric series
  x = 0.999860e0 * np.cos(el) - 0.025127e0 * np.cos(g - el) + 0.008374e0 * np.cos(g + el) + 0.000105e0 * np.cos(g + g + el) + 0.000063e0 * t * np.cos(g - el) + 0.000035e0 * np.cos(g + g - el) - 0.000026e0 * np.sin(g - el - j) - 0.000021e0 * t * np.cos(g + el) + 0.000018e0 * np.sin(2.e0 * g + el - 2.e0 * v) + 0.000017e0 * np.cos(c) - 0.000014e0 * np.cos(c - 2.e0 * el) + 0.000012e0 * np.cos(4.e0 * g + el - 8.e0 * m + 3.e0 * j) - 0.000012e0 * np.cos(4.e0 * g - el - 8.e0 * m + 3.e0 * j) - 0.000012e0 * np.cos(g + el - v) + 0.000011e0 * np.cos(2.e0 * g + el - 2.e0 * v) + 0.000011e0 * np.cos(2.e0 * g - el - 2.e0 * j)
  y = 0.917308e0 * np.sin(el) + 0.023053e0 * np.sin(g - el) + 0.007683e0 * np.sin(g + el) + 0.000097e0 * np.sin(g + g + el) - 0.000057e0 * t * np.sin(g - el) - 0.000032e0 * np.sin(g + g - el) - 0.000024e0 * np.cos(g - el - j) - 0.000019e0 * t * np.sin(g + el) - 0.000017e0 * np.cos(2.e0 * g + el - 2.e0 * v) + 0.000016e0 * np.sin(c) + 0.000013e0 * np.sin(c - 2.e0 * el) + 0.000011e0 * np.sin(4.e0 * g + el - 8.e0 * m + 3.e0 * j) + 0.000011e0 * np.sin(4.e0 * g - el - 8.e0 * m + 3.e0 * j) - 0.000011e0 * np.sin(g + el - v) + 0.000010e0 * np.sin(2.e0 * g + el - 2.e0 * v) - 0.000010e0 * np.sin(2.e0 * g - el - 2.e0 * j)
  z = 0.397825e0 * np.sin(el) + 0.009998e0 * np.sin(g - el) + 0.003332e0 * np.sin(g + el) + 0.000042e0 * np.sin(g + g + el) - 0.000025e0 * t * np.sin(g - el) - 0.000014e0 * np.sin(g + g - el) - 0.000010e0 * np.cos(g - el - j)
  
  #Precess_to new equator?
  if equinox is not None:   
     x, y, z = precess_xyz(x, y, z, 1950, equinox)
  
  xvel = -0.017200e0 * np.sin(el) - 0.000288e0 *  np.sin(g + el) - 0.000005e0 * np.sin(2.e0 * g + el) - 0.000004e0 * np.sin(c) + 0.000003e0 * np.sin(c - 2.e0 * el) + 0.000001e0 * t * np.sin(g + el) - 0.000001e0 * np.sin(2.e0 * g - el)
  yvel = 0.015780 * np.cos(el) + 0.000264 * np.cos(g + el) + 0.000005 * np.cos(2.e0 * g + el) + 0.000004 * np.cos(c) + 0.000003 * np.cos(c - 2.e0 * el) - 0.000001 * t * np.cos(g + el)
  zvel = 0.006843 * np.cos(el) + 0.000115 * np.cos(g + el) + 0.000002 * np.cos(2.e0 * g + el) + 0.000002 * np.cos(c) + 0.000001 * np.cos(c - 2.e0 * el)
  
  #Precess to new equator?
  if equinox is not None:   
    xvel, yvel, zvel = precess_xyz(xvel, yvel, zvel, 1950, equinox)
  
  return x, y, z, xvel, yvel, zvel
