# -*- coding: utf-8 -*-
"""
List of physical constants (SI and CGS)

Taken from NIST (http://physics.nist.gov/cuu/Constants/index.html)
"""
# SI-units
# value                    name                          unit          reference
#===============================================================================
cc     = 299792458.        # speed of light              m/s
Msol   = 1.98892e30        # solar mass                  kg
Rsol   = 6.955e8           # solar radius                m
Lsol   = 3.839e26          # solar luminosity            W
Tsol   = 5778.             # solar effective temperature K
Mearth = 5.9742e24         # earth mass                  kg            google
Rearth = 6371e3            # earth radius                m             wikipedia
Mjup   = 1.8986e27         # Jupiter mass                kg            wikipedia
Mlun   = 7.346e22          # Lunar mass                  kg            wikipedia
au     = 149598.e6         # astronomical unit           m
pc     = 3.08568025e+16    # parsec                      m
ly     = 9.4605284e+15     # light year                  m
hh     = 6.6260689633e-34  # Planck constant             J/s
hhbar  = 1.05457162853e-34 # reduced Planck constant     J/s
kB     = 1.380650424e-23   # Boltzmann constant          J/K
NA     = 6.0221417930e23   # Avogadro constant           1/mol
sigma  = 5.67040040e-8     # Stefan-Boltzmann constant   W/m2/K4
GG     = 6.6742867e-11     # gravitational constant      m3/kg/s2
RR     = 8.31447215        # (ideal) gas constant        J/K/mol
aa     = 7.5657e-16        # radiation constant          J/m3/K4
a0     = 52.9177e-12       # Bohr radius of hydrogen     m

# CGS-units
# value                        name                          unit           reference
#====================================================================================
cc_cgs     = 29979245800.      # speed of light              cm/s
Msol_cgs   = 1.98892e33        # solar mass                  g
Rsol_cgs   = 6.955e10          # solar radius                cm
Lsol_cgs   = 3.839e33          # solar luminosity            erg/s
Mearth_cgs = 5.9742e27         # earth mass                  g              google
Rearth_cgs = 6371e5            # earth radius                cm             wikipedia
Mjup_cgs   = 1.8986e30         # Jupiter mass                g              wikipedia
Mlun_cgs   = 7.346e25          # Lunar mass                  g              wikipedia
au_cgs     = 149598.e8         # astronomical unit           cm
pc_cgs     = 3.08568025e+18    # parsec                      cm
ly_cgs     = 9.4605284e+17     # light year                  cm
hh_cgs     = 6.6260689633e-27  # Planck constant             erg/s
hhbar_cgs  = 1.05457162853e-27 # reduced Planck constant     erg/s
kB_cgs     = 1.380650424e-16   # Boltzmann constant          erg/K
sigma_cgs  = 5.67040040e-5     # Stefan-Boltzmann constant   erg/cm2/s/K4
GG_cgs     = 6.67428e-8        # gravitational constant      cm3/g/s2
RR_cgs     = 8.31447215e7      # (ideal) gas constant        erg/K/mol
aa_cgs     = 7.5657e-15        # radiation constant          erg/cm2/K4
a0_cgs     = 52.9177e-10       # Bohr radius of hydrogen     cm

# other stuff
Mabs_sol = 4.83                # solar bolometric abs mag    mag
Mapp_sol = -26.74              # solar visual apparent mag   mag