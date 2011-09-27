# -*- coding: utf-8 -*-
"""
List of physical constants (SI, CGS and solar)

Many constants are taken from NIST (http://physics.nist.gov/cuu/Constants/index.html)

Others have their explicit reference listed.
"""
# SI-units
# value                    name                          unit          reference
#===============================================================================
cc     = 299792458.        # speed of light              m/s
Msol   = 1.988547e30       # solar mass                  kg            Harmanec & Prsa 2011
Rsol   = 6.95508e8         # solar radius                m             Harmanec & Prsa 2011
Lsol   = 3.846e26          # solar luminosity            W             Harmanec & Prsa 2011
Tsol   = 5779.5747         # solar effective temperature K             Harmanec & Prsa 2011
Mearth = 5.9742e24         # earth mass                  kg            google
Rearth = 6371e3            # earth radius                m             wikipedia
Mjup   = 1.8986e27         # Jupiter mass                kg            wikipedia
Mlun   = 7.346e22          # Lunar mass                  kg            wikipedia
au     = 149597870700.     # astronomical unit           m             Harmanec & Prsa 2011
pc     = 3.085677581503e+16# parsec                      m             Harmanec & Prsa 2011
ly     = 9.4605284e+15     # light year                  m
hh     = 6.6260689633e-34  # Planck constant             J/s
hhbar  = 1.05457162853e-34 # reduced Planck constant     J/s
kB     = 1.380650424e-23   # Boltzmann constant          J/K
NA     = 6.0221417930e23   # Avogadro constant           1/mol
sigma  = 5.67040040e-8     # Stefan-Boltzmann constant   W/m2/K4       Harmanec & Prsa 2011
GG     = 6.67384e-11       # gravitational constant      m3/kg/s2      Harmanec & Prsa 2011
GGMsol = 1.32712442099     # grav. constant x solar mass m3/s2         Harmanec & Prsa 2011
RR     = 8.31447215        # (ideal) gas constant        J/K/mol
aa     = 7.5657e-16        # radiation constant          J/m3/K4
a0     = 52.9177e-12       # Bohr radius of hydrogen     m

# CGS-units
# value                        name                          unit           reference
#====================================================================================
cc_cgs     = 29979245800.      # speed of light              cm/s
Msol_cgs   = 1.988547e33       # solar mass                  g              Harmanec & Prsa 2011
Rsol_cgs   = 6.95508e10        # solar radius                cm             Harmanec & Prsa 2011
Lsol_cgs   = 3.846e33          # solar luminosity            erg/s          Harmanec & Prsa 2011
Mearth_cgs = 5.9742e27         # earth mass                  g              google
Rearth_cgs = 6371e5            # earth radius                cm             wikipedia
Mjup_cgs   = 1.8986e30         # Jupiter mass                g              wikipedia
Mlun_cgs   = 7.346e25          # Lunar mass                  g              wikipedia
au_cgs     = 149597.870700e8   # astronomical unit           cm             Harmanec & Prsa 2011
pc_cgs     = 3.085677581503e+18# parsec                      cm             Harmanec & Prsa 2011
ly_cgs     = 9.4605284e+17     # light year                  cm
hh_cgs     = 6.6260689633e-27  # Planck constant             erg/s
hhbar_cgs  = 1.05457162853e-27 # reduced Planck constant     erg/s
kB_cgs     = 1.380650424e-16   # Boltzmann constant          erg/K
sigma_cgs  = 5.67040040e-5     # Stefan-Boltzmann constant   erg/cm2/s/K4   Harmanec & Prsa 2011
GG_cgs     = 6.67384e-8        # gravitational constant      cm3/g/s2       Harmanec & Prsa 2011
GGMsol_cgs = 1327124.42099     # grav. constant x solar mass cm3/s2         Harmanec & Prsa 2011
RR_cgs     = 8.31447215e7      # (ideal) gas constant        erg/K/mol
aa_cgs     = 7.5657e-15        # radiation constant          erg/cm2/K4
a0_cgs     = 52.9177e-10       # Bohr radius of hydrogen     cm

# solar units
# value                        name                          unit           reference
#====================================================================================
Msol_sol  = 1.                 # solar mass                  Msol           Harmanec & Prsa 2011
GG_sol    = 3.944620808727e-07 # gravitional constant        Rsol3/Msol/s2  Harmanec & Prsa 2011
GGMsol_sol= 3.944620719386e-27 # grav. constant x solar mass Rsol3/s2       Harmanec & Prsa 2011

# other stuff
Mabs_sol = 4.75                # solar bolometric abs mag    mag            Harmanec & Prsa 2011
Mapp_sol = -26.74              # solar visual apparent mag   mag