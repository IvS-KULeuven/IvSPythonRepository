"""
Access atomic data

Mendelev's Table::

    1  H  Hydrogen   11 Na Natrium    21 Sc Scandium  31 Ga Gallium    41 Nb Niobium    51 Sb Antimony
    2  He Helium     12 Mg Magnesium  22 Ti Titanium  32 Ge Germanium  42 Mo Molybdenum 52 Te Tellurium
    3  Li Lithium    13 Al Aluminium  23 V  Vanadium  33 As Arsenic    43 Tc Technetium 53 I  Iodine
    4  Be Beryllium  14 Si Silicon    24 Cr Chromium  34 Se Selenium   44 Ru Ruthenium  54 Xe Xenon
    5  B  Boron      15 P  Phosphorus 25 Mn Manganese 35 Br Bromine    45 Rh Rhodium    55 Cs Caesium
    6  C  Carbon     16 S  Sulfur     26 Fe Iron      36 Kr Krypton    46 Pd Palladium  56 Ba Barium
    7  N  Nitrogen   17 Cl Chlorine   28 Ni Nickel    37 Rb Rubidium   47 Ag Silver     57 La Lanthanum
    8  O  Oxygen     18 Ar Argon      27 Co Cobalt    38 Sr Strontium  48 Cd Cadmium    58 Ce Cerium
    9  F  Fluorine   19 K  Potassium  29 Cu Copper    39 Y  Yttrium    49 In Indium     59 Pr Praseodymium
    10 Ne Neon       20 Ca Calcium    30 Zn Zinc      40 Zr Zirconium  50 Sn Tin        60 Nd Neodymium

    61 Pm Promethium 71 Lu Lutetium   81 Tl Thallium  91 Pa Protactinium 101 Md Mendelevium    111 Rg Roentgenium
    62 Sm Samarium   72 Hf Hafnium    82 Pb Lead      92 U  Uranium      102 No Nobelium       112 UubUnunbium
    63 Eu Europium   73 Ta Tantalum   83 Bi Bismuth   93 Np Neptunium    103 Lr Lawrencium     113 UutUnuntrium
    64 Gd Gadolinium 74 W  Tungsten   84 Po Polonium  94 Pu Plutonium    104 Rf Rutherfordium  114 UuqUnunquadium
    65 Tb Terbium    75 Re Rhenium    85 At Astatine  95 Am Americium    105 Db Dubnium        115 UupUnunpentium
    66 Dy Dysprosium 76 Os Osmium     86 Rn Radon     96 Cm Curium       106 Sg Seaborgium     116 UuhUnunhexium
    67 Ho Holmium    77 Ir Iridium    87 Fr Francium  97 Bk Berkelium    107 Bh Bohrium        118 UuoUnunoctium
    68 Er Erbium     78 Pt Platinum   88 Ra Radium    98 Cf Californium  108 Hs Hassium
    69 Tm Thulium    79 Au Gold       89 Ac Actinium  99 Es Einsteinium  109 Mt Meitnerium
    70 Yb Ytterbium  80 Hg Mercury    90 Th Thorium   100 Fm Fermium     110 Ds Darmstadtium

End
"""
import os

import numpy as np
import pylab as pl
from numpy import inf

from ivs import config
from ivs.io import ascii

atomcode = ['X', 'H','He','Li','Be', 'B', 'C', 'N', 'O', 'F','Ne',
                'Na','Mg','Al','Si', 'P', 'S','Cl','Ar', 'K','Ca',
                'Sc','Ti', 'V','Cr','Mn','Fe','Ni','Co','Cu','Zn',
                'Ga','Ge','As','Se','Br','Kr','Rb','Sr', 'Y','Zr',
                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
                'Sb','Te', 'I','Xe','Cs','Ba','La','Ce','Pr','Nd',
                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
                'Lu','Hf','Ta', 'W','Re','Os','Ir','Pt','Au','Hg',
                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
                'Pa', 'U','Np','Pu','Am','Cu','Bk','Cf','Es','Fm',
                'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
                'Rg']
roman = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII']



stellar = os.path.join('linelists/mask')

#{ Main functions

def get_lines(teff,logg,z=0,atoms=None,ions=None,wrange=(-inf,inf),\
                blend=0.0):
    """
    Retrieve line transitions and strengths for a specific stellar type
    
    Selection wavelength range in angstrom.
    
    Ions should be a list of ions to include. This can either be a string or
    a number
    
    A lines is considerd a blend if the closest line is closer than C{blend} angstrom.
    
    Example usage:
    
    Retrieve all Silicon lines between 4500 and 4600 for a B1V star.
    
    >>> data = get_lines(20000,4.0,atoms=['Si'],wrange=(4500,4600))
    >>> p = pl.figure()
    >>> p = pl.vlines(data['wavelength'],1,1-data['depth'])
    
    See how the depth of the Halpha line varies wrt temperature:
    
    >>> teffs = range(5000,21000,1000) + range(22000,32000,2000) + range(30000,50000,50000)
    >>> depths = np.zeros((len(teffs),7))
    >>> for i,teff in enumerate(teffs):
    ...     data = get_lines(teff,5.0,ions=['HI'],wrange=(3800,7000))
    ...     depths[i] = data['depth']
    
    >>> p = pl.figure();p = pl.title('Depth of Balmer lines (Halpha-Heta)')
    >>> p = pl.plot(teffs,1-depths,'o-')
    >>> p = pl.xlabel('Effective temperature');p = pl.grid()
    
    """
    #-- get filepath
    filename = 'mask.%d.%02d.p%02d'%(int(teff),int(logg*10),int(z))
    filename = config.get_datafile(stellar,filename)
    #-- read in the data and extract relevant columns
    data = ascii.read2recarray(filename,skip_lines=1,dtype=[('wavelength','f8'),('ion','f8'),('depth','f8'),('c3','f8'),('c4','f8'),('c5','f8')])
    data = pl.mlab.rec_drop_fields(data,['c3','c4','c5'])
    data['wavelength'] *= 10.
    #-- remove blends
    if blend>0:
        blends_left = np.hstack([0,np.diff(data['wavelength'])])
        blends_right= np.hstack([np.diff(data['wavelength']),1e10])
        keep = (blends_left>blend) & (blends_right>blend)
        data = data[keep]
    
    #-- only keep those transitions within a certain wavelength range
    keep = (wrange[0]<=data['wavelength']) & (data['wavelength']<=wrange[1])
    data = data[keep]
    #-- only keep those transitions that belong to certain ions or atoms
    if atoms is not None or ions is not None:
        keep = np.array(np.zeros(len(data)),bool)
    else:
        keep = np.array(np.ones(len(data)),bool)
    if atoms is not None:
        #-- convert all atoms to their number and select the appropriate ones
        atoms = [(isinstance(atom,str) and atomcode.index(atom.title()) or atom) for atom in atoms]
        for atom in atoms:
            keep = keep | (np.abs(data['ion']-atom)<0.5)
    if ions is not None:
        #-- convert all ions to their number and select the appropriate ones
        ions = [(isinstance(ion,str) and name2ioncode(ion) or ion) for ion in ions]
        for ion in ions:
            keep = keep | (np.abs(data['ion']-ion)<0.005)
    
    return data[keep]


#}

#{ Convert ion codes to ion names and back

def ioncode2name(ioncode):
    """
    Convert 14.01 to SiII
    """
    atom = int(np.floor(ioncode))
    ion = int(np.round((ioncode-atom)*100))
    return atomcode[atom]+roman[ion]

def name2ioncode(name):
    """
    Convert SiII to 14.01
    """
    atomname,ionstage = splitname(name)
    return atomcode.index(atomname) + roman.index(ionstage)/100.

def splitname(name):
    """
    Split SiII into Si and II
    """
    atomname = ''
    ionstage = ''
    for char in name:
        if not atomname or char.islower():
            atomname+=char
        elif char.isupper():
            ionstage+=char
    return atomname,ionstage


#}


if __name__=="__main__":
    import doctest
    doctest.testmod()
    pl.show()
    