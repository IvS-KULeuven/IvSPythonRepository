"""
Retrieve linelists from VALD.

With this program, you can request VALD lines files for all available ions in a wavelength range from 1000 to 10500 Angstrom.
The requested linesfiles contain the following data:
Wavelength (Angstrom), atom number Z + ionisation (.0 = neutral, .1 = first ionisation, .2 = second ionisation, etc.), excitation potential (eV), oscillator strength (log gf)
There are two request options: 
Option 1: request one or multiple linelists of ions seperately in specific wavelength ranges.
Option 2: request a wavelength range with all known spectral lines in that wavelength range.
"""
import logging
import os
from ivs import config
from ivs.io import ascii
import numpy as np
#Option 1

logger = logging.getLogger("IVS.SPECTRA")

# With this definition, you can request linelists for each ion seperately within a specific wavelength range.
# elem = an array of ions e.g. ['CI','OII'], xmin and xmax: wavelength range in which the spectral lines are searched, outputdir = output directory chosen by the user.

def Single_VALD_files(elem=None,xmin=3200.,xmax=4800.,outputdir=None):
  """
  Request linelists for each ion seperately within a specific wavelength range.
  
  elem = an array of ions e.g. ['CI','OII'], xmin and xmax: wavelength range in which the spectral lines are searched, outputdir = output directory chosen by the user.
  
  If no elements are given, this function returns all of them.
  
  @param elem: list of ions
  @type elem: list of str
  """
  if elem is None:
    files = sorted(config.glob('VALD_individual','VALD_*.lijnen'))
    elem = [os.path.splitext(os.path.basename(ff))[0].split('_')[1] for ff in files]
  
  all_lines = []
  for i in range(len(elem)):
    print elem[i]
    filename = config.get_datafile('VALD_individual','VALD_' + elem[i] + '.lijnen')
    if not os.path.isfile(filename):
      logger.info('No data for element ' + str(elem[i]))
      return None
   
    newwav,newexc,newep,newgf = np.loadtxt(filename).T
    lines = np.rec.fromarrays([newwav,newexc,newep,newgf],names=['wavelength','ion','ep','gf'])
    keep = (xmin<=lines['wavelength']) & (lines['wavelength']<=xmax)
    if not hasattr(keep,'__iter__'):
      continue
    lines = lines[keep]
    if len(lines) and outputdir is not None:
      ascii.write_array(lines,outputdir + 'VALD_' + str(elem[i]) + '_' + str(xmin) + '_' + str(xmax) + '.dat',auto_width=True,header=True,formats=['%.3f','%.1f','%.3f','%.3f'])
    elif len(lines):
      all_lines.append(lines)
    else:
      logger.info('No lines of ' + str(elem[i]) + ' in given wavelength range')
  return np.hstack(all_lines)