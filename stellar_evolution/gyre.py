import logging
import numpy as np
from ivs.units import constants
from ivs.units import conversions
from ivs.aux import numpy_ext as ne

logger = logging.getLogger('STAR.GYRE')

def make_gyre_inputfile(**kwargs):
    """
    keyword arguments:
    model_name
    degree
    omega_min
    omega_max
    n_omega
    """
    control_file = 'gyre_nad.in'
    contents = """
&eqmodel
        file = '{model_name}'
        file_type = 'B3'
        G = {G:.10e}
/

&grid
/

&bvp
        l = {degree}
        scheme = 'MAGNUS_GL2'
/

&freqs
        omega_min = {omega_min:f}
        omega_max = {omega_max:f}
        n_omega = {n_omega:d}
        omega_grid = 'MIXED'
/

&output
      eigval_file = 'eigvals.h5'
      eigfunc_prefix = 'eigfunc-'

/
"""
    with open(control_file,'w') as ff:
        ff.write(contents.format(**kwargs))
    return control_file
        
 