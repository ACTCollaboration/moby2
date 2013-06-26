from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Some backwards compatibility stuff.
"""

import moby2

import numpy as np
import os

mcfg = moby2.util.get_user_config()

MOBY = mcfg.get_deep(('ACT_configs', 'MOBY_DIR'), '')
MOBYSHARED = mcfg.get_deep(('ACT_configs', 'MOBYSHARED_DIR'), '')

if MOBYSHARED == '':
    moby2.util.log.trace('moby', 1, 'MOBYSHARED_DIR not set; some ACT products not available')
if MOBY == '':
    moby2.util.log.trace('moby', 1, 'MOBY_DIR not set; some ACT products not available')


def get_ACT_shunt_R(array_name):
    """
    Expects array_name in ['mbac145',...]
    """
    AR = AR_NAMES[array_name].upper()  # e.g. ar1
    filename = os.path.join(MOBYSHARED, 'detectors', '2008',
                            '%s_shunt_data.txt' % AR)
    R = np.loadtxt(filename)
    R[np.isnan(R)] = 0.
    baddies = R <= 0.
    if AR=='AR2':
        R[baddies]=0.000797
    elif AR=='AR1':
        R[baddies]=0.0007601722
    elif AR=='AR3':
        col_mask = np.ones(baddies.shape, 'bool')
        col_mask[:,24:] = False
        R[baddies*~col_mask] = 0.000769
        R[baddies* col_mask] = 0.002063
    else:
        raise
    R[R<=0] = 0.
    return R

    


##
## Data :P
##

AR_NAMES = {
            'mbac145': 'AR1',
            'mbac215': 'AR2',
            'mbac280': 'AR3',
            'ar1': 'actpol1',
            'ar2': 'actpol2',
}

