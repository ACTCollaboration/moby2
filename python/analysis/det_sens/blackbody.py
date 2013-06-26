from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Low-level blackbody and Rayleigh-Jeans spectra.

In all of this, T is a temperature in Kelvin and f is a frequency in GHz.
"""

import numpy as np
import moby2.util.constants as const

c = const.SPEED_OF_LIGHT
k_pJ = const.BOLTZMANN_pJ
h_pJ = const.PLANCK_pJ
T_cmb = const.CMB_TEMPERATURE
pi = np.pi

GHz = 1e9

def blackbody_x(T, f):
    """
    Return dimensionless parameter
         x = h f / k T
    """
    return h_pJ / k_pJ * (f*GHz) / T

def rayleigh(T, f):
    """
    RJ spectral radiance, in (pW / m**2 / s / steradian) / (GHz).
    """
    return 2. * (f*GHz)**2 * k_pJ * T / c**2 * GHz

def blackbody(T, f):
    """
    Blackbody spectral radiance, in (pW / m**2 / s / steradian) / (GHz).
    """
    x = blackbody_x(T, f)
    return rayleigh(T, f) * x / (np.exp(x) - 1)

def spectrumToBrightness(S, f):
    """
    Given spectral radiance at f, returns the corresponding blackbody temperature.
    """
    x = np.log( 1. + ((c**2 / 2 / h_pJ / (f*GHz)**3) * S / GHz)**-1 )
    return h_pJ * (f*GHz) / k_pJ / x

def spectrumToRJ(S, f):
    """
    Given spectral radiance at f, returns the corresponding RJ temperature.
    """
    return 0.5 / (f*GHz/c)**2 / k_pJ * S / GHz

def blackbodyToRJ(T, f):
    return spectrumToRJ(blackbody(T, f), f)

def RJToBlackbody(T, f):
    return spectrumToBrightness(rayleigh(T, f), f)

def DCMB_factor(f):
    """
    Returns differential CMB temperature linearization factor X such
    that spectral radiance I is
      I = X * dT
    """
    x = blackbody_x(T_cmb, f)
    return GHz * k_pJ * \
           2 * (f*GHz / c)**2 * (x/2 / np.sinh(x/2))**2

def DCMB(dT, f):
    """
    Returns spectral radiance associated with source with temperature
    dT measured in differential CMB blackbody units.
    """
    return dT * DCMB_factor(f)

def spectrumToDCMB(S, f):
    """
    Given spectral radiance, return dT_CMB.
    """
    return S / DCMB_factor(f)
    
def RJToDCMB(f):
    """
    Returns the factor that converts RJ temperature to differential
    CMB temperature, at frequency f.
    """
    return spectrumToDCMB(rayleigh(1, f), f)
    
if __name__ == '__main__':
    print('Check inversions...')
    T = 30.
    f = 148.
    S = blackbody(T, f)
    print(T, '-->', spectrumToBrightness(S, f))

    S = rayleigh(T, f)
    print(T, '-->', spectrumToRJ(S, f))
