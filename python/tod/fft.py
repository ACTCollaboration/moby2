from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

def fft(data, dets=None, all_freqs=True):
    if dets is None:
        dets = np.arange(data.shape[0])
    return moby2.libactpol.fft_f_r(data, dets, all_freqs)

def ifft(fdata, dets=None, nsamps=None):
    if dets is None:
        dets = np.arange(fdata.shape[0])
    if nsamps is None:
        nsamps = fdata.shape[1]
    return moby2.libactpol.ifft_f_r(fdata, dets, nsamps)
