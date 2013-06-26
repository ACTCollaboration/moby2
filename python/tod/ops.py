from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Routines for simple processing of time-stream data arrays (float32
arrays with dimensions (n_det, n_sample)).

This is a good place from which to call optimized / parallelized /
weird custom C language code in libactpol.

It's best if such functions accept a tod=.. parameter, but also work
if only a data array is passed in (through data=...).  Then callers
don't have to load up a TOD to use optimized routines.
"""

import numpy as np
from moby2 import libactpol
import moby2

def detrend_tod(tod=None, dets=None, data=None):
    """
    Subtract a line through the first and last datapoints of each tod,
    preserving the data mean.

    Instead of a tod, you can pass in a data array through data=... .
    
    Returns (y0, y1), which are vectors of the values removed from the
    first and last samples, respectively.
    """
    if data is None:
        data = tod.data
    one_d = data.ndim == 1
    if one_d:
        data.shape = (1,-1)
    if dets is None:
        dets = np.ones(data.shape[0], 'bool')
    if np.asarray(dets).dtype == 'bool':
        dets = np.nonzero(dets)[0]
    y0, y1 = data[dets,0], data[dets,-1]
    slopes = (y1-y0) / (data.shape[1]-1)
    y0, y1 = -(y1-y0)/2, (y1-y0)/2
    x = np.arange(data.shape[1]).astype('float')
    for di, _y0, _slope in zip(dets, y0, slopes):
        data[di] -= x*_slope + _y0
    if tod is not None:
        tod.abuses += [{'name':'detrendData', 'type':'modeRemoval'}]
    if one_d:
        data.shape = -1
        return np.array((y0[0], y1[0]))
    return np.array((y0, y1))

def retrend_tod(trends, tod=None, dets=None, data=None):
    """
    Put the trends back in.
    """
    if data is None:
        data = tod.data
    one_d = data.ndim == 1
    if one_d:
        data = data.reshape(1,-1)  # better view for us
        trends = np.array(trends).reshape(1,2)  # sanity check?
    if dets is None:
        dets = np.ones(data.shape[0], 'bool')
    if np.asarray(dets).dtype == 'bool':
        dets = np.nonzero(dets)[0]
    y0, y1 = trends
    slopes = (y1-y0) / (data.shape[1]-1)
    y0, y1 = -(y1-y0)/2, (y1-y0)/2
    x = np.arange(data.shape[1]).astype('float')
    for di, _y0, _slope in zip(dets, y0, slopes):
        data[di] += x*_slope + _y0
    if tod is not None:
        tod.abuses += [{'name':'detrendData', 'type':'modeRemoval'}]

def remove_mean(tod=None, dets=None, data=None):
    """
    Perform
        data[dets] -= data[dets].mean(axis=1)

    If data is not provided, tod.data is used.
    """
    if data is None:
        data = tod.data
    if dets is None:
        dets = np.arange(data.shape[0])
    return libactpol.remove_mean(data, np.asarray(dets, 'int64'))

def remove_median(tod=None, dets=None, data=None):
    """
    Perform
        data[dets] -= data[dets].mean(axis=1)

    If data is not provided, tod.data is used.
    """
    if data is None:
        data = tod.data
    if dets is None:
        dets = np.arange(data.shape[0])
    data[dets] -= np.median(data[dets], axis=1)[:,np.newaxis]

def data_mean_axis0(data=None, dets=None):
    """
    Compute mean along axis 0.  To get a common mode or whatever.
    """
    if data is None:
        data = tod.data
    if dets is None:
        dets = np.arange(data.shape[0])
    dets = np.asarray(dets, 'int32', order='C')
    return libactpol.data_mean_axis0(data, dets)

def apply_calibration(data, dets, amps):
    """
    * data is float32 array with shape [n_det, n_samp]
    * dets is int* array/list of length n
    * amps is float* array/list with shape [n_samp]

    Perform
         data[dets,:] *= amps[dets][None,:]
    quickly.
    """
    dets = np.asarray(dets, 'int32', order='C')
    amps = np.asarray(amps, 'float32', order='C')
    if len(dets) != len(amps):
        raise ValueError("target dets not same length as amps.")
    return libactpol.apply_calibration(data, dets, amps)

def data_dot_modes(data, dets, modes, cuts=None):
    """
    Gets
         amps = dot(data[dets,d], modes)
    quickly.

    Result has size [len(dets), n_mode].
    """
    dets = np.asarray(dets, 'int32', order='C')
    modes = np.asarray(modes, 'float32', order='C')
    if cuts is not None:
        cuts = cuts._encode_c()
    return libactpol.data_dot_modes(data, dets, modes, cuts)

def remove_modes(data, dets, modes, amps):
    """
    Remove modes from a data array, with a different amplitude
    specified for each mode and detector.

    * data is float32 array with shape [n_det, n_samp]
    * dets is int* array/list of length n
    * modes is float* array with shape [n_mode,n_samp]
    * amps is float* array with shape [n_mode,n_det]

    Like:
       data[dets,:] -= (amps[:,None,:] * modes[:,:,None]).sum(axis=0)

    but fast.
    """
    dets = np.asarray(dets, 'int32', order='C')
    modes = np.asarray(modes, 'float32', order='C')
    amps = np.asarray(amps, 'float64', order='C')
    return libactpol.remove_modes(data, dets, modes, amps)

def remove_filter_gain(tod):
    """
    Apply the readout filter gain
    """
    cal_params = {'type': 'remove_readout_filter_gain'}
    cal_filtgain = moby2.scripting.get_calibration(cal_params, tod=tod)
    tod.data[cal_filtgain.det_uid] *= cal_filtgain.cal[:,np.newaxis] 
    return tod.info.mce_filter.gain()
