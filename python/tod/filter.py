from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

def apply_simple(tod_data, filt, detrend=False, retrend=False, dets=None):
    # Note sample time is not needed because time constants are not applied.
    one_dim = (tod_data.ndim == 1)
    if one_dim:
        assert(dets is None)
        tod_data = tod_data.reshape(1,-1)
        dets = np.array([0], dtype='int32')
    elif dets is None:
        dets = np.arange(tod_data.shape[0], dtype='int32')
    else:
        dets = np.asarray(dets, dtype='int32')
    moby2.libactpol.filter_tod_data(
        tod_data, dets, filt.astype('complex64'),
        None, 0., 0,
        int(detrend), int(retrend))
    if one_dim:
        return tod_data[0]
    return tod_data


def prefilter_tod(tod, time_constants=None,
                  deconvolve_readout=True,
                  detrend=False, retrend=False,
                  global_filt=None):
    """
    Deconvolve butterworth readout filter and detector time constants.
    Optionally apply some additional filter global_filt (which must be
    len(tod.data[-1])).
    """
    dets = tod.cuts.get_uncut()
    if global_filt is None:
        filt = np.ones(tod.data.shape[-1], 'complex128')
    else:
        filt = np.array(global_filt, dtype='complex128')
    if deconvolve_readout:
        filt /= MCEButterworth(tod)
    if time_constants is None:
        tau = np.zeros(len(dets))
    else:
        tau = np.asarray(time_constants)[dets]
    tau_do_inverse = True
    moby2.libactpol.filter_tod_data(
        tod.data,
        np.asarray(dets, dtype='int32'),
        filt.astype('complex64'),
        tau.astype('float32'),
        tod.ctime[1] - tod.ctime[0],
        int(tau_do_inverse), int(detrend), int(retrend))
    

def genFreqs( ndata, tSample ):
    """
    @brief Generate frequency vector of specified length and corresponding
    to a specified total time period.
    
    @param ndata    Number of elements in vector
    @param tSample  Sample time
    """
    dn = 2      # if you like the central frequency to be negative, change dn to 1
    return 1./(ndata*tSample) * np.hstack((np.arange(0., (ndata+dn)//2),
                                           np.arange(-(ndata+dn)//2+dn, 0)))


def MCEButterworth(tod):
    """
    @brief Return the butterworth filter evaluated at frequencies in the vector freq.
    """
    import moby2.util.mce as mce
    ntimes = len(tod.ctime)
    sampleTime = (tod.ctime[-1] - tod.ctime[0]) / (ntimes-1)
    f = genFreqs( ntimes, sampleTime )
    mce_filter = mce.MCEButterworth.from_runfile(tod.info.runfile)
    return mce_filter.transfer(f) / mce_filter.transfer(0)
