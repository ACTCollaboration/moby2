from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits       # requires pyfits or astropy

import numpy as np

class TODModeSet:
    """
    The TODModeSet combines three pieces of information:
    - det_uid, a (n_det,) array.
    - weights, an (n_det,n_modes) array.
    - modes, an (n_modes,n_samp) array.
    """
    def __init__(self, det_uid, shape=None, dtype=None):
        self.det_uid = det_uid
        if shape is not None:
            if len(shape) != 2:
                raise ValueError('Expected shape=(n_modes, n_samp)')
            self.modes = np.zeros(shape, dtype)
            self.weights = np.zeros((len(self.det_uid), self.modes.shape[0]))

    @classmethod
    def from_fits_file(cls, filename):
        def extract_table(sdb, keyfmt, dtype=None):
            count = 0
            while True:
                if (keyfmt % count) not in sdb.dtype.names:
                    break
                count += 1
            if dtype is None:
                dtype = sdb[keyfmt % 0].dtype
            output = np.zeros((count, len(sdb)), dtype)
            for i in range(count):
                output[i,:] = sdb[keyfmt%i]
            return output
        data1 = moby2.util.StructDB.from_fits_table(filename, index=1)
        data2 = moby2.util.StructDB.from_fits_table(filename, index=2)
        self = cls(det_uid=data1['det_uid'])
        self.weights = extract_table(data1, 'weight%i').transpose()
        self.modes = extract_table(data2, 'mode%i')
        return self

    def to_fits_file(self, filename=None):
        prihdr = fits.Header()
        n_modes, n_samp = self.modes.shape
        prihdr['n_modes'] = n_modes
        prihdu = fits.PrimaryHDU(header=prihdr)
        tb0 = moby2.util.StructDB.from_data(
            [('det_uid', self.det_uid)] + [
             ('weight%i'%i, self.weights[:,i]) for i in range(n_modes)]
            ).to_fits_table()
        tb1 = moby2.util.StructDB.from_data(
            [('mode%i'%i, self.modes[i]) for i in range(n_modes)]
            ).to_fits_table()
        hdulist = fits.HDUList([prihdu, tb0, tb1])
        if filename is not None:
            hdulist.writeto(filename, clobber=True)
        return hdulist

    @classmethod
    def from_hdf(cls, target):
        cls.check_class(target, 'tod_modeset', 1)
        self = cls(det_uid=target['det_uid'])
        self.weights = np.array(target['weights'])
        self.modes = np.array(target['modes'])
        return self

    def to_hdf(self, target):
        kw = {'compression': 'gzip'}
        target.create_dataset('det_uid', data=self.det_uid.astype('uint32'), **kw)
        target.create_dataset('weights', data=self.weights.astype('float32'), **kw)
        target.create_dataset('modes', data=self.modes.astype('float32'), **kw)
        cls.set_class(target, 'tod_modeset', 1)

    def get_tod(self, dets=None, dtype=None, mode_idx=None):
        """
        Return weights dot modes for the desired dets.
        """
        if dets is None:
            dets = list(range(0, self.weights.shape[0]))
        if mode_idx is None:
            mode_idx = list(range(0, len(self.modes)))
        if np.asarray(dets).ndim == 0:
            return np.dot(self.weights[dets,mode_idx], self.modes[mode_idx])
        output = np.empty((len(dets), len(self.modes[0])), dtype=dtype)
        for j,i in enumerate(dets):
            output[j,:] = np.dot(self.weights[i,mode_idx], self.modes[mode_idx])
        return output

    def remove_modes(self, target, dets=None):
        if dets is None:
            dets = range(0, self.weights.shape[0])
        amps = np.array(np.transpose(self.weights), order='C')
        if self.modes.dtype == np.float64:
            moby2.libactpol.remove_modes64(
                target, np.array(dets).astype('int32'), self.modes, amps)
        elif self.modes.dtype == np.float32:
            moby2.libactpol.remove_modes(
                target, np.array(dets).astype('int32'), self.modes, amps)
        else:
            raise ValueError('Fast mode removal only supported for '
                             'self.modes.dtype float32 and float64.')


class TODModeSetArchive(moby2.util.HDFArchive):
    _moby2_class_name = 'tod_flags_archive'
    _moby2_class_version = 0
    _moby2_read_func = TODModeSet.from_hdf
    _moby2_write_func = TODModeSet.to_hdf


class PCA:
    @classmethod
    def for_data(cls, data):
        """Compute the covariance matrix of data array, and diagonalize it.
        Return PCA object containing this information (in attributes
        cov, E, R).
        """
        # Compute covariance matrix.
        cov = np.cov(data)
        # Diagonalize it.  eigh returns eigenvalues and matrix such
        # that   cov =    R . diag(E) . R^T
        E, R = np.linalg.eig(cov)  # eigh nans sometimes...
        self = cls()
        self.cov, self.E, self.R = cov, E, R
        # Sort by eigenvalue.
        self.E[np.isnan(self.E)] = 0.
        idx = np.argsort(-self.E)
        self.E = self.E[idx]
        self.R = self.R[:,idx]
        return self
    
    def get_decomp(self, data, n_modes=None):
        """Return a TODModeSet based on the strongest n_modes eigenmodes of
        this PCA.  Note that the input data can have a different shape
        than the original data used to compute the covariance, as long
        as it has the same number of detectors.  The returned
        TODModeSet contains the modes and their projection weights
        onto the original input vectors.
        """
        if n_modes is None:
            n_modes = len(self.E)
        modes = TODModeSet(np.arange(data.shape[0]),
                           shape=(n_modes, data.shape[1]))
        R = self.R[:,:n_modes]
        e = self.E[:n_modes]**-.5
        modes.modes[:] = np.dot(R.transpose(), data)
        modes.weights[:] = R
        return modes

