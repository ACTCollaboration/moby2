from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np


class HWPModes(object):
    """
    Store mode coefficients for a bunch of detectors.
    
    self.det_uid, of shape (n_det), says what detectors they are.

    self.coeff has shape (n_det, mode_max*2).  For each detector the modes are
    ordered

       (1c, 1s, 2c, 2s, ..., Nc, Ns)

    where, e.g., 2c means cos(2*chi).  And N is self.mode_max.  
    """
    det_uid = None
    coeffs = None
    mode_max = None

    def extract(self, det_uid):
        suid = moby2.util.StructDB.from_data([('x', self.det_uid)])
        idx = suid.select_inner({'x': det_uid})
        s = idx>=0
        output = HWPModes()
        output.coeffs = self.coeffs[idx,:]
        output.coeffs[~s,:] = 0.
        output.det_uid = np.array(det_uid)
        output.mode_max = self.mode_max
        return output, s

    @classmethod
    def zeros(cls, det_uid, mode_max):
        self = cls()
        self.det_uid = np.array(det_uid)
        self.mode_max = mode_max
        self.coeffs = np.zeros((len(self.det_uid), self.mode_max*2))
        return self

    @classmethod
    def from_coeffs(cls, det_uid, coeffs):
        self = cls()
        self.det_uid = np.array(det_uid)
        self.coeffs = np.array(coeffs)
        self.mode_max = self.coeffs.shape[1]//2
        return self

    """
    FITS read/write.
    """

    @classmethod
    def from_fits_table(cls, filename, index=1, field_prefix=''):
        # No this is not very efficient.
        indata = moby2.util.StructDB.from_fits_table(filename, index=1)
        max_mode = 0
        coeffs = []
        while True:
            field = '%smode_%d' % (field_prefix, max_mode+1)
            if field+'_c' not in indata.dtype.names:
                break
            coeffs.append(indata[field+'_c'])
            coeffs.append(indata[field+'_s'])
            max_mode += 1
        self = cls()
        self.det_uid = indata['det_uid']
        self.mode_max = max_mode
        self.coeffs = np.transpose(coeffs)
        return self

    def to_fits_table(self, filename, **kwargs):
        data = [('det_uid', self.det_uid)]
        for i in range(self.mode_max):
            data.extend([('mode_%d_c' % (i+1), self.coeffs[:,i*2  ]),
                         ('mode_%d_s' % (i+1), self.coeffs[:,i*2+1])])
        return moby2.util.StructDB.from_data(data).to_fits_table(filename, **kwargs)

    """
    HDF read/write.
    """

    @classmethod
    def from_hdf(cls, target):
        det_uid, coeffs = [np.array(target[k]) for k in ['det_uid', 'coeffs']]
        self = cls()
        self.det_uid = det_uid
        self.coeffs = coeffs
        self.mode_max = self.coeffs.shape[1]//2
        return self

    def to_hdf(self, target):
        target.create_dataset('det_uid', data=self.det_uid.astype('int32'),
                              compression='gzip')
        target.create_dataset('coeffs', data=self.coeffs.astype('float32'),
                              compression='gzip')

    def get_reconstructor(self, chi):
        """
        Returns an object with sine and cosine vectors initialized so
        that individual detector A(chi) can be quickly constructed.
        If you want to jump straight to the answer, do

          achi = hwp_modes.get_projector(chi).get()
        """
        return HWPReconstructor(self, chi)

    # depot support.
    _depot_structure = '{class}/{tag}/{tod_name}_achi.fits'
    write_to_path = to_fits_table
    read_from_path = from_fits_table


class HWPModesArchive(moby2.util.HDFArchive):
    _moby2_write_func = HWPModes.to_hdf
    _moby2_read_func  = HWPModes.from_hdf


class HWPReconstructor:
    def __init__(self, modes, chi):
        self.hwp_modes = modes
        # Create cos and sin vectors. Shape (n_modes*2, n_chi).
        self.vects = np.empty((modes.mode_max*2, len(chi)))
        for k in range(modes.mode_max):
            self.vects[2*k  ] = np.cos(chi*(k+1))
            self.vects[2*k+1] = np.sin(chi*(k+1))

    def get_achi(self, dets=None, mask=None):
        if dets is None:
            dets = np.arange(len(self.hwp_modes.det_uid))
        if np.asarray(dets).ndim == 0:
            return self.get_achi([dets], mask=mask)[0]
        # Now you can assume that dets is an array.
        # coeffs has shape (n_dets, n_modes*2)
        if mask is None:
            return np.dot(self.hwp_modes.coeffs[dets], self.vects)
        else:
            return np.dot(self.hwp_modes.coeffs[dets], self.vects[:,mask])

class HWPProjector:
    def __init__(self, chi, n_modes):
        # Create cos and sin vectors. Shape (n_modes*2, n_chi).
        self.vects = np.empty((n_modes*2, len(chi)))
        for k in range(n_modes):
            self.vects[2*k  ] = np.cos(chi*(k+1))
            self.vects[2*k+1] = np.sin(chi*(k+1))
        self.M_inv = None

    def recompute_M_inv(self):
        self.M_inv = np.linalg.inv(np.tensordot(self.vects, self.vects, (1,1)))

    def get_coeffs_fast(self, data, nonzero_fraction=None):
        """
        This is fast because the mode-coupling matrix is not
        recomputed for each detectors individual mask; the solution is
        thus only approximate.  Rather, it's assumed that the user set
        data[i][~mask] = 0 and passed us nonzero_fraction (an array of
        shape (n_det), with values between 0 and 1) indicating the
        fraction of samples that were not masked.  This approximation
        degrades as the nonzero_fraction gets smaller.
        """
        if self.M_inv is None:
            self.recompute_M_inv()
        # (n_det, 2*n_modes) <- (n_det, n_chi) . (2*n_modes, n_chi)
        proj = np.tensordot(data, self.vects, (1,1))
        coeffs = np.dot(proj, self.M_inv)
        if nonzero_fraction is not None:
            coeffs /= nonzero_fraction[:,None]
        return coeffs

    def get_coeffs(self, data, dets=None, mask=None):
        """
        This fits the modes to the detector data.  Inputs are:
          data:    size (n_det, n_samp)
          dets:    optional, 1d; indexes into data[].
          mask:    optional, shape (n_samp); selects samples in each data[i].

        Returns: mode coefficients as an (n_det, n_mode) array.

        If data is 1-d, it's assumed to be a single detector
        time-stream, and the function returns 1-d coefficients with
        shape (n_mode).

        This routine is slower than get_coeffs_fast because it does an
        exact least-squares fit for the modes.
        """
        collapse = (data.ndim == 1)
        if collapse:
            # inflate...
            data = data.reshape(1,-1)
        if dets is None:
            dets = np.arange(data.shape[0])
        # Project.
        coeffs = np.empty((len(dets), len(self.vects)))
        # Now you can assume that dets is an array.
        if mask is None:
            vects = self.vects
        else:
            vects = self.vects[:,mask]

        I = np.linalg.inv(np.tensordot(vects, vects, (1,1)))

        for k in range(len(dets)):
            d = data[k]
            if mask is not None:
                d = d[mask]
            coeffs[k,:] = np.dot(vects, d)
        # Finish solution...
        coeffs = np.dot(I, coeffs)
        if collapse:
            return coeffs.reshape(-1)
        return coeffs
