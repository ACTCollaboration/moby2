from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np
import os

import moby2
from moby2.util import MobyDict

class _SimpleDetData(object):
    """
    Track some per-detector data.  Read and write it like an ACTDict.
    Record indices of referenced detectors in self.det_uid.
    """
    det_uid = None
    fields = []
    source_file = None
    columns_format_str = None
    header = None

    def __init__(self, det_uid=None):
        if det_uid is not None:
            det_uid = np.array(det_uid)
        for f in self.fields:
            setattr(self, f, None)
        self.det_uid = det_uid

    def copy(self, idx=None):
        """
        Return a copy of the data structure.
        """
        out = self.__class__()
        out.source_file = self.source_file
        out.columns_format_str = self.columns_format_str
        out.fields = [f for f in self.fields]
        if idx is None:
            out.det_uid = self.det_uid.copy()
            for f in out.fields:
                setattr(out, f, getattr(self, f).copy())
        else:
            out.det_uid = self.det_uid[idx].copy()
            for f in out.fields:
                setattr(out, f, getattr(self, f)[idx].copy())
        return out

    def write(self, filename):
        params = MobyDict()
        for f in self.fields + ['det_uid']:
            params[f] = getattr(self, f)
            if isinstance(params[f], np.ndarray):
                params[f] = list(params[f])
        return params.write_to_file(filename)

    def write_columns(self, filename, format_str=None, header=None):
        if header is None:
            header = self.header
        if header is None:
            header = ''
        if len(header) > 0 and header[-1] != '\n':
            header = header + '\n'
        if format_str is None:
            format_str = self.columns_format_str
        if isinstance(filename, basestring):
            fout = open(filename, 'w')
        else:
            fout = filename
        fout.write(header)
        for idx in range(len(self.det_uid)):
            fout.write(format_str.format(**self.get_row(idx)))

    def read(self, filename):
        # This probably should not be the default format.  If you want to
        # use the dict-style storage, load it using classmethod
        # .from_dict(...) for forwards compatibility.
        params = MobyDict.from_file(filename)
        if not 'det_uid' in params:
            # Translate from ACT-style
            if 'row' in params:
                params['det_uid'] = [r*32+c for r,c in zip(params['row'], params['col'])]
            elif 'rows' in params:
                params['det_uid'] = [r*32+c for r,c in zip(params['rows'], params['cols'])]
            else:
                raise RuntimeError('Failed to retrieve det_uid from file %s' % filename)
        for f in self.fields + ['det_uid']:
            # Convert to array on read.
            if f in params:
                if hasattr(params[f], '__len__') and \
                        not isinstance(params[f], basestring):
                    params[f] = np.array(params[f])
                setattr(self, f, params[f])
        self.source_file = filename

    @classmethod
    def from_columns_file(cls, filename, columns=None, fields=None,
                          casts={}):
        if fields is None:
            fields = ['det_uid'] + list(self.fields)
        if columns is None:
            columns = list(range(len(fields)))
        self = cls()
        casts_ = {'det_uid': int}
        casts_.update(casts)
        data = []
        for line in open(filename):
            w = line.split()
            if len(w) == 0 or w[0][0] == '#':
                continue
            data.append([casts_.get(f,float)(w[c]) for f,c in 
                         zip(fields, columns)])
        # Transpose and make vectors
        data = [np.array(x) for x in zip(*data)]
        for f, d in zip(fields, data):
            setattr(self, f, d)
        return self

    @classmethod
    def from_fits_table(cls, filename, index=1, boolify=None):
        fits_data = moby2.util.StructDB.from_fits_table(
            filename, index=index)
        self = cls()
        for name in fits_data.dtype.names:
            d = fits_data[name]
            if boolify is not None and name in boolify:
                d = (d != 0)
            setattr(self, name, d)
        return self

    @classmethod
    def from_dict(cls, dict_file):
        self = cls()
        self.read(dict_file)
        return self

    def get_row(self, idx=None, det_uid=None):
        """
        Return a dict of values for a single detector.
        """
        if idx is None:
            try:
                idx = list(self.det_uid).index(idx)
            except:
                return None
        output = {'det_uid': self.det_uid[idx]}
        for f in self.fields:
            output[f] = getattr(self, f)[idx]
            # Some numpy-derived scalars defy .format; naturalize them.
            if isinstance(output[f], np.bool_):
                output[f] = bool(output[f])
        return output

    def get_index(self, det_uid):
        """
        For each det_uid (which should be a list or 1-d array),
        determine index into the data sets contained in this object,
        based on self.det_uid.  Returns (mask, indices), each of which
        is the same length as det_uid.  mask is a boolean showing
        which det_uid were matched, indices are the corresponding
        index into self.<whatever>.  Non-matched indices are returned
        as 0.
        """
        _det_uid = list(self.det_uid)
        idx = np.zeros(len(det_uid), 'int')
        mask = np.zeros(len(det_uid), 'bool')
        for i,d in enumerate(det_uid):
            try:
                idx[i] = _det_uid.index(d)
            except ValueError:
                continue
            mask[i] = True
        return mask, idx

    def get_property(self, prop, det_uid=None, default=0):
        if det_uid is None: det_uid = self.det_uid
        if isinstance(prop, basestring):
            mask, out = self.get_property([prop], det_uid, default=default)
            return mask, out[0]
        if np.isscalar(det_uid):
            mask, out = self.get_property(prop, [det_uid], default=default)
            return mask[0], [o[0] for o in out]
        mask, idx = self.get_index(det_uid)
        out = []
        for p in prop:
            if hasattr(self, p):
                out.append(getattr(self, p)[idx])
                out[-1][~mask] = default
            else:
                raise ValueError("object does not have property %s" % p)
        return mask, out

    def set_property(self, prop, value):
        if isinstance(prop, basestring):
            prop = [prop]
            value = [value]
        for p,v in zip(prop,value):
            if len(v) != len(self.det_uid):
                raise ValueError("property %s has wrong number of elements" % p)
            if hasattr(self, p):
                setattr(self, p, v)
            else:
                raise ValueError("object does not have property %s" % p)
            
class DetectorList(_SimpleDetData):
    pass

class DarkCandidates(_SimpleDetData):
    fields = ['dark']

    @classmethod
    def from_dict(cls, filename):
        dark = _SimpleDetData.from_dict(filename)  # fields=[]
        self = cls()
        self.det_uid = dark.det_uid
        if hasattr(dark, 'dark'):
            self.dark = dark.dark
        else:
            self.dark = np.ones(self.det_uid.shape, 'bool')
        self.source_file = dark.source_file
        return self

    @classmethod
    def for_tod(cls, tod):
        if tod.info.array.startswith('mbac'):
            from moby2.instruments import mbac
            return cls.from_dict(mbac.get_ACT_dark_candidate_file(
                    tod.info.season, tod.info.array))
        # Implement ACT-pol...
        raise

class LiveCandidates(_SimpleDetData):
    fields = ['live']

    @classmethod
    def from_ACT_file(cls, filename):
        # Columns det, row, col.
        self = cls()
        self.det_uid = np.loadtxt(filename, unpack=1)[0]
        self.live = np.ones(self.det_uid.shape, 'bool')
        return self

class TimeConstants(_SimpleDetData):
    fields = ['tau']

    @classmethod
    def from_ACT_file(cls, filename, f3dB):
        # ACT-style time constants come in as a n_rows x n_cols array
        taus = np.loadtxt(filename)
        # But sometimes they are actually f3dB values...
        if f3dB:
            taus[taus!=0] = 1./(2 * np.pi * taus[taus!=0])
        # And usually 0 means they aren't valid.
        mask = (taus!=0)
        rows, cols = mask.nonzero()
        det_uid = np.array([r*32 + c for r,c in zip(rows, cols)])
        self = cls()
        self.det_uid = det_uid
        self.tau = taus[mask]
        self.source_file = filename
        return self

    # Depot support
    _depot_structure = '{class}/{tag}/{first_five}/{tod_name}.tau'

    def write_to_path(self, path, *args, **kwargs):
        return self.write(path)

    @staticmethod
    def read_from_path(path, *args, **kwargs):
        det_uid, tau = moby2.util.ascii.read_columns(path, [0,1])
        ok = (tau>0)
        self = TimeConstants(det_uid=det_uid[ok])
        self.tau = tau[ok]
        return self

class PolarizationAngles(_SimpleDetData):
    """
    Polarization angles by detector.  Attributes:

       pol_angle - polarization angle in radians, measured CCW from West
                   ("right" as projected onto the sky).

       det_uid   - corresponding detector UIDs 

    Since this is a SimpleDetData, all listed data are implicitly
    valid.
    """
    fields = ['pol_angle']

    def write_to_path(self, path):
        fout = open(path, 'w')
        for d, phi in zip(self.det_uid, self.pol_angle):
            fout.write('%4i %8.2f\n' % (d, f, phi * 180/np.pi))
        del fout

    # Depot support
    _depot_structure = 'PolAngles/{tag}.txt'

    def write_to_path(self, path, *args, **kwargs):
        return self.write(path)

    @staticmethod
    def read_from_path(path, format=None, columns=None, units=None,
                       convention=None, bad_pol_value=None):
        if format is None:
            # Same as write format
            det_uid, ok, phi = moby2.util.ascii.read_columns(path)
        elif format == 'ascii_with_mask':
            if columns is None:
                columns = [0,1,2]
            det_uid, phi, ok = moby2.util.ascii.read_columns(
                path, columns=columns)
            if bad_pol_value is None:
                # Unless the meaning of the flag is passed, assume 0 is bad
                ok = ok.astype('bool')
            else:
                ok = (ok != bad_pol_value)
        elif format == 'ascii':
            if columns is None:
                columns = [0,1]
            det_uid, phi = moby2.util.ascii.read_columns(
                path, columns=columns)
            if not bad_pol_value is None:
                ok = (phi != bad_pol_value)
            else:
                ok = np.ones(phi.shape)
        else:
            raise ValueError("uknown file format %s" % format)
            
        # Rescaling or reference change?
        if units in [None, 'degrees']:
            phi *= np.pi/180
        elif units == 'radians':
            pass
        else:
            raise ValueError("uknown units %s" % units)
        # Convention?
        if convention in [None, 'phi', 'convention0']:
            pass
        elif convention in ['gamma', 'convention1']:
            phi = np.pi/2 - phi
        else:
            raise ValueError("uknown angle convention %s" % convention)
        # Return
        self = PolarizationAngles(det_uid=det_uid[ok])
        self.pol_angle = phi[ok]
        return self

class Calibration(_SimpleDetData):
    fields = ['cal', 'calRMS']

    def __init__(self, det_uid=None):
        super(Calibration, self).__init__(det_uid=det_uid)
        if det_uid is not None:
            self.cal = np.ones(self.det_uid.shape)
            self.calRMS = np.zeros(self.det_uid.shape)

    # Depot support
    _depot_structure = '{class}/{tag}/{first_five}/{tod_name}.cal'
    def write_to_path(self, path, *args, **kwargs):
        return self.write(path)

    @staticmethod
    def read_from_path(path, *args, **kwargs):
        return Calibration.from_dict(path)

    def apply(self, tod=None, data=None, det_uid=None,
              zero_invalid=True):
        """
        Apply calibration to data, using multiplication.  User should
        pass in:
            data:     a float32 array with shape (ndet, nsamp)
            det_uid:  an integer array with shape (ndet)
        or
            tod:      a TOD object, with .data and .det_uid attributes
                      satisfying the above requirements.

        This function returns (mask, cal), where mask is a boolean
        array indicating which detectors had valid calibration numbers
        in this object, and cal is a float array of calibration
        numbers.  Both arrays are of size (ndet).

        If zero_invalid is True, then detectors with invalid
        calibration will be scaled by 0.  Otherwise, they will be left
        unscaled.
        """
        # Valid argument combos?
        if ((tod is None, data is None, det_uid is None) not in [
                (True, False, False),
                (False, True, True)]):
            raise ValueError("Specify data=...,det_uid=... OR tod=... ")
        if tod is not None:
            data, det_uid = tod.data, tod.det_uid
        mask, cal = self.get_property('cal', det_uid=det_uid)
        idx = np.arange(len(det_uid))
        if zero_invalid:
            # Set invalid dets to 0.
            cal[~mask] = 0.
        else:
            # Restrict to only detectors with valid calibration.
            idx = idx[mask]
        # Apply the calibration.  Quickly.
        moby2.libactpol.apply_calibration(data, idx.astype('int32'),
                                          cal[idx].astype('float32'))
        return mask, cal

    @classmethod
    def from_hdf(cls, target):
        # In our storage scheme, the target is a dataset...
        data = np.array(target)
        
        assert('det_uid' in data.dtype.names)
        self = cls(data['det_uid'])
        for k in self.fields:
            if k in data.dtype.names:
                setattr(self, k, data[k])
        return self

    def write_hdf(self, target):
        data = moby2.util.StructDB.from_data([
            ('det_uid', self.det_uid),
            ('cal', self.cal),
            ('calRMS', self.calRMS)])
        # Replace target (a group) with our dataset.
        name = str(target.name)
        hfile, parent = target.file, target.parent
        del hfile[name]
        data.to_hdf(hfile, dataset=name, compression='gzip')
        return hfile[name]

class RelCal(Calibration):
    fields = ['cal', "calRMS", "stable"]
        
class FlatField(RelCal):
    pass

class CalGC(RelCal):
    pass

class CalibrationArchive(moby2.util.HDFArchive):
    _moby2_class_name = 'tod_cal_archive'
    _moby2_class_version = 0
    _moby2_read_func = Calibration.from_hdf
    _moby2_write_func = Calibration.write_hdf
