from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

from moby2.tod import TODCuts, CutsVector
from moby2.util import HDFArchive

class TODFlags:
    """
    Manages the association of a set of flags to each sample of each
    detector in a TOD.  Helps with I/O of such information, and the
    conversion of different combinations of flags to TODCuts (which
    are effectively are carry a single flag bit).
    """

    def __init__(self, ndets=None, nsamps=None, det_uid=None,
                 sample_offset=None, flag_names=None):
        """To initialize the class, pass at least the det_uid vector or ndets
        (from which a simple det_uid vector will be generated).
        Parameters nsamps and sample_offset have the same meaning as
        in TODCuts.  flag_names can be used to initialize a set of
        base flags.
        """
        if ndets is None:
            ndets = len(det_uid)
        if det_uid is None:
            det_uid = np.arange(ndets)
        self.det_uid = det_uid
        self.nsamps = nsamps
        self.sample_offset = sample_offset
        self.flag_names = []
        self.base_flags = []
        self.flag_defs = {}
        if flag_names is not None:
            for n in flag_names:
                self.add_base_flag(n)

    @classmethod
    def for_cuts(cls, cuts):
        """Returns an empty TODFlags instance based on the TODCuts instance
        passed in.  This just sets det_uid, sample_offset, nsamps.
        The cuts themselves are not copied in.
        """
        return cls(nsamps=cuts.nsamps, det_uid=cuts.det_uid,
                   sample_offset=cuts.sample_offset)

    def add_base_flag(self, name, data=None, make_copy=True):
        """Create a new base_flag, with the given name.  If data is passed in,
        it should be a TODCuts object, and it will be used to
        initialize the flag.
        """
        assert name not in self.flag_names
        self.flag_defs[name] = ('base', len(self.flag_names))
        self.flag_names.append(name)

        if not make_copy:
            if data is None:
                make_copy = True
            else:
                # Check that cuts are compatible...
                if data.nsamps != self.nsamps or data.sample_offset != self.sample_offset:
                    print(' ...overriding no-copy request ("%s").' % name)
                    make_copy = True
        if make_copy:
            self.base_flags.append(TODCuts(
                det_uid=self.det_uid, nsamps=self.nsamps,
                sample_offset=self.sample_offset))
            if data is not None:
                self.base_flags[-1].merge_tod_cuts(data)
        else:
            self.base_flags.append(data)
        
    def add_derived_flag(self, name, cut_components, uncut_components=[]):
        """Create a new derived flag, called `name`.  The flag is derived by
        or-ing together the flags listed in cut_components, """
        assert name not in self.flag_names
        comp_i0 = [self.flag_names.index(c) for c in cut_components]
        comp_i1 = [self.flag_names.index(c) for c in uncut_components]
        self.flag_defs[name] = ('derived', comp_i0, comp_i1)
        self.flag_names.append(name)

    def get_cuts(self, name):
        if not name in self.flag_defs:
            raise ValueError("Flag '%s' not found in flag_defs: %s" % (
                name, list(self.flag_defs.keys())))
        flag_def = self.flag_defs[name]
        if flag_def[0] == 'base':
            return self.base_flags[self.flag_names.index(name)].copy()
        elif flag_def[0] == 'derived':
            i_pos, i_neg = flag_def[1:]
            cuts = TODCuts(det_uid=self.det_uid, nsamps=self.nsamps,
                           sample_offset=self.sample_offset)
            for i in i_pos:
                cuts.merge_tod_cuts(self.base_flags[i])
            for i in i_neg:
                cuts.merge_tod_cuts(self.base_flags[i].get_complement())
            return cuts

    """
    I/O.
    """

    @staticmethod
    def _explode_bits(bits, total_bits):
        """Create an array of uint8 of sufficient size to express `total_bits`
        bits.  Set bits to 1 according to the list of bits in `bits`.
        For example, if bits=[0,1,5,15] and total_bits=24, the
        function returns array([35,128,0]).
        """
        out = np.zeros((total_bits + 7) // 8, 'uint8')
        for b in bits:
            out[b >> 8] |= (1 << (b % 8))
        return out

    def _prepare_output(self):
        """Service function for output routines; recode internal state for
        easy storage.
        """
        # Use C routine to convert raw cuts vectors to packed form.
        X = moby2.libactpol.pack_flags(
            [x.cuts for x in self.base_flags], self.nsamps)
        stack_bounds, flag_stack, index_stack = X
        base_names = [HDFArchive.encode_23(n) for n in self.flag_names
                      if self.flag_defs[n][0] == 'base']
        output = [
            ('!sample_offset', self.sample_offset),
            ('!sample_count', self.nsamps),
            ('det_uid', self.det_uid),
            ('stack_bounds', stack_bounds),
            ('index_stack', index_stack),
            ('flag_stack', flag_stack),
            ('flag_names', np.array(base_names)),
            ]
        n_words = flag_stack.shape[1]
        dt, word_size = np.uint8, 8

        # Add on the derived field information.
        derived = []
        for n in self.flag_names:
            if self.flag_defs[n][0] != 'derived':
                continue
            flagval = np.zeros((2, n_words), dtype=dt)
            i_pos, i_neg = self.flag_defs[n][1:]
            for i in i_pos:
                flagval[0, i >> word_size] |= (1 << i)
            for i in i_neg:
                flagval[1, i >> word_size] |= (1 << i)
            derived.append((HDFArchive.encode_23(n), flagval))
        if len(derived) > 0:
            dnames, dmasks = list(map(np.array, list(zip(*derived))))
        else:
            dnames = np.zeros((0,), dtype='S0')
            dmasks = np.zeros((0,2,n_words), dtype=dt)  # if empty, give right shape.
        output.append(('derived_names', dnames))
        output.append(('derived_masks', dmasks))
        return output
        
    def write_hdf(self, target):
        """Encode this object into an HDF5 file.  The target must be an empty
        h5py HDF5 group.

        By design, this function does not support a string address
        (file + group name); use TODFlagsArchive for container
        management.
        """
        output = self._prepare_output()
        for k, v in output:
            if k[0] == '!':
                target.attrs[HDFArchive.encode_23(k[1:])] = v
            else:
                target.create_dataset(HDFArchive.encode_23(k), data=v,
                                      compression='gzip')
        HDFArchive.set_class(target, 'tod_flags', 1)

    @classmethod
    def from_hdf(cls, target):
        """Load an instance of this class from the HDF5 object pointed to by
        target.  This could be an open h5py.File, or a group within
        one.

        By design, this function does not support a string address
        (file + group name).  For such container management, use a
        TODFlagsArchive.
        """
        HDFArchive.check_class(target, 'tod_flags', 1)

        self = cls(det_uid=np.array(target['det_uid']),
                   nsamps=target.attrs['sample_count'],
                   sample_offset=target.attrs['sample_offset'])
        for name in target['flag_names']:
            name = HDFArchive.decode_23(name)  # py2/3
            self.add_base_flag(name)
        bit_masks = [self._explode_bits([b], len(self.base_flags))
                     for b in range(len(self.base_flags))]

        # Call C layer to get lists
        cvecs = moby2.libactpol.unpack_flags(
            np.array(target[b'stack_bounds'], 'uint32'),
            np.array(target[b'flag_stack'], 'uint8'),
            np.array(target[b'index_stack'], 'uint32'),
            target.attrs[b'sample_count'],
            len(target[b'flag_names']))

        for i in range(len(self.det_uid)):
            for b, m in enumerate(bit_masks):
                self.base_flags[b].cuts[i] = CutsVector(
                    cuts_in=cvecs[b][i], nsamps=self.nsamps)

        # The derived fields.
        for i,n in enumerate(target['derived_names']):
            pos = [HDFArchive.decode_23(_n) for _n,m in
                   zip(target[b'flag_names'], bit_masks)
                   if np.any(target[b'derived_masks'][i][0] & m)]
            neg = [HDFArchive.decode_23(_n) for _n,m in
                   zip(target[b'flag_names'], bit_masks)
                   if np.any(target[b'derived_masks'][i][1] & m)]
            self.add_derived_flag(HDFArchive.decode_23(n), pos, neg)

        return self
        
    def write_hdf_deprecated(self, hdf_file, group_name=None, clobber=False,
                  compression='gzip'):
        """Store the TODFlags flags object to the given hdf_file, which can be
        a filename, or a node in an open h5py.File instance.  Data are
        written to group `group_name`.  If group exists, function will
        fail unless clobber=True, in which case the existing group
        will be deleted.
        """
        import h5py

        open_file_in = not isinstance(hdf_file, basestring)
        if not open_file_in:
            hdf_file = h5py.File(hdf_file, 'a')
        if group_name is None:
            group_name = '/'
        if group_name in hdf_file:
            print('Found existing group %s' % group_name)
            group = hdf_file[group_name]
            if clobber:
                for k,v in output:
                    if k[0] == '!' and k[1:] in group.attrs:
                        del group.attrs[k[1:]]
                    elif k in group:
                        del group[k]
        else:
            group = hdf_file.create_group(group_name)

        self.write_hdf(group)

        if not open_file_in:
            # This is not reliable... safer to pass in an open h5py.File.
            hdf_file.close()

    @classmethod
    def from_hdf_deprecated(cls, hdf_file, group_name='/'):
        """Initialize a TODFlags flags object from the specified hdf_file
        filename.  Data are loaded from HDF5 group `group_name`.
        """
        import h5py
        open_file_in = not isinstance(hdf_file, basestring)
        if not open_file_in:
            hdf_file = h5py.File(hdf_file, 'r')
        group = hdf_file[group_name]

        self = cls.from_hdf(group)

        if not open_file_in:
            hdf_file.close()
        return self


class TODFlagsArchive(HDFArchive):
    _moby2_class_name = 'tod_flags_archive'
    _moby2_class_version = 0
    _moby2_read_func = TODFlags.from_hdf
    #_moby2_write_func = TODFlags.write_hdf

    def _moby2_write_func(self, item, dest):
        return item.write_hdf(dest)
