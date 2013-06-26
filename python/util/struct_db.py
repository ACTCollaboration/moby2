from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np
import os

from .ascii import read_columns
from . import moby_fits

class StructDB(np.ndarray):
    """
    Mini-database for storage of simple records.  For example,
    detector properties (that you want to look up by det_uid) or TOD
    properties (that you want to look up by TOD name).

    In addition to the normal interface of a numpy structured array,
    some database-like functions are supported, for finding single or
    multiple records based on selected keys.
    """

    def __new__(cls, shape, dtype=float, buffer=None, offset=0,
                strides=None, order=None):
        self = super(StructDB, cls).__new__(cls, shape, dtype, buffer,
                                            offset, strides, order)
        self.formats = {}
        return self

    def __array_finalize__(self, obj):
        if isinstance(obj, StructDB):
            self.formats = {}
            if hasattr(obj, 'formats'):
                for k in list(obj.formats.keys()):
                    if k in self.dtype.name:
                        self.formats[k] = obj.formats[k]

    def copy(self, order='C'):
        out = super(StructDB, self).copy(order=order)
        out.formats = self.formats.copy()
        return out

    @classmethod
    def from_data(cls, data, formats={}):
        """
        Construct a StructDB array from numpy arrays.  data can be one of:
        
        * A dictionary of key -> value associations
        * A list of tuples [(key1,value1), (key2,value2), ...]
        """
        if isinstance(data, dict):
            data = list(data.items())
        dtype = [(d[0],np.asarray(d[1]).dtype) for d in data]
        self = cls(shape=len(data[0][1]), dtype=dtype)
        for d in data:
            self[d[0]][:] = d[1]
            if len(d) > 2:
                self.formats[d[0]] = d[2]
        self.formats.update(formats)
        return self

    def get_row(self, i):
        return dict([(k,self[k][i]) for k in self.dtype.names])

    def index(self, props={}, fail=True):
        l_props = dict([(k,[v]) for k,v in list(props.items())])
        idx = self.select_inner(l_props)
        if idx[0] >= 0:
            return idx[0]
        if fail:
            raise ValueError("Could not match specification %s" % props)
        return -1

    """
    Database-like operations.
    """
    def select_outer(self, props={}, mask=False, fields_out=None):
        """
        Select a set of records based on outer product of the
        (key,values) pairs in the props dict.  "props" is a dictionary
        where each key is a field property and each value is a list of
        desired values for that property.  For example, if props is:

        {'row': [1,2,3], 'col': [2,4]}

        Then the function returns the (probably) 6 detectors that
        have 'row' property in [1,2,3] and 'col' property in [2,4].

        If a value is given as "None", it is ignored (i.e. no attempt
        is made to match that key).

        Returns the matching indices, or a boolean mask of size
        len(self) if mask==True.
        """
        # Create a mask to track the matches
        selection = np.ones(len(self), 'bool')
        # Reduce the mask with each property
        for k, values in list(props.items()):
            if values is None:
                continue
            if k in self.dtype.names:
                selection *= [v in values for v in self[k]]
            elif k == 'det_uid':
                selection[values] *= True
            else:
                print("Property %s not defined for array." % k)
                selection[:] = False
                break
        if fields_out is None:
            if mask:
                return selection
            return selection.nonzero()[0]
        if isinstance(fields_out, basestring):
            return self[fields_out][selection]
        else:
            return StructDB.from_data([(f, self[f][selection]) for f in fields_out])

    def select_inner(self, props, mask=False):
        """
        Select a set of records based on inner product of the
        (key,values) pairs in the props dict.  "props" is a dictionary
        where each key is a property and each value is a list of
        desired values for that property.  The value lists must all be
        the same length.  For example, if props is:

        {'row': [1,2,3], 'col': [2,4,5]}

        Then the function returns the index of 3 detectors with
        (row,col) in [(1,2),(2,4),(3,5)].

        If a value is given as "None", it is ignored (i.e. no attempt
        is made to match that key).

        When mask==False (the default), the array that is returned
        gives the det_uid of the detectors matching the requested
        properties, in the order that they were matched, with -1
        dropped in when a match could not be made.

        When mask==True, a boolean mask is returned that identifies
        the position of successful matches; the ordering is lost, and
        match failures are ignored.
        """
        # Get tuples for items we are trying to find
        keys = [k for k,v in list(props.items()) if v is not None]
        props = list(zip(*[props[k] for k in props]))
        # Init matched indices to -1.
        matched = np.zeros(len(props), int)-1
        for k in keys:
            if not k in self.dtype.names:
                print("Property %s not defined for array." % k)
                return matched
        # Find each.
        for ip, p in enumerate(props):
            s = np.ones(len(self), 'bool')
            for k, v in zip(keys, p):
                s *= self[k] == v
            i = s.nonzero()[0]
            if len(i) > 0:
                matched[ip] = i
        if mask:
            mask = np.zeros(self.ndets, 'bool')
            mask[matched[matched>=0]] = True
            return mask
        return matched

    @classmethod
    def from_column_file(cls, filename, field_map, dtypes={},
                         comments='#', skip=0):
        """
        Load data from columnated ascii file.  field_map a map from
        field name to column number (0-indexed), which can be a dict
        or a list of tuples.
        
        dtypes, comments, and skip are all passed to ascii.read_columns.

        E.g.:
          det_data = StructDB.from_column_file(
              'clusters.txt', [('id', 0),('RA', 1), ('Dec', 2)])
        """
        t_dtypes = {}
        for k,v in list(dtypes.items()):
            t_dtypes[field_map[k]] = v
        # Get column indices
        if isinstance(field_map, dict):
            field_map = list(field_map.items())
        names, t_columns = list(zip(*field_map))
        # Load data arrays
        data = read_columns(filename, t_columns, skip=skip, dtypes=t_dtypes,
                            comments=comments)
        # That's it
        return cls.from_data(list(zip(names, data)))

    def to_column_file(self, filename, fields=None, formats={},
                       header=True):
        if isinstance(filename, basestring):
            fout = open(filename, 'w')
        else:  # must be file-like, then...
            fout = filename

        if fields is None:
            fields = self.dtype.names
        formats_list = []
        for k in fields:
            formats_list.append(formats.get(k, self.formats.get(k)))
            if formats_list[-1] is None:
                formats_list[-1] = moby_fits.get_fmtcode(self[k])
        data_refs = [self[k] for k in fields]

        #fmt_str = ' '.join([all_formats[k] for k in fields]) + '\n'
        #for i in xrange(len(self)):
        #    fout.write(fmt_str % tuple([self[k][i] for k in fields]))
        #del fout
            
        # Write data
        for row in range(len(data_refs[0])):
            tokens = [fmt % d[row] for fmt,d in zip(formats_list, data_refs)]
            if header:
                # Pad, when possible, to keep header aligned with data
                header, didx = '#', 0
                for f, t in zip(fields, tokens):
                    didx += len(t)
                    header += ' ' * max(1,didx-len(header)-len(f)) + f
                    didx += 1
                fout.write(header + '\n')
                # Mark header as written
                header = False
            fout.write(' '.join(tokens) + '\n')
        del fout

    @classmethod
    def from_fits_table(cls, filename, index=1):
        ftr = moby_fits.MobyFitsTableReader(filename, index)
        self = cls.from_data(ftr.get_data())
        self.formats = ftr.get_formats()
        return self

    def to_fits_table(self, filename=None, fields=None, formats=None,
                      fits_formats=None, clobber=True):
        """
        Write the data to a FITS file, including printf formatting
        information for easy ascii dumping.

        'fields' is a list of fields to write; this will default to
        self.fields.  'fits_formats' and 'formats' also default to
        self.<that>, but if these formats are missing they will be
        guessed based on the data.

        This function returns the resulting fits BinTableHDU, so
        pass filename=None if the HDU is what you really want.
        """
        ftw = moby_fits.MobyFitsTableWriter(self, self.formats)
        hdu = ftw.get_hdu()
        if not filename is None:
            ftw.write(filename, clobber=clobber)
        return hdu

    @classmethod
    def from_hdf(cls, filename=None, dataset=None, group=None):
        import h5py
        if isinstance(filename, basestring):
            hfile = h5py.File(filename, 'r')
        else:
            hfile = filename
        assert (group is None or dataset is None)
        if group is not None:
            if isinstance(group, basestring):
                group = hfile[group]
            node = group
            items = list(group.items())
            try:
                dt = [(str(k),v.dtype) for (k,v) in items]
            except:
                raise RuntimeError(
                    'Could not handle %s:%s; perhaps choose a dataset from: %s' % 
                    (hfile.filename, dataset.name, [v.name for k,v in items]))
            self = cls(len(items[0][1]), dt)
            for k,v in items:
                self[str(k)] = v
        else:
            if dataset is None:
                dataset = '/'
            if isinstance(dataset, basestring):
                dataset = hfile[dataset]
            self = cls(dataset.shape, dataset.dtype)
            self[:] = dataset.value
            node = dataset
        self.formats = dict(node.attrs.get('_formats', []))
        if isinstance(filename, basestring):
            hfile.close()
        return self

    def to_hdf(self, filename, dataset=None, group=None, clobber=False,
               compression=None):
        """
        Write the StructDB array to the indicated HDF5 file.  If a
        dataset is named, the array is written as a single, structured
        dataset with that name.  If group is specified, the StructDB
        is written as multiple, simple datasets, with specified group
        name.  We hereby declare that the former approach should be
        taken whenever convenient.
        """
        import h5py
        if isinstance(filename, basestring):
            hfile = h5py.File(filename, 'a')
        else:
            hfile = filename
        assert(int(group is None) + int(dataset is None) == 1)
        dest = dataset
        if dataset is None:
            dest = group
        if dest in hfile:
            if clobber:
                del hfile[dest]
            else:
                raise RuntimeError("Location %s exists in HDF file.  Pass "
                                     "clobber=True to overwrite." % dest)
        if dataset is not None:
            node = hfile.create_dataset(dataset, data=self, compression=compression)
        if group is not None:
            node = hfile.create_group(group)
            for n in self.dtype.names:
                node.create_dataset(n, data=self[n], compression=compression)
        if len(self.formats):
            node.attrs['_formats'] = list(self.formats.items())
        if isinstance(filename, basestring):
            hfile.close()
