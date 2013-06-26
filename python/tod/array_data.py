from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
import numpy as np
import moby2

"""
Terminology: det_uid

The index that is meant to describe the absolute, invariant
enumeration of detectors in an array should be called det_uid, and NOT
det_idx or det_index or anything else.

The other things: det_idx, det_index, ... can be used as indicies into
a dimension of an array that represents different detectors.  For
example, tod.data might contain data for column 0 only.  Then when we
loop over those:

for det_idx in range(tod.data.shape[0]):
  print 'Time stream ', det_idx
  print ' The detector uid at this position is ', tod.det_uid[det_idx]
  print ' This corresponds to multiplexing column ', \
      my_det_data['col'][tod.det_uid[det_idx]]

We more or less keep things straight:
 Time stream 0
  The detector uid at this position is 0
  This corresponds to multiplexing column 0
 Time stream 1
  The detector uid at this position is 32
  This corresponds to multiplexing column 1
 ...

It gets confusing if det_idx and det_uid are not consistently named in
a way that distinguishes absolute detector enumeration from the thing
that indexes dimension 0 of tod.data.
"""


class ArrayData(dict):
    """
    Mini-database for storage of detector information.  This is most
    useful for storing constant properties of detectors.  For example,
    to associate each detector index with a multiplexing row and
    column.  Since pointing, time constants, etc., tend to change they
    should be stored elsewhere.

    Data are stored as dictionary items.  They should all be numpy

    arrays, and have a consistent length equal to self.ndets.  Then
    you can just grab the det_uid like:

      my_uid = (ddata['row'] == 15).nonzero()[0]

    Which is the same as

      my_uid = ddata.select_outer({'row': [15]})
    """

    # What do we have here?
    ndets = None

    def __init__(self, ndets):
        self.ndets = ndets

    def add_property(self, name, data):
        self[name] = np.array(data)

    def get_property(self, name, det_uid=None):
        if isinstance(name, basestring):
            return self.get_property([name], det_uid=det_uid)[0]
        if det_uid is None:
            output = [self[n].copy() for n in name]
        return [self[n][det_uid] for n in name]

    def get_info(self, det_uid):
        output = dict([(k,v[det_uid]) for k,v in list(self.items())])
        output['det_uid'] = det_uid
        return output

    def copy(self):
        new = self.__class__(self.ndets)
        for k,v in list(self.items()):
            new[k] = v.copy()
        return new
    
    """
    Database-like operations.
    """

    def select_outer(self, props={}, det_uid=None, mask=False):
        """
        Select a set of detectors based on outer product of the
        (key,values) pairs in the props dict.  "props" is a dictionary
        where each key is a detector property and each value is a list
        of desired values for that property.  For example, if props
        is:

        {'row': [1,2,3], 'col': [2,4]}

        Then the function returns the (probably) 9 detectors that
        have 'row' property in [1,2,3] and 'col' property in [2,4].

        If a value is given as "None", it is ignored (i.e. no attempt
        is made to match that key).

        Returns the matching indices, or a boolean mask if mask==True.

        If det_uid is provided, the mask or indices are with respect
        to the indices in det_uid, in the sense that
        my_det_uid[select_outer(..., det_uid=my_det_uid)] are the det_uid
        that match.
        """
        # Create a mask to track the matches
        selection = np.ones(self.ndets, 'bool')
        # Reduce the mask with each property
        for k, values in list(props.items()):
            if values is None:
                continue
            if k in self:
                try:
                    # values should be a list...
                    selection *= [v in values for v in self[k]]
                except TypeError:
                    # ...but allow the odd single item
                    selection *= [v==values for v in self[k]]
            elif k == 'det_uid':
                selection[values] *= True
            else:
                print("Property %s not defined for array." % k)
                selection[:] = False
                break
        if det_uid is not None:
            selection = selection[det_uid]
        if mask:
            return selection
        return selection.nonzero()[0]

    def select_inner(self, props, det_uid=None, mask=False):
        """
        Select a set of detectors based on inner product of the
        (key,values) pairs in the props dict.  "props" is a dictionary
        where each key is a detector property and each value is a list
        of desired values for that property.  The value lists must all
        be the same length.  For example, if props is:

        {'row': [1,2,3], 'col': [2,4,5]}

        Then the function returns the det_uid of 3 detectors with
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
        for k in keys:
            if not k in self:
                print("Property %s not defined for array." % k)
                return [None] * len(props)
        # Get tuples of all detectors
        if det_uid is None:
            data = list(zip(*[self[k] for k in keys]))
        else:
            data = list(zip(*[self[k][det_uid] for k in keys]))
        output = np.zeros(len(props), 'int')
        output -= 1
        for i,p in enumerate(props):
            try:
                output[i] = data.index(p)
            except ValueError:
                pass
        if mask:
            mask = np.zeros(self.ndets, 'bool')
            mask[output[output>=0]] = True
            return mask
        return output


    """
    Class methods for I/O from various formats.
    """

    @classmethod
    def from_dict(cls, data, exclude=[]):
        exclude_ = ['det_uid', 'inverted_column']
        exclude_.extend(exclude)
        det_uid = np.array(data['det_uid'])
        self = cls(det_uid.max()+1)
        for k, v in list(data.items()):
            if k not in exclude_:
                v = np.array(v)
                self.add_property(k, np.zeros(self.ndets, v.dtype))
                self[k][det_uid] = v
        if 'inverted_column' in list(data.keys()):
            self.add_property('inverted_column', data['inverted_column'])
        return self

    @classmethod
    def from_fits_table(cls, filename):
        # This eventually should just *be* a StructDB... that's where we head.
        mfr = moby2.util.StructDB.from_fits_table(filename)
        det_uid = mfr['det_uid']
        assert(np.all(det_uid == np.arange(len(det_uid))))
        self = cls(len(det_uid))
        for k in mfr.dtype.names:
            self.add_property(k, mfr[k])
        return self
        
    @classmethod
    def from_column_file(cls, filename, field_map, casts={}):
        """
        Load columnar data from a file, with each uncommented row
        identified as a new detector.

        field_map is a dict mapping a detector propertery to a column
        number.  "casts" are optional converters; e.g. {'hex_type':
        str}.  The default cast is int.

        E.g.:
          det_data = DetectorData.from_column_file(
                   'ar1.txt',
                   {'row': 1, 'col': 2', 'hex_name': 3},
                   casts={'hex_name': str})
        """
        keys = list(field_map.keys())
        cols = [field_map[k] for k in keys]
        casts = [casts.get(k, int) for k in keys]
        values = [[] for k in keys]

        for line in open(filename):
            w = line.split()
            if len(w) == 0 or w[0] == '#':
                continue
            for col, cast, value in zip(cols, casts, values):
                value.append(cast(w[col]))
        
        for k, v in zip(keys, values):
            self.add_property(k, v)

        return self

    def to_struct_db(self, det_uid=None):
        """
        I want this to be a StructDB so badly now.
        """
        data = [('det_uid', self['det_uid'])]
        for k in sorted(self.keys()):
            if k == 'det_uid':
                continue
            data.append((k, self[k]))
        db = moby2.util.StructDB.from_data(data)
        if det_uid is None:
            return db
        idx = db.select_inner({'det_uid': det_uid})
        if np.any(idx < 0):
            raise ValueError("ArrayData does not contain all "\
                "requested det_uid.")
        return db[idx]

            
