from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

"""
Generic database support, I guess.  Instrument-specific support should
be provided elsewhere (e.g. scripting).
"""

import numpy as np
from . import ctime
import os

class TODRecordList(list):
    """
    A list that can also project out certain attributes from its
    members, using, e.g., .get_fields('basename') or
    .get_fields(['min_RA_south', 'max_RA_south']).
    """
    def get_fields(self, fields, as_array=True):
        if isinstance(fields, basestring):
            values = [getattr(x, fields) for x in self]
            if not as_array:
                return values
            return np.array(values)
        return [self.get_fields(f) for f in fields]

def convert_ctimes(ctimes):
    """
    Process ctimes, converting string dates to floats.

    If ctimes is a list or tuple, each entry is converted and returned
    in the same kind of container.
    """
    is_tuple = isinstance(ctimes, tuple)
    is_list = isinstance(ctimes, list)
    if not is_tuple and not is_list:
        return convert_ctimes((ctimes,))[0]
    results = []
    for x in ctimes:
        if isinstance(x, basestring):
            # Parse as UTC date
            x = ctime.from_string(x)
        results.append(float(x))
    if is_list:
        return results
    return tuple(results)


class TODList(list):
    """
    @brief Modified list that implement methods to operate with TOD lists:
           Operator:
                        a + b:  Add to a all TODs in b not present in a.
                        a += b: Same as a = a + b.
                        a - b:  Removes from a all matching TODs in b.
                        a -= b: Same as a = a - b.
                        a * b:  Leaves only TODs in a and b.
                        a *= b: Same as a = a * b.
    """
    @classmethod
    def from_file(cls, filename):
        self = cls()
        f = open(filename)
        for l in f:
            if l.strip()[0] != "#":
                self.append(l.split()[0].strip())
        f.close()
        return self

    def to_file(self, filename, header = None):
        f = open(filename, "w")
        if header is not None: f.write(header)
        for obs in self: f.write("%s\n"%obs)
        f.close()

    def __add__(self, other):
        newlist = TODList()
        newlist.extend(self)
        newlist.extend(self._find_new(other)[0])
        return newlist

    def __iadd__(self, other):
        self.extend(self._find_new(other)[0])
        return self

    def __sub__(self, other):
        newlist = TODList()
        newlist.extend(self)
        old = self._find_new(other)[1]
        for obs in old:
            if obs in newlist: newlist.remove(obs)
        return newlist

    def __isub__(self, other):
        old = self._find_new(other)[1]
        for obs in old: 
            if obs in self: self.remove(obs)
        return self

    def __mul__(self, other):
        newlist = TODList()
        newlist.extend(self._find_new(other)[1])
        return newlist

    def __imul__(self, other):
        old = self._find_new(other)[1]
        self = TODList()
        self.extend(old)
        return self

    def _find_new(self,other):
        n_self = []
        for obs in self: 
            if os.path.basename(obs)[-4:] == '.zip':
                n_self.append(os.path.basename(obs)[:-4])
            else:
                n_self.append(os.path.basename(obs))
        n_other = []
        for obs in other: 
            if os.path.basename(obs)[-4:] == '.zip':
                n_other.append(os.path.basename(obs)[:-4])
            else:
                n_other.append(os.path.basename(obs))
        new = []
        old = []
        for name, obs in zip(n_other, other):
            if name in n_self: old.append(self[n_self.index(name)])
            else: new.append(obs)
        return (new,old)
