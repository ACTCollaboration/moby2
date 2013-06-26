from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

def get_selection(cut_params, *data):
    """
    Return a boolean mask that indicates which elements of [data,...]
    pass the cuts specificed in cut_params.

    Each of the *data arguments must be a dict or StructDB such that
    data[i][key] returns an array of the same length.

    The cut_params is a list containing tuples of the form
       (key, vals)
    
    The datasets passed in data are checked to see if they contain
    data under the key, and if they do the data are compared to vals.

    Comparison to vals is achieved differently depending on whether
    vals is an instance of list, tuple, or neither.

    If vals is a list, the cut will pass only when the data are
    members of the list.  E.g.
       ('array', ['AR1', 'AR2'])

    If vals is a tuple of length 2, the cut will pass only when the
    data are greateer than or equal to vals[0] but less than vals[1].
    If vals is a tuple of length 3, the cut will only pass when the
    data lie in the specified interval, modulo vals[2].  E.g.
       ('ctime', (1234560000, 123470000))
       ('hour_utc', (-1,11,24))
    """
    n = len(data[0])
    mask = np.ones(n, bool)
    for cut_p in cut_params:
        key, vals = cut_p
        for d in data:
            if isinstance(d, np.ndarray) and key in d.dtype.names:
                break
            elif key in d:
                break
        else:
            raise ValueError('Key %s not found in any provided '
                               'data sets.' % key)
        x = d[key]
        if isinstance(vals, list):
            mask *= [_x in vals for _x in x]
        elif isinstance(vals, tuple):
            if len(vals) == 3:
                _x = (x - vals[0]) % vals[2] + vals[0]  # in-branch
            elif len(vals) == 2:
                _x = x
            else:
                raise ValueError('Key %s constraint is a tuple of '
                                   'length %i; expected 2 or 3 items.' % 
                                   (key, len(vals)))
            mask *= (vals[0] <= _x) * (_x < vals[1])
        else:
            mask *= (x == vals)
    return mask

