from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np

def read_columns(filename, columns=None, skip=0, dtypes={},
                 comments='#'):
    """
    Read whitespace-delimited columnar data from a text file.

    Pass (0-indexed) list of columns if you only want some of them.

    Each column of data is converted to a numpy array.  To force a
    data type, pass it in the column's entry in the dtypes array.
    E.g. dtypes = {0: int, 1: 'float32'}

    Skips the first skip rows.  Ignores blank lines and lines starting
    with anything character in the comments string.

    Returns a tuple of arrays.
    """
    data = None
    if isinstance(filename,basestring):
        fin = open(filename)
    else:
        fin = filename
    for line in fin:
        row = line.split()
        if len(row) == 0 or row[0][0] in comments:
            continue
        if skip > 0:
            skip -= 1
            continue
        if data is None:
            if columns is None:
                nc = len(row)
                columns = [i for i in range(nc)]
            else:
                nc = len(columns)
            data = [[] for c in columns]
        for i,c in enumerate(columns):
            data[i].append(row[c])
    if data is None:
        return None
    for i,c in enumerate(columns):
        dtype = dtypes.get(c)
        if dtype is not None:
            data[i] = np.array(data[i], dtype=dtype)
        else:
            for dtype in [int,float,None]:
                try:
                    data[i] = np.array(data[i], dtype=dtype)
                except ValueError:
                    continue
                break
            else:
                raise ValueError("Could not convert column %i of %s" % \
                    (c, filename))
            
    return tuple(data)
