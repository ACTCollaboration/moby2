from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
rethread -- many numpy operations are not natively parallel, and thus
can be accelerated using python threading.

Stolen from http://www.scipy.org/Cookbook/Multithreading , by
AMArchibald.  Augmented to determine machine processor count.

"""

import sys, os, threading
from itertools import count

# Use OMP_NUM_THREADS, or CPU count, as default threads value.
DEFAULT_THREADS = os.getenv('OMP_NUM_THREADS')
if DEFAULT_THREADS is None:
    try:
        DEFAULT_THREADS = max([ int(line.split()[-1]) \
                                for line in open('/proc/cpuinfo') \
                                if line.startswith('processor') ]) + 1
    except:
        DEFAULT_THREADS = 4
else:
    DEFAULT_THREADS = int(DEFAULT_THREADS)


def foreach(f,l,threads=None,return_=True):
    """
    Apply f to each element of l, in parallel.
    """
    if threads is None:
        threads = DEFAULT_THREADS
    if threads>1:
        iteratorlock = threading.Lock()
        exceptions = []
        if return_:
            n = 0
            d = {}
            i = zip(count(),l.__iter__())
        else:
            i = l.__iter__()

        def runall():
            while True:
                iteratorlock.acquire()
                try:
                    try:
                        if exceptions:
                            return
                        v = next(i)
                    finally:
                        iteratorlock.release()
                except StopIteration:
                    return
                try:
                    if return_:
                        n,x = v
                        d[n] = f(x)
                    else:
                        f(v)
                except:
                    e = sys.exc_info()
                    iteratorlock.acquire()
                    try:
                        exceptions.append(e)
                    finally:
                        iteratorlock.release()
        
        threadlist = [threading.Thread(target=runall) for j in range(threads)]
        for t in threadlist:
            t.start()
        for t in threadlist:
            t.join()
        if exceptions:
            a, b, c = exceptions[0]
            raise a(b).with_traceback(c)
        if return_:
            r = list(d.items())
            r.sort()
            return [v for (n,v) in r]
    else:
        if return_:
            return [f(v) for v in l]
        else:
            for v in l:
                f(v)
            return

def parallel_map(f,l,threads=None):
    return foreach(f,l,threads=threads,return_=True)

# Numpy extensions
# -- these are things that numpy does not seem to do well on its own

def _1d_idx(data, mask):
    """
    Returns indices into 1st dimension of data, based on mask.
    data      an array
    mask      Either None, an array of bool with size data.shape[0],
              or an array of integers.
    """
    from numpy import arange
    if mask is None:
        return arange(data.shape[0])
    elif mask.dtype == 'bool':
        return nmask.nonzero()[0]
    else:
        return mask # must be indices...

def pnorm(data, mask=None, **kwargs):
    """
    Equivalent to
       (data[mask]**2).sum(axis=1)
    """
    from numpy import array, dot
    idx = _1d_idx(data, mask)
    return array(foreach((lambda i: dot(data[i],data[i])), idx, **kwargs))

def pdot(data, v, mask=None, **kwargs):
    """
    Equivalent to
        dot(data[mask], v)
    for 2-d data and 1-d v.
    """
    idx = _1d_idx(data, mask)
    from numpy import array, dot
    return array(foreach((lambda i: dot(data[i],v)), idx, **kwargs))

def pmax(data, mask=None, **kwargs):
    from numpy import array
    idx = _1d_idx(data, mask)
    return array(foreach((lambda i: data[i].max()), idx, **kwargs))

def pmin(data, mask=None, **kwargs):
    from numpy import array
    idx = _1d_idx(data, mask)
    return array(foreach((lambda i: data[i].min()), idx, **kwargs))
