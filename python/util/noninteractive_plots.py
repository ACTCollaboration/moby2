from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
from matplotlib import use as matplotlib_use
matplotlib_use('Agg')

def noshow(filename='_noshow.png'):
    import pylab
    pylab.savefig(filename)
    pylab.clf()
