from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
if os.getenv('DISPLAY') is None or os.getenv('DISPLAY')=='':
    #print 'Warning: defaulting to non-interactive matplotlib backend'
    from .noninteractive_plots import show
else:
    from pylab import show
