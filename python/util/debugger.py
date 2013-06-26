from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
#
# Based on
#  http://code.activestate.com/recipes/65287-automatically-start-the-debugger-on-an-exception/

import sys

def info(type, value, tb):
   if hasattr(sys, 'ps1') or not sys.stderr.isatty():
      # we are in interactive mode or we don't have a tty-like
      # device, so we call the default hook
      sys.__excepthook__(type, value, tb)
   else:
      import traceback, pdb
      # we are NOT in interactive mode, print the exception...
      traceback.print_exception(type, value, tb)
      print()
      # ...then start the debugger in post-mortem mode.
      pdb.pm()

def interactive(on=True):
    if on:
        sys.excepthook = info
    else:
        sys.excepthook = sys.__excepthook__
