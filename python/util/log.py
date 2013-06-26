from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
A logger, with levels and file output.  There is a global logger here.
There are module-level functions that mimic psLib logging.
"""

import sys
import time
import traceback

class Logger:
    verbosity = 1
    show_time = True
    log_file = None
    
    def __init__(self, verbosity=1, show_time=True, log_file=None,
                 hook_exceptions=False, append=False):
        self.verbosity = verbosity
        self.show_time = show_time
        self.set_file(log_file, append=append)
        if hook_exceptions:
            self.hook_exceptions()

    def hook_exceptions(self):
        ExceptionLogger.attach(self)

    def set_file(self, log_file=None, append=False):
        if log_file is not None:
            log_file = open(log_file, {True: 'a', False: 'w'}[append])
        self.log_file = log_file

    def set_level(self, level):
        if level is not None:
            self.verbosity = level

    def __call__(self, *args):
        return self.trace(0, ' '.join(map(str,args)))

    def trace(self, level, msg):
        if self.show_time:
            s = time.asctime() + ' '
        else:
            s = ''
        s = s + ' '*max(0,level)
        for line in msg.split('\n'):
            if level <= self.verbosity:
                print(s + line)
            if self.log_file is not None and level < 6:
                self.log_file.write(s + line + '\n')
                if level < self.verbosity:
                    self.log_file.flush()

    def trace_exception(self, stack):
        self.trace(0, '****************************')
        self.trace(0, 'Unhandled exception!  Trace:')
        for line in stack:
            self.trace(0, line)

    def init_from_params(self, params):
        self.set_level(params.get('verbosity'))
        self.set_file(params.get('log_file'), params.get("append", False))
        log_exceptions = params.get('log_exceptions', False)
        if log_exceptions:
            self.hook_exceptions()


class ExceptionLogger:
    """
    Replacement for sys.excepthook function that copies the traceback to a
    logger.
    """

    def __init__(self, parent):
        self.parent = parent

    def __call__(self, etype, value, tb):
        self.parent.trace_exception(traceback.format_exception(etype, value, tb))
        sys.except_hook = self.old_excepthook
        
    @classmethod
    def attach(cls, parent):
        self = cls(parent)
        self.old_excepthook = sys.excepthook
        sys.excepthook = self

    def remove(self):
        if self.old_excepthook is not None:
            if hasattr(sys, 'excepthook'):
                sys.excepthook = self.old_excepthook

    def __del__(self):
        self.remove()

#
# Global logger?
#

logger = Logger()

def traceSetLevel(unit, level):
    global logger
    if logger is None:
        logger = Logger()
    logger.set_level(level)

def trace(unit, level, msg):
    global logger
    if logger is None:
        logger = Logger()
    logger.trace(level, msg)

def set_logger(new_logger):
    global logger
    logger = new_logger

def get_logger():
    global logger
    return logger
