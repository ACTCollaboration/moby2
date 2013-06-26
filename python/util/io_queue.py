from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
The module faciliates the queueing of IO operations.  That means you
can have the next TOD loading, in the background, while you analyze
the first.

This will only work if the main blocking IO function doesn't hold the
global interpreter lock while active.  This is true of most python
native IO, and can be forced onto swigged functions where desired.
(See the wrapping of mobyLib.readTODThreadful, for swig example.)
"""

import threading

def _thread_func(target, callback, args, **kwargs):
    """
    Launches target(*args, **kwargs), and upon completion passes
    the result to callback(...).

    This helps us get around threading.Thread's unwillingness
    to store the return value of the target function.
    """
    output = target(*args, **kwargs)
    callback(output)

class IOQueuer:
    def __init__(self, target, target_kwargs={},
                 args=[], kwargs=[]):
        self.target, self.target_kwargs = target, target_kwargs
        self.args, self.kwargs = args, kwargs
        self.thread = None
        self.index = 0
        self.prefetched = None
        self.prefetch()
        
    def prefetch(self):
        """
        Begin the next IO job.
        """
        if self.index >= len(self.args):
            self.thread = None
            return
        # Construct kwargs
        kw = {}
        kw.update(self.kwargs[self.index])
        kw.update(self.target_kwargs)
        # Embed args
        args = [self.target, self.callback, self.args[self.index]]
        # Launch thread
        self.thread = threading.Thread(target=_thread_func,
                     args=args, kwargs=kw)
        self.index += 1
        self.thread.start()
        return self.thread
        
    def callback(self, output):
        self.prefetched = output

    def get_next(self, block=True):
        """
        Obtain the next IO job result.  If a prefetched result is
        available, it is returned immediately.

        If block=False and no prefetched result is available, then
        get_next will return None.  Otherwise it will block until the
        next result is available.
        """
        if self.thread is None:
            if self.prefetch() is None:
                raise StopIteration
        if not block and self.thread.isAlive():
            return None
        self.thread.join()
        tod = self.prefetched
        self.prefetched = None
        self.prefetch() # Start the next read
        return tod

class TODQueuer(IOQueuer):
    """
    Instantiate with a list of tod filenames, and a list of kwargs to
    pass to TOD.read with each filename.
    
    The IOQueuer will immediately begin loading the first TOD.  It can
    be retrieved through the get_next method.  Upon retrieval, the next
    TOD will begin queuing.
    """
    def __init__(self, filenames, kwarg_list):
        import moby2
        IOQueuer.__init__(self, moby2.TOD.from_dirfile, {},
                          [[f] for f in filenames],
                          kwarg_list)

