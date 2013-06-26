from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import os
import errno

def make_sure_place_exists(filename):
    path = os.path.split(filename)[0]
    if path == '':
        return
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

class OutputFiler:
    def __init__(self, prefix=None, create_missing=True, data={}):
        if prefix is None:
            prefix = './'
        self.prefix = prefix
        self.create_missing = create_missing
        self.data = data.copy()

    def __call__(self, *args, **kwargs):
        """
        Same as self.get_filename.
        """
        return self.get_filename(*args, **kwargs)

    def get_filename(self, filename, data=None, tod_info=None, create=None,
                     remove_existing=False):
        """
        Return the filename, taking into account this OutputManager's
        prefix, and including a call to filename.format(**stuff),
        where stuff includes the contents of "data" (a dict) and
        tod_info (a TODInfo), if provided.
        
        Note that filenames beginning with explicit root paths, such
        as '/', './', or '../' will not have the prefix applied.

        If create is True, or self.create_missing == True, the output
        path will be created if it does not alread exist.
        """
        if filename is None:
            return None
        if \
                (filename[:1] != '/') and \
                (filename[:2] != './') and \
                (filename[:3] != '../'):
            filename = self.prefix + filename
        # Create formatting data
        fmt_data = self.data.copy()
        if data is not None:
            fmt_data.update(data)
        if tod_info is not None:
            fmt_data.update(tod_info.get_dict())
        # Apply format; create dirs if necessary, and return.
        full_path = filename.format(**fmt_data)
        if create is None:
            create = self.create_missing
        if create:
            make_sure_place_exists(full_path)
        # Remove if exists?
        if remove_existing and os.path.exists(full_path):
            os.remove(full_path)

        return full_path

    def savefig(self, name, clf=True, figure=None, **kwargs):
        import pylab as pl
        if figure is None:
            figure = pl.gcf()
        figure.savefig(self.get_filename(name, **kwargs))
        if clf:
            figure.clf()
        
class MobyContext(object):
    outputm = None
    logger = None
    def get_outputm(self, outputm=None):
        if outputm is not None:
            return outputm
        if self.outputm is None:
            self.outputm = OutputFiler(prefix='./')
        return self.outputm
    def get_logger(self, logger=None):
        if logger is not None:
            return logger
        if self.logger is None:
            self.logger = moby2.util.log.get_logger()
        return self.logger
