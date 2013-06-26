from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
filebase.py

Objects that help find the locations of datafiles, given some more
basic piece of information (such as the TOD basename).

Note that the factory function

   moby2.scripting.products.get_filebase(**params)

is a good way to get a filebase from a scripting context.
"""

import os
import time

class MobyFilebase:
    """
    Abstract class from which other Filebases inherit.  But mostly
    just holds these comments...

    Inheritors should implement one or both of

        filename_from_name(self, name, single=False)
        filename_from_id(self, id_, single=False)
    """
    def get_full_path(self, name, multi=False):
        return self.filename_from_name(name, single=not multi)

format_finders = {
    'dirfile': lambda x: x,
    'dirfile.zip': lambda x: x+'.zip',
}


class FilesystemFilebase(MobyFilebase):
    def __init__(self, root=None, structure=None, formats=None):
        self.root = root
        self.structure = structure
        self.formats = formats
        if self.formats is None:
            self.formats = ['dirfile.zip', 'dirfile']
    def filename_from_name(self, name, single=False):
        if self.structure is None:
            fullname = os.path.join(self.root, name)
        elif self.structure == 'first5':
            fullname = os.path.join(self.root, name[:5], name)
        elif self.structure == 'yyyymmdd':
            nightname = time.strftime('%Y%m%d', time.gmtime(int(name[:10])))
            fullname = os.path.join(self.root, nightname, name)
        else:
            raise ValueError
        for format in self.formats:
            formatted_name = format_finders[format](fullname)
            if os.path.exists(formatted_name):
                if single:
                    return formatted_name
                return [formatted_name]
        if single:
            return None
        return []
    def filename_from_id(self, id_):
        return self.filename_from_name(id_)


class MultiFilebase(MobyFilebase):
    """
    A multi-drop extension.
    """
    def __init__(self, components=[]):
        self.components = []
        for c in components:
            self.add_filebase(c)
    def add_filebase(self, fb):
        self.components.append(fb)
    def filename_from_name(self, name, single=False):
        results = []
        for c in self.components:
            fn = c.filename_from_name(name, single=single)
            if single and fn is not None:
                return fn
            if not single:
                results.extend(fn)
        if single:
            return None
        return results

class NullFilebase(MobyFilebase):
    def __init__(self):
        pass
    def filename_from_name(self, name, single=False):
        if single:
            return None
        return []

