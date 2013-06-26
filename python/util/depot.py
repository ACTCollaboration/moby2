from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
The Depot is used to generate filenames, and provide a unified
interface for saving and loading analysis results on disk.

Once a storage scheme is standardized, then objects can be loaded
through a sequence such as

  depot = Depot(depot_path=...)
  cuts = depot.read_object(TODCuts, tag='new', tod=tod)

If an object supports methods read_from_path and write_to_path, it
should be compatible with the depot.  The depot read_object and
write_object will pass keyword arguments in "options" through to the
object's methods.  E.g.:

  depot.write_object(my_result, structure='{class}/{season}.txt',
                     season=2013,
                     options={'format': 'ascii'})

Classes intended for use with the depot can provide a default file
system structure through the _depot_structure attribute.  This will be
used to generate the path that is passed to read_from_path,
write_to_path.
"""

import moby2

#from moby2.instruments import actpol
import os
import time

from .decorators import deprecatedmethod

DEFAULT_NIGHT_OFFSET = -9*3600
DEFAULT_NIGHT_FORMAT = '%Y-%m-%d' 


DEFAULT_STRUCTURE_SIMPLE = '{class}/{tag}'
DEFAULT_STRUCTURE_PERTOD = '{class}/{tag}/{first_five}/{tod_name}'

class Depot:
    night_offset = DEFAULT_NIGHT_OFFSET
    night_format = DEFAULT_NIGHT_FORMAT
    class_aliases = None

    def __init__(self, depot_path=None, depot_env=None, structures={}):
        self.class_aliases = {}
        if depot_path is None and depot_env is None:
            raise ValueError("No depot or depot environment variable was specified.")
        if depot_path is None:
            depot_path = os.getenv(depot_env)
        self.depot_path = depot_path
        if not os.path.exists(self.depot_path):
            raise RuntimeError("No depot found at %s" % self.depot_path)
        self.structures = {}
        self.structures.update(structures)

    def get_codes(self, **kwargs):
        """
        Generate a dictionary of strings to use when building
        filenames.  Data to parse must be passed in as keyword
        arguments, e.g.

            depot.get_codes(cls=TimeConstants, array='ar1', tag='new')
        
        Recognized arguments:
            obj
            cls              (string or classobj; otherwise extracted from object)
            tag
            tod
            tod_info         (will be preferred to tod.info)
            ctime            (else: extracted from tod_info / tod.info)
            array            (else: extracted from tod_info / tod.info)
            instrument       (else: extracted from tod_info / tod.info)
            night_name       (else: computed from ctime, night_format, night_offset)
            night_format     (overrides default night_format)
            night_offset     (overrides default night_offset
        """
        data = {}
        if not kwargs.get('cls') is None:
            cls = kwargs['cls']
            #if isinstance(cls, DepotResult):
            #    # They passed us an instance; we can deal with that.
            #    cls_name = cls.__class__.__name__
            #el
            if isinstance(cls, basestring):
                # It's just the name of the class.
                cls_name = cls
            elif hasattr(cls, '__name__'):
                cls_name = cls.__name__
            else:
                cls_name = cls.__class__.__name__
        elif not kwargs.get('obj') is None:
            cls_name = kwargs['obj'].__class__.__name__
        else:
            cls_name = None
        if not cls_name is None:
            data['class'] = cls_name
        # TOD info?
        tod_info = kwargs.get('tod_info', None)
        tod = kwargs.get('tod', None)
        if tod_info is None:
            if isinstance(getattr(tod, 'info', None),
                          moby2.TODInfo):
                tod_info = tod.info
            elif isinstance(tod, basestring):
                tod_info = moby2.scripting.get_tod_info({"filename":tod})
        # Copy instrument, array, season from tod_info unless otherwise provided
        for k in ['instrument', 'array', 'season', 'ctime']:
            if k in kwargs:
                data[k] = kwargs[k]
            elif not tod_info is None:
                v = getattr(tod_info, k, None)
                if not v is None:
                    data[k] = v
        # And TOD name
        if 'tod_name' in kwargs:
            data['tod_name'] = kwargs['tod_name']
        elif not tod_info is None and \
                not getattr(tod_info, 'name', None) is None:
            data['tod_name'] = tod_info.name
        if 'tod_name' in data and data.get('ctime') is None:
            try:
                data['ctime'] = int(data['basename'].split('.')[0])
            except:
                pass
        # Derivatives of ctime
        if data.get('ctime') is not None:
            if data.get('night_name') is None:
                night_time = data['ctime'] + self.night_offset
                data['night_name'] = time.strftime(
                    self.night_format, time.gmtime(night_time))
            if data.get('first_five') is None:
                data['first_five'] = str(data['ctime'])[:5]
        if 'tag' in kwargs:
            data['tag'] = kwargs['tag']
        return data

    def get_structure(self, cls=None, obj=None, data={}):
        if cls is None and not obj is None:
            cls = obj.__class__
        if not cls is None:
            if isinstance(cls, basestring):
                cls_name = cls
            else:
                cls_name = cls.__name__
            if cls_name in self.structures:
                structure = self.structures[cls_name]
            elif hasattr(cls, '_depot_structure'):
                structure = cls._depot_structure
            else:
                if 'tod_name' in data:
                    structure = DEFAULT_STRUCTURE_PERTOD
                else:
                    structure = DEFAULT_STRUCTURE_SIMPLE
            return structure
        return None

    def get_relative_path(self, cls=None, tag=None, tod=None, tod_info=None,
                          structure=None, **kwargs):
        """
        Returns relative path, of the form <class>/<tag>/<stuff...>.
        """
        # Assemble dictionary of strings.
        data = self.get_codes(cls=cls, tag=tag, tod=tod, tod_info=tod_info,
                              **kwargs)
        # Identify a structure
        if structure is None:
            # Check this depot's overrides
            cls_name = data.get('class')
            if cls_name in self.structures:
                structure = self.structures[cls_name]
            else:
                structure = self.get_structure(cls=cls, obj=kwargs.get('obj'),
                                               data=data)
        filename = structure.format(**data)
        if not getattr(cls, '_depot_simple_suffix', None) is None:
            filename = filename + cls._depot_simple_suffix
        return filename

    def get_full_path(self, cls=None, tag=None, tod=None, tod_info=None,
                      structure=None, make_dirs=False, **kwargs):
        rel_path = self.get_relative_path(cls, tag, tod, tod_info, structure,
                                          **kwargs)
        full_path = os.path.join(self.depot_path, rel_path)
        if make_dirs:
            di, bn = os.path.split(full_path)
            if not os.path.exists(di):
                os.makedirs(di)
        return full_path

    def read_object(self, cls, **kwargs):
        filename = self.get_full_path(cls=cls, **kwargs)
        return cls.read_from_path(filename, **kwargs.get('options', {}))

    def write_object(self, obj, **kwargs):
        kwargs = kwargs.copy()
        if kwargs.get('make_dirs') is None:
            kwargs['make_dirs'] = True
        filename = self.get_full_path(obj=obj, **kwargs)
        return obj.write_to_path(filename, **kwargs.get('options', {}))
        
    @deprecatedmethod
    def read_result(self, cls, tag, tod=None, tod_info=None, format=None):
        return self.read_object(cls, tag=tag, tod=tod, tod_info=tod_info,
                                format=None)
    
    @deprecatedmethod
    def write_result(self, obj, tag, tod=None, tod_info=None, format=None,
                     force=False):
        return self.write_object(obj, tag=tag, tod=tod, tod_info=tod_info,
                                 format=format)

    @deprecatedmethod
    def check_result(self, cls, tag=None, tod=None, tod_info=None,
                     format=None):
        filename = self.get_full_path(cls=cls, tag=tag, tod=tod,
                                      tod_info=tod_info)
        return os.path.exists(filename)
        
