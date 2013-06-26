from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os

INCLUDE_KEY = '_include'

class MobyDict(dict):
    """
    Dangerous config file reader / writer.  It's like ACTDict but can
    handle declarations spread over multiple lines.

    This uses "eval" so be careful.  It also isn't very smart about
    comments, so don't try to get away with too much.

    When loading from file use the class method from_file:

       my_config = MobyDict.from_file('config.dict')

    You can also update a MobyDict from a file, much as one would
    dict.update(another_dict).
    """
    do_includes = False
    do_subs = False
    include_path = None
    key_order = None

    @classmethod
    def from_file(cls, filename, do_includes=False, do_subs=False,
                  include_path=None):
        self = cls()
        self.do_includes = do_includes
        self.do_subs = do_subs
        if include_path is None:
            include_path = ['','./']
        self.include_path = include_path
        # Go
        self.update_from_file(filename)
        return self

    def update_from_file( self, filename, do_includes=None, do_subs=None,
                          include_path=None,
                          include_depth=10, subs_depth=10):
        # Do recursion and string substitution?
        if do_includes is None:
            do_includes = self.do_includes
            if include_path is None:
                include_path = self.include_path 
        if do_subs is None:
            do_subs = self.do_subs
        # Load file.
        f = open( filename )
        line_num = 0
        exp = None
        for line in f:
            line_num += 1
            line = line.strip()
            s = line.split('#')
            if len(s) == 0 or len(s[0]) == 0:
                continue
            s = s[0]
            cont = s.strip().endswith('\\')
            if cont:
                s = s.split('\\')[0]
            if exp is None:
                exp = s
                exp_line = line_num
            else:
                exp = exp + s
            # Forced continuation?
            if cont:
                continue
            # Attempt evaluation
            if not '=' in exp:
                continue
            idx = exp.find('=')
            k, v = exp[:idx].strip(), exp[idx+1:].strip()
            try:
                v = eval(v)
            except NameError:
                print('name error in key %s, replacing nans...' % k)
                v = eval(v.replace('nan', 'float(\'nan\')'))
            except SyntaxError:
                # Any multi-line expression will get here at some point
                #print 'syntax error in key %s' % k
                continue
            # We key = value.  Any special handling?
            if do_includes and k == INCLUDE_KEY:
                # Find that file...
                sub_files = v
                if isinstance(sub_files, basestring):
                    sub_files = [sub_files]
                for sub_file in sub_files:
                    for prefix in self.include_path:
                        if prefix == '':
                            prefix = os.path.split(filename)[0]
                        childfile = os.path.join(prefix, sub_file)
                        if os.path.exists(childfile):
                            break
                    else:
                        raise ValueError("Could not include file %s using path %s" % \
                              (sub_file, self.include_path))
                    self.update_from_file(childfile, do_includes=True, do_subs=do_subs,
                                          include_depth=include_depth-1)
            else:
                self[k] = v
            exp = None
        if exp is not None:
            raise RuntimeError('Error while parsing %s -- incomplete or invalid eval block '\
                'starting on line %i' % (filename, exp_line))
        if do_subs:
            self.evaluate_subs(subs_depth)
        f.close()

    def evaluate_subs(self, subs_depth=10):
        for i in range(subs_depth):
            changes = 0
            # Expand strings
            for k, v in list(self.items()):
                if isinstance(v, str) and '%' in v:
                    try:
                        self[k] = v % self
                    except ValueError:
                        raise ValueError('Failed to expand value '\
                            '"%s" in file %s' % (v, filename))

                    if v != self[k]:
                        changes += 1
            if changes == 0:
                break
        
    def write( self, filename, mode = 'w' ):
        f = open( filename, mode )
        keys = list(self.keys())
        keys.sort()
        if self.key_order in [None, 'alphabetical']:
            first_keys = []
        else:
            first_keys = self.key_order
            for k in first_keys:
                if k in keys:
                    keys.remove(k)
        for key in first_keys:
            f.write( "%s = %s\n" % (key,repr(self[key])) )
        for key in keys:
            f.write( "%s = %s\n" % (key,repr(self[key])) )
        f.close()

    write_to_file = write

    def get_deep(self, keys, default=None):
        """
        Call as:
             self.get_deep(key1, key2, ..., keyn, default)

        Returns
             self[key1][key2]...[keyn]
        or the default if, at any point in the descent, a key is not found.
        """
        if isinstance(keys, basestring):
            raise ValueError("first argument should be tuplish, not stringy.")
        p = self
        for k in keys:
            if not k in p:
                return default
            p = p[k]
        return p
                
    def set_deep(self, keys, value):
        """
        Call as:
             self.set_deep(key1, key2, ..., keyn, value)

        Tries to set
             self[key1][key2]...[keyn] = value

        If at any point one of the keys is not found, a dictionary is
        created at that location.
        """
        if isinstance(keys, basestring):
            raise ValueError("first argument should be tuplish, not stringy.")
        depth = len(keys)-1
        p = self
        for k in keys[:-1]:
            if not k in p:
                p[k] = {}
            p = p[k]
        p[keys[-1]] = value
        
    """
    Depot support: this is most useful when subclassing MobyDict for
    object you want to save to / load from a Depot.  The read_result
    and write_result calls will work with or without the tod argument.

    e.g.
         class myWeirdNumbers(MobyDict):
             pass

         wn = myWeirdNumbers()
         wn['a'] = 1.98765
         depot = moby2.util.Depot()
         depot.write_result(wn, 'my_tag', tod)

         ...

         wn = depot.read_result(myWeirdNumbers, 'my_tag', tod)
    """

    @classmethod
    def depot_filename_from_path(cls, path):
        return '%s/%s.dict' % (path, cls.__name__)

    @classmethod
    def readFromPath(cls, path, *args):
        o = cls.from_file(cls.depot_filename_from_path(path))
        return o

    def writeToPath(self, path):
        self.write_to_file(self.depot_filename_from_path(path))

    read_from_path = readFromPath
    write_to_path = writeToPath
