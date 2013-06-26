from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""Functions and classes to assist with using HDF5 files for archives
of TOD meta-data.  The main class here is HDFArchive, which can, for
example, be subclassed as a TODFlagsArchive and used to store TODFlags
objects.

"""

# String handling in python2 / python3: HDF5 distinguishes between
# ascii and unicode strings, and so does python.  In python2 and
# python3, reading / writing b'this' and u'this' will interoperate
# consistently.  The problem arises when reading or writing 'this',
# i.e. an undecorated "string".  In python2 this is encoded into HDF
# in the same way as b'this'.  But in python3 it is equivalent to
# u'this'.
#
# Note that this seems only to affect data _values_.  It is not
# relevant for group and dataset names, for whatever reason...
#
# Because moby2 is primarily a legacy code associated with python2, we
# defer to data written by python2.  This means encoding string-like
# data as bytes whenever python2 would have encoded that data as
# bytes.  For the bas HDFArchive class, this only really affects the
# polymorphism management string _moby2_class_name.  Unfortunately,
# subclasses will need to be careful when encoding objects into HDF --
# they can use HDFArchive static methods encode_23 and decode_23 to
# auto-convert on write / read.


# @deprecated
def hdf_set_class(target, class_name, class_version):
    return HDFArchive.set_class(target, class_name, class_version)

# @deprecated
def hdf_get_class(target):
    return HDFArchive.get_class(target)


class HDFArchive:
    _moby2_class_name = 'generic_archive'
    _moby2_class_version = 0
    _moby2_read_func = None
    _moby2_write_func = None

    @staticmethod
    def decode_23(x):
        """
        Convert x to a str.  Used when loading stuff from hdf5, where
        python2/3 complications arise.  The argument will usually be a
        bytes() object, such as b'field_name', and in python2 that's the
        same as a str so it will be returned unmodified.  In python3 it
        will be decoded to unicode, which is the same as a string.
        """
        if isinstance(x, str):
            return x
        return x.decode()

    @staticmethod
    def encode_23(x):
        """
        Convert x to a bytes object.  Used when writing stuff to hdf5,
        where python2/3 complications arise.  The argument will usually be
        a string, and in p ython2 that's the same as str.  In python3 we
        encode it to utf8.
        """
        if isinstance(x, bytes):
            return x
        return x.encode('utf8')

    @staticmethod
    def set_class(target, class_name, class_version):
        # Store in the archive as bytes('class name'), for backwards compat.
        target.attrs['_moby2_class_name'] =  HDFArchive.encode_23(class_name)
        target.attrs['_moby2_class_version'] =  class_version

    @staticmethod
    def get_class(target):
        return (HDFArchive.decode_23(target.attrs['_moby2_class_name']),
                target.attrs['_moby2_class_version'])

    @staticmethod
    def check_class(target, class_name, class_version):
        c, v = HDFArchive.get_class(target)
        assert(c == class_name)
        assert(v == class_version)

    def __init__(self, filename, mode='r', group=None, do_init=False, check_class=True):

        """Open the specified HDF5 file for use as an HDFArchive.  If group is
        specified, write to that address in the file.

        The mode can be 'r', 'w', or 'a'.  Their meanings are:

        * 'r': Open existing file for reading.  Fails if file does not
          exist (or if specified group is not found inside).  This is
          intended to prevent the creation of empty files all over the
          place, when an incorrect file or group name is given.

        * 'a': Open a file for writing and reading.  If file does not
          exist, it is created.

        * 'w': Same as 'a', but deprecated, because typical usage of
          mode flags would have 'a' describe this behavior.

        do_init=True can be used to force the application of class
        name/version attributes onto the output target; use with
        caution.
        """
        import h5py, os

        if mode not in ['r', 'w', 'a']:
            raise ValueError("Expected mode r or a, or maybe w")

        self.mode = mode
        if mode is 'r':
            if not os.path.exists(filename):
                raise IOError("File does not exist: %s" % filename)
            self.hfile = h5py.File(filename, 'r')
            if group is None:
                self.target = self.hfile
            else:
                self.target = self.hfile[group]
        else:
            do_init = do_init or (not os.path.exists(filename))
            self.hfile = h5py.File(filename, 'a')
            if group is None:
                self.target = self.hfile
            elif group in self.hfile:
                self.target = self.hfile[group]
            else:
                self.target = self.hfile.create_group(group)
                do_init = True
        
        if do_init:
            HDFArchive.set_class(self.target, self._moby2_class_name,
                                 self._moby2_class_version)
        elif check_class:
            c, v = hdf_get_class(self.target)
            assert(c == self._moby2_class_name)
            assert(v == self._moby2_class_version)

    def close(self):
        if hasattr(self, 'hfile') and self.hfile is not None:
            self.hfile.close()
            self.hfile = None

    def __del__(self):
        self.close()

    # Context management.  HDF5 files like to be closed properly.

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    # Read / write / check individual items.

    def get_item(self, tod_id, read_func=None):
        if read_func is None:
            return self._moby2_read_func(self.target[tod_id])
        else:
            return read_func(self.target[tod_id])

    def set_item(self, tod_id, item, clobber=False, write_func=None):
        assert (self.mode != 'r')
        if tod_id in self.target:
            assert clobber
            print('Replacing %s/%s...' % (self.target.name, tod_id))
            del self.target[tod_id]
        dest = self.target.create_group(tod_id)
        if write_func is None:
            self._moby2_write_func(item, dest)
        else:
            write_func(item, dest)

    def del_item(self, tod_id):
        assert (self.mode != 'r')
        if not tod_id in self.target:
            return False
        del self.target[tod_id]
        return True

    def has_item(self, tod_id):
        return tod_id in self.target

    def keys(self):
        return list(self.target.keys())

    # For use in parallel processing.

    def join_from(self, file_list, delete_from_source=False,
                  delete_when_empty=True, clobber=True):
        """Open the files in file_list in sequence, and copy archived items
        into this archive.  If delete_from_source, the items are
        removed from the source archives.  If delete_when_empty, each
        source file is deleted once all items have been transferred.

        If clobber=False, then an error is raised if this archive
        already contains an item with an imported id.

        """
        # Note this function uses h5py primitives for the data copy,
        # rather than get_item and set_item.  This allows one to use
        # HDFArchive efficiently and generically to mess with any
        # subclass Archive.
        import h5py, os
        friend_mode = {True: 'r+', False: 'r'}[delete_from_source]
        for f in file_list:
            with h5py.File(f, mode=friend_mode) as afriend:
                for k in list(afriend.keys()):
                    if self.has_item(k):
                        if not clobber:
                            raise AssertionError("Key %s exists and clobber=False." % k)
                        self.del_item(k)
                    # This works whether afriend[k] is a group or a dataset.
                    afriend.copy(k, self.target)
                    if delete_from_source:
                        del afriend[k]
                remove_it = (delete_from_source and delete_when_empty and
                             len(list(afriend.keys())) == 0)
            if remove_it:
                print('Deleting empty source %s.' % f)
                os.remove(f)


