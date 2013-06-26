from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Wrapping functions for dirfile handling.
"""

import moby2
import moby2.libactpol as libactpol
from .stuff import aStruct
import numpy as np

import os, sys

class DirfileManager:
    """
    Container for C-level dirfile data.

    Use this to load fields from a dirfile into a more structured object.

    This object manages non-reference counted C variables.  So don't
    copy its internal data or you will segfault later.
    """

    dirfile = None
    info = None
    
    def __init__(self, filename=None, dirfile=None):
        self.dirfile = dirfile
        if filename is not None:
            self.open(filename)

    def __del__(self):
        self.close()

    def open(self, filename):
        """
        Initialize the DirfileManager from a dirfile.
        """
        self.dirfile = libactpol.dirfile_open(filename)
        if self.dirfile is None:
            raise IOError('Could not open %s as dirfile!\n' % filename)
        self.filename = filename

    def close(self):
        """
        Close the handle to the dirfile.
        """
        if self.dirfile is not None:
            libactpol.dirfile_close(self.dirfile)
            self.dirfile = None

    def get_frame_count(self, channel=None):
        """
        Count frames and samples in the dirfile. Returns an object
        with attributes:

          n_frames:   number of complete frames in the dirfile
          spf:        samples per frame in the given channel
          n_bonus:    number of samples in any final, incomplete frame
          n_samples:  total number of samples, including incomplete frame

        If channel is None, then only n_frames is valid.  The other
        three fields can be obtained (though at the cost of doing a
        trial load of the data) for a particular channel if needed.
        """
        n_frames, spf, n_bonus = libactpol.dirfile_get_frame_count(
            self.dirfile, channel)
        info = aStruct({'n_frames': n_frames,
                        'spf': spf,
                        'n_bonus': n_bonus,
                        'n_samples': n_frames*spf + n_bonus})
        if self.info is None:
            self.info = info
        return info

    def get_channel_info(self, channel):
        """
        Query the existence and samples-per-frame of a channel.
        Returns a structure with attributes:
        
          exists:     True or False depending on whether channel is defined.
          spf:        samples per frame for the channel.
        """
        exists, spf = libactpol.dirfile_get_channel_info(self.dirfile, channel)
        return aStruct({'exists': exists,
                        'spf': spf})

    def load_channel(self, channel, start=0, count=None, dtype=None, output=None):
        """
        Get data for a channel.  The "start" and "count" are in terms
        of samples, not frames.  dtype is a numpy-style type
        specifier; though I bet that only 32-bit ints and 32/64 bit
        floats are supported.
        """
        if count is None:
            if self.info is None:
                self.get_frame_count(channel)
            spf = self.get_channel_info(channel).spf
            count = self.info.n_samples * spf // self.info.spf
            count -= start
        return libactpol.dirfile_load_channel(self.dirfile, channel,
                                              start, count,
                                              dtype, output)
    
    def load_channels(self, channel_list, start=0, count=None,
                      dtype=None, output=None,
                      extractor=None):
        """
        Get data for multiple channel.  The underlying routine
        exploits openmp to parallize the decompression.  The
        extractor, if not None, will be inspected to extract some
        bitfield from the underlying data.  The channel_list should be
        a list of valid field names.  The data are returned as a 2d
        array.  Other arguments are as in .load_channels.
        """
        if len(channel_list) == 0:
            return np.array([], dtype=dtype)
        # Check for spf consistency.
        spf = self.get_channel_info(channel_list[0]).spf
        assert all([self.get_channel_info(c).spf == spf
                    for c in channel_list[1:]])  #spf check.
        if count is None:
            if self.info is None:
                self.get_frame_count(channel_list[0])
            
            count = self.info.n_samples * spf // self.info.spf
            count -= start

        return libactpol.dirfile_load_channels(self.dirfile,
                                               list(channel_list),
                                               start, count,
                                               dtype, output,
                                               extractor)
    
    def load_channels_resampled(self, channel_list, spf0=None,
                                start=0, count=None, 
                                dtype=None, output=None):
        """
        Like load_channels, except that all channels are resampled to
        match the spf0 specified.  This can be used to load several
        channels with different samples-per-frame, projecting them
        into a preferred spf.
        """
        if self.info is None:
            self.get_frame_count(channel_list[0])
        if spf0 is None:
            spf0 = self.info.spf
        if count is None:
            count = self.info.n_samples * spf0 // self.info.spf
            count -= start
            
        # Read all channels of the same SPF with a single call.
        spf_groups = {}
        for ci, c in enumerate(channel_list):
            spf = self.get_channel_info(c).spf
            if not spf in spf_groups:
                spf_groups[spf] = []
            spf_groups[spf].append((c,ci))
        if output is None:
            output = np.empty((len(channel_list), count), dtype)
        for spf, chans in list(spf_groups.items()):
            chan_names = [c[1] for c in chans]
            i_min, i_max = (start * spf // spf0), ((start+count) * spf // spf0)
            data = self.load_channels([c for c,_ in chans],
                                      i_min, i_max - i_min, dtype)
            i_map = np.arange(start, start+count) * spf // spf0 - i_min
            for (c,i),d in zip(chans, data):
                output[i,:] = d[i_map]
        return output
    
    def get_fields(self):
        """
        Read all fields from the format file.  Returns a list of
        tuples of the form (type, name, detail) where 'type' is the
        dirfile entry type (e.g. "raw", "bit", "lincom"), 'name' is
        the name of the field, and detail is some type-specific
        additional information (for "raw", it's the SPF, for "bit"
        it's the source field).
        """
        format_file = os.path.join(self.filename, 'format')
        if os.path.exists(format_file):
            fin = open(format_file)
        else:
            filenamez = self.filename
            if not os.path.exists(self.filename):
                filenamez = self.filename + '.zip'
            import zipfile
            zipf = zipfile.ZipFile(filenamez)
            fin = zipf.open('format')
        fields = []
        for line in fin:
            if type(line) != type(''):
                line = line.decode('utf-8')
            w = line.split()
            if len(w) == 0 or w[0][0] == '#': continue
            name = w[0]
            try:
                ftype = w[1].lower()
                if ftype == 'raw':
                    extras = w[3]
                elif ftype == 'bit':
                    extras = w[2]
                else:
                    extras = None
            except:
                ftype = 'garbled'
                extras = None
            fields.append((ftype, name, extras))
        return fields

    
