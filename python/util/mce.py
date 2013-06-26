from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
# Stolen from mce_script/trunk:r1096
#

import numpy
import sys
from os import stat


#
# MCE Flatfile handling
#


# This is what an MCE looks like.
MCE_RC = 4
MCE_COL = 8
MCE_DWORD = 4

# This is how fast it goes
MCE_FREQ = 50e6

# This block read maximum (bytes) is to keep memory usage reasonable.
MAX_READ_SIZE = int(1e9)

class HeaderFormat:
    """
    Contains description of MCE header content and structure.
    """
    def __init__(self):
        self.offsets = {
            'status': 0,
            'frame_counter': 1,
            'row_len': 2,
            'num_rows_reported': 3,
            'data_rate': 4,
            'address0_ctr': 5,
            'header_version': 6,
            'ramp_value': 7,
            'ramp_addr': 8,
            'num_rows': 9,
            'sync_box_num': 10,
            'runfile_id':  11,
            'userfield': 12,
            }
        self.header_size = 43
        self.footer_size = 1

class BitField(object):
    """
    Describes the truncation and packing of a signal into a carrier word.
    """
    def define(self, name, start, count, scale=1., signed=True):
        self.name = name
        self.start = start
        self.count = count
        self.scale = scale
        self.signed = signed
        self.unwrap_period = 2**(self.count)
        return self

    def extract(self, data, rescale=True, unwrap=False, **kwargs):
        """
        Extracts bit field from a numpy array of 32-bit signed integers.
        Assumes a two's complement architecture!
        """
        rescale = deprecate_arg(rescale, kwargs, 'rescale', 'do_scale')
        unwrap = deprecate_arg(unwrap, kwargs, 'unwrap', 'do_unwrap')
        if len(kwargs) > 0:
            raise TypeError("%s: got an expected keyword argument '%s'" % \
                                (sys._getframe().f_code.co_name, list(kwargs.keys())[0]))

        if self.signed:
            # Integer division preserves sign
            right = 32 - self.count
            left = right - self.start
            if left != 0:
                data = numpy.array(data).astype('int32') * 2**left
            if right != 0:
                data = numpy.array(data).astype('int32') / 2**right
        else:
            # For unsigned fields, bit operations should be used
            data = (data >> self.start) & ((1 << self.count)-1)
        if unwrap:
            data = unwrap_array(data, self.unwrap_period)
        if not rescale:
            return data
        return data.astype('float') * self.scale

deprecation_warnings = True

def deprecate_arg(new_val, kwargs, new_arg, old_arg):
    # Note this pops the bad value from kwargs
    if old_arg in kwargs:
        if deprecation_warnings:
            print('Use of argument "%s" is deprecated, the new word is "%s".' % \
                (old_arg, new_arg))
        return kwargs.pop(old_arg)
    return new_val


class DataMode(dict):
    """
    A DataMode consists of a set of BitFields describing signal
    packing into MCE data words.
    """
    def __init__(self):
        dict.__init__(self)
        self.fields = []
        self.raw = False
    def define(self, *args, **kargs):
        for a in args:
            self.fields.append(a.name)
            self[a.name] = a
        for k in list(kargs.keys()):
            if k == 'raw':
                self.raw = True
                self.raw_info = kargs[k]
        return self


#Define the MCE data modes

MCE_data_modes = { \
    '0': DataMode().define(BitField().define('error', 0, 32)),
    '1': DataMode().define(BitField().define('fb', 0, 32, 2.**-12)),
    '2': DataMode().define(BitField().define('fb_filt', 0, 32)),
    '3': DataMode().define(BitField().define('raw', 0, 32),
                           raw={'n_cols':8, 'offsettable': False}),
    '4': DataMode().define(BitField().define('fb', 14, 18),
                           BitField().define('error', 0, 14)),
    '9': DataMode().define(BitField().define('fb_filt', 8, 24, 2.**1),
                           BitField().define('fj', 0, 8)),
    '10': DataMode().define(BitField().define('fb_filt', 7, 25, 2.**3),
                            BitField().define('fj', 0, 7)),
    '11': DataMode().define(BitField().define('row', 3, 7, signed=False),
                            BitField().define('col', 0, 3, signed=False)),
    '12': DataMode().define(BitField().define('raw', 0, 32),
                            raw={'n_cols':1, 'offsettable': True}),
}


class MCEData:
    """
    Container for MCE data (single channel) and associated header and origin information.
    """
    def __init__(self):
        self.data = {}
        self.source = None
        self.start_frame = 0
        self.n_frames = 0
        self.header = None
        self.data_is_dict = True
        self.data = []
        self.channels = []

def _rangify(start, count, n, name='items'):
    """
    Interpret start as an index into n objects; interpret count as a
    number of objects starting at start.  If start is negative,
    correct it to be relative to n.  If count is negative, adjust it
    to be relative to n.
    """
    if start < 0:
        start = n + start
    if start > n:
        print('Warning: %s requested at %i, beyond available %s.' %\
            (name, start, name))
        start = n
    if count is None:
        count = n - start
    if count < 0:
        count = n - start + count
    if start + count > n:
        print('Warning: %i %s requested, exceeding available %s.' %\
            (count, name, name))
        count = n - start
    return start, count
        
    
class SmallMCEFile:
    """
    Facilitate the loading of (single channels from) raw MCE
    flat-files.  Extraction and rescaling of data content is performed
    automatically by default.

    After instantiation with a data filename, a call to Read() will
    return the detector data as an MCEData object.

    See code for 'Reset' method for list of useful attributes.
    """
    def __init__(self, filename=None, runfile=True, basic_info=True):
        """
        Create SmallMCEFile object and load description of the data
        from runfile and header.

        filename: path to MCE flatfile 

        runfile: if True (default), filename.run is used.  If False,
          no runfile is used.  Pass a string here to override the
          runfile filename.

        basic_info: if True (default), basic file information is
          loaded from runfile and frame header.
        """
        # Initialize basic parameters
        self.Reset()        
        self.filename = filename      # Full path to file

        # Set runfile name
        if (filename is not None) and (runfile == True):
            self.runfilename = filename+'.run'
        else:
            self.runfilename = runfile

        # If the time is right, compute content and packing
        if (filename is not None) and basic_info:
            self._GetPayloadInfo()
        if (self.runfilename != False) and basic_info:
            self._ReadRunfile()
            if filename is not None:
                self._GetContentInfo()

    def Reset(self):
        # Describe data sources
        self.filename = None
        self.runfilename = None

        # Describe payload frame (set through _GetPayloadInfo)
        self.n_ro = 0                 # Readout frame count
        self.size_ro = 0              # Readout data payload size, in dwords, per RC.
        self.n_rc = 0                 # Number of RC reporting
        self.rc_step = 0              # RC interleaving stride (due to MAS reordering)
        self.frame_bytes = 0          # Readout frame total size, in bytes.

        # Describe data content (set through _GetContentInfo)
        self.n_frames = 0             # Number of samples per detector
        self.n_rows = 0               # Number of rows stored by RC
        self.n_cols = 0               # Number of cols stored by RC
        self.data_mode = 0            # Data mode of RCs
        self.raw_data = False         # Data is raw mode column data
        self.divid  = 1               # Period (in frames) of CC read queries to
                                      #  RC (i.e. data_rate)
        self.freq = 0.                # Mean sampling frequency, in Hz.

        # Members for storing meta-data
        self.header = None            # Becomes dict of first frame header
        self.runfile = None           # Becomes MCERunfile for this data.

    def _rfMCEParam(self, card, param, array=False, check_sys=True):
        """
        Look up MCE 'card, param' value in runfile.
        """
        data = self.runfile.Item('HEADER', 'RB %s %s'%(card,param), type='int', \
                                 array=array)
        if data is None and check_sys:
            # On SCUBA2, some things are stored in sys only.  Blech.
            data = self.runfile.Item('HEADER', 'RB sys %s'%(param),
                                     type='int', array=array)
            if data is not None:
                data = data[0]
        return data

    def _GetRCAItem(self, param):
        """
        Gets 'rc? <param>' for each RC returning data, warns if the
        setting is not consistent across acq cards, and returns the
        value from the first card.
        """
        rcs = [i+1 for i,p in enumerate(self.header['_rc_present']) if p]
        vals = [ self._rfMCEParam('rc%i'%r, param) for r in rcs ]
        for r,v in zip(rcs[1:], vals[1:]):
            if v is None and vals[0] is not None:
                print('Warning: param \'%s\' not found on rc%i.' % \
                    (param, r))
                continue
            if vals[0] != v:
                print('Warning: param \'%s\' is not consistent accross RCs.' % \
                    (param))
                break
        return vals[0]

    def _GetContentInfo(self):
        """
        Using frame header and runfile, determines how the RC data are
        packed into the CC readout frames.

        Sets members n_cols, n_rows, divid, data_mode, n_frames.
        """
        if self.runfile is None:
            if self.runfilename == False:
                raise RuntimeError('Can\'t determine content params without runfile.')
            self._ReadRunfile()
        # In a pinch we could get these params from the runfile.
        if self.size_ro == 0:
            raise RuntimeError('Can\'t determine content params without data file.')
        # Switch on firmware revision to determine 'num_cols_reported' support
        fw_rev = self._GetRCAItem('fw_rev')
        if fw_rev >= 0x5000001:
            self.n_cols = self._GetRCAItem('num_cols_reported')
            self.n_rows = self._GetRCAItem('num_rows_reported')
        else:
            self.n_cols = MCE_COL
            self.n_rows = self._rfMCEParam('cc', 'num_rows_reported', array=False)
        self.divid = self._rfMCEParam('cc', 'data_rate', array=False)

        # Get data_mode information
        self.data_mode = self._GetRCAItem('data_mode')
        dm_data = MCE_data_modes.get('%i'%self.data_mode)
        if dm_data is None:
            dm_data = MCE_data_modes['0']

        # For 50 MHz modes, the data is entirely contiguous
        if dm_data.raw:
            self.raw_data = True
            self.n_rows = 1
            self.n_cols = dm_data.raw_info['n_cols']
            self.n_frames = self.n_ro * self.size_ro // self.n_cols
            self.freq = MCE_FREQ
            return
            
        # For rectangle modes, check RC/CC packing is reasonable
        count_rc = self.n_rows * self.n_cols
        count_cc = self.size_ro
        
        # Check 1: Warn if count_rc does not fit evenly into count_cc
        if count_cc % count_rc != 0:
            print('Warning: imperfect RC->CC frame packing (%i->%i).' % \
                (count_rc, count_cc))

        # Check 2: Warn if decimation/packing is such that samples are
        #     not evenly spaced in time.
        if count_rc != count_cc:
            if count_rc * self.divid != count_cc:
                print('Warning: bizarro uneven RC->CC frame packing.')
        
        # Determine the final data count, per channel.  Any times
        # that are not represented in all channels are lost.
        self.n_frames = (count_cc // count_rc) * self.n_ro

        # Store mean sampling frequency
        nr, rl, dr = [self._rfMCEParam('cc', s) for s in \
                          ['num_rows', 'row_len', 'data_rate']]
        self.freq = (MCE_FREQ / nr / rl / dr) * (count_cc / count_rc)


    def _GetPayloadInfo(self):
        """
        Determines payload parameters using the data header and file size.

        Sets members n_ro, n_rc, size_ro, frame_bytes, rc_step.
        """
        if self.header is None:
            self._ReadHeader()
        # Compute frame size from header data.
        self.n_rc = self.header['_rc_present'].count(True)
        # Payload size (per-RC) is the product of two numbers:
        #  mult1 'num_rows_reported'
        #  mult2 'num_cols_reported', or 8 in pre v5 firmware
        mult1 = self.header['num_rows_reported']
        mult2 = (self.header['status'] >> 16) & 0xf
        if mult2 == 0:
            mult2 = MCE_COL
        self.rc_step = mult2
        self.size_ro = mult1*mult2
        self.frame_bytes = \
            MCE_DWORD*(self.size_ro * self.n_rc + \
                           self.header['_header_size'] + self.header['_footer_size'])
        # Now stat the file to count the readout frames
        if self.filename is not None:
            # This conditional caginess is for subclassing to MCEBinaryData.
            file_size = stat(self.filename).st_size
            self.n_ro = file_size // self.frame_bytes
            if file_size % self.frame_bytes != 0:
                print('Warning: partial frame at end of file.')

    def _UpdateNFrames(self):
        # Partial GetInfo... no error checking.
        file_size = stat(self.filename).st_size
        self.n_ro = file_size // self.frame_bytes
        count_rc = self.n_rows * self.n_cols
        count_cc = self.size_ro
        self.n_frames = (count_cc // count_rc) * self.n_ro

    def _ReadHeader(self, offset=None, head_binary=None):
        """
        Read the frame header at file position 'offset' (bytes),
        determine its version, and store its data in self.header.
        """
        # It's a V6, or maybe a V7.
        format = HeaderFormat()
        if head_binary is None:
            if self.filename is None:
                raise RuntimeError('Can\'t read header without data file.')
            fin = open(self.filename)
            if offset is not None:
                fin.seek(offset)
            head_binary = numpy.fromfile(file=fin, dtype='<i4',
                                         count=format.header_size)
        # Lookup each offset and store
        self.header = {}
        for k in format.offsets:
            self.header[k] = head_binary[format.offsets[k]]
        # Provide some additional keys to help determine frame size
        self.header['_rc_present'] = [(self.header['status'] & (1 << 10+i))!=0 \
                                          for i in range(MCE_RC)]
        self.header['_header_size'] = format.header_size
        self.header['_footer_size'] = format.footer_size

    def _ReadRunfile(self):
        """
        Load the runfile data into self.runfile_data, using the filename in self.runfile.
        Returns None if object was initialized without runfile=False
        """
        if self.runfilename == False:
            return None
        self.runfile = MCERunfile(self.runfilename)
        return self.runfile

    def ReadRaw(self, count=None, start=0, raw_frames=False):
        """
        Load data as CC output frames.  Most users will prefer the
        Read() method, which decodes the data into detector channels.

        Returns a (frames x dets) array of integers.
        """
        if self.size_ro <= 0:
            self._GetPayloadInfo()
        # Do the count logic, warn user if something is amiss
        start, count = _rangify(start, count, self.n_ro, 'frames')
        # Check max frame size
        if count * self.frame_bytes > MAX_READ_SIZE:
            # Users: override this by changing the value of mce_data.MAX_READ_SIZE
            print('Warning: maximum read of %i bytes exceeded; limiting.' % \
                MAX_READ_SIZE)
            count = MAX_READ_SIZE // self.frame_bytes

        # Open, seek, read.
        f_dwords = self.frame_bytes // MCE_DWORD
        fin = open(self.filename)
        fin.seek(start*self.frame_bytes)
        a = numpy.fromfile(file=fin, dtype='<i4', count=count*f_dwords)
        n_frames = len(a) // f_dwords
        if len(a) != count*f_dwords:
            print('Warning: read problem, only %i of %i requested frames were read.'% \
                  (len(a)//f_dwords, count))
        a.shape = (n_frames, f_dwords)
        if raw_frames:
            # Return all data (i.e. including header and checksum)
            return a
        else:
            # Return the detector data only
            ho = self.header['_header_size']
            return a[:,ho:ho+self.size_ro*self.n_rc]

    def _NameChannels(self, row_col=False):
        """
        Determine MCE rows and columns of channels that are read out
        in this data file.  Return as list of (row, col) tuples.  For
        raw mode data, only a list of columns is returned.
        """
        if self.runfile is None:
            self._ReadRunfile()
        rc_p = self.header['_rc_present']
        rcs = [i for i,p in enumerate(rc_p) if p]

        # Is this raw data?  Special handling.
        dm_data = MCE_data_modes['%i'%self.data_mode]
        if dm_data.raw:
            if dm_data.raw_info['offsettable']:
                offsets = [ self._rfMCEParam('rc%i'%(rc+1), 'readout_col_index') \
                                for rc in rcs ]
            else:
                offsets = [ 0 for rc in rcs]
            return [i + o + r*MCE_COL for i in range(dm_data.raw_info['n_cols']) \
                        for r,o in zip(rcs, offsets)]

        row_index = [ self._rfMCEParam('rc%i'%(r+1), 'readout_row_index') \
                          for r in rcs ]
        col_index = [ self._rfMCEParam('rc%i'%(r+1), 'readout_col_index') \
                          for r in rcs ]
        for i in range(len(rcs)):
            if row_index[i] is None: row_index[i] = 0
            if col_index[i] is None: col_index[i] = 0
        col_index = [ c + r*MCE_COL for r,c in zip(rcs, col_index) ]

        if row_col:
            # Use the row-indexing provided by the first RC.
            rows = [i+row_index[0] for i in range(self.n_rows)]
            # Columns can be non-contiguous
            cols = [col_index[rc]+c for rc in range(self.n_rc)
                    for c in range(self.n_cols) ]
            return (rows, cols)
        else:
            # Assemble the final list by looping in the right order
            names = []
            for row in range(self.n_rows):
                for rc in range(self.n_rc):
                    r = row + row_index[rc]
                    for col in range(self.n_cols):
                        c = col + col_index[rc]
                        names.append((r, c))
        return names
                    
    def _ExtractRect(self, data_in, dtype='float'):
        """
        Given CC data frames, extract RC channel data assuming
        according to data content parameters.
        """
        # Input data should have dimensions (n_cc_frames x self.size_ro*self.n_rc)
        n_ro = data_in.shape[0]

        # Reshape data_in to (cc_frame, cc_row, cc_col) so we can work
        # with each RC's data one-by-one
        data_in.shape = (n_ro, -1, self.n_rc * self.rc_step)

        # Probably should leave the type the same, oops.
        if dtype is None:
            dtype = data_in.dtype

        # Short-hand some critical sizes and declare output data array
        f = self.n_cols*self.n_rows          # RC frame size
        p = self.size_ro // f                # CC/RC packing multiplier
        data = numpy.zeros((self.n_rows, self.n_rc, self.n_cols, n_ro * p),
                           dtype=dtype)

        # The only sane way to do this is one RC at a time
        for rci in range(self.n_rc):
            # Get data from this rc, reshape to (cc_frame, cc_idx)
            x = data_in[:,:,self.rc_step*rci:self.rc_step*(rci+1)].reshape(n_ro, -1)
            # Truncate partial data and reshape to RC frames
            x = x[:,0:f*p].reshape(-1, self.n_rows, self.n_cols)
            # Transpose to (rc_row, rc_col, rc_time) and store
            data[:,rci,:,:] = x.transpose((1,2,0))
        # Just return with one space, one time index.
        return data.reshape(self.n_rc*f, -1)


    def _ExtractRaw(self, data_in, n_cols=8):
        """
        Extract 50 MHz samples from raw frame data.
        """
        # In raw data modes, the RCs always return a perfect set of contiguous data.
        n_samp = data_in.shape[0] * data_in.shape[1] // self.n_rc // n_cols
        data = numpy.zeros((n_cols*self.n_rc, n_samp), dtype='int')

        # Reshape data_in to (cc_frame, cc_row, cc_col) so we can work
        # with each RC's data one-by-one
        data_in.shape = (-1, self.size_ro//self.rc_step, self.n_rc * self.rc_step)
        for rci in range(self.n_rc):
            # Get data from this rc as 1d array.
            x = data_in[:,:,self.rc_step*rci:self.rc_step*(rci+1)].reshape(-1)
            # Truncate partial data and reshape to (rc_sample, column)
            nf = n_cols * (x.shape[0] // n_cols)
            x = x[0:nf].reshape(-1, n_cols)
            # Transpose to (column, rc_sample) and store
            data[n_cols*rci:n_cols*(rci+1),:] = x.transpose()
        return data

    def Read(self, count=None, start=0,
             extract=True, rescale=True, data_mode=None,
             field=None, fields=None, row_col=False,
             raw_frames=False, cc_indices=False,
             n_frames=None, unfilter=False, unwrap=False,
             **kwargs):
        """
        Read MCE data, and optionally extract the MCE signals.

        count       Number of samples to read per channel (default=None,
                    which means all of them).  Negative numbers are taken
                    relative to the end of the file.
        start       Index of first sample to read (default=0).
        extract     If True, extract signal bit-fields using data_mode.  You
                    actually can't turn this off.
        rescale     If True, rescale the extracted bit-fields to match a
                    reference data_mode.
        data_mode   Overrides data_mode from runfile, or can provide data_mode
                    if no runfile is used.
        field       A single field to extract.  The output data will contain
                    an array containing the extracted field.  (If None, the
                    default field is used.)
        fields      A list of fields of interest to extract, or 'all' to get
                    all fields.  Output will contain a dictionary will all
                    extracted fields.  Takes precedence over field argument.
        row_col     If True, detector data is returned as a 3-D array with
                    indices (row, column, frame).
        raw_frames  If True, return a 2d array containing raw data (including
                    header and checksum), with indices (frame, index_in_frame).
        cc_indices  If True, count and start are interpreted as readout frame
                    indices and not sample indices.  Default is False.
        unfilter    If True, deconvolve the MCE low pass filter from filtered
                    data.  If set to 'DC', divide out the DC gain only.
        unwrap      If True, remove effects of digital windowing to restore
                    full dynamic range of signal.  (Only works if fields are
                    extracted.)
        """
        # Deprecated args...
        extract = deprecate_arg(extract, kwargs, 'extract', 'do_extract')
        rescale = deprecate_arg(rescale, kwargs, 'rescale', 'do_scale')
        unwrap = deprecate_arg(unwrap, kwargs, 'unwrap', 'do_unwrap')
        count = deprecate_arg(count, kwargs, 'count', 'n_frames')
        if len(kwargs) > 0:
            raise TypeError("%s: got an expected keyword argument '%s'" % \
                                (sys._getframe().f_code.co_name, list(kwargs.keys())[0]))
        
        # When raw_frames is passed, count and start are passed directly to ReadRaw.
        if raw_frames:
            return self.ReadRaw(count=count, start=start, raw_frames=True)

        # We can only do this if we have a runfile
        if self.n_frames == 0:
            self._GetContentInfo()

        # Allow data_mode override
        if data_mode is not None:
            self.data_mode = data_mode

        if cc_indices:
            start *= pack_factor
            if count is not None:
                count *= pack_factor

        # Decode start and count arguments
        start, count = _rangify(start, count, self.n_frames, 'samples')

        # Convert sample indices to readout frame indices
        if self.raw_data:
            # Raw data is contiguous and uninterrupted
            cc_start = start * self.n_cols // self.size_ro
            cc_count = ((count+start)*self.n_cols + self.size_ro-1) // \
                self.size_ro - cc_start
        else:
            # For packed data, trim excess frame words
            pack_factor = self.size_ro // (self.n_rows * self.n_cols)
            cc_start = start // pack_factor
            cc_count = (count + start + pack_factor-1) // pack_factor - cc_start

        # Get detector data as (n_ro x (size_ro*n_rc)) array
        data_in = self.ReadRaw(count=cc_count, start=cc_start)

        # Check data mode for processing instructions
        dm_data = MCE_data_modes.get('%i'%self.data_mode)
        if dm_data is None:
            print('Warning: unimplemented data mode %i, treating as 0.'%self.data_mode)
            dm_data = MCE_data_modes['0']

        # Handle data packing
        if dm_data.raw:
            # Raw mode data is automatically contiguous
            data = self._ExtractRaw(data_in, dm_data.raw_info['n_cols'])
            # Trim
            offset = start - cc_start * self.size_ro
            data = data[:, offset:offset+count]
        else:
            # Normal/rectangle data may be packed incommensurately.
            data = self._ExtractRect(data_in)
            # Trim to caller's spec
            offset = start - cc_start * pack_factor
            data = data[:, offset:offset+count]

        # Create the output object
        data_out = MCEData()
        data_out.source = self.filename
        data_out.n_frames = data.shape[1]
        data_out.header = self.header

        # Unravel the field= vs. fields=[...] logic
        if field is None:
            field = 'default'
        force_dict = (fields is not None)
        if fields is None:
            fields = [field]
        elif fields == 'all':
            fields = dm_data.fields
        for i,f in enumerate(fields):
            if f=='default':
                fields[i] = dm_data.fields[0]                
        data_out.data_is_dict = (len(fields) > 1 or force_dict)
        if data_out.data_is_dict:
            data_out.data = {}

        # Get the filter description, if we need it
        if unfilter and ('fb_filt' in fields):
            filt = MCEButterworth.from_runfile(self.runfile)

        # Extract each field and store
        for f in fields:
            # Use BitField.extract to get each field
            new_data = dm_data[f].extract(data, rescale=rescale, unwrap=unwrap)
            if row_col:
                new_data.shape = (self.n_rows, self.n_cols*self.n_rc, -1)
            # Filter?
            if f == 'fb_filt':
                if unfilter == 'DC':
                    new_data /= filt.gain()
                elif unfilter == True:
                    new_data = filt.apply_filter(new_data, inverse=True,
                                                 decimation=1./self.divid)
                elif unfilter == False:
                    pass
                else:
                    raise ValueError("unexpected value for unfilter= argument to MCEFile.Read")
            if data_out.data_is_dict:
                data_out.data[f] = new_data
            else:
                data_out.data = new_data

        data_out.channels = self._NameChannels(row_col=row_col)
        return data_out

# Let's just hope that your MCEFile is a Small One.
MCEFile = SmallMCEFile


#
# MCE Runfile handling
#

class BadRunfile(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MCERunfile:
    def __init__(self, filename=None):
        self.filename = filename
        self.data = {}
        if filename is not None:
            self.Read(filename)

    def Read(self, filename):
        if isinstance(filename, basestring):
            f = open(filename, "r")
        else:
            f = filename

        # Make sure lines are strings; this is a python2/3 tension.
        # This isn't a problem when opening a normal file in mode 'r'
        # but Zipfile seems to open for binary read by default.
        lines = f.readlines()
        if len(lines) and type(lines[0]) != type(''):
            lines = [line.decode('utf-8') for line in lines]

        block_name = None
        block_data = {}
        self.data = {}

        for l in lines:
            key, data = runfile_break(l)
            if key is None: continue

            if key[0] == '/':
                if (block_name != key[1:]):
                    raise BadRunfile('closing tag out of place')
                if data != '':
                    raise BadRunfile('closing tag carries data')
                self.data[block_name] = block_data
                block_name = None
                block_data = {}
            elif block_name is None:
                if data is None or data == '':
                    if key in self.data:
                        raise BadRunfile('duplicate block \'%s\''%key)
                    block_name = key
                else:
                    raise BadRunfile('key outside of block!')
            else:
                block_data[key] = data
        return self.data
    
    def Item(self, block, key, array=True, type='string'):
        if block not in self.data or key not in self.data[block]:
            return None
        data = self.data[block][key]
        if type=='float':
            f = [float(s) for s in data.split()]
            if not array and len(f) <= 1: return f[0]
            return f
        if type=='int':
            f = [int(s) for s in data.split()]
            if not array and len(f) <= 1: return f[0]
            return f
        if type!='string':
            print('Unknown type "%s", returning string.' % type)
        if array:
            return data.split()
        return data

    def Item2d(self, block, key_format, array=True, type='string',
               first = 0, count = None):
        done = False
        result = []
        row = first
        while not done:
            g = self.Item(block, key_format % row, array=array, type=type)
            if g is None:
                break
            result.append(g)
            row = row + 1
            if count is not None and row - first == count:
                break
        return result

    def Item2dRC(self, block, key_format, array=True, type='string',
                 first = 0, count = None, rc_count=4, rc_start=1):
        rc_data = []
        for i in range(4):
            d = self.Item2d(block, key_format%(i+1), array=array,
                            type=type, first=first, count=count)
            if d is None:
                return None
            for column in d:
                rc_data.append(column)
        return rc_data

    def ReadoutFilter(self):
        return MCEFilter.from_runfile(self)

    def __getitem__(self, key):
        return self.data[key]

    @classmethod
    def from_dirfile(cls, filename, prefix=''):
        """
        Find runfile in dirfile, load it, return it.
        """
        import glob, zipfile
        if zipfile.is_zipfile(filename):
            # Locate the runfile
            zf = zipfile.ZipFile(filename)
            infos = zf.infolist()
            for rfi in infos:
                if rfi.filename.startswith(prefix) and \
                        rfi.filename.endswith('.run'):
                    runfile = zf.open(rfi.filename)
                    break
            else:
                runfile = None
            # We found it; open as a file.
        else:
            # Well it must be an ordo then.
            runfile = glob.glob(filename + '/' + prefix + '*.run')
            if len(runfile) > 0:
                runfile = open(runfile[0])
            else:
                runfile = None
        # Did any of that work?
        if runfile is None:
            raise IOError('Could not find Extras/*.run; try runfile=False.')
        return MCERunfile(runfile)


def runfile_break(s):
    reform = ' '.join(s.split())
    words = reform.split('>')
    n_words = len(words)
    
    if n_words == 0 or words[0] == '':
        return None, None

    if words[0][0] == '#':
        return None, None

    if words[0][0] != '<':
        raise BadRunfile(s)
    
    key = words[0][1:]
    data = ' '.join(words[1:])

    return key, data


#
# Functions for processing squid data
#

# unwrap -- this is a low-level portable version.  We should have a
# second routine that interprets the MCEData structure directly, to
# lookup period from data_mode or whatever.

def unwrap_array(data, period, in_place=False):
    """
    Removes jumps (due to fixed bit-width windowing, or something)
    from a data array.  "data" should be an array of signed values,
    with possible jumps of size "period" in its right-most dimension.

    With in_place=True, the data is unwrapped in place, rather than
    creating a new buffer for the unwrapped data.
    """
    ddata = data[...,1:] - data[...,:-1]
    ups = (ddata >  period/2).astype('int').cumsum(axis=-1)
    dns = (ddata < -period/2).astype('int').cumsum(axis=-1)
    if not in_place:
        data = data.astype('float')
    data[...,1:] += float(period) * (dns - ups)
    return data

def unwrap(*args, **kwargs):
    """
    This in a alias for unwrap_array, which you should use now.
    """
    if deprecation_warnings:
        print('Use of "unwrap" function is deprecated, the new name '\
            ' is "unwrap_array".')
    return unwrap_array(*args, **kwargs)

#
# MCE low-pass filters
#

class MCEFilter:
    @staticmethod
    def from_runfile(runfile):
        """
        Return a filter object based on the runfile information.
        """
        # It's probably a Butterworth.
        return MCEButterworth.from_runfile(runfile)

class MCEButterworth(MCEFilter):
    def __init__(self, params, f_samp=None):
        """
        Initialize with a list of 6 parameters, corresponding to 4
        coefficients and two gain magnitudes.
        """
        self.params = params
        self.f_samp = f_samp
    
    def transfer(self, f, f_samp=None, power=False):
        """
        Return filter transfer function at frequencies f.
        
        f is the array of frequencies at which to evaluate the response.

        f_samp is the sampling frequency.

        Setting power=True will return the power window function
        (square of the modulus of the transfer function).
        """
        if f_samp is None:
            f_samp = self.f_samp
        f = f / f_samp
        K = 1./2**14
        scalars = [K, K, K, K, 1., 1.]
        b11, b12, b21, b22, k1, k2 = [s*p for s,p in zip(scalars, self.params)]
        z = numpy.exp(-2j*numpy.pi*f)
        H = (1. + z)**4 / (1. - b11*z + b12*z**2) / (1. - b21*z + b22*z**2)
        H /= 2**(k1+k2)
        if power:
            return abs(H)**2
        return H

    def spectrum(self, *args, **kwargs):
        print('*** please use "transfer" method instead of "spectrum" method.')
        return self.transfer(*args, **kwargs)

    def gain(self):
        """
        Estimate the DC gain of the filter.
        """
        return self.transfer(0).real

    def f3dB(self, cutoff=0.5, f_samp=None):
        """
        Estimate the frequency at which the filter attenuates half of
        the signal power (relative to DC).
        """
        from scipy.optimize import fmin
        if f_samp is None:
            f_samp = self.f_samp
        g0 = self.gain()
        def _spec(x):
            if x.ndim > 0: x = x[0]
            if x > 0.5:
                x = 0.5 - x #flip
            if x < 0:
                return (1.-x) * g0
            return abs(cutoff - abs(self.transfer(x)/g0)**2)
        return fmin(_spec,0.1,disp=0)[0] * f_samp

    # Filter application

    def apply_filter(self, data, decimation=1., inverse=False,
                     gain0=None):
        """
        Apply or de-apply filter to the last dimension of array
        "data", using Fourier representation.

        decimation      ratio of data frequency (e.g. 400 Hz) to the
                        internal sampling frequency (e.g. 15151 Hz).

        inverse         apply inverse of filter (deconvolve its
                        effects)

        gain0           by default the application / deapplication
                        will include the DC gain of the filter.  If
                        instead you want the DC gain to be one, pass
                        gain0=1.
        """
        n = data.shape[-1]
        freqs = numpy.arange(float(n))/n
        freqs[int((n+1)/2):] -= 1.
        spec = self.transfer(freqs, f_samp=1./decimation)
        if gain0 is not None:
            spec *= gain0 / self.gain()
        if inverse:
            spec = 1./spec
        return numpy.fft.ifft(numpy.fft.fft(data)*spec).real

    def apply_filter_fir(self, data, truncate=False,
                         stages=None):
        """
        Apply filter to data by applying the discrete-time filter.

        truncate        If true, intermediate calculations are
                        truncated as they would be in the MCE's fixed
                        point implementation.  This allows for complete
                        simulation of digital artifacts.
        """
        import scipy.signal as scs
        # Special hack
        n = data.shape[-1]
        b = [1., 2., 1.]
        # First filter
        if stages is None or 0 in stages:
            a = [1., -self.params[0]/2.**14, self.params[1]/2.**14]
            data = scs.lfilter(b, a, data) / 2**self.params[5]
            if truncate:
                data = numpy.floor(data)
        # Second filter
        if stages is None or 1 in stages:
            a = [1., -self.params[2]/2.**14, self.params[3]/2.**14]
            data = scs.lfilter(b, a, data) / 2**self.params[4]
            if truncate:
                data = numpy.floor(data)
        return data

    # Static methods for convenient construction

    @classmethod
    def from_params(cls, ftype, fparams, f_samp=None):
        params = None
        if ftype is None or ftype == 0 or ftype == 1:
            # Classic filter
            params = [32092, 15750, 31238, 14895, 0, 11]
        elif ftype == 2:
            # That other hard-coded filter
            params = [32295, 15915, 32568, 16188, 3, 14]
        elif ftype == 255:
            # Parametrized
            params = fparams
        # Did this all work out?
        if params is None or len(params) != 6:
            raise ValueError("Invalid filter parameters for ftype='%i'" %\
                ftype)
        return cls(params, f_samp)

    @classmethod
    def from_runfile(cls, runfile):
        """
        Parses an MCERunfile for filter parameters and returns an
        MCEFilter.
        """
        if isinstance(runfile, str):
            # I guess it's a filename
            runfile = MCERunfile(runfile)
        # Preferred readout card
        rc = runfile.Item('FRAMEACQ', 'RC')
        if rc is None:
            rc = 'rc1'
        else:
            rc = 'rc' + rc[0]
        # Maybe there is a filter description
        ftype = runfile.Item('HEADER', 'RB %s fltr_type' % rc,
                             type='int', array=False)
        fparams = runfile.Item('HEADER', 'RB %s fltr_coeff' % rc, type='int')
        # Store the mux frequency too
        nr, rl = [runfile.Item('HEADER', 'RB cc %s' % s, type='int', array=False)
                  for s in ['num_rows', 'row_len']]
        # 2007: load from sys.
        if nr is None:
            nr, rl = [runfile.Item('HEADER', 'RB sys %s' % s, type='int', array=True)[0]
                      for s in ['num_rows', 'row_len']]
        f_samp = MCE_FREQ / nr / rl
        # That should be enough
        return cls.from_params(ftype, fparams, f_samp)

