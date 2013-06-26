from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
import numpy as np

# Local
import moby2
import moby2.libactpol as libactpol
from moby2.util import DirfileManager, mce, log

from .array_data import ArrayData


DET_DTYPE = 'float32'

pointing_fields = [
    {'name': 'az',  # some 2007 tods don't have Enc_Az_Deg_Astro
     'field': 'Enc_Az_Deg',
     'dtype': np.float64,
     'scale': np.pi / 180, 
     'shift': np.pi,
     },
    {'name': 'alt',
     'field': 'Enc_El_Deg',
     'dtype': np.float64,
     'scale': np.pi / 180, 
     },
    {'name': 'ctime',
     'field': 'C_Time',
     'dtype': np.float64,
     },
    {'name': 'enc_flags',
     'field': 'enc_flags',
     'dtype': np.uint16,
     },
    {'name': 'data_mask',
     'field': 'tes_flags',
     'dtype': np.int32,
     },
]


log_minor = 6
log_major = 1
def trace(level, msg):
    logger = log.get_logger()
    logger.trace(level, msg)

class TOD:

    # Detector data attribues
    data = None             #: The detector data, a numpy array of shape (n_det,n_samp)
    det_uid = None          #: Detector uids; integer array of shape (n_dept)
    nsamps = None           #: n_samp

    #Pointing data attributes.
    az = None               #: Boresight azimuth in radians (n_samp)
    alt = None              #: Boresight altitude in radians (n_samp)
    ctime = None            #: Time, as unix timestamp (n_samp)
    enc_flags = None        #: Pointing validity indicator (n_samp)

    # Associated information
    cuts = None             #: A preferred TODCuts object for the TOD
    info = None             #: A TODInfo object providing misc. TOD metadata
    abuses = None           #: A place to log processing operations

    # Constructors
    def __init__(self, n_samples, n_dets=None,
                 create_data=True,
                 create_pointing=True,
                 info=None):
        """
        The default constructor returns an empty TOD object.  The data
        and pointing vectors are initialized to zero by default.
        """
        if not n_dets is None:
            self.det_uid = np.arange(n_dets)
            if create_data:
                self.data = np.zeros((n_dets, n_samples), DET_DTYPE)
        self.nsamps = n_samples
        if create_pointing:
            for f in pointing_fields:
                setattr(self, f['name'], np.zeros(n_samples, f['dtype']))
        if info is None:
            info = TODInfo()
        self.info = info
        self.abuses = []

    # Copy constructor
    def copy(self,
             copy_data=True,
             copy_pointing=True,
             copy_id=True,
             copy_cuts=True,
             copy_info=True,
             resample=1,
             resample_filter=True,
             resample_offset=None
             ):
        """
        Return a copy of the TOD.  A complete copy is made of all
        detector and pointing data unless the corresponding copy_*
        field is passed in.

        This routine can also resample a TOD.  The TOD is resampled at
        every Nth index, where N is a (possibly non-integral) value
        given by the ``resample`` parameter.  If resample_filter=True,
        then a low-pass filter is used to resample the detector
        data... but that only currently works if resample is a power
        of 2.
        """
        if resample == 1:
            nsamps = self.nsamps
            sample_idx = None
        else:
            if resample_offset is None:
                resample_offset = 0
            nsamps = int(self.nsamps / resample)
            sample_idx = np.arange(nsamps) * resample + resample_offset
            sample_idx = sample_idx.astype('int')
        out = self.__class__(nsamps)

        if copy_data:
            if sample_idx is None:
                out.data = self.data.copy()
            else:
                if resample_filter == True:
                    # Use the simple cos roll-off
                    resamp_iters = np.log(resample) / np.log(2)
                    if abs(resamp_iters - round(resamp_iters)) > 1e-9:
                        raise ValueError("Filtered downsampling only supports powers of 2.")
                    resamp_iters = int(round(resamp_iters))
                    out.data = downsample_filter_simple(self.data, resamp_iters,
                                                        offset=resample_offset)
                else:
                    print('No filtering on downsample!')
                    out.data = np.empty((self.data.shape[0], nsamps), dtype=DET_DTYPE)
                    out.data[:] = self.data[:,sample_idx]
        if self.det_uid is None:
            out.det_uid = None
        else:
            out.det_uid = self.det_uid.copy()

        if copy_pointing:
            if sample_idx is None:
                out.az = self.az.copy()
                out.alt = self.alt.copy()
                out.ctime = self.ctime.copy()
                out.enc_flags = self.enc_flags.copy()
            else:
                out.az = self.az[sample_idx]
                out.alt = self.alt[sample_idx]
                out.ctime = self.ctime[sample_idx]
                out.enc_flags = self.enc_flags[sample_idx]
        if copy_cuts and self.cuts is not None:
            out.cuts = self.cuts.copy(resample=resample)
        if copy_info:
            out.info = self.info.copy(resample=resample)

        return out


    def get_hk(self, fields, filename=None, start=None, count=None,
               dtype='float32', fix_gaps=False):
        """
        Load housekeeping (or other) fields from a dirfile.  Assuming
        that the TOD was loaded using TOD.for_dirfile, this function
        will look to the same dirfile and load the data corresponding
        to the same sample indices.

        "fields" can be a string, or a list of string.  The function
        will return an array or a list of arrays, respectively.

        If you need to override the dirfile name, pass "filename".  If
        you want different starting or stopping indices, use "start"
        and "count".  To manipulate the output data type, set dtype
        (our wrapping of getdata probably only supports float32,
        float64, int32, int64, and maybe some unsigned ints).

        Set fix_gaps=True to attempt patching of non-400 Hz ACTPol HK
        data in cases where the encoder tick dropped out for a while.
        Missing samples will be interpolated.
        """
        if filename is None:
            filename = self.info.filename

        dirfile = DirfileManager(filename)
        frame_info = dirfile.get_frame_count('cpu_s')
        spf0 = frame_info.spf

        if start is not None or count is not None:
            #Ignore this TOD's start and count completely
            nsamps = frame_info.n_samples
            if start is None:
                start = 0
            elif start < 0:
                start = max(0, nsamps + start)
            else:
                start = min(start, nsamps)
            if count is None:
                count = nsamps - start
            elif count < 0:
                count = max(0, nsamps + count - start)
            else:
                count = min(count, nsamps - start)
        else:
            start = self.info.sample_index
            count = self.nsamps

        # Generate sample indices at TOD detector data resolution
        idx0 = (np.arange(count) * self.info.downsample_level).astype('int') + start
        if fix_gaps:
            # You need to see the enc_flags vector from the start.
            enc_flags = dirfile.load_channel(
                'enc_flags', 0, start+count, np.uint16)
            enc_mask = (enc_flags & 0x100)==0
            # Make map from our indices into recorded data (except for
            # 400 Hz channels, which have been patched
            idx_hk = np.zeros(enc_flags.shape, int)
            idx_hk[enc_mask] = np.arange(enc_mask.sum())
            # Now trim the HK lookup and mask to match the target.
            idx_hk = idx_hk[start:]
            enc_mask = enc_mask[start:]

        single_field = isinstance(fields, basestring)
        if single_field:
            fields = [fields]
        data = []
        for f in fields:
            spf = dirfile.get_frame_count(f).spf
            if spf == 0:
                print('field not recognized: %s' % f)
                data.append([None])
            elif not fix_gaps:
                idx = idx0 * spf // spf0
                d = dirfile.load_channel(f, idx[0], idx[-1]+1-idx[0], dtype=dtype)
                data.append(d[idx-idx[0]])
            else: # i.e., fix_gaps
                idx_chan = idx_hk[enc_mask]*spf//spf0
                lo, n = idx_chan[0], idx_chan[-1]-idx_chan[0]+1
                d_in = dirfile.load_channel(f, lo, n,
                                            dtype=dtype)
                d_out = np.zeros(len(enc_mask), 'float32')
                d_out[enc_mask] = d_in[idx_chan-idx_chan[0]]
                # Fill gaps
                cutv = moby2.tod.CutsVector.from_mask(~enc_mask)
                noise = np.zeros((~enc_mask).sum(), 'float32')
                libactpol.fill_cuts(d_out, cutv, noise,
                                    10, False)
                data.append(d_out.astype(dtype))
        if single_field:
            return data[0]
        return data
        

    #
    # I/O
    #

    @classmethod
    def init_from_dirfile(
                     cls, filename=None,
                     start=None, end=None,
                     channel_names=[],
                     hk_spec=[],
                     extractor=None,
                     read_data = True
                     ):
        """
        Create a TOD object with data taken from a (probably
        ACT-style) dirfile.  Since we use libactpol calls, the dirfile
        can be slimmed and zipped.
        """
        # Check it.
        trace(log_minor, 'Loading from %s' % filename)
        if not os.path.exists(filename):
            raise RuntimeError("File %s not found." % filename)

        # Get a handle to the dirfile.
        dirfile = DirfileManager(filename)
        
        # Check that desired detector channels exist.
        trace(log_minor, 'Getting channel info for %i channels' %
              len(channel_names))
        channels_ok = [dirfile.get_channel_info(n).exists for n in channel_names]
        channels_ok = np.nonzero(channels_ok)[0]
        if len(channels_ok) == 0 and len(channel_names) > 0:
            raise IOError("none of the requested tesdata fields were found")
        ## Get frame_count somehow
        frame_count = None
        for count_channel in [channel_names[i] for i in channels_ok] \
                + list(channel_names)[:1] + ['row_len']:
            trace(log_minor, 'Attempting to get frame count from %s...' %
                  count_channel)
            try:
                frame_count = dirfile.get_frame_count(count_channel)
                break
            except:
                continue
        if frame_count is None:
            raise ValueError("Could not determine framecount.")

        n_samples = frame_count.n_samples
        ## Handle start and end arguments
        if start is None:
            start = 0
        if start < 0:
            start = n_samples + start
        start = max(0, start)
        if end is None:
            end = n_samples
        if end < 0:
            end = n_samples + end
        n_samples = min(n_samples, end - start)

        # Create a TOD that can hold those channels
        trace(log_minor, 'Creating TOD object (%i,%s)' %
              (len(channels_ok), n_samples))
        self = cls(n_dets=len(channels_ok), n_samples=n_samples,
                   create_data=False, create_pointing=False)

        trace(log_minor, 'Loading data, samples %i to %i' % 
              (start, start+n_samples))

        # Try to load detector data
        if (len(channels_ok) > 0) and read_data:
            chans = channel_names[channels_ok]
            read_dtype = 'uint32'
            if extractor is not None:
                read_dtype = 'float32'
            self.data = dirfile.load_channels(
                list(channel_names[channels_ok]),
                start, n_samples,
                read_dtype, extractor=extractor)

        # Try to load pointing data
        for info in hk_spec:
            if not dirfile.get_channel_info(info['field']).exists:
                trace(log_major, "Failed to load pointing field '%s'" % 
                      info['field'])
                continue
            trace(log_minor, '...loading %s' % info['field'])
            data = dirfile.load_channel(info['field'], start, n_samples,
                                        info['dtype'], None)
            setattr(self, info['name'], data * info.get('scale', 1) + info.get('shift', 0))

        self.det_uid = channels_ok
        self.info = None

        return self

    @classmethod
    def from_dirfile(cls, *args, **kwargs):
        print('TOD.from_dirfile is (all of a sudden) deprecated.')
        print('Use scripting.get_tod as an alternative.')
        if len(args)>0:
            assert(len(args)==1 and not 'filename' in kwargs)
            kw = {'filename': args[0]}
            kw.update(kwargs)
        else:
            kw = kwargs
        from moby2.instruments import actpol
        return actpol.get_tod(**kw)
            


class TODInfo:
    """
    Container for information that describes the origins and basic
    features of the TOD.  Not all fields will be defined in all cases.
    """

    tod_id = None
    name = None
    basename = None
    filename = None
    runfile = None
    sample_index = 0
    downsample_level = 1
    ctime = None
    season = None
    array = None
    array_data = None
    instrument = None

    def copy(self, resample=1):
        out = self.__class__()
        # Only directly assign immutable items
        for attr in ['name', 'basename', 'filename','runfile',
                     'ctime', 'season', 'array', 'instrument']:
            setattr(out, attr, getattr(self, attr))
        # Then the weird ones
        out.array_data = self.array_data.copy()
        # And the resampling
        out.sample_index = int(self.sample_index / resample)
        out.downsample_level = self.downsample_level * resample
        out.det_uid = self.det_uid.copy()
        return out

    def get_dict(self, det_uid=None):
        """
        Get a dictionary of simple data; convenient for use in
        string formatting.
        """
        odict = {}
        for attr in ['tod_id', 'name', 'basename', 'filename',
                     'ctime', 'season', 'array', 'instrument',
                     'day_str', 'time_str']:
            odict[attr] = getattr(self, attr)
        if det_uid is not None:
            odict['det_uid'] = det_uid
            if not self.array_data is None:
                odict.update(self.array_data.get_info(det_uid))
        return odict
            

def downsample_filter_simple(data_in, n_iter=1, offset=0):
    """
    Do nearest-neighbor remixing to downsample data_in (..., nsamps)
    by a power of 2.
    """
    if n_iter <= 0:
        return None
    ns_in = data_in.shape[-1]
    ns_out = ns_in // 2
    dims = data_in.shape[:-1] + (ns_out,)
    data_out = np.empty(dims, dtype=data_in.dtype)
    # Central sample
    data_out[...,:] = data_in[...,offset:ns_out*2+offset:2] * 2
    # To the left (all output samples except maybe the first)
    l_start = 1-offset
    l_count = ns_out - l_start
    data_out[...,l_start:] += data_in[...,(1-offset):2*l_count:2]
    # To the right (all output samples except maybe the last)
    # depending on 2*ns_out+offset <= ns_in
    r_count = (ns_in - offset) // 2
    data_out[...,:r_count] += data_in[...,offset+1::2]
    # Normalization...
    data_out[...,:] /= 4
    if l_start > 0:
        data_out[...,0] *= 4./3
    if r_count < ns_out:
        data_out[...,-1] *= 4./3
    if n_iter <= 1:
        return data_out
    # Destroy intermediate storage, and iterate
    data_in = data_out
    return downsample_filter_simple(data_in, n_iter-1, offset)
