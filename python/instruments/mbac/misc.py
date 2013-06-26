from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import os
import numpy as np

from moby2.util import mce

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

ENC_FLAGS = {
    'ENCODER_REPEATED':  0x0001,
    'EXCESSIVE_LOST':    0x0004,
    'ERROR_READING_ENC': 0x0008,
    'SYNC_FAKED':        0x0010,
    'NO_READINGS':       0x0100,
    'BIT_MASK':          0xffff,
}

def get_tod(
    filename=None,
    start=None, end=None,
    det_uid=None,
    rows=None, cols=None, block_read=True,
    read_data=True,
    extract=True,
    data_mode=None,
    runfile=True,
    fix_sign=False,
    cuts=None,
    repair_pointing=False):
    """
    Populate a TOD object from an ACT-style dirfile.
    """
    # Check it.
    if not os.path.exists(filename):
        raise RuntimeError("File %s not found." % filename)

    # Generate channel names
    tod_info = moby2.scripting.get_tod_info({"filename":filename})

    ## If cuts are passed in, use them to set det_uid, start, end.
    if not cuts is None:
        if det_uid is None:
            det_uid = cuts.get_uncut(det_uid=True)
        if start is None:
            start = cuts.sample_offset
        if end is None:
            end = cuts.sample_offset + cuts.nsamps

    ## What detectors do you want?
    if det_uid is None:
        if block_read or rows is None or cols is None:
            det_uid = tod_info.array_data.select_outer(
                {'row': rows, 'col': cols})
        else:
            det_uid = tod_info.array_data.select_inner(
                {'row': rows, 'col': cols})
    ## Make sure det_uid is an array of ints.
    det_uid = np.asarray(det_uid, 'int')
    ## What are their names?
    channels = tod_info.array_data['df_name'][det_uid]

    # Load the runfile?
    if runfile is None or runfile is True:
        runfile = mce.MCERunfile.from_dirfile(filename, prefix='Extras/')
    elif runfile is False:
        runfile = None
    elif isinstance(runfile, basestring):
        runfile = mce.MCERunfile(runfile)
    else:
        pass # Treat it like an MCERunfile object.

    # Interpret runfile
    if runfile is not None and data_mode is None:
        data_mode = runfile.Item('HEADER', 'RB rc1 data_mode',
                                 type='int', array=False)
    if extract in [None, True]:
        extract = 'fb_filt'
    if extract:
        extractors = mce.MCE_data_modes[str(data_mode)]
        if not extract in extractors:
            raise ValueError("These data do not have field '%s'. Available fields: %s" % (
                extract, list(extractors.keys())))
        extractor_basic = extractors[extract]
        extractor = moby2.tod.accelerated_MCE_extractors.get(
            str(data_mode), {}).get(extract, None)
        if extractor is None:
            def f(din, dout):
                dout[:] = extractor_basic.extract(din)
            extractor = f
    else:
        extractor = None

    self = moby2.TOD.init_from_dirfile(
        filename,
        start=start,
        end=end,
        channel_names=channels,
        hk_spec=pointing_fields,
        extractor=extractor,
        read_data=read_data)

    # Update self.det_uid to the real det_uid.
    self.det_uid = det_uid[self.det_uid]

    # Convert our data_mask flags to a boolean mask
    if hasattr(self, 'data_mask'):
        self.data_mask = (self.data_mask == 0)
    else:
        self.data_mask = np.ones(self.nsamps, 'bool')

    # Guess the TODInfo
    self.info = tod_info
    if start is None:
        start = 0
    self.info.sample_index = start
    self.info.det_uid = self.det_uid

    if runfile is not None:
        if 'HEADER' in runfile.data:
            self.info.mce_filter = \
                mce.MCEButterworth.from_runfile(runfile)
        self.info.runfile = runfile

    # Assign some cuts?
    if self.det_uid is not None:
        if cuts is None:
            self.cuts = moby2.TODCuts.for_tod(self)
        else:
            # Note this fails if you override start / end to extend
            # beyond cuts validity.
            self.cuts = cuts.copy(det_uid=self.det_uid).extract(
                start, end-start)

    # Sign correction?
    if fix_sign and read_data:
        cal = tod_info.array_data.get_property('optical_sign', self.det_uid)
        moby2.tod.apply_calibration(self.data,
                                    np.arange(len(cal)),
                                    cal)

    # Fix glitches in scan
    if repair_pointing:
        max_gap_size = 5
        bit_mask = ENC_FLAGS['BIT_MASK'] ^ (
            ENC_FLAGS['ENCODER_REPEATED'] | ENC_FLAGS['SYNC_FAKED'])
        # Use enc_flags == 1 to identify short stretches of invalid encoders.
        usable = ((self.enc_flags & bit_mask) == 0)
        dups = ((self.enc_flags & ENC_FLAGS['ENCODER_REPEATED']) != 0)
        edges = np.hstack([True, usable*~dups, True])
        edges = (edges[1:] != edges[:-1]).nonzero()[0].reshape(-1,2)
        for e0, e1 in edges:
            if e1 - e0 > max_gap_size:
                usable[e0+1:e1] = False
                continue
            x = np.arange(-1, e1-e0+1)
            for v in [self.az, self.alt]:
                v[e0:e1] = np.polyval(np.polyfit(x[[0,-1]], v[[e0-1,e1]], 1), x[1:-1])
        self.pointing_mask = usable
        # Also do this, I guess.
        moby2.tod.repair_pointing(self)
    else:
        # Otherwise, set pointing_mask to distrust all flags except
        # "SYNC_FAKED", which are probably fine.
        bit_mask = ENC_FLAGS['BIT_MASK'] ^ ENC_FLAGS['SYNC_FAKED']
        self.pointing_mask = ((self.enc_flags & bit_mask) == 0)

    return self
