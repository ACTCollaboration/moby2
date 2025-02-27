from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
import numpy as np

import moby2
import moby2.libactpol as libactpol

import moby2.util.log as psLib

CUTS_DTYPE = 'int32'

class CutsVector(np.ndarray):
    """
    Represents the cut samples in a single time-stream.  This class
    extends a numpy array, enforcing a particular shape and data type.

    self.nsamps gives the number of samples in the data to which the
    cuts apply.

    self[:,:] is a (ncuts, 2) array of indices.  self[i,0] and
    self[i,1] are the index of the beginning and end of a cut region.
    """
    def __new__(cls, cuts_in=None, nsamps=None, ncuts=0, always_cut=False):
        if cuts_in is not None:
            ncuts = len(cuts_in)
        self = np.ndarray.__new__(cls, shape=(ncuts,2), dtype=CUTS_DTYPE)
        self.nsamps = nsamps
        if cuts_in is None:
            self[:] = 0
        elif len(cuts_in) > 0:
            self[:] = cuts_in
        return self

    def get_mask(self, nsamps=None, invert=False):
        """
        Return an array of booleans of length self.nsamps, with cut
        samples flagged as "True".  Pass invert=True to get not(this).

        The complementary routine is the classmethod from_mask.
        """
        if nsamps is None:
            nsamps = self.nsamps
        return libactpol.mask_from_cuts(self, nsamps, invert)

    def get_complement(self):
        """
        Get a cuts vector that is the complement of this one, in the
        sense that all uncut samples are cut and vice-versa.
        """
        if self.shape[0] == 0:
            return self.__class__([[0,self.nsamps]])
        # Naive complement includes possible trivial regions at start and end of vector.
        complement = np.empty((self.shape[0]+1, 2), CUTS_DTYPE)
        complement[ 0,:] = (0, self[0,0])
        complement[-1,:] = (self[-1,1], self.nsamps)
        complement[1:-1] = self.ravel()[1:-1].reshape(-1,2)
        # Trim trivial regions
        start_i = int(complement[0,0] == complement[0,1])
        end_i = complement.shape[0] - int(complement[-1,0] == complement[-1,1])
        return self.__class__(complement[start_i:end_i], nsamps=self.nsamps)

    def get_collapsed(self):
        """
        Return a CutsVector based on this one, with trivial regions
        and abutting indices removed.
        """
        if len(self) == 0:
            return self
        # Collapse abutting
        keepers = np.ones(self.shape, 'bool')
        gap_size = self[1:,0] - self[:-1,1]
        keepers[ 1:,0] *= (gap_size > 0)
        keepers[:-1,1] *= (gap_size > 0)
        # Remove trivial
        collapsed = self.ravel()[keepers.ravel()].reshape(-1,2)
        collapsed[:,0][collapsed[:,0]<0] = 0
        if self.nsamps is not None:
            collapsed[:,1][collapsed[:,1]>self.nsamps] = self.nsamps
        non_trivial = collapsed[:,0] < collapsed[:,1]
        return self.__class__(collapsed[non_trivial,:], self.nsamps)

    def get_resampled(self, resample=1):
        return self.__class__((self / resample).astype(CUTS_DTYPE), self.nsamps)

    def get_buffered(self, left, right=None):
        """
        Pad each cut region.  If only left=(int) is passed, the cuts
        are buffered on both sides.
        """
        if len(self) == 0:
            return self.__class__(nsamps=self.nsamps)
        if right is None:
            right = left
        new_vect = self + [[-left, right]]
        return self.__class__(new_vect, self.nsamps).get_collapsed()

    def _encode_c(self):
        assert (self.dtype == CUTS_DTYPE) and \
            (self.flags['C_CONTIGUOUS']) and \
            (self.ndim == 2) and \
            (self.shape[1] == 2)
        return self

    @classmethod
    def from_mask(cls, mask):
        mask = np.asarray(mask, dtype='bool', order='C')
        cuts = libactpol.cuts_from_mask(mask)
        return cls(cuts_in=cuts, nsamps=mask.shape[-1])

    @classmethod
    def new_always_cut(cls, nsamps):
        return cls([[0,nsamps]], nsamps)


def _encode_cuts_single(x, new=False):
    if len(x) == 0:
        return _new_cuts()
    y = np.asarray(x, dtype=CUTS_DTYPE, order='C')
    if new and (y is x):
        y = np.array(x)
    if not y.ndim == 2 and y.shape[1] == 2:
        raise ValueError("Cuts array does not have required shape")
    return y

def _new_cuts():
    return np.empty((0,2), dtype=CUTS_DTYPE)

# Decorator to iterate operations on all dets.
def iterate_dets(f):
    f0 = f
    def iterated(self, det, *args, **kwargs):
        if not hasattr(det, '__len__') or  \
                (hasattr(det, 'ndim') and det.ndim == 0):
            return f0(self, det, *args, **kwargs)
        return [f0(self, d, *args, **kwargs) for d in det]
    # Help, I need somebody's help.
    iterated.__doc__ =  f.__doc__
    return iterated


class TODCuts(object):
    cuts = None
    nsamps = None
    sample_offset = None
    det_uid = None

    _depot_structure = '{class}/{tag}/{first_five}/{tod_name}.cuts'

    def __init__(self, ndets=None, nsamps=None, det_uid=None,
                 sample_offset=None):
        if ndets is None:
            ndets = len(det_uid)
        if det_uid is None:
            det_uid = np.arange(ndets)
#        self.cuts = [_new_cuts() for i in range(ndets)]
        self.cuts = [CutsVector([], nsamps) for i in range(ndets)]
        self.nsamps = nsamps
        self.sample_offset = sample_offset
        self.det_uid = np.array(det_uid)

    def copy(self, resample=1, det_uid=None, cut_missing=False):
        """
        copy constructor... can do a bit more than copy.  Passing resample
        """
        nsamps = int(self.nsamps * resample)
        if self.sample_offset is None:
            sample_offset = None
        else:
            sample_offset = int(self.sample_offset * resample)
        if det_uid is None:
            # Keep the same det_uid.
            det_uid = self.det_uid
            det_idx = np.arange(len(det_uid))
        else:
            det_idx = np.zeros(len(det_uid), int) - 1
            for i,u in enumerate(det_uid):
                try:
                    det_idx[i] = list(self.det_uid).index(u)
                except ValueError:
                    if not cut_missing:
                        raise RuntimeError("If you want to promote the cuts to include "
                                             "absent detectors, you should specify "
                                             "cut_missing=True")
        out = self.__class__(det_uid=det_uid, nsamps=nsamps,
                             sample_offset=sample_offset)
        
        out.cuts = []
        for i in det_idx:
            if i < 0:
                out.cuts.append(CutsVector.new_always_cut(out.nsamps))
            else:
                out.cuts.append(self.cuts[i].get_resampled(resample))
        return out

    def extract(self, sample_offset, nsamps, cut_missing=False):
        """
        Extract a TODCuts from this one to match a particular range of
        indices.  For example, if you loaded a TOD using

           tod = moby2.TOD.from_dirfile(filename, start=START, end=START+1000)

        Then to resize and re-index your cuts object you should call:

           tod_cuts = tod_cuts0.extract(tod.info.sample_index, tod.nsamps)
           
        Note that the "sample_offset" is relative to the initial load
        of the cuts, not to the most recent "extract" call.  So a call
        sequence like this:

           cuts1 = cuts0.extract(1000, 216244)
           cuts2 = cuts1.extract(1000, 216244)

        will actually leave you with the same exact thing in cuts1 and
        cuts2.
        """
        # How far forward to start the extraction
        delta = sample_offset - self.sample_offset
        # Make sure extraction is in bounds, or at least cut oob samples.
        pad_left = pad_right = np.zeros((0,2),int)
        if delta < 0 or sample_offset+nsamps > self.sample_offset+self.nsamps:
            if not cut_missing:
                raise ValueError("CutsVector superset not permitted "\
                    "unless cut_missing=True.")
            if delta < 0:
                pad_left = [(0,-delta)]
            if sample_offset+nsamps > self.sample_offset+self.nsamps:
                pad_right = [(self.nsamps + self.sample_offset - sample_offset,
                             nsamps)]
        cvects = []
        for c in self.cuts:
            c1 = np.vstack((pad_left, c-delta, pad_right))
            ce = CutsVector(cuts_in=c1, nsamps=nsamps).get_collapsed()
            cvects.append(ce)
        out = self.__class__(det_uid=self.det_uid, nsamps=nsamps,
                             sample_offset=sample_offset)
        out.cuts = cvects
        return out

    def extract_for_tod(self, tod):
        """Return a version of these cuts that has been aligned to the
        sample_offset, nsamps, and det_uid of tod.  The cuts returned
        from this call should go through fill_cuts without problems.
        """
        return self.extract(tod.cuts.sample_offset,
                            tod.cuts.nsamps).copy(det_uid=tod.det_uid)

    def get_mask(self):
        """
        Return boolean mask identifying which detectors are not completely cut.
        """
        return np.array([not self.is_always_cut(i) for i in range(len(self.cuts))])

    def get_cut(self, det_uid=False):
        if det_uid:
            return self.det_uid[self.get_cut(det_uid=False)]
        return (~self.get_mask()).nonzero()[0]

    def get_uncut(self, det_uid=False):
        if det_uid:
            return self.det_uid[self.get_uncut(det_uid=False)]
        return (self.get_mask()).nonzero()[0]

    @iterate_dets
    def is_always_cut(self, det):
        if self.nsamps is None:
            raise ValueError('nsamps not set for cuts object.')
        c = self.cuts[det]
        if len(c) == 0 or len(c) > 1:
            return False
        return c[0,0] == 0 and c[0,1] >= self.nsamps

    @iterate_dets
    def set_always_cut(self, det_idx):
        if self.nsamps is None:
            raise ValueError('nsamps not set for cuts object.')
        self.cuts[det_idx] = CutsVector.new_always_cut(self.nsamps)
    
    @iterate_dets
    def add_cuts(self, det_index, new_cuts, mask=False):
        """
        Merge new_cuts (a CutsVector, equivalent array, or boolean
        mask) into the cuts for channel det_index.
        """
        if mask:
            new_cuts = CutsVector.from_mask(new_cuts)
        new_cuts = _encode_cuts_single(new_cuts)
        new_cuts = libactpol.merge_cuts([self.cuts[det_index], new_cuts])
        self.cuts[det_index] = CutsVector(new_cuts, self.nsamps)

    def merge_tod_cuts(self, tod_cuts, cut_missing=False):
        # Match the det_uid of self and tod_cuts
        tod_cuts = tod_cuts.extract(self.sample_offset, self.nsamps,
                                    cut_missing=cut_missing)
        assert(self.sample_offset == tod_cuts.sample_offset)
        if len(tod_cuts.det_uid) == 0:
            return
        new_uid = list(tod_cuts.det_uid)
        for i0, det_uid in enumerate(self.det_uid):
            if det_uid in new_uid:
                i1 = new_uid.index(det_uid)
                self.add_cuts(i0, tod_cuts.cuts[i1])
            elif cut_missing:
                self.set_always_cut(i0)

    def get_complement(self):
        out = self.__class__(det_uid=self.det_uid, nsamps=self.nsamps,
                             sample_offset=self.sample_offset)
        out.cuts = []
        for c in self.cuts:
            out.cuts.append(c.get_complement())
        return out
        
    def buffer(self, nbuf):
        for i, c in enumerate(self.cuts):
            self.cuts[i] = c.get_buffered(nbuf)

    def _encode_c(self, dets=None):
        cuts = [c._encode_c() for c in self.cuts]
        if dets is not None:
            if dets.dtype == 'bool':
                dets = dets.nonzero()[0]
            cuts = [cuts[i] for i in dets]
        return cuts

    """
    Reductions
    """

    def downsample(self, ndown):
        new_self = self.__class__(det_uid=self.det_uid.copy(), nsamps=self.nsamps)
        new_self.cuts = [c.resample(2**down) for c in self.cuts]
        return new_self


    """
    I/O
    """

    def write_act_cuts_file(self, filename):
        nr, nc = 33, 32
        fout = open(filename, 'w')
        fout.write('%i %i\n\n' % (nr, nc))
        for c in range(nc):
            for r in range(nr):
                uid = r*nc + c
                try:
                    i = list(self.det_uid).index(uid)
                except:
                    continue
#        for i, det_uid in enumerate(self.det_uid):
#            c, r = (det_uid % nc), (det_uid // nc)
                if len(self.cuts[i])==0:
                    continue
                fout.write('r%02ic%02i: ' % (r, c))
                # Orig moby segfaults if no trailing space at EOL
                for cut in self.cuts[i]:
                    fout.write('(%i,%i) ' % tuple(cut))
                fout.write('\n')

    def write_actpol_cuts_file(self, filename):
        fout = open(filename, 'w')
        fout.write('format = \'TODCuts\'\n'
                   'format_version = 2\n'
                   'n_det = %i\n' % (self.det_uid.max()+1))
        if self.nsamps is not None:
            fout.write('n_samp = %i\n' % self.nsamps)
        if self.sample_offset is not None:
            fout.write('samp_offset = %i\n' % self.sample_offset)
        fout.write('END\n')

        for i,uid in enumerate(self.det_uid):
            fout.write('%4i:' % (uid))
            if len(self.cuts[i])!=0:
                for cut in self.cuts[i]:
                    fout.write(' (%i,%i)' % tuple(cut))
            fout.write('\n')

    write = write_actpol_cuts_file

    write_to_path = write

    @classmethod
    def from_act_cuts_file(cls, filename, nsamps=None):
        """
        Load an ACT-style cuts file, for MBAC.
        """
        if isinstance(filename, basestring):
            filename = open(filename)
        # Improvise on lack of nsamps...
        if nsamps is None:
            nsamps = 2**31-1
        nrow, ncol = list(map(int, filename.readline().split()))
        self = cls(nrow*ncol, nsamps=nsamps)
        for line in filename:
            # process lines like: "r06c00: (204292,204354) (205402,205478)"
            w = line.split(':')
            if len(w) < 2:
                continue
            r, c = int(w[0][1:3]), int(w[0][4:6])
            pairs = w[1].split()
            pairs = [list(map(int, p[1:-1].split(','))) for p in pairs]
            self.cuts[r*ncol+c] = CutsVector(pairs, nsamps)
        return self

    @classmethod
    def from_actpol_cuts_file(cls, filename):
        """
        Load an ACTpol style cuts file.
        """
        if isinstance(filename, basestring):
            filename = open(filename)
        
        # Interpret header as key = value pairs.
        header = {}
        for line in filename:
            tokens = line.split('=')
            if len(tokens) == 0 or tokens[0][0] == '#':
                continue
            if tokens[0].strip() == 'END':
                break
            if len(tokens) != 2:
                raise IOError("Invalid cuts file format")
            key, value = tokens[0].strip(), tokens[1].strip()
            value = eval(value,{"__builtins__":None},{})
            header[key] = value
        else:
            raise IOError("Invalid cuts file format, header END not found.")

        self = cls(header['n_det'],
                   nsamps=header.get('n_samp'),
                   sample_offset=header.get('samp_offset'))

        if header['format_version'] in [1,2]:
            for line in filename:
                # process lines like: "   6: (204292,204354) (205402,205478)"
                # or: "   6 r06c00: (204292,204354) (205402,205478)"
                w = line.split(':')
                if len(w) < 2:
                    continue
                det_uid = int(w[0].strip().split()[0])
                pairs = w[1].split()
                pairs = [list(map(int, p[1:-1].split(','))) for p in pairs]
                self.cuts[det_uid] = CutsVector(pairs, self.nsamps)
        else:
            raise ValueError("Unknown format_version=%s" % repr(
                header['format_version']))
        return self

    @classmethod
    def read_from_path(cls, filename, tod=None):
        return cls.from_actpol_cuts_file(filename)

    @classmethod
    def for_tod(cls, tod, assign=True):
        self = cls(nsamps=tod.nsamps, det_uid=tod.det_uid,
                   sample_offset=tod.info.sample_index)
        if assign:
            tod.cuts = self
        return self


def get_constant_det_cuts(tod):
    cuts = TODCuts.for_tod(tod, assign=False)
    for i,d in enumerate(tod.data):
        if d.max() == d.min():
            cuts.set_always_cut(i)
    return cuts


def get_mce_cuts(tod):
    cuts = TODCuts.for_tod(tod, assign=False)
    for d in tod.det_uid:
        cuts.add_cuts(d, ~tod.data_mask, mask=True)
    return cuts
        
    

# def get_multiple_sources_cuts(tod, sources,
#                                map_pix=None,
#                                radius=10./60, offset=(0.,0.)):
#     """
#     Return TODCuts for tod that excise a circle around each source provided in a sources list.
    
#     Sources should be a list of (ra, dec) in radians
#     All other parameters have
#     units of degrees.
#     """
#     # Create a map centered on tod
#     wand = moby2.pointing.ArrayWand.for_tod(tod, coords='ra_dec')
#     ra, dec = wand.get_coords(tod.fplane)
#     ra = np.rad2deg(ra)
#     dec = np.rad2deg(dec)
#     if map_pix is None: map_pix = 0.5 / 60
#     maskMap = moby2.mapping.fits_map.spaceMap.simpleMap(
#         (-ra.min(), -ra.max()), (dec.min(), dec.max()), (map_pix,map_pix), dtype='float32', wtype='int32')
#     x0, y0 = offset
#     maskMap.x = maskMap.x - x0
#     maskMap.y = maskMap.y - x0

#     for source in sources:
#         r, d = np.rad2deg(source[1]), np.rad2deg(source[2])
#         my_mask = ((maskMap.x-(-r))**2 + (maskMap.y-d)**2 < radius**2)
#         maskMap.data[my_mask] = 1.
#     ny, nx = maskMap.data.shape
#     # psLib.trace('moby', 3, 'get_source_cuts map has %i x %i pixels' % (ny,nx))
#     # psLib.trace('moby', 3, ' map mask fraction is %f' % (float(maskMap.data.sum())/(nx*ny)))
#     # Project into TOD
#     gridding = moby2.pointing.GridPixelization.forFitsMap(maskMap)
#     proj = moby2.pointing.WandProjector(wand, gridding, tod.fplane)
#     # Warning this creates TOD-sized data...
#     mask_vect = proj.deproject_from_map(maskMap.data)
#     # Get cuts that is the same dets as focal plane entries.
#     cuts = moby2.TODCuts(nsamps=tod.nsamps, det_uid=tod.fplane.det_uid,
#                          sample_offset=tod.info.sample_index)
#     for i,d in enumerate(mask_vect):
#         cuts.add_cuts(i, d, mask=True)
#     return cuts


def get_source_cuts(tod, ra, dec,
                    map_size=None, map_pix=None,
                    radius=10./60, offset=(0., 0.)):
    """
    Return TODCuts for tod that excise a circle around (ra, dec).
    
    ra, dec should be provided in radians.  All other parameters have
    units of degrees.
    """
    if map_size is None:
        map_size = radius * 2
    if map_pix is None:
        map_pix = .003 * map_size
    # Create a small map centered on source
    g, dg = map_size, map_pix
    maskMap = moby2.mapping.fits_map.spaceMap.simpleMap(
        (-g, g), (-g, g), (dg,dg), dtype='float32', wtype='int32')
    x0, y0 = offset
    dx, dy = maskMap.x - x0, maskMap.y - y0
    my_mask = (dx**2 + dy**2 < radius**2)
    maskMap.data[my_mask] = 1.
    ny, nx = maskMap.data.shape
    psLib.trace('moby', 3, 'get_source_cuts map has %i x %i pixels' % (ny,nx))
    psLib.trace('moby', 3, ' centered at %f %f' % (x0, y0))
    psLib.trace('moby', 3, ' map mask fraction is %f' % (float(my_mask.sum())/(nx*ny)))
    # Project into TOD
    wand = moby2.pointing.ArrayWand.for_tod_source_coords(
        tod, ref_coord=(ra, dec), scan_coords=True)

    gridding = moby2.pointing.GridPixelization.forFitsMap(maskMap)
    proj = moby2.pointing.WandProjector(wand, gridding, tod.fplane)
    # Warning this creates TOD-sized data...
    mask_vect = proj.deproject_from_map(maskMap.data)
    # Get cuts that is the same dets as focal plane entries.
    cuts = moby2.TODCuts(nsamps=tod.nsamps, det_uid=tod.fplane.det_uid,
                         sample_offset=tod.info.sample_index)
    for i,d in enumerate(mask_vect):
        cuts.add_cuts(i, d, mask=True)
    return cuts


def get_turnaround_dets(pos_cuts, az, az_lim):
    """For each detector in pos_cuts, determine if the detector is inside
    the pos_cuts mask at turn-around.  Turn-around is defined as times
    when az is within az_lim of its minimum or maximum values (az and
    az_lim must both have the same units).
    """
    az_lo, az_hi = az.min(), az.max()
    turnaround_mask = ((az - az_lo > az_lim) + (az_hi - az > az_lim))
    det_ok = np.empty(len(pos_cuts.cuts))
    for i,c in enumerate(pos_cuts.cuts):
        det_ok[i] = ((c.get_mask()*turnaround_mask).sum() == 0)
    return det_ok


def get_scanning_cuts(az, t, buffer_time=0.3):
    """Based on az, determine the typical scan speed.  Cut portions where
    scan is significantly slower than that.  Buffer such cuts by the
    buffer_time, in seconds.
    """
    daz, dt = np.diff(az), np.diff(t)
    tstep = np.median(dt)
    vaz = daz/dt
    vtyp = np.median(abs(vaz))
    s_slow = abs(vaz) < vtyp / 3.
    cv = CutsVector.from_mask(s_slow)
    # Smooth any any minor glitches.
    cv = cv.get_complement().get_buffered(5).get_collapsed().get_complement()
    # Buffer the no-scan regions.
    cv = cv.get_buffered(max(buffer_time / tstep, 1))
    return cv


def find_glitches(data, dets, real_filter, nsig, maxGlitch, minSep):
    return libactpol.get_glitch_cuts(
        data,
        np.asarray(dets, dtype='int32'),
        np.asarray(real_filter, dtype='float32'),
        nsig, maxGlitch, minSep)


def get_glitch_cuts(data=None, dets=None, tod=None, params={}):
    # Get defaults from somewhere...
    gp = { 'nSig': 10.0, 'tGlitch' : 0.002, 'minSeparation': 30,
           'maxGlitch': 50000, 'highPassFc': 5.0, 'buffer': 6 }
    gp.update(params)
    #gp = actDict.ACTDict()
    #gp.read_from_file(constants.CUT_PARAMS[tod.info.season][tod.info.array])
    #gp = gp['glitchParams']
    #gp.update(params)
    if data is None:
        data = tod.data
    if dets is None:
        dets = np.arange(tod.data.shape[0])

    filtvec = moby2.tod.filters.sine2highPass(tod, fc=gp['highPassFc']) * \
        moby2.tod.filters.gaussianFilter(tod, timeSigma=gp['tGlitch'])

    glitch_cuts = libactpol.get_glitch_cuts(
        data,
        np.asarray(dets, dtype='int32'),
        np.asarray(filtvec, dtype='float32'),
        gp['nSig'], gp['maxGlitch'], gp['minSeparation'])
    
    if tod is not None:
        tod.abuses += [{'name' : 'filterAndCutGlitches', 'type' : 'cut',
                        'glitchParams': gp}]
        cuts = TODCuts.for_tod(tod, assign=False)
    else:
        cuts = TODCuts(data.shape[0], data.shape[1])

    for d, c in zip(dets, glitch_cuts):
        cuts.cuts[d] = CutsVector(c, cuts.nsamps)
    
    cuts.buffer(gp['buffer'])
    return cuts


def fill_cuts(tod=None, cuts=None, data=None,
              neighborhood=40, do_all=False, filterScaling=1.0,
              extrapolate=False,
              sample_index=0,
              no_noise = False):
    """Replace certain samples of a TOD with a straight line plus some
    white noise.  The samples to fill are specified by a TODCuts
    object, which is passed through the cuts= argument, or obtained
    from tod.cuts.  The data array to fill is passed through thte
    data= argument, or else tod.data is used.

    The cut samples will be replaced with a straight line and white
    noise, with the linear fit and RMS taken from a "neighborhood" of
    samples on either side of the cut region.  The white noise can be
    modulated by the filterScaling argument, or turned off entirely by
    passing no_noise=True.

    The extrapolate argument affects how cut regions at the beginning
    or end of a timestream are handled; False means that the linear
    fill will have its slope forced to zero.
    """
    psLib.trace("moby",4,"cuts.fill_cuts begins...")
    # Just passing in tod is enough.
    if cuts is None:
        cuts = tod.cuts
    if data is None:
        data = tod.data
        si = tod.info.sample_index
    else:
        si = sample_index
    if hasattr(tod, 'info'):
        tod_name = tod.info.name
    else:
        tod_name = 'data'

    nsamps = data.shape[-1]

    # Check alignment between cuts and data.
    if si != cuts.sample_offset:
        print("TOD was loaded from sample %i but cuts have sample offset %i. The first samples will be cut" %(si, cuts.sample_offset))
    else:
        assert((data.shape[-1] == cuts.nsamps) | (data.shape[-1] == cuts.nsamps+1))

    # Prevent det_uid mismatch.
    assert(data.shape[0] == len(cuts.det_uid))
    if data is None:
        assert(np.all(tod.det_uid == cuts.det_uid))

    if do_all:
        det_list = range(data.shape[0])
    else:
        det_list = cuts.get_uncut()

    for deti in det_list:
        mask = cuts.cuts[deti].get_mask()
        offset = cuts.sample_offset - si
        total_mask = np.ones(nsamps, dtype=bool)
        total_mask[max(0,offset):min(nsamps,cuts.nsamps+offset)] = \
            mask[max(0,-offset):min(nsamps-offset,cuts.nsamps)]
        cuts_list = CutsVector.from_mask(total_mask)
        if no_noise: 
            noise = None
        else:
            ncut = np.dot([-1,1], cuts_list.sum(axis=0))
            noise = np.random.normal(size=ncut).astype('float32') * filterScaling
        libactpol.fill_cuts(data[deti], cuts_list, noise,
                            neighborhood, extrapolate)

    if tod is not None:
        tod.abuses += [{ 'name':'fillCuts', 'neighborhood':neighborhood }]
    psLib.trace("moby",4,"cuts.fill_cuts exits.")
        

def fill_cuts_old(tod=None, cuts=None, data=None,
              neighborhood=40, do_all=False, filterScaling=1.0):
    """
    @brief fill regions of cut data with synthetic data
    @param tod TOD to fill
    @param cuts TODCuts to use (defaults to tod.cuts)
    @param data data array (defaults to tod.dta)
    @param neighborhood int: number of data points two either side of a cut
           region to consider when estimating baseline and rms.
    """
    # Just passing in tod is enough.
    if cuts is None:
        cuts = tod.cuts
    if data is None:
        data = tod.data
    if hasattr(tod, 'info'):
        tod_name = tod.info.name
    else:
        tod_name = 'data'

    if do_all:
        det_list = range(data.shape[0])
    else:
        det_list = cuts.get_uncut()

    psLib.trace("moby",3,"cuts.fill_cuts begins...")
    nsamps = data.shape[1]
    for deti in det_list:
        for cfirst, clast in cuts.cuts[deti]:
            if cfirst >= nsamps:
                break
            # This .trace call is not what makes the cuts slow:
            #psLib.trace("moby",6,"Cut between %d (%f sec), %d (%f sec)." % \
            #                ( cfirst, tod.ctime[cfirst], clast,
            #                  tod.ctime[clast-1] ) )
            # Left and right side sample windows and interior region.
            xStart = np.arange(max(0, cfirst-neighborhood), cfirst)
            xStop = np.arange(clast, min(clast+neighborhood, nsamps))
            xGap = np.arange(cfirst, min(nsamps, clast + 1))
            extr = (len(xStart)==0) or (len(xStop)==0)
            x = np.concatenate((xStart, xStop))
            if len(x) == 0 or len(xGap) == 0:
                c = c.__next__
                continue
            # Linear component, if both ends anchored.
            y = data[deti,x]
            if extr:
                m = 0.0
                b = y.mean()
            else:
                m, b = np.polyfit(x, y, 1)
            resid = y - b - x*m
            rms = resid.std()
            if rms == 0:
                psLib.trace("moby", 4, "fillCuts: Gap RMS 0 (%s: det %d, i0 %d, i1 %d)" % \
                        (tod_name, deti, cfirst, clast))
            yGap = b + xGap*m
            if filterScaling != 0.0 and rms > 0:
                yGap += np.random.normal(loc=0, scale=rms, size=len(xGap))
            data[deti,xGap] = yGap
    psLib.trace("moby",3,"cuts.fill_cuts ends.")
    if tod is not None:
        tod.abuses += [{ 'name':'fillCuts', 'neighborhood':neighborhood }]

