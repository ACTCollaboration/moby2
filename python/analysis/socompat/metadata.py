#
# Export ACT metadata to SO formats.
#

from sotodlib import core
from sotodlib.io import metadata

from .tod import actpol_load_observation

import so3g
import moby2

import numpy as np
import h5py
import re

__all__ = [
    'register_loaders',
    'extract_basename',
]

# Note some registration code runs on import; see bottom of source.

class ActTodCal(metadata.simple.PerDetectorHdf5):
    """ACT AbsCal are stored per-TOD, per-frequency band in HDF5 datasets.

    """
    @classmethod
    def _prefilter_data(cls, data_in):
        # Translations from ACT metadata names to new better names.
        return super()._prefilter_data(data_in, key_map={
            'tod_id': 'obs:obs_id',
            'band_id': 'dets:band'})


def get_abscal_proddb(filename, dataset):
    """Populate a ProdDb from the given HDF5 file, at the given dataset or
    group address.

    """

    fin = h5py.File(filename, mode='r')
    data = fin[dataset].value

    all_obs_id = sorted([str(x, 'ascii') for x in
                         list(set(list(data['tod_id'])))])

    scheme = core.metadata.ManifestScheme()\
             .add_exact_match('obs:obs_id')\
             .add_data_field('loader')\
             .add_data_field('dataset')
    man = core.metadata.ManifestDb(scheme=scheme)
    for obs_id in all_obs_id:
        man.add_entry({'obs:obs_id': obs_id,
                       'dataset': dataset,
                       'loader': 'actpol_abscal'},
                      filename,
                      commit=False)
    man.conn.commit()
    return man


def extract_basename(text):
    """Find a TOD name in text (ctime.ctime.ar1) and return it, along
    with each of the 3 components.

    """
    m = re.search('([0-9]+)\.([0-9]+)\.(ar[123456789])', text)
    return [m.group(i) for i in [0, 1, 2, 3]]


def load_act_cuts(filename, pa=None):
    """Read ACT cuts file and return an AxisManager containing the flags.

    The pa is needed to define the dets axis properly, but if passed
    as None then the code will try to determine it by searching for a
    standard ACT TOD basename in filename.
n
    """
    if pa is None:
        _, _, _, ar = extract_basename(filename)
        pa = 'pa' + ar[-1]

    cuts_in = moby2.tod.TODCuts.from_actpol_cuts_file(filename)
    det_names = ['%s_%04i' % (pa, u) for u in cuts_in.det_uid]

    flags = []
    for d, c in zip(det_names, cuts_in.cuts):
        i = so3g.proj.Ranges.from_array(c + cuts_in.sample_offset, cuts_in.nsamps)
        flags.append(i)
    flags = so3g.proj.RangesMatrix(flags)

    aman = core.AxisManager(
        core.LabelAxis('dets', det_names),
        core.OffsetAxis('samps', cuts_in.nsamps, cuts_in.sample_offset),
    )

    aman.wrap('flags', flags, [(0, 'dets'), (1, 'samps')])
    return aman


def load_act_cal(filename, pa=None):
    """Read ACT relcal file and return an AxisManager containing the
    info.

    The pa is needed to define the dets axis properly, but if passed
    as None then the code will try to determine it by searching for a
    standard ACT TOD basename in filename.

    """
    if pa is None:
        _, _, _, ar = extract_basename(filename)
        pa = 'pa' + ar[-1]

    cal_in = moby2.Calibration.from_dict(filename)
    det_names = ['%s_%04i' % (pa, u) for u in cal_in.det_uid]

    aman = core.AxisManager(
        core.LabelAxis('dets', det_names),
    )

    aman.wrap('cal', cal_in.cal, [(0, 'dets')])
    return aman


def load_detoffsets_file(position, polarization, pa=None):
    """Read ACT detector offsets + polarization angles files and return an
    AxisManager.

    """
    def guess_pa(text):
        # Find 'pa?' or 'ar?' in text (perhaps a filename) and return it.
        m = re.search('pa[123456789]', text)
        if m is not None:
            return m.group(0)
        m = re.search('ar[123456789]', text)
        if m is not None:
            return 'pa' + m.group(0)[-1]
        raise ValueError('Could not determine "pa" based on filename; '
                         'pass it in yourself.')

    if pa is None:
        pa = guess_pa(position)

    # Load position data...
    fp = moby2.scripting.products.get_detector_offsets(
        {'format': 'ascii',
         'columns': [('det_uid', 0), ('x', 2), ('y', 4), ('mask', 1)],
         'match_bad': 0,
         'filename': position})
    fp = fp.subset(fp.mask)

    det_names = ['%s_%04i' % (pa, u) for u in fp.det_uid]
    aman = core.AxisManager(core.LabelAxis('dets', det_names))
    aman.wrap('act_x', fp.x, [(0, 'dets')])
    aman.wrap('act_y', fp.y, [(0, 'dets')])
    aman.fp = fp

    if polarization is None:
        aman.wrap('act_pol', np.zeros(len(fp.x)), [(0, 'dets')])
    else:
        cols = moby2.util.ascii.read_columns(polarization)
        mask = cols[1] >= 0

        det_names = ['%s_%04i' % (pa, u) for u in cols[0][mask]]
        pol_ang = cols[1][mask]
        amanp = core.AxisManager(core.LabelAxis('dets', det_names))
        amanp.wrap('act_pol', pol_ang, [(0, 'dets')])
        aman.merge(amanp)

    # The coordinate convention for ACT and SO is ever-so-slightly
    # different.  Here's the libactpol construction:
    #    Quaternion_r3(q, -pol_angle);
    #    Quaternion_r2_mul(focalplane_x, q);
    #    Quaternion_r1_mul(-focalplane_y, q);

    # Note focalplane_{x,y} correspond to moby2.FocalPlane.{y,x}.
    # Also pol_angle is FocalPlane.phi - pi/2.
    #
    # This will differ from SO construction by O(x^3 + y^4).
    DEG = np.pi/180
    q = (so3g.proj.quat.euler(0, -aman['act_x']) *
         so3g.proj.quat.euler(1,  aman['act_y']) *
         so3g.proj.quat.euler(2, np.pi/2 - aman['act_pol']*DEG))
    aman.wrap('q', np.array(q), [(0, 'dets')])

    # The SO projection coordinates xi and eta:
    theta, phi, gamma_minus_phi = so3g.proj.quat.decompose_iso(q)
    gamma = gamma_minus_phi + phi
    gamma[(theta == 0) * (phi == 0)] = 0.

    aman.wrap('xi',  theta * np.cos(phi), [(0, 'dets')])
    aman.wrap('eta', theta * np.sin(phi), [(0, 'dets')])
    aman.wrap('gamma', gamma, [(0, 'dets')])

    return aman


# Support the SuperLoader api...

class ActLoader:
    def __init__(self, detdb=None, **kw):
        self.detdb = detdb

    def batch_from_loadspec(self, index_lines):
        return [self.from_loadspec(line) for line in index_lines]


class ActCutsLoader(ActLoader):
    def from_loadspec(self, index_line):
        return load_act_cuts(index_line['filename'])


class ActCalLoader(ActLoader):
    def from_loadspec(self, index_line):
        return load_act_cal(index_line['filename'])


class ActAbsCalLoader(ActLoader):
    def from_loadspec(self, index_line):
        return ActTodCal.from_loadspec(index_line, detdb=self.detdb)


class ActDetOfsLoader(ActLoader):
    def from_loadspec(self, index_line):
        return load_detoffsets_file(index_line['filename'],
                                    index_line.get('pol_filename'),
                                    pa=index_line['dets:array_name'])


class ActPointOfsLoader(ActLoader):
    def from_loadspec(self, index_line):
        return metadata.simple.PerDetectorHdf5.from_loadspec(
            index_line, detdb=self.detdb)


def register_loaders():
    metadata.loader.REGISTRY.update({
        'actpol_cuts': ActCutsLoader,
        'actpol_cal': ActCalLoader,
        'actpol_abscal': ActAbsCalLoader,
        'actpol_detofs': ActDetOfsLoader,
        'actpol_pointofs': ActPointOfsLoader,
    })
    from sotodlib.data import load
    load.OBSLOADER_REGISTRY['actpol_moby2'] = actpol_load_observation
