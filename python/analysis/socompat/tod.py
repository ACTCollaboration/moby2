import moby2
from sotodlib import core

import numpy as np

__all__ = ['get_tod']


def get_tod(params, aman=None):
    """params can be the usual kind of dictionary that you would pass to
    moby2.scripting.get_tod, but I bet you'd rather just pass in a
    string with the basename in it.  And that will work too.

    Returns pair (axisman, tod) containing the newly populated
    AxisManager and the husk of the moby2 TOD, with the latter of
    potential use through its .info and .get_hk attributes, for
    example.

    """
    if isinstance(params, str):
        params = {'filename': params}
    else:
        params = params.copy()

    if aman is not None:
        dets = aman._axes['dets'].vals
        det_uid = [int(x.split('_')[1]) for x in dets]
        if 'det_uid' not in params:
            params['det_uid'] = det_uid

    tod = moby2.scripting.get_tod(params)
    count = len(tod.ctime)
    if tod.data is None:
        # Hot-wire empty signal for no-dets case.
        tod.data = np.zeros((0, count), 'float32')

    pa = 'pa' + tod.info.array[-1]
    axes = {
        'dets':
        core.LabelAxis('dets', ['%s_%04i' % (pa, d) for d in tod.det_uid]),
        'samps':
        core.OffsetAxis('samps', count, params.get('start', 0),
                        tod.info.basename),
    }

    data = core.AxisManager(axes['dets'], axes['samps'])
    if params.get('read_data', True):
        data.wrap('signal', tod.data, [(0, 'dets'), (1, 'samps')])
    else:
        data.wrap('signal', None)
    data.wrap('timestamps', tod.ctime, [(0, 'samps')])

    del tod.data, tod.ctime

    # Include a cal to remove readout filter gain -- not needed if you
    # deconvolve the readout filter directly.
    data.wrap('readout_filter_cal',
              np.ones(data.dets.count) / tod.info.mce_filter.gain(),
              [(0, 'dets')])

    # Describe the readout filter in sotodlib language.
    iir = mce_to_iir(tod.info.mce_filter.params, tod.info.mce_filter.f_samp)
    iir_params = core.AxisManager()
    for k, v in iir.items():
        iir_params.wrap(k, v)
    data.wrap('iir_params', iir_params)

    flags = core.AxisManager(axes['samps'])
    flags.wrap('enc_flags', tod.enc_flags, [(0, 'samps')])
    flags.wrap('pointing_mask', tod.pointing_mask, [(0, 'samps')])
    data.wrap('flags', flags)

    del tod.enc_flags, tod.pointing_mask

    bore = core.AxisManager(axes['samps'])
    bore.wrap('az', (tod.az + np.pi) % (2 * np.pi) - np.pi,
              [(0, 'samps')])
    bore.wrap('el', tod.alt, [(0, 'samps')])
    bore.wrap('roll', np.zeros(tod.alt.shape, tod.alt.dtype),
              [(0, 'samps')])
    data.wrap('boresight', bore)

    del tod.az, tod.alt

    # Let's try jamming in the array data too...
    array_data = core.AxisManager(axes['dets'])
    ad = tod.info.array_data
    idx = ad.select_inner({'det_uid': tod.det_uid})
    for k, v in ad.items():
        if k == 'det_uid':
            continue
        array_data.wrap(k, v[idx], [(0, 'dets')])
    data.wrap('array_data', array_data)

    if aman is not None:
        data = aman.merge(data)

    return (data, tod)


def mce_to_iir(p, f_mux):
    """Converts MCE filter params p and muxing frequency f_mux to an
    iir_params block suitable for use with sotodlib
    (tod_ops.filters.iir_filter).

    Note the effective "DC gain" is included in these parameters and
    does not need to be carried separately.

    Returns:
      dict with keys 'a', 'b', and 'fscale'.

    """
    K = 1./2**14
    scalars = [K, K, K, K, 1., 1.]
    b11, b12, b21, b22, k1, k2 = [s*p for s,p in zip(scalars, p)]
    a = np.array([1., - b11 - b21, b11*b21 + b12 + b22, - b11 * b22 - b21 * b12, b12*b22])
    b = np.array([1., 4, 6, 4, 1]) / 2**(k1+k2)
    fscale = 1./f_mux
    return {'a': a, 'b': b, 'fscale': fscale}


def actpol_load_observation(db, obs_id, dets=None, samples=None,
                            no_signal=None, prefix=None, **kwargs):
    """Load observation -- keep this compatible with sotodlib.Context
    get_obs interface."""
    if len({k: v for k, v in kwargs.items() if v is not None}):
        raise RuntimeError(f'This function does not understand args: {kwargs}')

    aman = None
    if dets is not None:
        # Then it is a list of dets.  Restrict the list to only things
        # in our detset(s).
        detsets = db.get_detsets(obs_id)
        assert(len(detsets) == 1)
        dets_ref = db.get_dets(detsets[0])
        dets = [d for d in dets if d in dets_ref]
        # Pass it in through an aman.
        aman = core.AxisManager(core.LabelAxis('dets', dets))

    #To support moby2 filebases, the default prefix shall be '' --
    #i.e. the user must either record clean basenames or full paths
    #for this system or pass in a prefix explicitly.
    if prefix is None:
        prefix = ''

    start, end = 0, None
    if samples is not None:
        start, end = samples

    # Use the db to get the filename.
    files_by_detset = db.get_files(obs_id, prefix=prefix)
    assert(len(files_by_detset) == 1)
    for detset, fileidx in files_by_detset.items():
        assert(len(fileidx) == 1)
        filename, sample_start, sample_stop = fileidx[0]
        aman, tod = get_tod({'filename': filename,
                             'repair_pointing': True,
                             'read_data': (no_signal is not True),
                             'start': start,
                             'end': end}, aman=aman)

    return aman

