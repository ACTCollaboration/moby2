from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import pylab as pl

default_beam_image_opts = {
    'zoom': 6.,
    'center': True,
    'scale': 'log',
    'rescale': True,
    'size_inches': (6., 5.),
    'imshow_args': {},
}

def plot_beam_image(beam_map, params={}, beam_obs=None,
                    title=None,
                    filename=None, clf=True):
    popts = default_beam_image_opts.copy()
    popts.update(params)
    pl.gcf().set_size_inches(*popts['size_inches'])
    
    if beam_obs is not None:
        sig = beam_map.data - beam_obs['gaussian']['level']
        if popts.get('zero_radius'):
            r0 = popts['zero_radius']
            x0, y0 = beam_obs['gaussian']['center']
            r = beam_map.radii((x0,y0)) * 60.
            s = (beam_map.data!=0)*(r0 <= r)*(r < r0+1.)
            if s.sum() > 0:
                sig -= sig[s].mean()
    else:
        sig = beam_map.data

    if popts['rescale'] is True:
        sig /= beam_obs['gaussian']['amp']
    elif popts['rescale'] not in [False, None]:
        sig *= popts['rescale']
    if popts['scale'] in ['log', 'dB']:
        sig = pl.log10(abs(sig) + 1e-6*sig.max())
        if popts['scale'] == 'dB':
            sig = sig * 10.

    vmin, vmax = popts.get('vmin'), popts.get('vmax')
    if popts.get('std_dev_cut') is not None:
        xc, yc, r = popts.get('std_dev_cut_measure', (0,0,5./60))
        s = beam_map.circleMask(r, xc, yc)
        scut_lo, scut_hi = popts['std_dev_cut']
        m0, dm  = sig[s].mean(), sig[s].std()
        if vmin is None:
            vmin = max(m0 + scut_lo*dm, sig.min())
        if vmax is None:
            vmax = min(m0 + scut_hi*dm, sig.max())

    beam_map.imshow(data=sig, units='arcmin', vmin=vmin, vmax=vmax,
                    **params.get('imshow_args',{}))

    # Center and zoom
    if popts['center'] in [None,False]:
        x0, y0 = 0., 0.
    elif popts['center'] is True:
        x0, y0 = beam_obs['gaussian']['center']
        x0, y0 = x0*60, y0*60
    else:
        x0, y0 = popts['center']
    beam_map.zoom(popts['zoom'], x0, y0)

    if title is not None:
        pl.title(title)
    pl.xlabel('X (arcmin)')
    pl.ylabel('Y (arcmin)')
    pl.colorbar()
    if filename is not None:
        pl.savefig(filename)
        if clf:
            pl.clf()
    
