from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
from . import util

from . import beam_obs

import numpy as np
import os

fitsMap = moby2.mapping.fits_map.fitsMap
DEG = np.pi/180

#
# Elliptical even polynomial model for fitting peak
#
def peak_model(p, phi, r):
    a0, a2, a4, Q, U = p
    r = r * (1 + Q*np.cos(2*phi) + U*np.sin(2*phi))
    return  a0 + a2*r**2 + a4*r**4

def peak_chi2(p, x, y, z):
    phi, r = np.arctan2(y, x), (x**2 + y**2)**.5
    z0 = peak_model(p, phi, r)
    return ((z-z0)**2).sum()


#
# Cone model for fitting elliptical FWHM.
#

def cone_model(p, phi, dz):
    I,Q,U,m = p
    r0 = (I + Q * np.cos(2*phi) + U * np.sin(2*phi)) * (1. - m*dz)
    x0, y0 = r0*np.cos(phi), r0*np.sin(phi)
    return x0, y0

def cone_chi2(p, x, y, dz):
    phi = np.arctan2(y, x)
    x0, y0 = cone_model(p, phi, dz)
    return ((x-x0)**2 + (y-y0)**2).sum()


#
# Trumpet model for fitting elliptical wing.
#

def trumpet_model(p, phi, r):
    Q,U,a,base = p
    r0 = (1. + Q * np.cos(2*phi) + U * np.sin(2*phi)) * r
    return a*r0**-3 + base

def trumpet_chi2(p, phi, r, z):
    z0 = trumpet_model(p, phi, r)
    return ((z-z0)**2).sum()


def process_map(B, F, params={}, freq=None, context=None):
    def get_param(key, default):
        if key in params:
            if isinstance(params[key], dict) and freq in params[key]:
                return params[key][freq]
            return params[key]
        return default

    if context is None:
        context = moby2.scripting.MobyContext()
    fprint = context.get_logger()
    outputm = context.get_outputm()

    if freq is None:
        freq = params.get('force_frequency')
    if freq is None:
        try:
            freq = util.extract_freq(F)
        except:
            print('Warning, could not determine frequency band from filename.')
            freq = 'f000'

    fprint('map_name = %s' % B)
    fprint('map_file = %s' % F)
    fprint('freq_band = %s' % freq)

    row = []
    row.append(('name', B))
    row.append(('band', freq))
    m = fitsMap(F)

    if get_param('recenter', True):
        bol = beam_obs.BeamObs(map=m)
        bol.fitGaussian()
        x0, y0 = [bol['gaussian'][k] for k in ['dx', 'dy']]
        m.x -= x0
        m.y -= y0
    r = m.radii()
    R0, R1, R2 = [_x/60. for _x in
                  get_param('mask_radii_arcmin', (4., 5., 7.))]
    s0, s1, s2 = (r<R0), (r<R1), (r<R2)
    z0 = m.data[s1*~s0].mean()
    m.data -= z0
    dx, dy = m.pixn.truePitch()
    peak0 = m.data[s0].max()
    fprint('Rough peak estimate: ', peak0)
    s00 = (m.data / peak0 > .9) * s0
    p0 = (peak0, -peak0*10, 0., 0., 0.)
    fprint(p0)
    x, y, z = (m.x+m.y*0)[s00], (m.x*0+m.y)[s00], (m.data[s00])
    p1 = moby2.util.fitting.multi_fmin(peak_chi2, p0, args=(x, y, z),
                                       free=[True,True,False,False,False], disp=0)
    p2 = moby2.util.fitting.multi_fmin(peak_chi2, p1, args=(x, y, z),
                                       free=[True,False,False,True,True], disp=0)
    p3 = moby2.util.fitting.multi_fmin(peak_chi2, p2, args=(x, y, z), disp=0)
    pX = p3
    fprint('Refined: ', pX[0])
    fprint('Params: ', pX)
    
    import pylab as pl

    pl.subplot(211)
    pl.title(B)
    _y = m.data[s2]
    pl.scatter(r[s2]*60., _y, s=1, alpha=.4)
    z0 = peak_model(pX, np.arctan2(y, x), r[s00])
    pl.xlabel('Radius (arcmin)')
    pl.ylabel('Amplitude')
    pl.ylim(_y.min(), _y.max())
    pl.subplot(212)
    pl.scatter(r[s00]*60., (z-z0)/pX[0], marker='x')
    pl.xlabel('Radius (arcmin)')
    pl.ylabel('Residuals (peak fractional)')
    outputm.savefig('%s_reduced_00.png' % B)

    peak = pX[0]
    omega = (m.data[s0].sum() * abs(dx*dy) * (np.pi/180)**2 * 1e9 / peak)
    fprint('Solid angle (quick): %.4f nsr' % omega)
    row.extend(list(zip(['quick_peak', 'quick_omega'], [peak, omega])))
    
    # Wing fit and boost the solid angle.

    wing_radii = get_param('wing_fit_radii_arcmin', (3., 5.))
    s = (r >= wing_radii[0]/60.) * (r < wing_radii[1]/60.)
    z = m.data / peak

    # Binned data.
    dr = .2/60
    rbins = np.arange(wing_radii[0]/60., wing_radii[1]/60.+dr, dr)
    counts,_ = np.histogram(r.ravel(), rbins)
    b_y,_ = np.histogram(r.ravel(), rbins, weights=z.ravel())
    b_y /= counts
    b_r = 0.5*(rbins[1:]+rbins[:-1])

    # Unbinned data, for fitting ellipse too.
    r_, phi_, z_ = r[s], np.arctan2(m.y, m.x)[s], z[s]

    if not get_param('wing_fit_ellipticity', True):
        fprint('*** warning: fitting to binned data.')
        p0 = [0., 0., .1, b_y[-1]]
        p1 = moby2.util.fitting.multi_fmin(trumpet_chi2, p0, args=(0*b_r, b_r*60., b_y), disp=0,
                                           free=[False, False,True, True])
        pwing = p1
    else:
        p0 = [0., 0., .1, .01]
        p1 = moby2.util.fitting.multi_fmin(trumpet_chi2, p0, args=(phi_, r_*60., z_), disp=0,
                                           free=[False, False,True, True])
        p2 = moby2.util.fitting.multi_fmin(trumpet_chi2, p1, args=(phi_, r_*60., z_), disp=0,
                                           free=[True, True,True, True])
        pwing = p2

    # Talk about it...
    fprint('Wing fit.')
    fprint('Fit parameters: ', pwing)
    a_wing, base_wing = pwing[2]/60.**3, pwing[3]
    r_wing = wing_radii[0]/60.
    s_inside = (r<r_wing)
    if (abs(base_wing) > .01):
        fprint('*** warning: base_wing is %.5f ***' % base_wing)
    omega_inside = (((z[s_inside]-base_wing).sum() * abs(dx*dy) * (np.pi/180)**2 * 1e9) / 
                    (1. - base_wing))
    phii = np.arange(0., 360, .01) * np.pi/180
    phi_correction = np.mean((1. + pwing[0]*np.cos(phii) + pwing[1]*np.sin(phii))**-3)
    fprint('ellipticity correction is %.5f' % phi_correction)
    if (abs(phi_correction-1) > .01):
        fprint('*** warning: ellipticity correction is %.5f ***' % phi_correction)
    omega_outside = 2.*np.pi*a_wing/r_wing * (np.pi/180)**2 * 1e9 * phi_correction / (1.-base_wing)
    fprint('Wing (r=%.2f arcmin) fit solid angle: %.2f + %.2f = %.2f' % (
            r_wing*60, omega_inside, omega_outside, omega_inside + omega_outside))
    fprint('')
    row.extend(list(zip(['om_inner', 'om_outer', 'om_total'],
                   [omega_inside, omega_outside, omega_inside + omega_outside])))

    # Plot all that.
    ax = pl.subplot(211)
    for h, fmt, val in [
            (.9, 'inner: %.1f', omega_inside),
            (.8, ' wing: %.1f', omega_outside),
            (.7, 'total: %.1f', omega_inside+omega_outside)
    ]:
        ax.text(.95, h, fmt%val, ha='right', transform=ax.transAxes)
    pl.title(B)
    pl.scatter(r_*60, z_, s=10, alpha=.4)
    rr = np.arange(r_.min()*.9, r_.max()*1.1, .001)
    pl.plot(rr*60, trumpet_model(pwing, rr*0, rr*60.), lw=2, c='r')
    pl.scatter(b_r*60., b_y, marker='x', s=30, color='green')
    pl.xlim(0, rr.max()*60)
    pl.subplot(212)
    pl.scatter(r_*60, z_ - trumpet_model(pwing, phi_, r_*60),
               s=10, alpha=.4)
    pl.xlim(0, rr.max()*60)
    outputm.savefig('%s_reduced_10.png' % B)

    # This is pretty cool.  Fit a cone to points near the Half-Max.
    fprint()
    s = (abs(m.data / peak - 0.5) < .1) * s1
    fprint('Points in FWHM fit:', s.sum())
    p0 = (r[s].mean(), 0., 0., 1./r[s].mean())
    x, y, dz = (m.x+m.y*0)[s], (m.x*0+m.y)[s], (m.data / peak - 0.5)[s]
    p1 = moby2.util.fitting.multi_fmin(cone_chi2, p0, args=(x, y, dz), disp=0)
    fprint('Fit results: ', p1)
    I, P = p1[0] * 60., (p1[1]**2 + p1[2]**2)**.5 * 60.
    fprint(P/I)
    PHI = 0.5*np.arctan2(p1[2], p1[1]) / DEG
    GAM = PHI % 180 - 90.
    fprint('I.e., Mean FWHM = %.4f arcmin' % (2*I))
    fprint('Major axis FWHM = %.4f arcmin' % (2*I+2*P))
    fprint('Minor axis FWHM = %.4f arcmin' % (2*I-2*P))
    fprint('Major axis angle = %.4f degrees CCW from N' % (GAM))
    row.extend(list(zip(['fwhm', 'fwhm_maj', 'fwhm_min', 'fwhm_ang'],
                   [2*I, 2*I+2*P, 2*I-2*P, GAM])))

    pl.title(B)
    for _s, col in [(dz>0, '+'), (dz<=0, 'x')]:
        pl.scatter(x[_s], y[_s], marker=col)
    for i in range(len(x)):
        _x0, _y0 = cone_model(p1, np.arctan2(y[i], x[i]), dz[i])
        pl.plot([x[i],_x0],[y[i],_y0], c='green')
    phi = np.arange(0., 361.) * np.pi/180
    for d in np.linspace(-.1, .1, 5):
        _x, _y = cone_model(p1, phi, d)
        if (d==0):
            _r = (_x**2+_y**2)**.5
        pl.plot(_x, _y, c='black', ls={True: 'dotted', False: 'solid'}[d!=0])
    for r,p,lw in [(I+P, PHI, 3), (I-P, PHI+90, 1)]:
        _x, _y = r/60*np.cos(p*DEG), r/60*np.sin(p*DEG)
        pl.plot([-_x,_x], [-_y,_y], color='green', lw=lw)
        pl.gca().text(_x/2, _y/2, '%.3f\'' % (r*2))
    pl.xlim(-.03, .03)
    pl.gca().set_aspect('equal', 'datalim')
    outputm.savefig('%s_reduced_20.png' % B)

    fprint()
    return row


def summary_op(params, data, cat, mask=None, verbose=False):
    """Process the solid angle fit info in data, according to the
    parameters in params, using the matched obs_catalog in cat.
    Optionally apply mask to the data, too.

    The params should be a simple dictionary with entry 'select'::

      {'select': [...]}

    The list is passed to util.get_obs_mask.  Results are
    automatically grouped by band and array (pa).

    """
    def vprint(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    s0, n = util.get_obs_mask(cat, params['select'])
    if mask is not None:
        s0 *= mask
    pas = sorted(list(set(cat['pa'][s0])))
    bands = sorted(list(set(data['band'][s0])))
    grouped = {}
    for pa in pas:
        grouped[pa] = {}
        for band in bands:
            grouped[pa][band] = {}
            s = s0 * (cat['pa'] == pa) * (data['band'] == band)
            vprint(pa, band, s.sum())
            if s.sum() == 0: continue
            for k in ['om_inner', 'om_outer', 'om_total']:
                y = data[k][s]
                vprint ('  %-10s  %8.3f  %8.3f' % (k, y.mean(), y.std()))
                grouped[pa][band][k] = (y.mean(), y.std())
            # FWHM business...
            phi = data['fwhm_ang'][s] * DEG
            P = data['fwhm_maj'][s] - data['fwhm_min'][s]
            Q = P * np.cos(2*phi)
            U = P * np.sin(2*phi)
            subd = {'fwhm': data['fwhm'][s],
                    'fwhm_Q': Q,
                    'fwhm_U': U}
            for k, y in subd.items():
                vprint ('  %-10s  %8.3f  %8.3f' % (k, y.mean(), y.std()))
                grouped[pa][band][k] = (y.mean(), y.std())
            I, Q, U = [grouped[pa][band][k][0] for k in ['fwhm', 'fwhm_Q', 'fwhm_U']]
            P = np.sqrt(Q**2+U**2)
            vprint ('  %-10s  %8.3f' % ('fwhm_maj', I + P))
            vprint ('  %-10s  %8.3f' % ('fwhm_min', I - P))
            vprint ('  %-10s  %8.3f' % ('fwhm_ang', np.arctan2(U, Q)/2/DEG))

    return grouped


def driver(args):
    from argparse import ArgumentParser
    o = ArgumentParser()
    o.add_argument('param_file')
    o.add_argument('map_name', nargs='*', help="Maps to purge (if --purge).")
    o.add_argument('-r', '--refit', action='store_true')
    o.add_argument('-u', '--update', action='store_true')
    o.add_argument('--purge', action='store_true')
    args = o.parse_args(args)

    params = moby2.util.MobyDict.from_file(args.param_file)

    outputm = moby2.scripting.OutputFiler(
        prefix=params.get_deep(('output', 'prefix'), './'))

    context = moby2.scripting.output.MobyContext()
    context.outputm = outputm
    logger = context.get_logger()
    logger.set_file(outputm.get_filename('output.txt'), append=False)
    logger.show_time = False
    
    # Line-by-line results file
    ofile = outputm.get_filename('table.fits')

    if not os.path.exists(ofile):
        args.refit = True

    if args.refit:
        # Get map list.
        basenames, filenames = util.get_map_list(params['source_maps'])

        rows = []
        for B,F in zip(basenames, filenames):
            row = process_map(B, F, params['analysis'],
                              context=context)
            rows.append(row)
        
        # Bundle it.
        header = list(zip(*rows[0]))[0]
        columns = list(map(np.array, list(zip(*[list(zip(*row))[1] for row in rows]))))
        odb = moby2.util.StructDB.from_data(list(zip(header, columns)))

        logger('\n\nWriting results to %s' % ofile)
        odb.to_fits_table(ofile)

    if args.purge:
        # Remove specified outputs from the table.
        data = moby2.util.StructDB.from_fits_table(ofile)
        mask = np.ones(len(data), bool)
        for m in args.map_name:
            mask *= (data['name'] != m)
        print('Keeping %i of %i ...' % (mask.sum(), len(mask)))
        data = data[mask]
        data.to_fits_table(ofile)
        return

    # Summary operations.
    summary_ops = params.get_deep(('output', 'summaries'))
    if summary_ops is None or len(summary_ops) == 0:
        print('No summary operations defined in param file.')
    else:
        print('Reading %s for summary operations...' % ofile)
        data = moby2.util.StructDB.from_fits_table(ofile)
        cat = moby2.scripting.get_obs_catalog()
        basename = [x.split('_')[0] for x in data['name']]
        idx = cat.select_inner({'tod_name': basename})
        cat = cat[idx]

        for i, summary_params in enumerate(summary_ops):
            print('Processing summary block %i...' % i)
            results = summary_op(summary_params, data, cat, verbose=True)
