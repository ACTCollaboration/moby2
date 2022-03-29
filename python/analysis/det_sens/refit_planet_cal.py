from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

"""This module performs similar things to fit_planet_cal but is more
flexible in the data selection and parameter couplings.  For example,
you can load two seasons of data and fit a single loading slope but
two different amplitudes (?).
"""

import moby2
import numpy as np
import pylab as pl

from . import planet_cal
from moby2.analysis import beam_ana

def chi2(p, *args):
    b, m = p[:-1], p[-1]
    chi2 = 0
    for i, _b in enumerate(b):
        chi2 += ((args[i*2+1] - args[i*2]*m - _b)**2).sum()
    return chi2


def main(args):
    import os
    from argparse import ArgumentParser
    o = ArgumentParser()
    o.add_argument('param_file')
    o.add_argument('--write-abscal', action='store_true')
    args = o.parse_args(args)

    cal_cfg = moby2.util.MobyDict.from_file(args.param_file)
    model_cand_template = 'models_{tag}_cand.fits'.format(**cal_cfg)
    model_final_template = 'models_{tag}.fits'.format(**cal_cfg)

    if args.write_abscal:
        if os.path.exists(model_final_template):
            if os.path.exists(model_cand_template):
                o.error('Both %s and %s found.  Promote your candidate or '
                        'delete it.' %
                        (model_final_template, model_cand_template))
            else:
                write_abscal(cal_cfg, [model_final_template])
                o.exit()
        else:
            o.error('Finalized model file %s not found... promote '
                    'your template?' % model_final_template)

    cal_subsets = cal_cfg['cal_configs']

    models = []
    for name, pa, fcode, omega_in, t_ranges in cal_subsets:
        om = moby2.scripting.OutputFiler(prefix=cal_cfg['plot_prefix'].format(pa=pa, tag=cal_cfg['tag']))
        dats = []
        for scode, omega_per, load_cut, t_lo, t_hi in t_ranges:
            tid, x, y = moby2.util.ascii.read_columns(
                cal_cfg['data_src'] + '%s_%s_%s_accepted.txt' % (pa, scode, fcode))
            t = np.array([float(_t[:10]) for _t in tid])
            s = (t_lo <= t) * (t < t_hi) * (x < load_cut)
            om_recal = omega_per / omega_in
            dats.extend([x[s], np.log(y[s] * om_recal)])
        p0 = [.1] * (len(dats)//2) + [-.001]  #.1, .1, -.001
        print(p0)
        p0, dats = tuple(p0), tuple(dats)
        p1 = moby2.util.fitting.multi_fmin(chi2, p0, args=dats, free=[-1])
        p1 = moby2.util.fitting.multi_fmin(chi2, p1, args=dats, free=range(len(p0)-1))
        p1 = moby2.util.fitting.multi_fmin(chi2, p1, args=dats, free=range(len(p0)))
        print(p1)
        m = p1[-1]
        ax0 = pl.gca()
        ax1 = pl.axes((0.65,0.2,0.2,0.2))
        bins = np.arange(-.205, .206, .01)
        for i,(_b, ri) in enumerate(zip(p1[:-1], t_ranges)):
            x_, y_ = dats[i*2:][:2]
            r = y_ - _b - m*x_
            xx = np.linspace(x_.min(), x_.max(), 20)
            a, = ax0.plot(xx, np.exp(_b + m*xx), label='%s scatter=%.4f' % (ri[0], r.std()))
            a = ax0.scatter(x_, np.exp(y_), color=a._color, s=10)
            ax1.hist(r, bins=bins, alpha=.6)
        ax0.set_ylim(0)
        ax0.legend(loc='lower left')
        ax0.set_title('%s %s %s' % (name, pa, fcode))
        ax0.set_xlabel('PWV/sin(alt) [mm]')
        ax0.set_ylabel('Cal [pW/K]')

        om.savefig('%s_%s.png'% (fcode, name.replace('/','-')))
        for (scode, omega_per, load_cut, t_lo, t_hi), b in zip(t_ranges, p1):
            print(scode, fcode, t_lo, t_hi, b, m)
            models.append((pa, scode, fcode, t_lo, t_hi, b, m))

    cols = list(map(np.array, list(zip(*models))))
    moby2.util.StructDB.from_data(
        list(zip(['pa', 'scode', 'fcode', 't_lo', 't_hi', 'b', 'm'], cols))).\
        to_fits_table(model_cand_template)


def write_abscal(cal_cfg, model_files=[]):
    # Master catalog.
    print('Loading master catalog...')
    cat = moby2.scripting.get_obs_catalog()

    print('... loaded %i items.' % len(cat))
    mask = np.ones(len(cat), bool)
    for t in cal_cfg.get('blacklist', []):
        mask = mask * (cat['tod_name'] != t)
    cat = cat[mask]
    print('... and %i items remain after processing blacklist.' % len(cat))

    # Additional recal?
    recals = cal_cfg['recal_instructions']  # a list of (type, h5filename, h5datasetname).
    recals = [
        (r[0], moby2.util.StructDB.from_hdf(r[1], dataset=r[2])) for r in recals]

    # Loop through our models.
    from glob import glob
    output_cols = [[],[],[]]

    for mf in model_files:
        print('Loading model file %s' % mf)
        models = moby2.util.StructDB.from_fits_table(mf)
        for ri in range(len(models)):
            row = models[ri]
            print(row)
            #if row['pa'] != 'pa3': continue
            s = ((cat['pa'] == row['pa']) *
                 (cat['scode'] == row['scode']) *
                 (row['t_lo'] <= cat['ctime']) *
                 (cat['ctime'] < row['t_hi']))
            x = cat[s]['loading']
            patch_load = cal_cfg.get('patch_bad_loading', 1.5)
            if patch_load is None or patch_load < 0:
                print('  Dropping %i of %i that have loading < 0.' % ((x<0).sum(), len(x)))
                s *= (cat['loading'] >= 0)
                x = cat[s]['loading']
            else:
                print('  Patching %i of %i that have loading < 0.' % ((x<0).sum(), len(x)))
                x[x<0] = patch_load
            y = np.exp(-(row['b'] + row['m'] * x))
            for rdb_type, rdb in recals:
                tfcodes = np.array(['%s_%s' % (t, row['fcode'])
                                    for t in cat[s]['tod_name']])
                #idx = rdb.select_inner({'tod_id': tfcodes})
                idx = rdb.select_inner({'tod_id': cat[s]['tod_name'],
                                        'band_id': [row['fcode'] for i in s.nonzero()[0]]})
                n_found = (idx>=0).sum()
                print('  Recal found for %i of %i' % (n_found, len(idx)))
                if rdb_type != 'optional' and n_found < len(idx):
                    raise RuntimeError("Missing cals for non-optional recal input.")
                y[idx>=0] *= rdb['cal'][idx[idx>=0]]
            new_cols = [list(cat[s]['tod_name']),
                        [row['fcode'] for i in s.nonzero()[0]],
                        list(1e6 * y)]
            for c,n in zip(output_cols, new_cols):
                c.extend(n)

    TAG = cal_cfg['tag']
    hfile, tfile = 'abscal_%s.h5' % (TAG), 'abscal_%s.txt' % (TAG)
    tfile = 'abscal_%s.txt' % (TAG)

    print('Writing to %s and %s.' % (hfile, tfile))
    db = moby2.util.StructDB.from_data(list(zip(['tod_id', 'band_id', 'cal'], output_cols)))
    db = db[db.argsort()]
    db.to_hdf(hfile, dataset='abscal', clobber=True)
    db.to_column_file(tfile)
