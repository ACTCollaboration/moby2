from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Load detector sensitivities for some planet observations and analyze
detector and array sensitivity over time, loading, etc.
"""

import moby2
import numpy as np
import pylab as pl
import os

class ArrayDataArchive:
    def __init__(self):
        self.adata = {}
    def get(self, info):
        season = info.get('season')
        if season is None:
            season = '20' + info['scode'][1:]
        array = info.get('array')
        if array is None:
            array = 'ar' + info['pa'][-1:]
        if not season in self.adata:
            self.adata[season] = {}
        if not array in self.adata[season]:
            self.adata[season][array] = moby2.scripting.get_array_data(
                {'season': season, 'array_name': array.upper()})
        return self.adata[season][array]

def match_info(info, db):
    for rules, result in db:
        for k, v in rules:
            if info[k] != v:
                break
        else:
            return result
    print('Failed to match ', info, db)
    return None

class GroupMasker:
    def __init__(self, group_keys, **data):
        self.keys = group_keys
        self.data = {}
        for k in self.keys:
            self.data[k] = np.asarray(data[k])
        self.vals = []
        self.n_items = len(list(data.values())[0])
        for k in self.keys:
            self.vals.append(sorted(list(set(self.data[k]))))
        self.i = 0
    def get_mask(self, i=None):
        if i is None:
            i = self.i
            self.i += 1
        j = len(self.keys)-1
        idx = []
        while j >= 0:
            if len(self.vals[j]) == 1:
                idx.insert(0,0)
            else:
                idx.insert(0, i % len(self.vals[j]))
                i //= len(self.vals[j])
            j -= 1
        if i > 0:
            return None, None
        mask = np.ones(self.n_items, bool)
        info = {}
        for k,v,i in zip(self.keys, self.vals, idx):
            info[k] = v[i]
            mask *= (self.data[k] == info[k])
        return info, mask

def main(args_in=[]):

    import optparse
    o = optparse.OptionParser()
    o.add_option('--loading-range',type=float,nargs=2)
    opts, args = o.parse_args(args_in)

    cfg = moby2.util.MobyDict.from_file(args[0])

    if cfg['array_sens'].get('tod_list', {}).\
            get('source', 'obs_catalog') == 'obs_catalog':
        blist = None
    else:
        blist = moby2.scripting.get_tod_list(cfg['array_sens']['tod_list'],
                                             cmdline_args=args[1:])

    planet_db = moby2.scripting.get_obs_catalog(cfg['cal_noise'].get('obs_db'))
    print('Input database has %i items.' % len(planet_db))

    if blist is not None:
        idx = planet_db.select_inner({'tod_name': blist})
        planet_db = planet_db[idx]
        print('Intersection with tod_list leaves %i items.' % len(planet_db))

    # Reduce using properties from DB.
    mask0 = moby2.scripting.get_selection(
        cfg['cal_noise']['obs_select'], planet_db)
    planet_db = planet_db[mask0]
    print('Selection with obs_select leaves %i items.' % len(planet_db))

    # Final basename list?
    blist = planet_db['tod_name']

    # Get 'planet_day' for plotting.
    t0 = [moby2.util.ctime.from_string('%s-01-01 0:0:0' % s)
          for s in planet_db['season']]
    planet_day = (planet_db['ctime'] - t0) / 86400.

    # 
    # We are left with
    #
    MAX_UID = 2048
    all_uid = np.arange(MAX_UID)
    all_sens = np.zeros((MAX_UID, len(blist)))

    ifiler = moby2.scripting.OutputFiler(prefix=cfg['cal_noise']['sens_prefix'])
    for bi, b in enumerate(blist):
        ofile = ifiler.get_filename(b+'.txt')
        try:
            u, s = moby2.util.ascii.read_columns(ofile)
        except:
            continue
        all_sens[u,bi] = s

    if not np.any(all_sens):
        print('No calibrated data found... did you run get_fpfit_planetcal or ...')
        o.exit(10)

    adata_src = ArrayDataArchive()

    ofiler = moby2.scripting.OutputFiler(prefix=cfg['array_sens']['prefix'])

    group_by = cfg['array_sens']['group_by'] # e.g. ['scode', 'pa']
    grp_data = dict([(k, planet_db[k]) for k in group_by])
    grouper = GroupMasker(group_by, **grp_data)
    
    while True:
        info, mask = grouper.get_mask()
        if info is None:
            break
        if mask.sum() < 0:
            continue
        tag = '_'.join([info[k] for k in group_by])
        TAG = tag.upper()
        print('Selecting', tag)
        print(' count: ', mask.sum())
        if mask.sum() == 0:
            print()
            continue

        adata = adata_src.get(info)
        n_uid = adata['det_uid'].max()+1

        #print 'Typical PWV:', np.median(planet_db[kept]['pwv'])
        #print 'Typical loading:', np.median(planet_db[kept]['loading'])
        demands = cfg['array_sens']['demands']
        red_sens = np.zeros(n_uid)
        # Loop over detectors
        for i in range(n_uid):
            x = all_sens[i,:] * mask
            _s = (x!=0)
            if _s.sum() >= demands['min_count']:
                red_sens[i] = np.median(x[_s])

        absurd_level = demands['min_det_noise']
        all_freqs = sorted(list(set(adata['fcode']).difference({'f000'})))

        # Total band sensitivity vs. loading?
        for f,p in zip(all_freqs, 'bgrk'):
            arr_sens = np.zeros(all_sens.shape[1])
            x = all_sens[:n_uid,:].copy()
            x_s = (x > absurd_level)*mask[None,:]*(adata['fcode']==f)[:,None]
            x[~x_s] = 1.
            if np.all(~x_s):
                print(' ...no survivors at %s.  Next.' % f)
                continue
            fig, axes = pl.subplots(3, 1, figsize=(6., 7),
                                    sharex=True,sharey=False, squeeze=True)
            arr_sens = (x**-2 * x_s).sum(axis=0)  #1/(uK/rtHz)^2
            arr_sens[arr_sens!=0] = arr_sens[arr_sens!=0] **-.5 / 2**.5 # uK rtsec
            min_dets = cfg['array_sens'].get('min_dets', 50)
            det_count = x_s.sum(axis=0)
            det_count_mask = (det_count >= min_dets)
            print(' surviving n_dets > %i: %i' % (min_dets, det_count_mask.sum()))
            arr_sens[~det_count_mask] = 0.
            a_s = (arr_sens!=0)
            n_eff = match_info(
                info, cfg['array_sens']['standard_det_count'])[f]
            pl.sca(axes[0])
            scargs = {'c': planet_day[a_s],
                      'alpha': .7, 's': 30}
            ylargs = {'fontsize':10}
            pl.title('%s - %s' % (TAG, f))
            pl.scatter(planet_db['loading'][a_s], det_count[a_s], **scargs)
            pl.axhline(n_eff, ls='dotted', color='k')
            pl.colorbar()
            pl.ylabel('N_det', **ylargs)
            pl.sca(axes[1])
            pl.scatter(planet_db['loading'][a_s], arr_sens[a_s], **scargs)
            pl.ylabel('Array sens [uK rtsec]', **ylargs)
            pl.colorbar()
            pl.sca(axes[2])
            sens_eff = arr_sens[a_s] * det_count[a_s]**.5 / n_eff**.5
            pl.scatter(planet_db['loading'][a_s], sens_eff, **scargs)
            pl.xlabel('PWV/sin(alt) [mm]')
            pl.ylabel('Sens @%i dets [uK rtsec]' % n_eff, **ylargs)
            pl.colorbar()
            pl.subplots_adjust(hspace=.1)
            for ax in axes:
                span = ax.get_ylim()
                span = span[1] - span[0]
                delta = span / 4
                distances = []
                for gmod in [1,2,5]:
                    delta1 = gmod * 10**np.round(np.log10(delta/gmod))
                    distances.append((abs(delta1/delta-1), delta1))
                _, delta = min(distances)
                ax.yaxis.set_major_locator(pl.MultipleLocator(delta))
            ofiler.savefig('sens_%s_30_%s_sens_load.png' % (tag, f))

        # Save that?
        s = (red_sens!=0)
        ssave = s*(red_sens > absurd_level)
        if ssave.sum() < 10:
            print(' Only %i reduced sensitivities, skipping.' % ssave.sum())
            continue
        moby2.util.StructDB.from_data(
            [('det_uid', all_uid[:n_uid][ssave]),
             ('sens', red_sens[:n_uid][ssave])],
            formats={'sens': '%9.2f'}).to_fits_table(
                ofiler.get_filename('sens_%s.fits' % tag))

        pl.figure(figsize=(6., 5))

        # Plot them.
        # histogram...
        max_sens = min(np.median(red_sens[s])*5, red_sens[s].max())
        bins = np.arange(0, max_sens, max_sens/40)
        for f,p in zip(all_freqs, 'bgrk'):
            fs = (adata['fcode']==f)
            if not fs.any(): continue
            x = red_sens[fs*s]
            pl.hist(x, bins=bins,
                     label= '%s' % (f), alpha=.8, color=p)
            x = x[x>absurd_level]
            print('  array sens @%s: ' % f, (x**-2).sum()**-.5 / 2**.5, 'uK rtsec')

        pl.legend(loc='upper right')
        pl.xlabel('Sens at 20 Hz (uK/rtHz)')
        pl.title(TAG)
        ofiler.savefig('sens_%s_00_hist.png' % tag)

        # Cumulative detector weight.
        for f,p in zip(all_freqs, 'bgrk'):
            fs = (adata['fcode']==f) * s
            if not fs.any(): continue
            x = red_sens[fs]
            xs = np.sort(x)
            w = xs**-2
            n = np.cumsum(w)**-.5
            pl.semilogy(np.arange(len(n))+1, n,
                     label='%s' % (f), alpha=.8, color=p)
            pl.axhline(n[-1], label='%.1f uK/rtHz' % (n[-1]),
                       color=p, ls='dashed')

        pl.legend(loc='upper right')
        pl.xlabel('N detectors contributing')
        pl.ylabel('Accumulated sensitivity (uK/rtHz)')
        pl.title(TAG)
        ofiler.savefig('sens_%s_02_cumul.png' % tag)


        # Sensitivity by position.
        #fig, axes = pl.subplots(len(all_freqs),1,figsize=(6., 4*len(all_freqs)+1),
        #                        sharex=True,sharey=True, squeeze=False)
        #axes = axes.ravel()
        siz, spa = 30, .015
        for spi, (f,p) in enumerate(zip(all_freqs, 'bgrk')):
            pl.figure(figsize=(5.5,4.5))
            #pl.sca(axes[spi])
            fs = (adata['fcode']==f) * s
            if not fs.any(): continue
            z = red_sens[fs]
            x, y = adata['sky_x'] + spa*(adata['pol_family']=='B'), adata['sky_y']
            max_sens = min(np.median(red_sens[fs])*3, red_sens[fs].max())
            pl.scatter(x[fs], y[fs], s=siz, c=z, vmin=0, vmax=max_sens,
                       cmap='viridis')
            #label='%i GHz' % f)
            #pl.legend(loc='upper right')
            pl.colorbar()
            pl.title('%s - %s' % (TAG, f))
            ofiler.savefig('sens_%s_01_%s_pos.png' % (tag, f))

        print()
