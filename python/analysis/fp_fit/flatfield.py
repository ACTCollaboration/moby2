from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import os
import numpy as np
import pylab as pl

def load_reduced_fpfit(template, tod_id):
    """
    Load the reduced fp_fit data for this observation.
    """
    fn = template.format(tod_id=tod_id)
    if not os.path.exists(fn):
        return None
    if fn.endswith('.txt'):
        return moby2.util.StructDB.from_column_file(
            fn, [('det_uid',0), ('peak_dac',1), ('SN',2), ('tau',3)])
    return moby2.util.StructDB.from_fits_table(fn)

def load_cal(tod_name, cfg):
    try:
        cal = moby2.scripting.get_calibration(
            {'type': 'depot_cal',
             'tag': cfg['local_cal_tag'],
             'depot': cfg['local_cal_depot']},
            tod=tod_name)
        return cal
    except:
        pass
    print('Making cal result for %s...' % tod_name)
    tod = moby2.scripting.get_tod({'filename': tod_name, 'det_uid': [0], 'end': 400,
                                   'use_filebase': True})
    ndets = tod.info.array_data.ndets
    try:
        cal = moby2.scripting.get_calibration(
            cfg['local_cal'],
            tod=tod, det_uid=np.arange(ndets))
    except:
        print('Warning:  FAILED TO MAKE CAL for %s' % tod_name)
        return
    depot = moby2.scripting.get_depot({'path': cfg['local_cal_depot']})
    depot.write_object(cal, tag=cfg['local_cal_tag'], tod=tod_name)
    return cal

def limit_mask(vals, med_mul=3):
    med = np.median(vals[vals!=0])
    med_lim = min(med * med_mul, max(vals) + med * 0.2)
    return [(0, med_lim), (vals <= med_lim)]

def get_fiducial_main(args_in):
    from optparse import OptionParser as o
    o = o()
    o.add_option('-k','--key', default='get_fid')
    opts, args = o.parse_args(args_in)
    cfg = moby2.util.MobyDict.from_file(args[0])[opts.key]

    ulist = moby2.util.ascii.read_columns(cfg['tod_list'])[0]
    if cfg.get('blacklist') != None:
        ulist = sorted(list(set(ulist).difference(
            moby2.util.ascii.read_columns(cfg['blacklist'])[0])))

    # Loop through TODs; get the fp_fit and the calibration.  Reduce
    # them, and compute a few stats for cutting downstream.

    rows = []
    rms_all = []
    for b in ulist:
        print(b, end=' ')
        fp = load_reduced_fpfit(cfg['fpfit_reduced'], b)
        if fp is None:
            print(' ... no fp')
            continue
        if len(fp) == 0:
            print(' ... no survivors in fp')
            continue
        cal = load_cal(b, cfg)
        if cal is None or len(cal.cal) == 0:
            print(' ... no calibrated dets.')
            continue
        cal_mask, calv = cal.get_property('cal', det_uid=fp['det_uid'])
        cal_mask *= (calv != 0)
        fp_mask = (fp['peak_dac'] != 0)
        amps = fp['peak_dac'] * calv * cal_mask
        mask = amps!=0
        # At this point, we are keeping it.
        if (amps[mask]<0).sum() > mask.sum() / 10.:
            print('... lots of dets with negative amplitude... apply '
                  'optical cal?')
        else:
            print(' ... ok')
        rows.append([b, fp['det_uid'], amps, cal_mask, fp_mask])

    print()
    print('Working with %i observations.' % len(rows))
    print()

    rms_cut = cfg.get('rms_cut', 1e3)

    # Second pass -- do frequency bands individually.  Get the
    # frequency band assignment.
    tod_info = moby2.scripting.products.get_tod_info_new(
        {'tod_id': rows[0][0]})
    adata = moby2.scripting.products.get_array_data_new(
        {'tod_id': rows[0][0]})
    pretitle = '{array} {season}'.format(**tod_info)

    # Permit output manager to refer to tod_info items...
    outputm = moby2.scripting.OutputFiler(prefix=cfg.get('output_prefix','./'),
                                          data=tod_info)

    nom_freqs = adata['nom_freq']
    all_nom_freqs = sorted([x for x in set(nom_freqs) if x!=0])

    n_det = len(nom_freqs)
    all_amps = np.zeros((len(rows), n_det))
    inputs_ok = np.zeros((len(rows), n_det, 2), bool)

    for nom_freq in all_nom_freqs:
        fcode = 'f%03i' % nom_freq
        print()
        print('Frequency band: ', fcode)
        print()
        print('  TOD_band                        n_det        amp     rms')
        
        stats = []
        for (name, uids, amps, cal_mask, fp_mask), out, inp in \
            zip(rows, all_amps, inputs_ok):
            s0 = (nom_freqs[uids] == nom_freq) # channel select
            inp[uids[s0],0] = fp_mask[s0]
            inp[uids[s0],1] = cal_mask[s0]
            s = s0 * cal_mask * fp_mask
            if not s.any(): continue
            amps = amps[s]
            # Outlier rejection...
            bounds = np.percentile(amps[amps!=0], [10,90])
            s2 = (bounds[0] <= amps) * (amps <= bounds[1])
            a0, da =  np.median(amps[s2]), np.std(amps[s2])
            r = da/a0
            print('  %s_%s %6i   %8.3f %8.3f' % (name, fcode, (amps!=0).sum(), a0, r))
            if (r > rms_cut):
                amps *= 0
                print(' ... cut for high rms.')
            if np.isnan(r) or r < 0:
                print(' ... forcing fraction rms = 1.')
                r = 1.
            out[uids[s]] = amps
            stats.append((a0, da, r, ))

        # Plot showing the rmses
        rms_all = np.array(stats)[:,2]
        print(rms_all)
        pl.scatter(np.arange(len(rms_all)), np.log10(rms_all))
        pl.axhline(np.log10(rms_cut), ls='dotted')
        pl.title('Use this to tune the rms_cut.')
        pl.xlabel('Observation index')
        pl.ylabel('log10(rms) of the obs')
        outputm.savefig('10_rms_hist_%s.png' % fcode)

    # Now... compute median and scatter of each detector's valid obs.
    fits = np.zeros((n_det, 3))
    for i, fit_out in enumerate(fits):
        amps = all_amps[:,i]
        s = amps!=0
        fit_out[0] = s.sum() # count
        if fit_out[0] == 0:
            continue
        fit_out[1:3] = (amps[s].std(), amps[s].mean())

    ok = (fits[:,2]!=0)
    fits = moby2.util.StructDB.from_data([
        ('det_uid', np.arange(n_det)),
        ('n_obs', fits[:,0].astype(int)),
        ('fp_fails', (inputs_ok[:,:,1] * ~inputs_ok[:,:,0]).sum(axis=0)),
        ('cal_fails', (inputs_ok[:,:,0] * ~inputs_ok[:,:,1]).sum(axis=0)),
        ('scat',  fits[:,1]),
        ('amp',   fits[:,2]),
        ('frms',  (fits[:,1]*ok) / (fits[:,2] + ~ok)),
        ('flat_ok', ok),
        ('flat_recal', np.zeros(len(ok))),
        ('is_fid', np.zeros(len(ok),bool))],
        formats={
            'scat': '%7.5f', 'amp': '%8.3e', 'frms': '%7.4f',
            'flat_recal': '%7.5f'})

    # Enforce a minimum attendance to qualify for the flat field.
    min_n_obs = cfg.get('min_n_obs', 5)
    fits['flat_ok'] *= (fits['n_obs'] >= min_n_obs)

    # Sub-enforce a minimum fractional attendance to quality for fiducial.
    min_frac_obs = cfg.get('min_frac_obs', 0.8)
    fits['is_fid'][:] = fits['flat_ok'] * (fits['n_obs'] >= min_frac_obs * max(fits['n_obs']))

    # Load the selection limits.
    #
    sel_lims = {}
    for entry in cfg['fiducial_selection']:
        if 'nom_freq' in entry:
            sel_lims[entry['nom_freq']] = entry
        else:
            sel_lims['_default'] = entry

    for nf in all_nom_freqs:
        s = (nom_freqs == nf)*fits['flat_ok']
        print('Considering %i for band %s' % (s.sum(), nf))
        qual = fits['frms'][s]
        pl.hist(qual, alpha=.4, bins=np.arange(.01, .5, .01),
                label='F=%3i GHz' % nf)
        lims = sel_lims.get(nf, sel_lims.get('_default', {}))
        for x in lims.get('frms', []):
            pl.axvline(x, ls='dotted')

    pl.xlabel('Scatter parameter')
    pl.ylabel('Number of dets')
    pl.title(pretitle)
    pl.legend(loc='upper right')
    outputm.savefig('20.png')

    #
    # Scatter plot of fractional RMS vs. amplitude.
    #
    for nf in all_nom_freqs:
        s = nom_freqs == nf
        s1 = (~fits['is_fid'])[s]
        amp, frms = fits['amp'][s], fits['frms'][s]
        xlim, xmask = limit_mask(amp)
        ylim, ymask = limit_mask(frms)
        pl.scatter(amp[~s1], frms[~s1], alpha=.4, s=4, label=
                   'showing %i/%i points' % ((xmask*ymask).sum(), len(amp)))
        pl.scatter(amp[s1], frms[s1], alpha=.3, s=2, label=
                   '(%i are ineligible; low n_obs)' % s1.sum())
        lims = sel_lims.get(nf, sel_lims.get('_default', {}))
        for x in lims.get('amp', []):
            pl.axvline(x, ls='dotted')
        for y in lims.get('frms', []):
            pl.axhline(y, ls='dotted')
        pl.xlabel('Mean amp')
        pl.ylabel('Fractional rms')
        pl.xlim(*xlim)
        pl.ylim(*ylim)
        pl.legend(loc='best', fontsize=10)
        pl.title(pretitle + ' %3ighz - Fiducial det selection' % nf)
        outputm.savefig('21_f%03d.png' % nf)

    # In each band: mark the fiducial detectors, and compute flatfield.
    print('Summary:')
    for nf in all_nom_freqs:
        s = (nom_freqs == nf) * fits['flat_ok']
        is_fid = s * fits['is_fid']
        lims = sel_lims.get(nf, sel_lims.get('_default', {}))
        amp, frms = fits['amp'], fits['frms']
        for k in ['amp', 'frms']:
            L, x = lims[k], fits[k]
            is_fid *= (L[0] <= x) * (x < L[1])
        # amp0: the amplitude we have decided is representative, and
        # will serve as a target for the flatfield.
        amp0 = fits['amp'][is_fid].mean()
        # recal is the factor that takes all amps to amp0.
        recal = (amp0 * s) / (fits['amp'] + (~s))
        fits['flat_recal'][s] = recal[s]
        fits['is_fid'][s] = is_fid[s]
        print(' Frequency: ', nf)
        for k,v in [('detectors', (nom_freqs==nf).sum()),
                    ('flatfielded', s.sum()),
                    ('stable', fits['is_fid'][s].sum())]:
            print('   %-20s:  %6i' % (k,v))
    print()

    x0 = adata['sky_x']
    y0 = adata['sky_y']
    sep = 1./60
    x0 = x0 + (adata['pol_family']=='A')*sep

    for nf in all_nom_freqs:
        s = (nom_freqs == nf) * fits['flat_ok']
        pl.scatter(x0[s], y0[s], c=1./fits['flat_recal'][s],
                   vmin=0, vmax=2, s=30, label='flat field')
        pl.colorbar()
        not_fid = s * (~fits['is_fid'])
        if np.any(not_fid):
            pl.scatter(x0[not_fid], y0[not_fid], marker='x', color='k',
                       s=20, label='non fiducial')
        pl.legend(loc='best', fontsize=10)
        pl.title(pretitle + ' %3ighz - Relative efficiency' % nf)
        outputm.savefig('30_f%03d_flatfield.png' % nf)

    # Write out full db...
    fits.to_column_file(outputm.get_filename('stats.txt'))

    # Write out simple flatfield + fiducial flag.
    det_uid = np.arange(n_det)
    s = fits['flat_ok']
    db = moby2.util.StructDB.from_data([
            ('det_uid', det_uid[s].astype('int16')),
            ('cal', fits['flat_recal'][s].astype('float32')),
            ('stable', fits['is_fid'][s].astype('uint8'))])
    db.to_fits_table(outputm.get_filename('flatfield.fits'))

    if cfg.get('common_mode_check') is None:
        return

    ### The common mode check might be broken now, as stuff above this
    ### point has been updated significantly.
    ###

    print()
    print('Common mode check...')
    print('# TOD                    ', end=' ')
    for nf in all_nom_freqs:
        print('n%03i %6s %6s' % (nf, 'mean', 'rms'), end=' ')
    print()

    cmcfg = cfg['common_mode_check']
    depot = cmcfg['planet_depot']
    for b in tod_names:
        tod = moby2.scripting.get_tod({'filename': b, 'det_uid': [0],
                                       'end': 400, 'use_filebase': True})
        ndets = tod.info.array_data.ndets
        _uid = np.arange(ndets)
        cal0 = moby2.scripting.get_calibration(
            cmcfg['cal_in'], tod=tod, det_uid=_uid)
        cal1 = moby2.scripting.get_calibration(
            cmcfg['cal_out'], tod=tod, det_uid=_uid)
        fn = depot + cmcfg['cm_file'].format(tod_name=b, first_five=b[:5])
        modes = moby2.tod.TODModeSet.from_fits_file(fn)
        n_mo = len(modes.modes)
        idx = db.select_inner({'det_uid': modes.det_uid})
        s = (idx>=0)*(cal0.cal!=0)*(cal1.cal!=0) * (db['stable'][idx]==1)
        # We expect cal1 / cal0 to be the same as the flatfield value.
        print(b, end=' ') 
        for nf in all_nom_freqs:
            s1 = (nom_freqs[idx] == nf) * s
            r = cal0.cal[s1] / cal1.cal[s1] * db['cal'][idx[s1]]
            print('%4i %6.3f %6.3f' % (len(r), r.mean(), r.std()), end=' ')
        print()

