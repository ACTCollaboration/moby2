from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

PLANETCAL_FORMATS = {
    'U_dac': '%11.2f',
    'U_pW': '%11.7f',
    'dac_to_uK': '%12.2f',
    'pW_to_uK': '%12.5e',
}


def get_fpfit_planetcal(args):
    """
    Convert fp_fit to Uranus peak amplitude in pW.
    """

    import moby2
    from moby2.analysis import fp_fit
    from moby2.analysis import det_sens

    import os
    import numpy as np

    import optparse as o
    o = o.OptionParser()
    opts, args = o.parse_args(args)

    cfg = moby2.util.MobyDict.from_file(args[0])
    blist = moby2.scripting.get_tod_list(cfg['tod_list'],
                                         cmdline_args=args[1:])

    mfb = moby2.scripting.get_filebase()

    for b in blist:
        print(b)
        fpf = cfg['fp_fits']['pattern'].format(tod_name=b)
        if not os.path.exists(fpf):
            print('... no fpfit.')
            continue
        if fpf.endswith('.fits') or fpf.endswith('.fits.gz'):
            fp = fp_fit.FPFitFile.from_fits_table(fpf)
        else:
            fp = fp_fit.FPFitFile.from_columns_file(fpf)

        fn = mfb.get_full_path(b)
        if not os.path.exists(fn):
            print('... no TOD.')
            continue
        tod = moby2.scripting.get_tod({'filename': fn, 'end': 400, 'read_data': False})

        ofiler = moby2.scripting.OutputFiler(prefix=cfg['output']['prefix'],
                                             data=tod.info.get_dict())
        ofile = ofiler.get_filename('pcal.fits')
        if cfg['output'].get('skip_existing', False) and os.path.exists(ofile):
            print('... output exists, skipping.')
            continue

        # Get the calibration applied in fp_fit, and our desired calibration to pW.
        s0 = fp.ok.astype(bool)
        cal1 = moby2.scripting.get_calibration(cfg['fp_fits']['cal'], tod=tod)
        cal2 = moby2.scripting.get_calibration(cfg['output']['cal'], tod=tod)
        s1, calv1 = cal1.get_property('cal', det_uid=fp.det_uid)
        s2, calv2 = cal2.get_property('cal', det_uid=fp.det_uid)

        s = (s0*s1*s2*(calv1!=0)*(calv2!=0))
        print(' valid fp output: %i' % s0.sum())
        print(' valid fp cal: %i' % (s0*s1).sum())
        print(' valid pW cal: %i' % (s0*s2).sum())
        print(' total out:    %i' % s.sum())

        if s.sum() == 0:
            print(' skipping!')
            continue

        uid = fp.det_uid[s]
        h = fp.h[s]
        h_dac = h / calv1[s]
        h_cal = h_dac * calv2[s]

        f = tod.info.array_data['nom_freq'][fp.det_uid[s]]
        ar = tod.info.array
        t = tod.info.ctime

        h_uK = h_cal*0

        # Include a temperature model?
        source_name = cfg.get('source_name')
        print('Source is %s' % (source_name))
        if source_name is not None:
            Tmodel = det_sens.planet_cal.get_brightness_model(source_name)
            if Tmodel is None:
                raise ValueError("No model for calibration source %s" % source_name)
            source_omega = det_sens.planet_cal.get_planet_solid_angle(
                    source_name, t)
            nom_freqs = sorted(list(set(f)))
            for nom in nom_freqs:
                s = (f==nom)*(f!=0)
                if not np.any(s):
                    continue
                beam_omega = cfg['beam_solid_angle_sr'][ar][nom]
                band_center = cfg['band_centers'][ar][nom]
                print(' For f_nom=%.1f / f_cen=%.1f' % (nom, band_center))
                print('  Beam omega = %.1f nsr' % (beam_omega*1e9))
                print('  Planet omega = %.5f nsr' % (source_omega*1e9))
                source_T = Tmodel.get_uK_minus_background(band_center)
                print('  Planet T = %.1f K_cmb' % (source_T/1e6))
                peak_uK = source_T * source_omega / beam_omega
                print('  Peak brightness = %.6f K_cmb' % (peak_uK/1e6))
                h_uK[s] = peak_uK

        print()
        
        output = moby2.util.StructDB.from_data(
            [('det_uid', uid),
             ('U_dac', h_dac),
             ('U_pW', h_cal),
             ('dac_to_uK', h_uK/h_dac),
             ('pW_to_uK', h_uK/h_cal),
             ],
            formats=PLANETCAL_FORMATS)

        output.to_fits_table(ofile)

        # Write to a cal archive?
        cal_archive_file = cfg['output'].get('hdf_archive')
        if cal_archive_file is not None:
            with moby2.detectors.CalibrationArchive(cal_archive_file, 'w') as carc:
                calo = moby2.detectors.Calibration(output['det_uid'])
                calo.cal[:] = output['pW_to_uK']
                carc.set_item(b, calo, clobber=True)

        # Plot... ?
        pprefix = cfg['output'].get('plot_prefix')
        if pprefix is None:
            continue

        import pylab as pl
        pfiler = moby2.scripting.OutputFiler(prefix=pprefix,
                                             data=tod.info.get_dict())
        fig, sps = pl.subplots(2,1)
        for spi, x in enumerate([abs(h_cal), abs(h_uK/h_cal)/1e6]):
            x = abs(x)  # is this the right sign convention?
            bins = np.linspace(0, min(x.max(), np.median(x)*5), 30)
            for nom,col in zip(nom_freqs, 'kbgr'):
                label='f=%iGHz' % nom
                s = (f==nom)
                pl.sca(sps[spi])
                pl.hist(x[s], bins=bins, color=col,
                        label=label, alpha=.8)
                pl.legend(loc='upper right')
        pl.xlabel('Calibration factors (K/pW)')
        pl.sca(sps[0])
        pl.xlabel('%s Peak Height (pW)' % source_name.capitalize())
        pl.title('%s - %s' % (b, moby2.util.ctime.to_string(tod.info.ctime)))
        pl.tight_layout()
        pfiler.savefig('cal.png')


def get_map_planetcal(args):
    """This is similar to get_fpfit_planetcal, in that it produces
    per-detector absolute calibration data for a single planet TOD.
    However, instead of using per-detector amplitude fits, it combines
    with the amplitude of the planet in an all-detector planet map
    with the (atmospheric) calibration info that were used to make
    that map.

    """
    import moby2
    from moby2.analysis import det_sens
    from moby2.analysis import beam_ana

    import os
    import numpy as np

    import optparse as o
    o = o.OptionParser()
    opts, args = o.parse_args(args)

    cfg = moby2.util.MobyDict.from_file(args[0])
    blist = moby2.scripting.get_tod_list(cfg['source_tods'],
                                         cmdline_args=args[1:])

    source_name = cfg['planet_maps']['source_name']
    Tmodel = det_sens.planet_cal.get_brightness_model(source_name)
    if Tmodel is None:
        raise ValueError("No model for calibration source %s" % source_name)

    for b in blist:
        print(b)
        ti = moby2.scripting.get_tod_info({'filename': b})
        adata = moby2.scripting.get_array_data({}, tod_info=ti)
        # Get cuts
        cuts = moby2.scripting.get_cuts(cfg['planet_maps']['cuts'], tod=b)
        cal = moby2.scripting.get_calibration(cfg['planet_maps']['calibration'], tod=b)
        cal_mask, calv = cal.get_property('cal', det_uid=adata['det_uid'])
        calv[~cal_mask] = 0.
        cal_mask[:] = False
        cal_mask[cuts.get_uncut(det_uid=True)] = True
        cal_mask *= (calv!=0)
        calv = calv * cal_mask
        print(b, cal_mask.sum())
        peaks_cal = calv * 0
        peaks_uK = calv * 0

        source_omega = det_sens.planet_cal.get_planet_solid_angle(
            source_name, ti.ctime)

        freq_codes = sorted([('f%03i' % f, f) for f in list(set(adata['nom_freq'])) if f != 0])
        for fc, f in freq_codes:
            map_file = cfg['planet_maps']['map_file'].format(
                tod_id=b, tod_name=b, freq=fc)
            print(map_file, os.path.exists(map_file))
            bo = beam_ana.BeamObs(mapfile=map_file)
            bo.fit()
            peaks_cal[(adata['nom_freq'] == f)] = bo['peak']['amp']

            ar = ti.array
            beam_omega = cfg['beam_solid_angle_sr'][ar][f]
            band_center = cfg['band_centers'][ar][f]
            source_T = Tmodel.get_uK_minus_background(band_center)
            print('  Planet T = %.1f K_cmb' % (source_T/1e6))
            peak_uK = source_T * source_omega / beam_omega
            print('  Peak brightness = %.6f K_cmb' % (peak_uK/1e6))
            peaks_uK[adata['nom_freq'] == f] = peak_uK
            
        # The cal to "pW" includes the filter gain removal, but we do
        # not want that in our DAC to PW conversion.  Same goes for
        # optical_sign; we want that as part of the pW calibration
        # rather than the pW->uK numbers.
        tod = moby2.scripting.get_tod({'filename': b, 'use_filebase': True,
                                       'end': 100, 'det_uid': [0]})
        cal_G = moby2.scripting.get_calibration(
            [{'type': 'remove_readout_filter_gain'},
             {'type': 'array_data_parameter', 'parameter': 'optical_sign'}], tod=tod,
            det_uid=adata['det_uid'])
        print(cal_G.det_uid)

        #print list(peaks_uK[peaks_cal!=0] / peaks_cal[peaks_cal!=0] )
        ## det_uid  U_dac        U_pW    dac_to_uK     pW_to_uK
        #0       14.31  -0.0026499      2190.60 -1.18282e+07
        s = (peaks_cal != 0)*(peaks_uK != 0)*(calv != 0)
        cal_G = cal_G.cal[s]
        peaks_dac = peaks_cal[s] / (calv[s] / abs(cal_G))
        output = moby2.util.StructDB.from_data(
            [('det_uid', adata['det_uid'][s]),
             ('U_dac', peaks_dac),
             ('U_pW', abs(peaks_cal[s])),
             ('dac_to_uK', peaks_uK[s] / peaks_dac),
             ('pW_to_uK', peaks_uK[s] / abs(peaks_cal[s]))],
            formats=PLANETCAL_FORMATS)

        ofiler = moby2.scripting.OutputFiler(prefix=cfg['output']['prefix'],
                                             data=tod.info.get_dict())
        ofile = ofiler.get_filename('pcal.fits')
        output.to_fits_table(ofile)
        
def get_cal_noise(args):
    """
    Combine noise spectra with an absolute calibration to get
    calibrated detector noise levels in some band(s).
    """

    import moby2
    import numpy as np
    import os

    from moby2.analysis import det_sens

    import optparse
    o = optparse.OptionParser()
    opts, args = o.parse_args(args)

    cfg = moby2.util.MobyDict.from_file(args[0])
    blist = moby2.scripting.get_tod_list(cfg['cal_noise']['tod_list'],
                                         cmdline_args=args[1:])
    ofiler = moby2.scripting.OutputFiler(prefix=cfg['cal_noise']['sens_prefix'])

    cal_file_key = cfg['cal_noise'].get('planet_cal_key', 'pW_to_uK')
    if cal_file_key == 'uK_to_uK':
        print('Not using planet cal; assumes data are in uK already.')
        print()

    for b in blist:
        print(b)
        # Get calibration
        if cal_file_key == 'uK_to_uK':
            cal = None
        else:
            cal_file = cfg['cal_noise']['planet_cal'].format(name=b)
            if not os.path.exists(cal_file):
                print('Not found: ', cal_file)
                continue
            cal = moby2.util.StructDB.from_fits_table(cal_file)
        # Get spectrum
        spec_file = cfg['cal_noise']['tod_spectra'].format(name=b)
        if not os.path.exists(spec_file):
            print('Not found: ', spec_file)
            continue
        spec = det_sens.TODSpectrum.from_file(spec_file)
        # Frequency band.
        band_data = []
        for band_name, band_lo, band_hi in cfg['cal_noise']['freq_bands']:
            fs = (band_lo <= spec.f) * (spec.f < band_hi)
            noise = (spec.spec[:,fs]**2).mean(axis=1)**.5
            band_data.append((band_name, noise))
        band_names, noises = list(zip(*band_data))
        if len(spec.meta['det_uid']) == 0:
            print('No dets: ', spec_file)
            continue
        max_uid = spec.meta['det_uid'].max() + 1
        if cal is not None:
            max_uid = max(max_uid, cal['det_uid'].max()+1)
        all_uid = np.arange(0, int(max_uid))
        sens = np.zeros((len(noises), len(all_uid)))
        sens[:,spec.meta['det_uid']] = noises
        if cal is not None:
            sens[:,cal['det_uid']] *= abs(cal[cal_file_key])
        _s = (sens[0]!=0)
        formats = {}
        for n in band_names:
            formats[n] = '%11.1f'

        ofile_fits = ofiler.get_filename(b+'.fits')
        ofile_txt = ofiler.get_filename(b+'.txt')
        if _s.sum() == 0:
            print('No calibrated results.')
            for f in [ofile_fits, ofile_txt]:
                if os.path.exists(f):
                    print('(removing %s)' % f)
                    os.remove(f)
                continue
        data = [('det_uid', all_uid[_s])]
        for name, sen in zip(band_names, sens):
            data.append((name, sen[_s]))
        odb = moby2.util.StructDB.from_data(data, formats=formats)
        odb.to_fits_table(ofile_fits)
        odb.to_column_file(ofile_txt)
