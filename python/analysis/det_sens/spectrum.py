from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np
import os

class TimeConstantTransfer:
    def __init__(self, time_constants):
        self.tau = time_constants
    def get_transfers(self, f, index=None, power=True):
        # Positive sign in denominator sends this clockwise in C as f
        # increases; that's in the same direction as the MCE filter.
        tr = 1. / (1. + 2.j*np.pi*self.tau[:,None]*f)
        if power:
            return abs(tr)**2
        return tr

class MCETransfer:
    def __init__(self, mcef, n_det):
        self.mcef = mcef
        self.n_det = n_det
    def get_transfers(self, f, index=None, power=True):
        output = np.empty((self.n_det, len(f)), np.complex128)
        output[:,:] = self.mcef.transfer(f) / self.mcef.transfer(0)
        if power:
            return abs(output)**2
        return output

class FourierBinner:
    """
    Re-bin 1d power spectrum (or other spectra) with logarithmic bin
    spacing.
    """
    def __init__(self, n, n_fft=None, f_samp=1., res=200, f_min=None):
        """
        Create a binner for rebinning vectors of length n.  FFTs will
        automatically truncate to the largest power of 2 that fits
        into n, unless n_fft is specified.  The frequency vector is
        generated based on sampling frequency f_samp.
        self.set_res(res) will be called to set the logarithmic
        binning.
        """
        if n_fft is None:
            n_fft = int(2**np.floor(np.log2(n)))
        self.n, self.n_fft = n, n_fft
        self.f = np.fft.fftfreq(n_fft) * f_samp
        self.set_res(res)

    def set_bins(self, bins):
        self.bins = np.array(bins)
        self.counts, _ = np.histogram(self.f, bins=self.bins)
        self.f_bin, _ = np.histogram(self.f, bins=self.bins,
                                     weights=self.f)
        self.mask = self.counts>0
        self.f_bin = self.f_bin[self.mask] / self.counts[self.mask]
        self.idx = np.digitize(self.f, self.bins).astype('int32') - 1
        self.idx[self.idx>=len(self.bins)-1] = -1
        
    def set_res(self, n_bin, f_min=.01, f_max=None):
        """
        Create n_bin logarithmically spaced bins, between f_min and
        f_max (which defaults to the Nyquist frequency).
        """
        if f_max is None:
            f_max = self.f.max()
        bins = 10**np.linspace(np.log10(f_min), np.log10(f_max), n_bin)
        self.set_bins(bins)

    def rebin(self, spectra):
        """Rebin spectra, shape (n_det, n_freq) or (n_freq).  The input array
        can be float32 or complex64.  Returns re-binned spectra, which
        will be float64 or complex128, depending on the input.

        """
        cplx = np.iscomplexobj(spectra)
        if cplx:
            spectra = np.asarray(spectra, 'complex64')
        else:
            spectra = np.asarray(spectra, 'float32')
        if spectra.ndim < 2:
            return self.rebin_accel(spectra[None,:])[0]
        n_spec = spectra.shape[0]
        n_bin = len(self.counts)
        if cplx:
            output = np.zeros((n_spec, n_bin), 'complex')
            moby2.libactpol.bin_data_complex(self.idx, spectra, output)
        else:
            output = np.zeros((n_spec, n_bin), 'float')
            moby2.libactpol.bin_data(self.idx, spectra, output)
        return output[:,self.mask] / self.counts[self.mask]

    rebin_accel = rebin


class TODTransformer:
    """
    The main driver class for getting TOD noise power spectra.
    """
    def __init__(self, params):
        self.p = params
        self.log = moby2.util.log.get_logger()
        self.hwp_angles = None
        self.hwp_reconstructor = None

    def init_from_data(self, data, f_samp=None, cuts=None, mask=None):
        self.transfers = []
        self.f_samp = f_samp
        self.cuts = cuts
        if mask is None:
            mask = np.ones(data.shape[-2], bool)
        self.mask = mask
        self.frb = FourierBinner(data.shape[-1], f_samp=self.f_samp,
                                 res=self.p.get('res', 200))
        n_fft = self.frb.n_fft
        # Including factor of 2 here because all power must come from
        # interval between 0 and the Nyquist frequency.  Output
        # spectra are then per Hz.
        self.norm = (n_fft*self.f_samp/2)**-.5
        
    def init_from_tod(self, tod_file):
        """
        Note that this currently requires you to have initialized the
        object with params that include cuts, calibration, and time
        constants.
        """
        # Load cuts first
        self.log.trace(1, 'Loading cuts')
        cuts = moby2.scripting.get_cuts(self.p['cuts'], tod=tod_file)
        if self.p['source_cuts'] is not None:
            cuts_p = moby2.scripting.get_cuts(self.p['source_cuts'],
                                              tod=tod_file)
            cuts.merge_tod_cuts(cuts_p)
        # Choose a good length
        fft_size = self.p.get('fft_size', 2**16)
        n_fft = min(int(2**np.floor(np.log2(cuts.nsamps))), fft_size)
        self.log.trace(2, 'Will process %i samples' % n_fft)
        # Load TOD.
        start = cuts.sample_offset
        end = start + n_fft
        self.log.trace(1, 'Loading tod %s' % tod_file)
        tod = moby2.scripting.get_tod({
                'filename': tod_file,
                'cuts': cuts,
                'start': start,
                'end': end})
        # Store the adjusted cuts that are loaded into tod.  But kill
        # sample_index, or fill_cuts will be messed up.
        self.cuts = tod.cuts.copy()
        self.cuts.sample_offset = 0
        # Load and apply calibration.
        self.log.trace(1, 'Loading calibration...')
        calv = moby2.scripting.products.get_calibration(
            self.p['calibration'],
            tod.info, det_uid=tod.det_uid)
        ok, recal = calv.get_property('cal', det_uid=tod.det_uid)
        tod.data *= (recal*ok)[:,None]
        ok *= (recal!=0)
        self.cal = recal
        # Time constants
        self.transfers = []
        has_time_constants = (self.p['time_constants'] is not None)
        if has_time_constants:
            self.log.trace(1, 'Loading time constants...')
            tau = moby2.scripting.get_time_constants(
                self.p['time_constants'],
                tod_info=tod.info)
            ok_tau, tau = tau.get_property('tau', det_uid=tod.det_uid)
            ok = ok*ok_tau
            self.tau = tau
            self.transfers.append(TimeConstantTransfer(self.tau))
        else:
            self.log.trace(1, 'No time constants.')
            self.tau = None
        # Also the MCE readout filter.
        self.log.trace(1, 'Setting up transfer functions')
        readout_filter = self.p.get('readout_filter', 'mce')
        if readout_filter in [True, 'mce']: 
            mcef = tod.info.mce_filter
            self.transfers.append(MCETransfer(mcef, len(tod.det_uid)))
        elif readout_filter in [False, None]:
            pass
        else:
            raise ValueError('Invalid readout_filter = %s' % readout_filter)

        # hwp stuff
        if 'hwp_angles' in self.p:
            self.hwp_angles = moby2.scripting.get_hwp_angles(
                self.p['hwp_angles'], tod=tod)
            self.hwp_angles *= np.pi/180 # convert to radians.
        if 'hwp_modes' in self.p:
            from moby2.analysis import hwp
            achi_params = self.p['hwp_modes']
            depot = moby2.scripting.get_depot(achi_params.get('depot'))
            achi_filename = depot.get_full_path(hwp.HWPModes, tag=achi_params.get('tag'),
                                                tod=tod)
            modes = hwp.HWPModes.from_fits_table(achi_filename)
            modes, modes_ok = modes.extract(det_uid=tod.det_uid)
            # here we go.
            self.hwp_reconstructor = modes.get_reconstructor(self.hwp_angles)
            ok *= modes_ok
                        
        # Update cuts with the cal/time constant mask?
        self.mask = ok
        if not np.all(self.mask):
            tod.data[~self.mask,:] = 0.

        # Get sample rate -- note median is not robust and we should not trust
        # a simple mean.
        dt = np.diff(tod.ctime)
        dt_median = np.median(dt[::100])
        dt_cut = (dt_median*.5 < dt) * (dt < dt_median*1.5)
        # FFT prepare...
        self.f_samp = 1/np.mean(dt[dt_cut])
        self.log.trace(2, 'Sampling rate is %.5f' % self.f_samp)
        assert (self.f_samp > 50) and (self.f_samp < 1000)
        # Factor of 2 puts units to 1/rt(Hz).
        self.norm = (n_fft*self.f_samp/2)**-.5
        self.frb = FourierBinner(tod.nsamps, n_fft=n_fft, f_samp=self.f_samp,
                                 res=self.p.get('res', 200))
        return tod

    def preprocess_data(self, data):
        # Remove common mode to fill cuts; then reestimate common mode:
        self.log.trace(1, 'Pre-processing (%i,%i)' % tuple(data.shape))
        if self.cuts is not None:
            data -= data.mean(axis=1)[:,None]
            self.log.trace(1, 'Filling cuts to compute common mode...')
            data_copy = data.copy()
            moby2.tod.fill_cuts(data=data, cuts=self.cuts, extrapolate=False)
            cm = np.dot(self.mask, data) / self.mask.sum()
            self.log.trace(1, 'Restoring data and removing common mode...')
            data[:] = data_copy - cm
            self.log.trace(1, 'Filling cuts for real.')
            moby2.tod.fill_cuts(data=data, cuts=self.cuts, extrapolate=False)
            data += cm
        if self.hwp_reconstructor is not None:
            for i in self.mask.nonzero()[0]:
                data[i] -= self.hwp_reconstructor.get_achi(i) * self.cal[i]
        data[~self.mask] = 0.
        data -= data.mean(axis=1)[:,None]
        if self.p.get('get_spec', {}).get('detrend'):
            moby2.tod.detrend_tod(data=data)

        if self.p.get('get_spec', {}).get('test_rms'):
            data[:] = np.random.normal(scale=self.p['get_spec']['test_rms'],
                                       size=data.shape)

    def get_spectrum(self, data):
        # Window all FFTs.
        n_samps = data.shape[-1]
        window = np.sin(np.linspace(0, 1, n_samps)*np.pi).astype('float32')
        renorm = self.norm * (window**2).mean()**-.5  # Use renorm to undo window effects.

        # Get spectra and rebin power
        self.log.trace(1, 'FFT...')
        #spectra_w = np.fft.fft(data*window)
        spectra_w = moby2.libactpol.fft_f_r(data*window, np.arange(len(data)), 1)
        self.log.trace(1, 'Rebinning.')
        spec_w = self.frb.rebin_accel(abs(spectra_w)**2) * renorm**2

        freq = self.frb.f_bin
        # Deconvolve various things.
        for t in self.transfers:
            spec_w /= t.get_transfers(freq)

        speco = TODSpectrum(f=freq.astype('float32'),
                            spec=(spec_w**.5).astype('float32'),
                            counts=self.frb.counts[self.frb.mask])
        return speco

    def get_cross(self, data, idx=None):
        # Window all FFTs.
        n_samps = data.shape[-1]
        window = np.sin(np.linspace(0, 1, n_samps)*np.pi).astype('float32')
        renorm = self.norm * (window**2).mean()**-.5  # Use renorm to undo window effects.

        if idx is None:
            idx = range(0, data.shape[0])

        # Get spectra and rebin power
        self.log.trace(1, 'FFT...')
        spectra_w = np.fft.fft(data[idx,:]*window)
        self.log.trace(1, 'Rebinning.')
        n_det = len(idx)
        freq = self.frb.f_bin
        xspecs = np.zeros((n_det,n_det,len(freq)), np.complex128)
        for i in range(n_det):
            print(i)
            xc = spectra_w[i].conj()*spectra_w[i:]
            xspecs[i,i:] = (self.frb.rebin_accel(xc.real) + 
                            1.j * self.frb.rebin_accel(xc.imag))
            xspecs[i,:i] = xspecs[:i,i].conj()

        # Deconvolve various things.
        xspecs *= renorm**2
        for t in self.transfers:
            trans = t.get_transfers(freq, power=False)
            xspecs /= trans[idx,None,:].conj()
            xspecs /= trans[None,idx,:]

        return freq, xspecs


class TODSpectrum:
    """
    Conventions: we will store the spectrum amplitude, not the
    spectrum power.  For TODs, we will use meta['name'] for the
    tod_name/basename, meta['det_uid'] for the detector ids,
    meta['cal'] for the calibration vector used on the TOD and
    meta['cal_units'] has a string describing the calibration units.
    (The spectrum itself should be in cal_units/sqrt(Hz).)
    """
    def __init__(self, f, spec, **kw):
        self.f = f
        self.spec = spec
        self.meta = kw

    def to_file(self, filename):
        np.savez(filename,
                 freq=self.f, spec=self.spec, meta=self.meta)

    @classmethod
    def from_file(cls, filename):
        data = np.load(filename)
        self = cls(data['freq'], data['spec'], **data['meta'].tolist())
        return self
        
 
def get_tod_spec_main(args):
    from optparse import OptionParser
    o = OptionParser()
    o.add_option('-i', '--interactive-debugger', action='store_true')
    o.add_option('--no-plots', action='store_true')
    o.add_option('--force', default=False, action='store_true')

    opts, args = o.parse_args(args)

    if opts.interactive_debugger:
        moby2.util.debugger.interactive(True)

    cfg_filename, args = args[0], args[1:]
    config = moby2.util.MobyDict.from_file(cfg_filename)
    
    tod_list = moby2.scripting.get_tod_list(config['get_spec']['tod_list'],
                                            cmdline_args=args)

    stats = {
        'failed': 0,
        'skipped': 0,
        'ok': 0,
    }

    for basename in tod_list:
        try:
            result = get_tod_spec_one(config, basename, no_plots=opts.no_plots, force=opts.force)
            result = {True: 'ok', False: 'skipped'}[result]
        except Exception as e:
            print(e)
            result = 'failed'
        stats[result] += 1

    print('\nSummary:')
    for (k, v) in stats.items():
        print('  %-10s:  %4i' % (k, v))

def get_tod_spec_one(config, basename, no_plots=None, force=None):

    bn = basename
    fb = moby2.scripting.products.get_filebase()
    fn = fb.filename_from_name(bn)
    print('tod_name:', bn)
    print('filename:', fn[0])

    outputd = moby2.scripting.OutputFiler(
        prefix=config['get_spec']['data_prefix'], data={'basename': bn})
    ofile = outputd.get_filename('.spec.npz')
    print('check ofile', ofile)

    if config['get_spec']['skip_existing'] and os.path.exists(ofile) \
       and not force:
        print('  ... output exists; skipping.')
        return False

    ttrans = TODTransformer(config)
    tod = ttrans.init_from_tod(fn[0])

    ttrans.preprocess_data(tod.data)

    spec = ttrans.get_spectrum(tod.data)

    # Save spectra.
    units = config['get_spec']['cal_units']
    spec.spec = spec.spec[ttrans.mask]
    spec.meta.update({'name': bn,
                      'det_uid': tod.det_uid[ttrans.mask],
                      'cal': ttrans.cal[ttrans.mask],
                      'cal_units': units})
    ttrans.log.trace(0, 'Writing to %s.' % ofile)
    spec.to_file(ofile)

    # Plot all
    ppref = config['get_spec'].get('plot_prefix')
    if no_plots:
        ppref = None
    if ppref is None:
        return True

    ttrans.log.trace(0, 'Plotting.')
    outputp = moby2.scripting.OutputFiler(
        prefix=config['get_spec']['plot_prefix'], data={'basename': bn})
    
    import pylab as pl

    # Some stuff is band-dependent.
    freqs = tod.info.array_data['nom_freq']
    freq_list = sorted(list(set(freqs)))
    colors = dict(list(zip(freq_list, 'kbgr')))
    legends = [(f, '%3i GHz' % f, pl.matplotlib.patches.Patch(color=colors[f]))
               for f in freq_list]

    pl.rcParams.update({'legend.fontsize': 10})
    transparent = max(.05, 3./ttrans.mask.sum())
    tstep = 10
    tvec = np.arange((tod.data.shape[1]+tstep-1)//tstep) / (ttrans.f_samp/tstep)
    to_leg = []
    for f in freq_list:
        s = (freqs[spec.meta['det_uid']] == f)
        if s.sum() == 0:
            continue
        e = tod.data[ttrans.mask,::tstep][s,:]
        pl.plot(tvec, e.transpose(), alpha=transparent, c=colors[f])
        to_leg.append(f)

    labs, arts = list(zip(*[x[1:] for x in legends if x[0] in to_leg]))
    pl.legend(arts, labs, loc='upper right')

    cm = np.dot(ttrans.mask, tod.data) / ttrans.mask.sum()
    pl.plot(tvec, cm[::tstep], lw=2, color='k')
    pl.xlabel('Time [s]')
    pl.ylabel('Signal [%s]' % (spec.meta['cal_units']))
    pl.title(bn)
    outputp.savefig('00_tod.png')

    # Plot it.
    to_leg = []
    for f in freq_list:
        s = (freqs[spec.meta['det_uid']] == f)
        if s.sum() == 0:
            continue
        pl.loglog(spec.f, spec.spec[s,:].transpose(),
                  alpha=transparent, c=colors[f])
        to_leg.append(f)

    labs, arts = list(zip(*[x[1:] for x in legends if x[0] in to_leg]))
    pl.legend(arts, labs, loc='upper right')
    pl.title(spec.meta['name'])
    pl.xlim(1e-2, spec.f.max())
    pl.xlabel('f [Hz]')
    pl.ylabel('Spectral density (%s/rtHz)' % spec.meta['cal_units'])

    outputp.savefig('10_spectra.png')

    # Show ~white noise by frequency band.
    white_band = config['get_spec'].get('white_noise_band', (10,20))

    return True
