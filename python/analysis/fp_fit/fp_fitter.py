#!/usr/bin/python
from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

"""
fpFitter.py : focal plane fitter

Tries to locate detectors in the focal plane, determines other stuff along the way.
"""

import moby2

from matplotlib import pyplot as plt
import numpy as np
import os, sys, time

# Local imports
from .fp_file import FPFitFile
from .fp_map import FPMap
from .fp_source import FitParams, BeamModel

DEG = np.pi/180
ARCMIN = DEG/60


# Fit validation

def check_convergence(p_in, p_out, tols, abs_flags):
    """
    All four input arguments are lists or tuples of the same length.
    Compares values in 'p_in' against those in 'p_out', and returns
    False if the difference is greater than the corresponding value in
    'tols'.  If the value in 'abs_flags' is True, the values are
    compared in absolutely rather than fractionally.
    """
    for i, (pi, po, t, a) in enumerate(zip(list(p_in), list(p_out), tols, abs_flags)):
        if a:
            if abs(pi - po) > t:
                #print 'busted %i (%f -> %f       t= %f'%(i, pi, po, t)
                return False
        else:
            d = abs(pi+po)/2
            if d < 1e-7:
                d = 1e-7
            if abs(pi - po) / d > t:
                #print 'busted %i (%f -> %f d= %f t= %f'%(i, pi, po, d, t)
                return False
    return True

def cut_check(p, cut_rules):
    """
    Returns False if beam params in 'p' fail the cut tests in 'rules'.
    """    
    keys = ['x0', 'y0', 'h', 'tau', 'w', 'sn']
    scalars = {'x0': np.pi/180, 'y0': np.pi/180}
    limk = ['cut_x0', 'cut_y0', 'cut_h', 'cut_tau', 'cut_w', 'cut_sn']
    for k, l in zip(keys, limk):
        if cut_rules.get(l) is None:
            continue
        lo, hi = cut_rules[l]
        scale = scalars.get(k, 1)
        if p[k] > hi*scale or p[k] < lo*scale or not np.isfinite(p[k]):
            message = 'rejecting %f based on %s cuts'%(p[k],k)
            return (False, message)
    return (True, 'Ok')
    

def apply_cuts(opts, fpf, tod_info=None, logger=None):
    fpf_out = FPFitFile(fpf.det_uid)
    if '_execcfg' in opts:
        tod_id = moby2.scripting.products.get_tod_id(tod_info=tod_info)
        ic = moby2.scripting.execcfg.InputChooser()
        opts.update(ic.get_config(opts['_execcfg'], tod_id=tod_id))

    do_reset = opts.get('reset')
    logger.trace(1, 'Resetting cuts, by request.')
    for i in range(len(fpf.det_uid)):
        b = fpf.get_row(i)
        if b['ok'] or do_reset:
            b['ok'], msg = cut_check(b, opts)
            if not b['ok'] and logger != None:
                logger.trace(2, 'det_uid=%i: %s' % (fpf.det_uid[i], msg))
        fpf_out.update_row(i, b)
    return fpf_out

def simpler_upper_bound(n):
    best_upper = 2*n
    for i0 in [1,3,5,7,9,11,15,21,22,25,27,33]:
        n2 = 1.*n/i0
        nu = 2**int(np.ceil(np.log2(n2)))
        best_upper = min(best_upper, nu*i0)
    return best_upper

def simpler_lower_bound(n):
    # This is for numpy... return largest number m <= n that
    # contains factors of 2 and 3 only.
    n_better = []
    n3 = 0
    while True:
        n2 = n // 3**n3
        if n2 <= 0:
            break
        n2 = np.floor(np.log(n2) / np.log(2))
        n_better.append(2**n2 * 3**n3)
        n3 += 1
    if len(n_better) == 0:
        return n
    return max(n_better)

def apply_friend_cuts(opts, fpf, logger=None):
    """
    Validate detectors based on nearby detectors being at some
    expected distance.  This is most useful after a coarse fit, which
    is noisy and has crazy outliers.
    """
    counts = np.zeros(fpf.det_uid.shape, int)
    r0, r1 = opts['ring_radii_arcmin']  # e.g. (1.7, 2.4)
    x, y= fpf.x0/DEG, fpf.y0/DEG
    for i in range(len(fpf.det_uid)):
        r = ((x - x[i])**2 + (y - y[i])**2)**.5 * 60
        s = (r0 <= r) * (r < r1) * fpf.ok
        counts[i] = s.sum()
    fpf_out = fpf.copy()
    fpf_out.ok *= counts >= opts['min_friend_count']
    return fpf_out

def apply_node_cuts(opts, fpf, logger=None):
    """Mask detectors if too many solutions are bunched in a small
    location.  This usually indicates some sort of readout glitch
    which has fooled the coarse fitter.
    """
    counts = np.zeros(fpf.det_uid.shape, int)
    r0 = opts['distance']
    x, y= fpf.x0/DEG, fpf.y0/DEG
    for i in range(len(fpf.det_uid)):
        r = ((x - x[i])**2 + (y - y[i])**2)**.5 * 60
        s = (r < r0) * fpf.ok
        counts[i] = s.sum()
    fpf_out = fpf.copy()
    fpf_out.ok *= (counts < opts['count'])
    return fpf_out
    

class FPScriptable(object):
    """
    Provide canned handling for an output manager and logger.  Objects
    can extend this to provide standardized handling of outputm and
    logger arguments.
    """
    def __init__(self, outputm=None, logger=None):
        if outputm == None:
            outputm = moby2.scripting.OutputFiler()
        self.outputm = outputm
        if logger == None:
            logger = moby2.util.log.get_logger()
        self.logger = logger
        
    def trace(self, level, msg):
        self.logger.trace(level, msg)


class CoarsePeakExtractor(FPScriptable):
    """
    Gets detector position estimate based on position of peak signal
    in time ordered data.
    """
    def __init__(self, opts={}, outputm=None, logger=None):
        FPScriptable.__init__(self, outputm=outputm, logger=logger)
        self.opts = opts

    def fit_tod(self, fptod, opts={}):
        # Band pass filter?
        _opts = self.opts.copy()
        _opts.update(opts)
        opts = _opts

        plot_prefix = opts.get('plot_prefix')

        tod = fptod.tod
        filt_opts = opts.get('filter_band_pass')
        if filt_opts != None:
            mean, sigma = filt_opts
            filt = moby2.tod.filters2.gaussianBand(fc=mean, df=sigma)
            sample_time = (tod.ctime[-1] - tod.ctime[0]) /\
                (tod.nsamps-1)
            # We're using fft3w so you need to pad the FFT for it to not be absurdly slow.
            n1 = simpler_upper_bound(tod.data.shape[1])
            g1 = np.zeros(n1)
            #g1[:len(g)] = g - g.mean()
            fw = filt.getExecFor(dt=sample_time, n=len(g1))
            #fw = filt.getExecFor(dt=sample_time, n=tod.data.shape[1])

        x, y = fptod.source_coords
        results = FPFitResultList()
        for i in range(len(tod.data)):
            uid = tod.det_uid[i]
            self.trace(1, 'Coarse fitting det_uid=%i' % uid)
            g = tod.data[i].copy()
            if filt_opts != None:
                self.trace(2, 'Filtering det_uid=%i (n=%i for %i)' % (uid, n1, len(g)))
                g1[:len(g)] = g
                g1[len(g):] = g.mean()
                g[:] = fw.applyTime(g1)[:len(g)]

            # Apodize.
            n_apod = 800
            g[:n_apod]  *= np.arange(n_apod) / float(n_apod)
            g[-n_apod:] *= (1 - np.arange(n_apod) / float(n_apod))
            peak = g.argmax()
            h = (g.max() - g.min()) / g.std()
            result = {'index': peak, 'dx': x[peak], 'dy': y[peak],
                      'h': h}
            #          'h': (g.max() - g.min())/2}
            # Reject peaks occurring too close to beginning or end of time stream.
            if (peak < len(g)*.05) or (peak > len(g)*.95):
                results.append(FPFitResult(ok=False, message='Coarse fit failed.',
                                           det_uid=uid, data=(result, {}, {})))
            else:
                results.append(
                    FPFitResult(ok=True, message='Coarse fit ok.',
                                det_uid=uid, data=(result, {}, {})))
        
            # Plot that?
            if plot_prefix != None:
                self.trace(2, 'Plotting det_uid=%i' % uid)
                self.plot(result, g, self.outputm.get_filename(
                        plot_prefix + '{det_uid:04d}.png',
                        data={'det_uid': uid}))
            self.trace(2, 'Finished processing for det_uid=%i' % uid)

        return results

    def go(self):
        if not np.any(self.mask):
            return False
        self.peak_index = self.z[self.mask].argmax()
        return True

    def plot(self, result, z, filename=None, title=''):
        peak_index = result['index']
        lo, hi = max(0, peak_index-100), min(len(z), peak_index+100)
        plt.subplot(211)
        plt.title(title)
        plt.plot(z)
        plt.axvline(peak_index, alpha=.2, color='red', lw=5)
        plt.subplot(212)
        plt.plot(z[lo:hi])
        x0, y0 = [result[k]/DEG for k in ['dx', 'dy']]
        plt.title('Peak near X,Y = {:.3f},{:.3f} deg'.format(x0, y0))
        if filename is not None:
            plt.savefig(filename)
            plt.clf()
    
    def get_params(self):
        return {'x0': self.x[self.mask][self.peak_index],
                'y0': self.y[self.mask][self.peak_index],
                'h':  self.z[self.mask][self.peak_index],}
        

class FPModelFitter(FPScriptable):
    def __init__(self, opts={}, outputm=None, logger=None):
        FPScriptable.__init__(self, outputm=outputm, logger=logger)
        self.opts = opts.copy()
        self.model = BeamModel.get_model(opts)
        self.n_calls = 0

    def _get_updated_opts(self, more_opts={}):
        opts = self.opts.copy()
        opts.update(more_opts)
        return opts

    def __call__(self, p, filt_opts, x, y, z=None):
        self.n_calls += 1
        if isinstance(p, FitParams):
            p = p.tuple()
            tau = filt_opts.get('tau', 0.)
        else:
            p, tau = p[:-1], p[-1]

        z0 = self.model(x, y, p)
        # Filter it?
        if 'filt' in filt_opts:
            f = filt_opts['time_const_freq']
            filt = filt_opts['filt'] / (1. + 2j*np.pi*abs(tau)*f)
            z0f = np.fft.fft(z0)
            z0 = np.fft.ifft(z0f * filt).real

        if z is None:
            return z0
        return z - z0

    def fit_model_tod(self, fptod, init_fits, opts={}, idx=None):
        """
        Call self.fit_model on each detector in fptod (an FPTOD
        object) that has a valid fit in init_fits (an FPFitFile
        object).  You can just do a subset of the dets by passing in
        the indices into fptod.tod.data in parameter idx.
        """
        opts = self._get_updated_opts(opts)

        if idx is None:
            idx = range(len(fptod.tod.data))
        
        fit_ok, (dx, dy) = init_fits.get_property(['x0', 'y0'],
                                                  det_uid=fptod.tod.det_uid)
        plot_data={'basename': fptod.tod.info.name}
        readout_filter = opts.get('readout_filter')
        if readout_filter == True:
            readout_filter = fptod.readout_filter
        sample_time = opts.get('sample_time')
        if sample_time is None:
            sample_time = fptod.sample_time
        results = FPFitResultList()
        for i in idx:
            uid = fptod.tod.det_uid[i]
            if not fit_ok[i]:
                results.append(FPFitResult(
                        ok=False, det_uid=uid,
                        message='Rejected for lack of input fit.'))
                continue
            init_fit = {'dx': dx[i], 'dy': dy[i]}
            s1, z1 = fptod.extract_data(i, dx[i], dy[i])
            x1, y1 = [_x[s1] for _x in fptod.source_coords]
            plot_data['det_uid'] = uid
            result = self.fit_model(
                x1, y1, z1, opts=opts, init_fit=init_fit,
                readout_filter=readout_filter, plot_data=plot_data,
                sample_time=sample_time,
                msg_prefix='det_uid=%4i: ' % uid)
            result.det_uid = uid
            results.append(result)
        return results

    def fit_model(self, x, y, data, mask=None, opts={}, init_fit=None,
                  readout_filter=None, sample_time=None, plot_data={},
                  fft_trim=True, msg_prefix=''):
        """
        Fit a beam model to timestream data.  The coordinates of the
        source should be passed in as 'x' and 'y', and the detector
        timestream in 'data'.  These are all 1-d arrays of the same
        length.  To fit some subset of these data, pass in a boolean
        array 'mask' of the same length.
        """
        opts = self._get_updated_opts(opts)

        if fft_trim:
            # Trim vectors at the ends so that fft factorizes.  This
            # is a huge speed-up...
            if mask is None:
                n_fft = len(x)
            else:
                n_fft = mask.sum()
            n_better = simpler_lower_bound(n_fft)
            dn = int(n_fft - n_better)
            if dn > 0:
                self.trace(2, msg_prefix + 
                           'Shortening time stream from %i to %i for FFT' % (
                               n_fft, n_better))
                if mask is None:
                    mask = np.ones(n_fft, bool)
                i = mask.nonzero()[0]
                n_left = dn//2
                n_right = (dn+1)//2
                mask[i[:n_left]] = False
                mask[i[-n_right:]] = False
                assert (mask.sum() == n_better)

        if not mask is None:
            x, y, data = x[mask], y[mask], data[mask]
        
        if len(x) < 100:
            message = 'starting point leaves only %i samples under mask' % len(x)
            self.trace(1, msg_prefix + message)
            return FPFitResult(ok=False, message=message)

        # Prepare time constant fitting and other Fourier stuff
        fit_tau = opts.get('fit_time_const')
        filt_opts = {'time_const': True,
                     'time_const_freq': np.fft.fftfreq(len(x)) / sample_time,
                     'filt': np.ones(len(x))}
        if readout_filter is not None:
            filt_opts['filt'] = filt_opts['filt'] * \
                readout_filter.transfer(filt_opts['time_const_freq']) / \
                readout_filter.transfer(0.)

        # Get a simple tuple from the input beam parameters
        beam = self.model.init_params(init_fit, opts=opts)
        p = beam.tuple()
        # Append time constant according to opts.
        tau_source = opts.get('init_time_const', 1e-6)
        if tau_source == 'input_fit':
            tau = init_fit['tau']
        else:
            tau = float(tau_source)
        # That's how the log drivers learn not to nan.
        tau = min(1e-9, tau)
        p = p + (tau,)

        # Get convergence tolerances
        tols = [opts['tolerance'][k] for k in 
                ['x0_atol', 'y0_atol', 'h_ftol', 'sig_ftol', 'base_ftol', 'tau_ftol']]
        abs_tol = [1, 1, 0, 0, 0, 0]
             
        # Define groups of parameters to fit simultaneously
        flag_sets = self.model.get_fit_flag_sets()
        flag_sets = [list(f) + [fit_tau and (i!=0)] for i,f in enumerate(flag_sets)]

        # Fit
        converged = False
        max_iters = 10 * len(flag_sets)
        max_time = 4.
        t_start = time.time()
        iter_i = 0
        while not converged and iter_i < max_iters and time.time() - t_start < max_time:
            p_0 = p
            f = flag_sets[iter_i % len(flag_sets)]
            # leastsq works faster than fmin for whatever reason.
            self.n_calls = 0
            fit_out = moby2.util.fitting.multi_leastsq(
                self, p, (filt_opts, x, y, data),
                free=f, full_output=1)
            p = fit_out[0]        # xopt
            iter_i += 1
            if iter_i < len(flag_sets):
                # Don't assess convergence until all flags have been tested.
                continue
            ok1 = (fit_out[4] in [1,2,3,4])  # ier
            ok2 = check_convergence(p_0, p, tols, abs_tol)
            converged = ok1 and ok2
            if converged:
                # Reject way-off fits.
                dist = ((p[0] - init_fit['dx'])**2 + (p[1] - init_fit['dy'])**2)**.5
                if dist > 10 * ARCMIN:
                    converged = False
                    self.trace(2,msg_prefix + 'Rejecting solution due to '
                               '>10 arcmin position change.')
                    break

        if converged:
            self.trace(1,msg_prefix + 'Converged after %i iters (t=%.1f)' % 
                       (iter_i, time.time() - t_start))
        else:
            self.trace(1,msg_prefix + 'No convergence after %i tries! (t=%.1f)' %
                       (iter_i, time.time() - t_start))
        
        # Return the usual fit stuff
        elapsed = time.time() - t_start

        beam.store(p[:-1])
        self.trace(2, msg_prefix + 'result: %s' % beam.pretty())

        if 'filt' in filt_opts:
            # Store absolute value of time constant.
            filt_opts['tau'] = abs(p[-1])

        if not opts.get('plot_prefix') is None:
            ofile = self.outputm.get_filename(
                    (opts['plot_prefix'] + '{det_uid:04d}.png'),
                    data=plot_data)
            self.trace(3, msg_prefix + 'plotting to %s' % ofile)
            self.plot_model(beam, filt_opts, x, y, data,
                            sample_time=sample_time,
                            filename=ofile,
                            title_text='Model fit {basename} - uid={det_uid:04d}'.\
                                format(**plot_data))

        if converged:
            message = 'Fit successful.'
        else:
            message = 'Failed to converge.'
        # You need to store some S/N results too.
        zm = self(beam, filt_opts, x, y, None)
        weight = (zm - zm.min()) / (zm.max() - zm.min())
        resid = (zm - data).std()
        wstd, rstd = (data*weight).std(), ((zm-data)*weight).std()
        if (wstd==0):
            merit = 0.
        else:
            merit = 1./(1.+(rstd/wstd)**2)  #Map into [0,1]
        quality = {
            'merit': merit,
            'rms_resid': resid,
            'S/N': (zm.max() - zm.min()) / resid
            }
        return FPFitResult(ok=converged, message=message, data=(beam, filt_opts, quality))

    def plot_model(self, beam, filt_opts, x, y, z=None, sample_time=None, 
                   filename=None, title_text=None):
        z_model = self(beam, filt_opts, x, y, None)
        fig = plt.gcf()

        tscale = sample_time
        tunit = 'seconds'
        ascale = 180.*60/np.pi
        unit = 'arcmin'

        x = ascale*(x - beam['dx'])
        t = tscale*np.arange(len(x))
        gg, bb = z, z_model

        arcmin_lims = self.opts.get('arcmin_lims', (2., 1.))
        XLIM = (-arcmin_lims[0], arcmin_lims[0])
        XMUL = arcmin_lims[1]

        plt.subplot(221)
        plt.plot(t, gg, 'k')
        plt.plot(t, bb, 'b')
        plt.subplot(223)
        plt.plot(t, gg-bb, 'r')
        plt.xlabel('Time (%s)'%tunit)
        plt.subplot(222)
        plt.plot(x, gg, 'k')
        plt.plot(x, bb, 'b')
        plt.xlim(*XLIM)
        plt.gca().xaxis.set_major_locator(plt.MultipleLocator(XMUL))
        plt.subplot(224)
        plt.plot(x, gg-bb, 'r')
        plt.xlabel('X (%s)'%unit)
        plt.xlim(*XLIM)
        plt.gca().xaxis.set_major_locator(plt.MultipleLocator(XMUL))

        if title_text != None:
            plt.figtext(0.5, 0.94, title_text, ha='center', va='bottom')

        if not filename is None:
            odir = os.path.split(filename)[0]
            if odir != '' and not os.path.exists(odir):
                os.makedirs(odir)
            plt.savefig(filename)
            plt.clf()

        return fig

class FPTOD(FPScriptable):
    # Manage a TOD for the purposes of focal plane fitting.

    def __init__(self, tod=None, outputm=None, logger=None):
        FPScriptable.__init__(self, outputm=outputm, logger=logger)
        self.set_tod(tod)

    def set_tod(self, tod):
        self.tod = tod
        self.sample_time = (tod.ctime[-1] - tod.ctime[0]) / (tod.nsamps-1)
        if hasattr(tod.info, 'mce_filter'):
            self.readout_filter = tod.info.mce_filter
        else:
            self.readout_filter = None

    def preprocess(self, opts={}):
        """
        Condition the TOD according to parameters in opts dictionary.
        Parameters recognized:

        scale: (float) Rescale data by factor.
        
        glitch_cuts: (dict) Cut and fill the data using standard moby2
          glitch cuts with these parameters.

        detrend: (bool) Remove a trend for each detector.

        And much, much more.
        """
        tod = self.tod

        if opts.get('remove_readout_filter'):
            tod.data /= tod.info.runfile.ReadoutFilter().gain()

        if opts.get('scale'):
            tod.data *= opts['scale']

        if opts.get('glitch_cuts'):
            cuts = moby2.scripting.products.get_cuts(opts['glitch_cuts'],
                                                     tod=tod)
            tod.cuts = cuts.copy(det_uid=tod.det_uid)
            moby2.tod.fill_cuts(tod)

        if hasattr(tod, 'data_mask'):
            i = (~tod.data_mask).nonzero()[0]
            self.trace(2, 'Filling masked samples (%i).' % len(i))
            for d in tod.data:
                d[i] = d[i-1]

        self.trace(3, 'Removing mean')
        moby2.tod.ops.remove_mean(tod=tod)
        if opts.get('detrend'):
            self.trace(3, 'Detrending')
            moby2.tod.ops.detrend_tod(tod=tod)

        if opts.get('common_mode') is not None:
            self.trace(3, 'Removing common mode')
            cm_src = opts['common_mode']['source']
            if cm_src[0] == 'file':
                self.trace(3, 'Common mode source is file %s' % cm_src[1])
                cm = np.loadtxt(cm_src[1])
                start, n = tod.info.sample_index, tod.nsamps
                cm = cm[start:n+start]
                cm -= cm.mean()
            elif cm_src[0] == 'compute':
                # You would be crazy to do this without cuts.
                cuts0 = moby2.scripting.products.get_cuts(opts['glitch_cuts'],
                                                          tod=tod)
                cuts1 = moby2.scripting.products.get_cuts(opts['position_cuts'],
                                                          tod=tod)
                cuts0.merge_tod_cuts(cuts1)
                cuts = cuts0.copy(det_uid=tod.det_uid)
                data_copy = data.copy()
                moby2.tod.fill_cuts(data=data_copy,cuts=cuts)
                cm = np.zeros(data.shape[1])
                for i in cuts.get_uncut():
                    cm[:] += data_copy[i] - data_copy[i].mean()
                cm -= cm.mean()
                del data_copy
                
            else:
                raise ValueError("Unknown common_mode source specifier '%s'"\
                    % cm_src[0])
            if opts['common_mode'].get('fit', True):
                amps = np.zeros(len(tod.det_uid))
                for i in range(len(tod.det_uid)):
                    amps[i] = np.dot(cm, tod.data[i]) / np.dot(cm,cm)
            else:
                amps = np.ones(len(tod.det_uid))
            for i in range(len(tod.det_uid)):
                tod.data[i] -= amps[i]*cm
            if opts['common_mode'].get('plot'):
                plt.plot(cm)
                ofile = self.outputm.get_filename(opts['common_mode']['plot'])
                plt.savefig(ofile)
                plt.clf()

        if opts.get('poly_order') is not None:
            self.trace(3, 'Removing polys')
            order = opts.get('poly_order', 10)

            # Orthogonal polynomial basis
            vects = np.empty((order+1, tod.nsamps))
            all_idx = np.arange(tod.nsamps).astype('float32') / tod.nsamps
            for i in range(vects.shape[0]):
                vects[i] = all_idx**i
                for j in range(0, i):
                    vects[i] -= np.dot(vects[i], vects[j]) * vects[j]
                vects[i] /= (vects[i]**2).sum()**.5
            # Fit and remove.
            decimation = opts.get('poly_decimation', 10)
            for i in range(tod.data.shape[0]):
                amps = np.dot(vects[:,::decimation], tod.data[i,::decimation])
                tod.data[i] -= np.dot(amps*decimation, vects)

        if opts.get('fix_sign'):
            if opts.get('guess_sign'):
                raise RuntimeError("You cannot pass fix_sign and guess_sign; pick one.")
            self.trace(3, 'Sign-correcting using ArrayData')
            assert('optical_sign' in tod.info.array_data)
            recal = tod.info.array_data.get_property('optical_sign',
                                                     det_uid=tod.det_uid)
            moby2.tod.apply_calibration(tod.data, np.arange(len(recal)), recal)
            
        if opts.get('guess_sign'):
            self.trace(3, 'Sign-correcting using timestream median.')
            for i in range(tod.data.shape[0]):
                mx, mn = tod.data[i].max(), tod.data[i].min()
                md = np.median(tod.data[i,::100])
                if (md - mn) > (mx - md):
                    tod.data[i] *= -1

        plot_format = opts.get('plot_prefix')
        if plot_format is not None:
            self.trace(3, 'Plotting conditioned TODs')
            decimation = len(tod.data[0]) // 1000
            n_plots = 0
            for i in range(len(tod.det_uid)):
                info = tod.info.get_dict(det_uid=tod.det_uid[i])
                plt.title('{name}  r{row:02d}c{col:02d} ({det_uid:d})'.\
                             format(**info))
                y = tod.data[i]
                y = y[:len(y)//decimation * decimation].reshape(-1, decimation)
                ylo, yhi = y.min(axis=1), y.max(axis=1)
                x = np.arange(len(y))*decimation
                plt.fill_between(x, ylo, yhi, color='blue')
                self.outputm.savefig(plot_format+'{det_uid:04d}.png', data=info)
                if opts.get('max_plots', 0) > 0:
                    n_plots += 1
                    if n_plots >= opts['max_plots']:
                        break

    def compute_source_coords(self, source_name=None):
        """
        Return focal plane coordinates of the target.
        """
        # It's rather important that tod.fplane be trivial.
        if source_name is None:
            fplane = moby2.pointing.FocalPlane([0.],[0.])
            if hasattr(self.tod, 'fplane'):
                self.trace(1,'Replacing tod.fplane.')
            self.tod.fplane = fplane
            sources = moby2.ephem.get_sources_in_patch(tod=self.tod)
            if len(sources) == 1:
                source_name = sources[0][0]
            else:
                raise RuntimeError("Multiple sources identified, pick one: %s" % \
                    str(sources))

        t, az, alt = self.tod.ctime, self.tod.az, self.tod.alt
        ra0, dec0 = moby2.ephem.get_source_coords(source_name, t.mean())
        theta, phi = moby2.pointing.get_coords_inverse(
            t, az, alt, t*0+ra0, t*0+dec0)
        self.source_coords = theta, -phi
        return source_name, self.source_coords

    def get_offset_coords(self, xi, eta):
        """
        Using beam parameters in 'beam', convert from planet
        coordinates in the focal plane to planet coordinates relative
        to the beam model.
        """
        # FIXME: this should be a rotation, not a flat-sky translation.
        dx = self.source_coords[0] - xi
        dy = self.source_coords[1] - eta
        return dx, dy
        
    def get_proximity_mask(self, dx, dy, opts={}):
        mask_radius_arcmin = opts.get('mask_radius', 0.)
        x, y = self.get_offset_coords(dx, dy)
        if mask_radius_arcmin > 0.:
            rM = (x**2 + y**2)**.5 < mask_radius_arcmin *np.pi/180/60
        else:
            rM = np.ones((x+y).shape, bool)
        # Convert to a CutsVector
        segments = moby2.tod.CutsVector.from_mask(rM)
        # Buffer it, and return the mask.  This used to be hard-coded for 400 Hz samples.
        lag = opts.get('lag_margin_samples', opts.get('lag_margin', 0)/self.sample_time)
        return segments.get_buffered(int(lag)).get_mask()

    def get_scan_mask(self, opts={}):
        # Left and right-going masks
        az = self.tod.az
        dirs = opts.get('scan_direction', 'both')
        if dirs == 'both':
            return np.ones(az.shape, 'bool')
        if dirs == 'left':
            return np.gradient(az) < 0
        if dirs == 'right':
            return np.gradient(az) > 0
        raise ValueError("scan_direction=%s is invalid" % dirs)

    def extract_data(self, i, dx, dy, opts={}):
        """
        Get subset of time-stream data, for detector i, for when
        source is near position (-dx, -dy).  So if (dx, dy) is a good
        guess for detector i's offset, you'll get the signal you want.

        Returns (mask, signal) where signal = tod.data[i,mask].
        Import keys in opts:

          'mask_radius': radius of cut-out region, in arcminutes.

          'destripe': boolean indicating whether to remove baseline
            levels from each segment.
        """
        do_remove = opts.get('destripe', True)
        r0 = opts.get('mask_radius', 10.)
        r1 = r0*1.5

        x, y = self.get_offset_coords(dx, dy)
        rM = (x**2 + y**2)**.5
        s0 = (rM <= r0*ARCMIN)
        if not do_remove:
            return s0, self.tod.data[i][s0]

        s1 = (rM <= r1*ARCMIN)
        edges = np.hstack([False, s1, False])
        edges = (edges[1:] != edges[:-1]).nonzero()[0].reshape(-1,2)
        output = np.empty(s1.sum())
        offset = 0
        source = self.tod.data[i]
        for e0,e1 in edges:
            x_ = x[e0:e1]
            s = ~s0[e0:e1]
            if s.sum() < 4:
                s0[e0:e1] = False
                continue
            order = 1
            if not s[0] or not s[-1]:
                order = 0
            p = np.polyfit(x_[s], source[e0:e1][s], order)
            output[offset:offset+e1-e0] = np.polyval(p, x_)
            offset += e1-e0
        return s0, source[s0] - output[s0[s1]]

     
class FPFitResult:
    """
    Container for a single detector's fit results, as returned by
    FPModelFitter.fit_model, and similar fitters.
    """
    def __init__(self, ok=None, det_uid=None, message=None, data=None):
        self.ok = ok
        self.det_uid = det_uid
        self.message = message
        self.data = data
   
class FPFitResultList(list):
    """
    Container for multiple FPFitResult objects, with methods to
    combine results for analysis or output.
    """
    def get_fp_file(self):
        if len(self) == 0:
            return FPFitFile(det_uid=[])
        fits_sorted = sorted([(x.det_uid, x) for x in self
                              if x.data != None ])
        if len(fits_sorted) == 0:
            return FPFitFile(det_uid=[])
        rows = [(x.det_uid, x.ok,
                 x.data[0]['dx'], x.data[0]['dy'],
                 x.data[0]['h'], x.data[2].get('w',0),
                 x.data[1].get('tau',0), x.data[2].get('S/N',0))
                for _,x in fits_sorted]
        uid, ok, x, y, h, w, tau, sn = list(map(np.array, list(zip(*rows))))
        fpf = FPFitFile(det_uid=uid)
        fpf.ok[:] = ok
        fpf.x0[:] = x
        fpf.y0[:] = y
        fpf.tau[:] = tau
        fpf.h[:] = h
        fpf.w[:] = w
        fpf.sn[:] = sn
        return fpf
    
    def get_struct_db(self):
        # Combine all data usefully.
        fields = [('det_uid', int, '%4i'),
                  ('ok', np.int16, '%1i'),
                  ('dx', float, '%11.4e'),
                  ('dy', float, '%11.4e'),
                  ('h', float, '%10.3e'),
                  ('w', float, '%10.3e'),
                  ('tau', float, '%10.3e'),
                  ('base', float, '%12.5e'),
                  ('S/N', float, '%8.2f'),
                  ('merit', float, '%8.3f'),
                  ('rms_resid', float, '%10.3e'),
                   ]
        n = len(self)
        data = {}
        for k, t, _ in fields:
            data[k] = np.zeros(n, t)
        data['det_uid'] = [r.det_uid for r in self]
        for i,r in enumerate(self):
            if r.data == None:
                continue
            data['ok'][i] = int(r.ok)
            for k in ['dx', 'dy', 'h', 'w', 'base']:
                data[k][i] = r.data[0][k]
            for k in ['S/N', 'merit', 'rms_resid']:
                data[k][i] = r.data[2][k]
            data['tau'][i] = r.data[1]['tau']
        formats = dict([(t,f) for t,_,f in fields if f != None])
        return moby2.util.StructDB.from_data([(k,data[k]) for k,_,_ in fields],
                                             formats=formats)
