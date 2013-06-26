#!/usr/bin/python
from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

import sys
import numpy as np

import moby2

trace = moby2.util.log.logger.trace

# transitional...
_fp_formats = {
    'det_uid': '%4d',
    'ok': '%1d',
    'x0': '%9.6f',
    'x0_err': '%9.6f',
    'y0': '%9.6f',
    'y0_err': '%9.6f',
    'tau': '%8.5f',
    'tau_err': '%8.5f',
    'h': '%.4e',
    'w': '%9.6f',
    'sn': '%9.1f',
    'base': '%.5e',
    'n_obs': '%3d',
    }
_fp_fields = ['ok', 'x0', 'x0_err', 'y0', 'y0_err', 'tau', 'tau_err',
              'h', 'w', 'sn', 'base', 'n_obs']
_fp_columns_format_str = ' '.join(['{%s:%s}'%(k, _fp_formats[k][1:])
                                   for k in _fp_fields]) + '\n'

class FPFitFile(moby2.detectors._SimpleDetData):
    fields = _fp_fields
    dtypes = {'ok': bool, 'n_obs': int}
    columns_format_str = _fp_columns_format_str
    xcfs = '{det_uid:4d} {ok:1d} '\
        '{x0:9.6f} {x0_err:9.6f} {y0:9.6f} {y0_err:9.6f} '\
        '{tau:8.5f} {tau_err:8.5f} '\
        '{h:.4e} {w:9.6f} {sn:9.1f} {n_obs:3d}\n'
    header = '# det_uid ok x0   x0_err    y0        y0_err    '\
        'tau      tau_err h           w        sn         n_obs'

    def __init__(self, det_uid=None):
        if det_uid is not None:
            self.det_uid = np.array(det_uid, dtype='int64')
            n = len(det_uid)
            for f in self.fields:
                setattr(self, f, np.zeros(n, self.dtypes.get(f, 'float64')))

    def __repr__(self):
        name = repr(self.__class__)
        return '%s with %i det_uid for fields ' % (name, len(self.det_uid)) + \
            ','.join(self.fields)

    def update_row(self, row, data):
        for k in self.fields:
            if k in data:
                getattr(self, k)[row] = data[k]

    @classmethod
    def from_columns_file(cls, filename):
        data = np.loadtxt(filename, unpack=1)
        det_uid = data[0].astype('int')
        self = cls(det_uid)
        self.ok = data[1].astype('int').astype('bool')
        if len(data[2:]) == 11:
            self.x0, self.x0_err, self.y0, self.y0_err, self.tau, self.tau_err, self.h, self.w, self.sn, self.base, self.n_obs = data[2:]
        elif len(data[2:-1]) == 9:
            self.x0, self.x0_err, self.y0, self.y0_err, self.tau, self.tau_err, self.h, self.w, self.sn = data[2:-1]
            self.base = 0 * self.w
        elif len(data[2:-1]) == 8:
            self.x0, self.x0_err, self.y0, self.y0_err, self.tau, self.tau_err, self.h, self.sn = data[2:-1]
            self.w = 0 * self.x0
            self.base = 0 * self.x0
        elif len(data[2:-1]) == 4:
            self.x0, self.x0_err, self.y0, self.y0_err = data[2:-1]
            self.base = 0
        else:
            raise ValueError("Strange number of columns in %s" % filename)
        self.n_obs = data[-1].astype('int')
        return self

    @classmethod
    def from_file(cls, filename):
        if filename.endswith('fits') or filename.endswith('fits.gz'):
            return cls.from_fits_table(filename)
        return cls.from_columns_file(filename)

    # This supercedes _SimpleDetData.write
    def write(self, filename, format=None):
        if format is None:
            if filename.endswith('fits') or filename.endswith('fits.gz'):
                format = 'fits'
            else:
                format = 'txt'
        data = [('det_uid', self.det_uid)]
        for k in self.fields:
            v = getattr(self, k)
            if v.dtype == bool:
                v = v.astype('int8')
            data.append((k, v))
        odb = moby2.util.StructDB.from_data(data,formats=_fp_formats)
        if format == 'fits':
            odb.to_fits_table(filename)
        elif format == 'txt':
            odb.to_column_file(filename)
        else:
            raise ValueError("Unknown format request, %s." % format)

    def write_reduced(self, filename, scale_amp=1.):
        format = 'txt'
        if filename.endswith('.fits') or filename.endswith('.fits.gz'):
            format = 'fits'
        s = self.ok.astype(bool)
        # det_uid peak_DAC SN    tau
        data = [('det_uid', self.det_uid[s]),
                ('peak_dac', self.h[s] * scale_amp),
                ('time_const', self.tau[s]),
                ('sn', self.sn[s]),
                ]
        odb = moby2.util.StructDB.from_data(
            data, formats={'peak_dac': '%12.3f',
                           'time_const': '%12.5f',
                           'sn': '%12.3f'})
        if format == 'txt':
            odb.to_column_file(filename)
        elif format == 'fits':
            odb.to_fits_table(filename)

    @classmethod
    def from_focal_plane(cls, fp):
        """
        Initialize from a FocalPlane object.
        """
        self = cls(fp.det_uid)
        self.x0 = fp.x.copy()
        self.y0 = fp.y.copy()
        self.ok = fp.mask.copy()
        zeros = np.zeros(self.ok.shape)
        self.tau, self.h, self.w = zeros.copy(), zeros.copy(), zeros.copy()
        self.base = zeros
        return self

    @classmethod
    def combine_fits(cls, fits, template=None, params={}):
        """
        Combine fits by shifting each one to match a template, and
        averaging the good fits for each detector.

        If a template is not provided, match to the first one.
        """
        trace(1, 'Fitting and averaging %i fits' % len(fits))
        if template is None:
            template = fits[0]
        # Start by shifting each fit to match the template.
        orig_fits, fits = fits, []
        fitter = FPTemplateFitter()
        fitter.set_template(template)
        fit_params = {'shift': True,
                      'rotation': False}
        fit_params.update(params)
        fit_results = [None for fi in range(len(orig_fits))]
        for fi,f0 in enumerate(orig_fits):
            if f0.ok.sum() < params.get('min_dets', 50):
                trace(2, 'Discarding fit with only %i good fits' % f0.ok.sum())
                continue
            ok, result = fitter.fit(f0, fit_params)
            if not ok:
                trace(2, 'Discarding fit due to failed template match')
                continue
            f1 = f0.copy()
            f1.x0 += result[0]
            f1.y0 += result[1]
            fits.append(f1)
            fit_results[fi] = result

        trace(1, 'Cut %i of %i fits (increase verbosity to see why).' % \
                  (len(orig_fits) - len(fits), len(orig_fits)))
        if len(fits) == 0:
            return None, None

        print([len(f.det_uid) for f in fits])
        n_det_uid = max([f.det_uid.max() for f in fits]) + 1
        output = cls(np.arange(n_det_uid))
        output.ok[:] = False

        ARCMIN = np.pi/180/60
        trace(1, 'Combining data for %i detectors' % n_det_uid)
        for uid in output.det_uid:
            ok = np.array([f.get_property('ok', det_uid=uid)[1]
                           for f in fits])
            x, y, tau = np.transpose([f.get_property(['x0','y0','tau'], det_uid=uid)[1]
                                      for f in fits])
            for _x in [x, y, tau]:
                # Yes, this happens...
                ok *= ~np.isnan(_x) * ~np.isinf(_x)
            x, y, tau = [_x[ok] for _x in [x,y,tau]]
            if ok.sum() < params.get('min_obs', 1):
                trace(2, 'Discarding det_uid=%i due to only %i contributors'
                      % (uid, ok.sum()))
                continue
            # Majority rules.
            x0, y0 = np.median(x), np.median(y)
            for iteration in [0,1,2]:
                d0 = ((x - x0)**2 + (y-y0)**2)**.5
                s0 = d0 < params.get('max_separation', 1)*ARCMIN
                if s0.sum() == 0:
                    break
                x0, y0 = x[s0].mean(), y[s0].mean()
            if s0.sum() <= 0:
                trace(2, 'Discarding det_uid=%i due to only %i items in '\
                          ' combination' % (uid, s0.sum()))
                continue
            vals = {
                'x0': x0, 'y0': y0,
                'x0_err': x[s0].std(),
                'y0_err': y[s0].std(),
                'tau': tau[s0].mean(),
                'tau_err': tau[s0].std(),
                'n_obs': s0.sum(),
                'ok': s0.sum() >= params.get('min_obs', 1) }
            output.update_row(uid, vals)
            trace(2, 'Result for det_uid=%i' % uid)
            for k in ['x0', 'y0', 'tau']:
                trace(2, ' %s = %10.5f +- %10.5f' %  (k, vals[k], vals[k+'_err']))
        return output, fit_results
        
    def plot_positions(self, filename, auto_zoom=True, params={},
                       title='', fig=None):
        import pylab as pl
        if fig is None:
            pl.figure()
            pl.gcf().set_size_inches(6., 6.)
        else:
            pl.figure(fig.number)

        s = self.ok
        if s.sum() == 0:
            pl.title(title + ' - no good fits')
            pl.savefig(filename)
            pl.clf()
        units = params.get('units', 'deg')
        scale = {'rad': 1., 'deg': 180/np.pi, 'arcmin': 60*180/np.pi}[units]
        x, y = self.x0[s]*scale, self.y0[s]*scale
        x0, y0 = np.median(x), np.median(y)
        r = ((x-x0)**2 + (y-y0)**2)**.5
        window = np.median(r)*3
        inside = r < params.get('zoom', scale*window)
        pl.scatter(x, y, alpha=0.5)
        if params.get('limits') is None:
            if np.any(inside):
                for vect,limiter in [(x,pl.xlim), (y,pl.ylim)]:
                    lo, hi = limiter()
                    lo = min(lo, vect[inside].min())
                    hi = max(hi, vect[inside].max())
                limiter(lo, hi)
        else:
            xlims, ylims = params['limits']
            pl.xlim(*xlims), pl.ylim(*ylims)
        pl.title(title + ' - %i dets outside window' % (~inside).sum())
        pl.xlabel('X (%s)' % units)
        pl.ylabel('Y (%s)' % units)

        def smart_locate(ax, n_max, bases=[1,2,5]):
            x0, x1 = ax.get_view_interval()
            if x1 == x0:
                return
            delta = (x1-x0) / (n_max-1)
            # Find smallest base and p such delta < base*10^p
            log_spacing = min([
                    np.ceil(np.log10(delta) - np.log10(b)) + np.log10(b)
                    for b in bases])
            loc = pl.MultipleLocator(10**log_spacing)
            ax.set_major_locator(loc)

        smart_locate(pl.gca().xaxis, 6)
        smart_locate(pl.gca().yaxis, 9)
        pl.savefig(filename)
        pl.clf()
        pl.figure()

    def plot_rowcol_summaries(self, filename, array_data):
        import pylab as pl
        def x_eyes(bads=None):
            # Mark bad fits with an x.
            if bads is None:
                bads = ~s
            pl.scatter(cols[bads], rows[bads], marker='x', edgecolor='gray')
        def limit_args(data, kw={}):
            lo, hi = data.min(), data.max()
            if s.sum() > 1:
                lo, hi = data[s].min(), data[s].max()
            if hi == lo:
                hi = lo + 1
            kw.update({'vmin': lo, 'vmax': hi})
            return kw
        def bin(data, dtype='float'):
            out = np.zeros((n_rows, n_cols), dtype)
            out[rows, cols] = data
            return out
        def imshow_reformat():
            # Tighten boundaries, add labels...
            pl.xlabel('Column')
            pl.ylabel('Row')
            pl.xlim(-0.5, n_cols-0.5)
            pl.ylim(-0.5, n_rows-0.5)

        s = self.ok
        rows, cols = array_data.get_property(['row', 'col'], det_uid=self.det_uid)
        n_rows, n_cols = rows.max()+1, cols.max()+1

        # Init plotting
        pl.figure()
        pl.gcf().set_size_inches(6., 6.)
        pl.subplots_adjust(left=.1, right=.95, top=.95, bottom=.1,
                           hspace=.2, wspace=.3)
        title_fs = 12
        # Time constants...
        #

        pl.subplot(2,2,1)
        z = self.tau * 1e3
        pl.imshow(bin(z), interpolation='nearest', **limit_args(z))
        pl.colorbar()
        x_eyes()
        pl.title('Time constants (ms)', fontsize=title_fs)
        imshow_reformat()

        pl.subplot(2,2,2)
        z = self.tau_err * 1e3
        pl.imshow(bin(z), interpolation='nearest', **limit_args(z))
        pl.colorbar()
        x_eyes()
        pl.title('Time constant errors (ms)', fontsize=title_fs)
        imshow_reformat()

        if self.ok.sum() > 10:
            pl.subplot(2,2,3)
            pl.hist(self.tau[self.ok]*1e3, bins=20) #min(20,self.ok.sum()//10)
            pl.xlabel('Time constant (ms)')
            pl.ylabel('N_dets')
            pl.subplot(2,2,4)
            pl.hist(self.tau_err[self.ok]*1e3, bins=self.ok.sum()//10)
            pl.xlabel('Time constant errors (ms)')
            pl.ylabel('N_dets')

        pl.savefig(filename+'time_const.png')
        pl.clf()
 
        # Positions and stuff
        #
        for i in [0,1]:
            pl.subplot(2,2,1+i)
            z = {0: self.x0_err, 1:self.y0_err}[i]
            z = z * 180*3600/np.pi # to arcseconds
            pl.imshow(bin(z), interpolation='nearest', **limit_args(z))
            pl.colorbar()
            x_eyes()
            imshow_reformat()
            pl.title('%s position RMS' % {0: 'X', 1: 'Y'}[i],
                     fontsize=title_fs)

        pl.subplot(2,2,3)
        z = self.n_obs
        pl.imshow(bin(z), interpolation='nearest')
        pl.colorbar()
        imshow_reformat()
        pl.title('N_obs', fontsize=title_fs)

        pl.savefig(filename+'positions.png')
        pl.clf()

        # Destroy our subplot adjustments
        pl.figure()


class FPTemplateFitter:
    """
    Class for shift/rotate/shearing a template FPFitFile to match a
    target FPFitFile.

    After initializing, set the template to use:

       fitter = FPTemplateFitter()
       fitter.set_template(my_template_fp)
       ok, params = fitter.fit(my_target_fp)

    Those params are stored internally, so you can get the model FP:

       model_for_target = fitter.get_modeled(my_target_fp)
    """

    param_names = ['dx', 'dy', 'theta', 'scale', 'shear_theta', 'shear_scale']

    formats = {'dx': '%9.6f',
               'dy': '%9.6f',
               'scale': '%11.4e',
               'n_dets': '%4i',
               'theta': '%9.6f',
               'shear_scale': '%11.4e',
               'shear_theta': '%9.6f',
               }

    @classmethod
    def from_params(cls, opts, tod_info=None):
        if '_execcfg' in opts:
            tod_id = moby2.scripting.products.get_tod_id(tod_info=tod_info)
            ic = moby2.scripting.execcfg.InputChooser()
            opts1 = ic.get_config(opts['_execcfg'], tod_id=tod_id)
            for k,v in list(opts1.items()):
                if not k in opts:
                    opts[k] = v

        if 'depot' in opts:
            depot = moby2.scripting.get_depot(opts['depot'])
            if not 'structure' in opts:
                opts['structure'] = '{tag}'
            filename = depot.get_full_path(**opts)
        else:
            filename = opts['filename']
        trace(2, 'Loading as template: %s' % filename)

        load_args = opts['column_def']

        pos_data = moby2.util.StructDB.from_column_file(filename, load_args)

        r = opts.get('template_rescale', (1.,1.))
        if 'ok' in pos_data.dtype.names:
            mask = (pos_data['ok'].astype(int) != 0)
        else:
            mask = np.ones(pos_data['x'].shape, bool)

        template_fits = FPFitFile(det_uid=pos_data['det_uid'][mask])
        template_fits.x0[:] = pos_data['x'][mask] * r[0]
        template_fits.y0[:] = pos_data['y'][mask] * r[1]
        template_fits.ok[:] = True
        self = cls()
        self.set_template(template_fits)
        return self

    def set_template(self, template):
        self.template = template
        self.pivot = self.template.x0[self.template.ok].mean(), \
            self.template.y0[self.template.ok].mean()

    @staticmethod
    def _rotate(theta, x, y):
        c, s = np.cos(theta), np.sin(theta)
        return x*c - y*s, y*c + x*s
    
    def model(self, params, x=None, y=None):
        """
        Shift, rotate, shear the current template according to params
        dict.  Return the resulting offsets (x, y).
        """
        dx, dy, theta, scale, sh_theta, sh_scale = params
        scale, sh_scale = np.exp(scale), np.exp(sh_scale)
        # Shift away array center and rescale
        if x is None:
            tp = self.template
            x, y = tp.x0, tp.y0
        out_x, out_y = scale*(x - self.pivot[0]), scale*(y - self.pivot[1])
        # Shear
        out_x, out_y = self._rotate(+sh_theta, out_x, out_y)
        out_x *= sh_scale
        out_x, out_y = self._rotate(-sh_theta, out_x, out_y)
        # Rotate
        out_x, out_y = self._rotate(theta, out_x, out_y)
        # Restore array center and apply additional shift.
        return out_x + self.pivot[0] - dx, out_y + self.pivot[1] - dy

    def model_inverse(self, params, out_x, out_y):
        """
        Inverse of self.model.  Keep it up to date!
        """
        dx, dy, theta, scale, sh_theta, sh_scale = params
        scale, sh_scale = np.exp(scale), np.exp(sh_scale)
        # Remove additional shift.
        x, y = out_x - self.pivot[0] + dx, out_y - self.pivot[1] + dy
        # Unrotate
        x, y = self._rotate(-theta, x, y)
        # Unshear
        x, y = self._rotate(+sh_theta, x, y)
        x /= sh_scale
        x, y = self._rotate(-sh_theta, x, y)
        x, y = x/scale + self.pivot[0], y/scale + self.pivot[1]
        return x, y

    def fit(self, fp, params, trace_level=0):
        """
        Fit positions to a template, which is also an FPFitFile but
        may represent different det_uid.  'params' should be a dict
        like this one:

          params = {
            'shift': True,
            'rotation': True,
            'scale': True,
            'shear': True,
            }
        
        Returns (ok, params).  The fitted_template has the same
        det_uid as self.
        """
        template = self.template
        # Get mask of items that are ok in both the template and fits
        fp_ok = fp.ok.astype('bool').copy()
        _, temp_ok = template.get_property('ok', fp.det_uid)
        fp_ok *= temp_ok

        # Get the template and fits positions for those ok items
        _, x0 = template.get_property('x0', fp.det_uid[fp_ok])
        _, y0 = template.get_property('y0', fp.det_uid[fp_ok])
        x1, y1 = fp.x0[fp_ok], fp.y0[fp_ok]

        self.A = x0,y0
        self.B = x1,y1

        # Identify parameters we want to vary
        free_params = [params.get('shift', True)]*2
        free_params.append(params.get('rotation', True))
        free_params.append(params.get('scale', False))
        free_params.extend([params.get('shear', False)]*2)

        if fp.ok.sum() == 0:
            trace(trace_level+0, 'No items for template fit')
            self.result = False, [0. for f in free_params]
            return self.result

        trace(trace_level+0, 'Fitting template using %i items' % fp_ok.sum())
        
        # Start fit with shift based on mean displacement
        params0 = [x1.mean()-self.pivot[0], y1.mean()-self.pivot[1],
                   0., 0., 0., 0.]

        trace(trace_level+1, 'Starting parameters: %s' % str(params0))

        trace(trace_level+1, 'Free parameters: %s' % str(free_params))
            
        def fit_chi2(params):
            x_model, y_model = self.model(params, x0, y0)
            var = (x1 - x_model)**2 + (y1 - y_model)**2
            #return var.sum()
            # Attenuate contribution of outliers?  Not clear this works...
            mvar = np.median(var)
            var_roll = var * (10*mvar / (10*mvar + var))
            return var_roll.sum()

        # Minimize... start with position or all is lost.
        params1 = params0
        for iters in [0,1]:
            for free_mask in [
                # Fit position only...
                [True , True , False, False, False, False],
                # Fit rotation and scale
                [False, False, True , True , False, False],
                # Fit skew
                [False, False, False, False, True , True ],
                # Fit skew and position
                [True , True , False, False, True , True ],
                # Let everything float
                [True , True , True , True , True , True ]]:
                free = np.array(free_params) * free_mask
                if free.sum() > 0:
                    params1 = moby2.util.fitting.multi_fmin(
                        fit_chi2, params1, free=free, disp=0,
                        xtol=1e-6, ftol=1e-6)
                    trace(trace_level+2, 'params snapshot: %s' % str(params1))

        trace(trace_level+1, 'Final parameters: %s' % str(params1))

        self.result = True, params1
        return self.result

    def check_result(self, opts):
        """
        Check self.result against ranges passed in by user.  User
        passes in a dict with keys like "<name>_range", where <name>
        is one of self.param_names.  The values are the range (lo, hi) of
        acceptable values.  If any range checks fail, the function
        returns false.
        """
        ok, params = self.result
        if not ok:
            return False
        for k, v in zip(self.param_names, params):
            k = '%s_range' % k
            if not k in opts: continue
            if not ((opts[k][0] <= v) and (v < opts[k][1])):
                return False
        return True

    def get_modeled(self, det_uid=None):
        """
        Return a FPFitFile with the modeled detector positions.  Pass
        in the desired det_uid, or the template det_uid will be
        used.
        """
        if det_uid is None:
            det_uid = self.det_uid
        matched = FPFitFile(det_uid=det_uid)
        _, ok = self.template.get_property('ok', matched.det_uid)
        _, x0 = self.template.get_property('x0', matched.det_uid)
        _, y0 = self.template.get_property('y0', matched.det_uid)
        matched.ok = ok

        params = self.result[1]
        matched.x0, matched.y0 = self.model(params, x0, y0)
        return matched

    def make_plots(self, fp, modeled, plot_prefix='./',
                   title=None):
        """
        Show fit quality in a few plots.
        """
        import pylab as pl
        def sane_axes():
            fig.gca().xaxis.set_major_locator(pl.MaxNLocator(4))
            fig.gca().yaxis.set_major_locator(pl.MaxNLocator(5))
            fig.gca().set_aspect('equal', 'datalim')

        DEG = 180./np.pi
        fig = pl.figure()
        fig.set_size_inches(8., 4.)
        pl.subplots_adjust(left=.1, right=.98, top=.85, bottom=.1,
                           hspace=.2, wspace=.3)
        pl.subplot(121)
        tp = self.template
        s, x, y = tp.ok, tp.x0, tp.y0
        pl.scatter(x[s], y[s], marker='o', s=4, alpha=.5)
        pl.xlabel('X')
        pl.ylabel('Y')
        pl.title('Input template')
        sane_axes()

        # The model positions
        pl.subplot(122)
        s, x, y = modeled.ok, modeled.x0 * DEG, modeled.y0 * DEG
        pl.scatter(x[s], y[s], alpha=.2)
        # And the fit positions
        s, x, y = fp.ok, fp.x0 * DEG, fp.y0 * DEG
        pl.scatter(x[s], y[s], marker='x')
        # Now connect them with lines...
        u = fp.det_uid[s]
        ok1, (x1, y1) = modeled.get_property(['x0','y0'], det_uid=u)
        x, y = x[s], y[s]
        for i in ok1.nonzero()[0]:
            pl.plot([x1[i]*DEG, x[i]], [y1[i]*DEG, y[i]], color='k', alpha=.4)
        pl.xlabel('X (deg)')
        pl.ylabel('Y (deg)')
        pl.title('Fitted result')
        sane_axes()

        if title != None:
            pl.figtext(0.5, 0.93, title, va='bottom', ha='center')

        pl.savefig(plot_prefix + 'fit.png')
        pl.figure() # destroy our settings...

    def old_make_plots(self, fp, modeled, plot_prefix='./',
                   title=None):
        """
        Show fit quality in a few plots.
        """
        import pylab as pl
        DEG = 180./np.pi
        pl.figure()
        pl.gcf().set_size_inches(6., 6.)
        pl.subplots_adjust(left=.15, right=.95, top=.90, bottom=.1,
                           hspace=.2, wspace=.3)
        tp = self.template
        s, x, y = tp.ok, tp.x0, tp.y0
        pl.scatter(x[s], y[s], marker='x')
        pl.savefig(plot_prefix + '0template.png')
        pl.clf()

        s, x, y = modeled.ok, modeled.x0 * DEG, modeled.y0 * DEG
        pl.scatter(x[s], y[s], alpha=.2)
        pl.xlabel('X (deg)')
        pl.ylabel('Y (deg)')
        pl.savefig(plot_prefix + '1model.png')
        pl.clf()

        # The model positions
        s, x, y = modeled.ok, modeled.x0 * DEG, modeled.y0 * DEG
        pl.scatter(x[s], y[s], alpha=.2)
        # And the fit positions
        s, x, y = fp.ok, fp.x0 * DEG, fp.y0 * DEG
        pl.scatter(x[s], y[s], marker='x')
        # Now connect them with lines...
        u = fp.det_uid[s]
        ok1, (x1, y1) = modeled.get_property(['x0','y0'], det_uid=u)
        x, y = x[s], y[s]
        for i in ok1.nonzero()[0]:
            pl.plot([x1[i]*DEG, x[i]], [y1[i]*DEG, y[i]], color='k', alpha=.4)
        pl.xlabel('X (deg)')
        pl.ylabel('Y (deg)')
        if title is not None:
            pl.title(title)
        pl.savefig(plot_prefix + '2fit.png')
        pl.figure() # destroy our settings...


    # Formatted output...

    def get_ascii(self, names=None, params=None):
        if names is None:
            names = self.param_names
        if params is None:
            params = self.result[1]
        idx = [self.param_names.index(f) for f in names]
        text = [ self.formats.get(n, '%11.4e') % params[i]
                 for n,i in zip(names,idx) ]
        return ' '.join(text)

    @staticmethod
    def write_fit_list(filename, keys, fits, format=None):
        if format == 'fits':
            columns = list(zip(*[f.result[1] for f in fits]))
            col_defs = ([('id', keys), ('ok', [int(f.result[0]) for f in fits])] + 
                        list(zip(fits[0].param_names, columns)))
            db_out = moby2.util.StructDB.from_data(
                col_defs, formats=fits[0].formats)
            db_out.to_fits_table(filename)
        else:
            if isinstance(filename, basestring):
                filename = open(filename, 'w')
            names = fits[0].param_names
            filename.write('# %s\n' % ' '.join(names))
            for key, fit in zip(keys, fits):
                text = fit.get_ascii(names=names)
                filename.write('%s %s\n' % (key, text))
    

        
