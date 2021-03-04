from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
This file contains a driver program for fitting an "opacity" model to
measured planet amplitudes to produce a calibration model (TODAbsCal).
"""

import moby2
import numpy as np
import pylab as pl

from . import planet_cal
from moby2.analysis import beam_ana

DEG = np.pi/180

class RecalModel:
    """
    Holds a simple recalibration model, currently only smart enough to
    undo atmospheric opacity (and similar effects).
    """
    @classmethod
    def write(cls, filename, data):
        output = moby2.util.MobyDict()
        output.update(data)
        output.write(filename)

    @classmethod
    def from_file(cls, filename):
        self = cls()
        self.data = moby2.util.MobyDict.from_file(filename)
        self.cal0 = self.data['cal0_pW_per_K'] * 1e-6  # pW per uK.
        self.tau_w = self.data['tau_w']
        return self

    def get_cal(self, loading, alt_deg=None):
        if alt_deg is not None:
            # assume 'loading' PWV and compute airmass.
            loading = loading / np.sin(alt_deg*np.pi/180)
        return self.cal0**-1 * np.exp(loading*self.tau_w)

    def get_cal_for_catalog(self, cat, columns=['tod_name', 'loading']):
        cal = self.get_cal(cat[columns[1]])
        odb = moby2.util.StructDB.from_data([(columns[0], cat[columns[0]]),
                                             ('cal_uK', cal)],
                                            formats={'cal_uK': '%12.5e'})
        return odb


def wlinfit(x, y, w):
    v = np.array([np.ones(x.shape), x])
    M = np.dot(v*w, np.transpose(v))
    b = np.dot(v*w, y)
    return np.dot(np.linalg.inv(M), b)[::-1]

def get_catalog(cat_opts, catalogs):
    key = id(cat_opts)
    if not key in catalogs:
        catalogs[key] = moby2.scripting.get_obs_catalog(cat_opts)
    return catalogs[key]

def update_cascade(defaults, params):
    output = defaults.copy()
    for k in list(params.keys()):
        if k in defaults and isinstance(defaults[k], dict):
            if params.get(k, {}).get('_reset'):
                output[k] = params[k]
            else:
                output[k] = update_cascade(defaults[k], params[k])
        else:
            output[k] = params[k]
    return output
        
def load_amplitudes(params, info_dict):
    """Read planet peak heights from some kind of file.  The ``info_dict``
    will be used to format the filename according to what data subset
    is being analyzed.  The ``params`` argument is a dict that looks
    like one of these:

      {'type': 'quick_pick',
       'filename': 'quick_beam.pik',
       'field': 'gaussian_amp',
       'masks': ...,      #optional
       'also_load': ...,  #optional
      }

      {'type': 'column_file',
       'filename': 'peaks_{fcode}.txt',
       'columns': [0,1],  #optional
      }

      {'type': 'solid_angle',
       'filename': 'solid_angle/table.fits'
      }

    """
    atype = params['type']
    also = {}
    if atype == 'quick_pik':
        filename = params['filename'].format(**info_dict)
        bol = beam_ana.BeamObsList(filename)
        # Mask?
        mask = np.ones(len(bol), bool)
        for m in params['masks']:
            mask *= bol.get(*m).astype('bool')
        bol = bol.select(mask)
        basenames = bol.get('init', 'basename')
        amps = bol.get(*params['field'])
        if 'also_load' in params:
            also = dict([(k2,bol.get(k1,k2)) for k1,k2 in params['also_load']])
            
    elif atype == 'column_file':
        filename = params['filename'].format(**info_dict)
        columns = params.get('columns', [0,1])
        basenames, amps = moby2.util.ascii.read_columns(filename, columns)

    elif atype == 'get_solid_angle':
        # Output from the get_solid_angle script.
        filename = params['filename'].format(**info_dict)
        db = moby2.util.StructDB.from_fits_table(filename)
        basenames, amps = db['name'], db['quick_peak']

    # Frequency band sanity.
    if 'freq_code' in info_dict:
        extracted = []
        for n in basenames:
            n, f = beam_ana.util.extract_basename(n), beam_ana.util.extract_freq(n)
            extracted.append((n,f))
        basenames, tod_freqs = list(map(np.array, list(zip(*extracted))))
        s = (tod_freqs == info_dict['freq_code'])
        basenames, amps = basenames[s], amps[s]
        for k in list(also.keys()):
            also[k] = also[k][s]
#    else:
#        if len(set(tod_freqs)) > 1:
#            print 'You selected maps covering multiple frequency bands:', set(tod_freqs)
#            print "Specify your choice of band, e.g. 'info': {'freq_code': 'f150'}."
#            raise RuntimeError
#

    return basenames, amps, also
   
def get_source_temperature(params, info_dict):
    ttype = params['type']
    if ttype == 'constant':
        return params['T']
    # Is this a planet name?
    Tmodel = planet_cal.get_brightness_model(ttype)
    if Tmodel is not None:
        freq = params['frequency_GHz']
        T_uK = Tmodel.get_uK_minus_background(freq)
        return T_uK / 1e6
    print('Failed to decode source: ', params)

def get_beam_solid_angle(params, info_dict):
    assert params['type'] == 'constant'
    return params['solid_angle_sr']

def fmt_err(x, dx, fmt='%6.3f'):
    return fmt%x + ' +/- ' + fmt%dx


# The factoring of the original binary such that do_cal_job is called,
# repeatedly, by main, was not done in a thoughtful way...  feel free
# to empower do_cal_job so it has more potential additional
# usefulness.

def do_cal_job(cfg, job_name, catalogs=None):
    if catalogs is None:
        catalogs = {}
    print()
    print('--------')
    print(job_name)
    print()
    cat_dict = get_catalog(cfg['info_catalog'], catalogs)
    info_dict = cfg['info']
    info_dict['job'] = job_name
    outputm = moby2.scripting.OutputFiler(
        prefix=cfg['output_prefix'].format(**info_dict))
    # Load / compute amplitudes.
    tod_names, amp, also = load_amplitudes(cfg['amplitudes'], info_dict)
    
    subcat = cat_dict.select_inner({'tod_name': tod_names})
    assert(np.all(subcat >= 0))
    cat_dict = cat_dict[subcat]
    print('Loaded %i amplitudes' % len(tod_names))

    # Compute the ephemeris and calibrated amp/K for all these obs...
    job_data = {
        'T_SOURCE': get_source_temperature(cfg['source'], info_dict),
        'OMEGA_BEAM': get_beam_solid_angle(cfg['beam'], info_dict),
        }

    if 'solid_angle' in cfg['source']:
        job_data['OMEGA_SOURCE'] = 0*cat_dict['ctime'] + cfg['source']['solid_angle']
    else:
        job_data['OMEGA_SOURCE'] = planet_cal.get_planet_solid_angle(
            cfg['source']['type'], cat_dict['ctime'])

    T_DILUTED = (job_data['T_SOURCE'] * job_data['OMEGA_SOURCE'] /
                 job_data['OMEGA_BEAM'])
    # Naive conversions
    cals = amp / T_DILUTED

    # Allow job to select certain observations from the list.
    include = moby2.scripting.get_selection(
        cfg['include'], cat_dict, {'amp': cals}, also)

    # Sub-select based on a white list?
    if cfg.get('whitelist'):
        whitelist = get_catalog(cfg['whitelist'], catalogs)
        idx = whitelist.select_inner({'tod_name': tod_names})
        s = (idx>=0)
        print('Whitelist keeps %i of %i otherwise included obs.' % (
            (include*s).sum(), include.sum()))
        include *= s

    cat_dict = cat_dict[include]
    tod_names, cals = tod_names[include], cals[include]
    T_DILUTED = T_DILUTED[include]
    for k,v in list(also.items()):
        also[k] = v[include]

    print('Restricted to %i observations' % len(tod_names))
    print('  -> T_diluted: (%s) K' % fmt_err(
        T_DILUTED.mean(), T_DILUTED.std(), '%.4f'))

    # The cuts, the cuts.
    mask = moby2.scripting.get_selection(
        cfg['cuts'], cat_dict, {'amp': cals}, also)

    pl.gcf().set_size_inches(6.5, 5)
    print('Kept %i out of %i' % (mask.sum(), len(mask)))
    print('  -> T_diluted: (%s) K' % fmt_err(
        T_DILUTED[mask].mean(), T_DILUTED[mask].std(), '%.4f'))
    if np.any(cals[mask]<=0):
        mask *= (cals>0)
        print('Killed non-positive amplitudes; %i items remain.' % mask.sum())
    x, y = cat_dict['loading'], np.log(cals*mask + (~mask).astype(int))
    if 0: #np.all(noise[s] != 0):
        w = noise[s]**-2
        w0 = np.median(w)
        markersize = w0/w*30
        alpha = .8 * (w/w.max())**.5
        n = np.linspace(0,1.,100)
        c = np.zeros((4,len(w)))
        c[0:2], c[2], c[3] = .5, 1, alpha
        c = c.transpose()
        alpha = None
    else:
        w = 1
        markersize = 10
        alpha = .5
        c = 'b'

    fit_type = cfg.get('fit_type', 'linear')
    if fit_type == 'linear':
        p = wlinfit(x[mask], y[mask], w)
    elif fit_type == 'constant':
        print('Fitting flat line!')
        p = [0., np.log(np.exp(y[mask]).mean())]

    scatter = (y[mask] - np.polyval(p, x[mask])).std()
    print('Log-linear fit:')
    print('  p = %e %e' % tuple(p))
    pl.scatter(x[mask], np.exp(y)[mask], c=c, s=markersize)
    if (~mask).sum() > 0:
        pl.scatter(x[~mask], cals[~mask], alpha=.2, marker='x', c='red')
    for s, filename in [(~mask, outputm.get_filename('rejected.txt')),
                        (mask, outputm.get_filename('accepted.txt'))]:
        fout = open(filename ,'w')
        for i in s.nonzero()[0]:
            fout.write('%s %8.3f %11.4e\n' %
                       (tod_names[i], x[i], np.exp(y[i])))
        fout.close()
    _xx = np.arange(0, x.max()*1.1, .1)
    pl.plot(_xx, np.exp(np.polyval(p, _xx)))
    pl.xlim(0, x.max()*1.05)
    pl.ylim(0, np.exp(y[mask].max())*1.1)
    def text(loc, text):
        xl, yl = pl.xlim(), pl.ylim()
        wx, wl = xl[1] - xl[0], yl[1] - yl[0]
        _x, _y = xl[1]-wl*.05, np.exp(p[1])/2+loc*wl*.2
        pl.text(_x, _y, text, ha='right', fontsize=10)
        print('  '+ text)
    #pl.text(x.max(), p[1], "Slope = %.4f / mm" % p[0], ha='right')
    pl.title('Fit for %10s (%i of %i points)' % (
            job_name, mask.sum(), len(mask)))
    text(-.0, "slope = %.4f / mm" % p[0])
    T0 = np.exp(np.polyval(p, 0.))
    Ttyp = np.exp(np.polyval(p, 0.5/np.sin(50*DEG)))
    text(-.2, "peak(loading=0) = %.5e" % T0)
    text(-.4, "peak(pwv=.5,alt=50) = %.5e" % Ttyp)
    text(-.6, "scatter = %.3f" % scatter)
    pl.ylabel('pW / [delta K CMB]')
    pl.xlabel('PWV / sin(alt) [mm]')
    # Inset plot with scatter histogram...
    ax = pl.axes([.2, .2, .22, .22])
    ax.patch.set_alpha(0.7)
    bins = np.arange(-0.21, 0.215, .02)
    delt = np.clip(y[mask] - np.polyval(p, x[mask]), bins[0]+1e-5, bins[-1]-1e-5)
    ax.hist(delt, bins=bins)
    for _x in [bins[1], bins[-2]]: pl.axvline(_x, c='k', ls='dotted')
    outputm.savefig('cal.png')
    model = moby2.util.MobyDict()
    model.update(job_data)
    model.update({'dataset': job_name,
                  #'T_SOURCE': T,
                  #'OMEGA_BEAM': cfg['OMEGA_BEAM'],
                  'cal0_pW_per_K': np.exp(p[1]),
                  'tau_w': -p[0],
                  'scatter': scatter})
    RecalModel.write(outputm.get_filename('model.par'), model)

    # Bootstrap!
    nbs = 1000
    bs = []
    no = mask.sum()
    for i in range(nbs):
        datai = np.floor(np.random.uniform(size=no)*no).astype('int')
        x1, y1 = x[mask][datai], y[mask][datai]
        pbs = np.polyfit(x1, y1, 1)
        bs.append(pbs)
    cov = np.cov(np.transpose(bs))
    # Slope error?
    print('Bootstrap:')
    print('  amp error:     %.4e' % cov[0,0]**.5)
    print('  slope error:   %.4e' % cov[1,1]**.5)
    # Eval...
    for _x in cfg.get('bootstrap_eval', []):
        _y0 = np.exp(np.polyval(p, _x))
        _y_all = np.array([np.polyval(_p, _x) for _p in bs])
        print(' @%.3f:  %.4e' % (_x, np.exp(_y_all).std() / _y0))


def main(args_in):
    import optparse as o
    o = o.OptionParser()
    opts, args = o.parse_args(args_in)

    if len(args) < 1:
        o.error('Expected: a single configuration file name.')

    cfg = moby2.util.MobyDict.from_file(args.pop(0))

    fcfg = cfg['fit_opacity']
    job_names, _ = list(zip(*fcfg['jobs']))
    job_cfgs = dict(fcfg['jobs'])

    if len(args) == 0:
        jobs = job_names
    else:
        jobs = args

    catalogs = {}
    for job_name in jobs:
        if not job_name in job_cfgs:
            o.error('job name "%s" not defined in config file.' % job_name)
        cfg = fcfg['defaults'].copy()
        cfg = update_cascade(fcfg['defaults'], job_cfgs[job_name])
        do_cal_job(cfg, job_name, catalogs)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

