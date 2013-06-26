from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Simple common mode removal.

Identify light / dark dets and the common mode of each.
"""

import pickle
import time

import numpy as np

# Dodge immediate import of pylab; let user pick a back end.
pl = None

def import_pylab():
    global pl
    if pl is None:
        import pylab
        pl = pylab


from moby2.util import rethread
import moby2.util.log as psLib

def trace(level, text):
    psLib.trace('moby', level, text)

class StatModeSet:
    default_rules = {
        #Pathology masking
        'patho_margin': 50,
        'patho_jump': 5e6,
        #Light mode assignment
        'light_noise_sigma': 1.,
        'light_proj_sigma': 1.,
        #Dark mode assignment
        'dark_lightness_sigma': 2.,
        'dark_darkness_sigma': 2.,
        #Noise cut, pW
        'noise_cut': 0.005,
        }
    
    # archive support -- store enough information to be able
    #  use the results and to run the plot* methods.
    depot_keys = ['rows', 'cols', 'shape',
                  'projs', 'norms', 'assignments', 'modes', 'noise',
                  'cal', 'darkCand', 'lightCand',
                  'masks', 'rules', 'lightLimits']

    def __init__(self, data, cal=None, lightCand=None, darkCand=None,
                 title='', shape=None, rows=None, cols=None, rules=None):
        """
        To do things properly for a tod, use the forTOD classmethod.
        
        data        array with shape (n_det, n_sample)
        cal         float vector containing to-pW calibration factors;
                    used to determine light mode membership (detectors
                    with cal==0 will be ignored in light mode assignment)
        lightCand   boolean vector of length n_det indicating dets to use
                    for light mode estimation.
        darkCand    boolean vector of length n_det indicating dets to use
                    for dark mode estimation.
        title       title string for plots
        shape       (nrow, ncol) for plotting matrix of assignments
        rows, cols  vectors of length n_det giving row and col of each det.
        """
        import_pylab()
        self.data = data
        self.cal = cal
        self.darkCand = darkCand
        self.lightCand = lightCand
        self.masks = {}
        self.title = title
        self.shape = shape
        self.rows = rows
        self.cols = cols
        self.rules = {}
        self.rules.update(self.default_rules)
        if rules is not None:
            self.rules.update(rules)

    @classmethod
    def forTOD(cls, tod, cal, lightCand=None, darkCand=None, rules=None):
        """
        Create statMode to process given TOD and calibration.
        
        tod        TOD.
        cal        calibration factors (length tod.data.shape[0]), with 0s
                   for uncalibrated detectors.

        darkCand   boolean vector of detectors to use for dark mode estimation.
                   Pass None to use default of (tod.rows==33).
        lightCand  boolean vector of detectors to use for light mode estimation.
                   Pass None to use default to (cal != 0).
        """
        ndet = tod.data.shape[0]
        nr, nc = 33, 32
        rows = np.arange(ndet) // nc
        cols = np.arange(ndet) % nc
        if darkCand is None:
            darkCand = (rows==nr-1)
        if lightCand is None:
            lightCand = (cal != 0)*(~darkCand)
        return cls(tod.data, cal=cal, lightCand=lightCand, darkCand=darkCand,
                   title=tod.info.name, shape=(nr,nc), rows=rows, cols=cols,
                   rules=rules)
        
    def decodeFiducialDets(self, fd):
        """
        Use dictionary 'fd' to select live mode candidates.
        """
        dets = self.lightCand
        if 'rows' in fd:
            dets *= np.array([r in fd['rows'] for r in self.rows])
        if 'cols' in fd:
            dets *= np.array([r in fd['cols'] for r in self.cols])
        if 'dets' in fd:
            rc = list(zip(self.rows, self.cols))
            dets *= np.array([_rc in rc for _rc in fd['dets']])
        self.lightCand = dets

    def go(self):
        if 1: #try:
            for m in ['calcModes', 'getProjs', 'getModeStats', 'assignDets', 'cutByNoise']:
                trace(3, 'doing %s' % m)
                getattr(self, m)()
            trace(3, 'done')
        else: #except:
            return False, 'darn.'
        return True, 'phew.'

    def pathologyMask(self):
        n = self.rules['patho_margin']
        #amp = amax(abs(self.data[:,n:] - self.data[:,:-n]), axis=1)
        amp = np.array(rethread.foreach((lambda y: (abs(y[n:] - y[:-n])).max()), self.data))
        return amp < self.rules['patho_jump']

    def calcModes(self):
        """
        Simplified mode identification -- user passes in light and dark suggestions.
        """
        def mean_axis0(data, mask=None):
            if mask is None: mask = np.ones(data.shape[0], dtype='bool')
            m = np.zeros(data.shape[1:],dtype='float')
            for i in mask.nonzero()[0]:
                m += data[i]
            return m / np.sum(mask)
        # Mask stupid dets and compute common mode
        data = self.data
        pm = self.pathologyMask()
        zm = (rethread.pmax(data) - rethread.pmin(data)) != 0
        self.masks['path'] = pm
        self.masks['zm'] = zm
        # Common mode, normalized
        cm_dets = pm*zm*self.lightCand
        self.masks['cm'] = cm_dets
        cm0 = mean_axis0(data, cm_dets)
        cm0_N = (cm0**2).sum()**.5
        cm0 /= cm0_N
        # Dark mode, normalized
        dm_dets = pm*zm*self.darkCand
        self.masks['dm'] = dm_dets
        if dm_dets.sum() > 0:
            dm0 = mean_axis0(data, dm_dets)
            dm0_N = (dm0**2).sum()**.5
            dm0 /= dm0_N
        else:
            dm0 = 0 * cm0
            dm0_N = 1.
        # Dot
        cm1 = cm0 - dm0*np.dot(cm0,dm0)
        cm1 /= (cm1**2).sum()**.5
        self.modes = {'light': cm0, \
                      'dark': dm0, \
                      'light1': cm1/(cm1**2).sum()**.5 }
        self.mode_norms = {'light': cm0_N, 'dark': dm0_N }

    def getProjs(self):
        """
        Compute projection of all modes in self.modes against all dets in tod
        (or optionally just those in dets).
        """
        data = self.data
        norms = rethread.pnorm(data)**0.5
        norms[norms==0.] = 1.
        projs = {}
        for k, m in zip(list(self.modes.keys()), list(self.modes.values())):
            projs[k] = rethread.pdot(data, m) / norms
        self.projs, self.norms = projs, norms
        return projs, norms

    def getModeStats(self):
        """
        Estimate mode cut parameters based on mode projections and
        candidate dets.
        """
        # Basic non-pathological dets
        mask0 = self.masks['zm'] * self.masks['path']
        # Recalibrate from projection to power in pW*n
        pc = self.projs['light'] * abs(self.cal) * self.norms
        # Rest of signal is noise, in pW*n
        ns = (1.-self.projs['light']**2)**0.5 * abs(self.cal) * self.norms
        ns0 = np.median(ns[mask0*self.lightCand])
        # Negative light projections are discouraged, as are outliers.
        s1 = (pc > 0)*(ns > 0)*(ns < ns0*10)*mask0*self.lightCand
        # Make a box around the light detectors
        #dn = (ns[s1]**2).mean()**0.5
        n, dn = ns[s1].mean(), ns[s1].std()
        p, dp = pc[s1].mean(), pc[s1].std()
        self.lightLimits = {'n': (n, dn), 'p': (p, dp)}

        # Dark stats, no calibration information
        pd, pl = self.projs['dark'], self.projs['light1']
        m = mask0*self.darkCand
        self.darkLimits = {
            'd': (1 - pd[m]**2).mean(),
            'l': (pl[m]**2).mean(),
        }
        return self.lightLimits, self.darkLimits
        
    def assignDets(self, ):
        (n, dn), (p, dp) = self.lightLimits['n'], self.lightLimits['p']
        # Basic non-pathological dets
        mask0 = self.masks['zm'] * self.masks['path']
        # Recalibrate from projection to power in pW*n
        pc = self.projs['light'] * abs(self.cal) * self.norms
        # Rest of signal is noise, in pW*n
        ns = (1.-self.projs['light']**2)**0.5 * abs(self.cal) * self.norms
        p_sigma, n_sigma = self.rules['light_proj_sigma'], \
                           self.rules['light_noise_sigma']
        lightOk = (abs(pc-p) < dp*p_sigma)*((ns-n) < dn*n_sigma) #abs(ns-n)<dn*n_sigma)
        self.masks['lightOk'] = lightOk

        # NEW VERSION          VVVVVV
        # Dark modes are uncalibrated, so we consider the normalized
        # power levels in the light1 and dark modes, and the residual.
        # These satisfy PL + PD + PN = 1.
        pd, pl = self.projs['dark'], self.projs['light1']
        m = mask0*self.darkCand
        # Minimum darkness:
        d_sig = self.rules['dark_darkness_sigma']
        d_lim = 1 - d_sig * self.darkLimits['d']
        #d_lim = 1 - (1-darkness[m]).mean() * d_sig
        # Maximum lightness:
        l_sig = self.rules['dark_lightness_sigma']
        #l_lim = lightness[m].mean() * l_sig
        l_lim = l_sig * self.darkLimits['l']
        # Cut
        darkOk = (pd**2 > d_lim)*(pl**2 < l_lim)*(pd>0)
        self.masks['darkOk'] = darkOk
        # Assignments are somewhat more strict than the Ok masks.
        self.assignments = {
            'light': lightOk*mask0,
            'dark': darkOk*~lightOk*mask0,
            }
        return self.assignments
    
    def cutByNoise(self):
        """
        Call after modes have been removed to compute residual noise
        and make further cuts to live assignment.
        """
        # Compute RMS based on norms and projections
        self.noise = (1-self.projs['dark']**2-self.projs['light1']**2)**0.5
        self.noise *= self.norms*abs(self.cal)/self.data.shape[1]**0.5
        self.masks['noiseCut'] = (self.noise < self.rules['noise_cut'])
        self.assignments['light'] *= self.masks['noiseCut']

    def remove(self, data=None):
        """
        Remove dark and light1 modes from data (or self.data).  Uses
        self.norms and self.projs.
        """
        if data is None:
            data = self.data
        ndet = data.shape[0]
        dets = [i for i in range(ndet)]
        for k in ['light1', 'dark']:
            mode = self.modes[k].astype('float32')
            weight = (self.norms*self.projs[k]).astype('float32')
            #xobyLib.removeModeDetsWeighted(data, dets, mode, weight)
            data -= weight[:,None] * mode
            #data[:,:] -= (self.norms*self.projs[k]).reshape(-1, 1)*m

    def recalibrate(self, apply=False):
        """
        Use light mode detectors to estimate true common mode strength
        in pW.  Use it to obtain calibration for all dets based on light1
        mode projections and norms.
        """
        #Estimate calibrated common mode strength, with RMS error
        fd = self.assignments['light']
        cal, alpha, norms = abs(self.cal), self.projs['light'], self.norms
        z = (cal * alpha * norms)[fd]
        pW, dpW = z.mean(), z.std()
        #Obtain calibration for all detectors based on this
        recal = pW / alpha / norms
        recal[(alpha*norms==0)] = 0
        self.recal = {
            'dets': fd,
            'cal': recal,
            'pW': (pW, dpW)
            }
        # Apply?
        if apply:
            sel = (self.recal['cal'] != 0)*(~self.lightCand)
            self.cal[sel] = self.recal['cal'][sel]
            self.assignDets()
        return self.recal

    def _guessDims(self):
        shape, rows, cols = self.shape, self.rows, self.cols
        if shape is None:
            shape = (33, n//33)
        if rows is None or cols is None:
            rows = np.array([i//shape[1] for i in range(n)])
            cols = np.array([i%shape[1] for i in range(n)])
        return shape, rows, cols
    
    def _prePlot(self, clf):
        if clf:
            pl.gcf().set_size_inches(5, 5)
            pl.clf()

    def _postPlot(self, title, filename, xlabel=None, ylabel=None):
        if title is None:
            pl.title(self.title)
        else:
            pl.title('%s - %s' % (self.title, title))
        if xlabel is not None: pl.xlabel(xlabel)
        if ylabel is not None: pl.ylabel(ylabel)
        if filename is not None: pl.savefig(filename)

    def plotModes(self, title=None, filename=None, clf=True):
        """
        Plot the computed light and dark common modes.
        """
        def downSample(data, points=2000):
            d = data.shape[-1] // points
            if d <= 1: return data
            return data[np.arange(0, data.shape[-1], d)]
        self._prePlot(clf)
        for name, color in [('dark', 'k'), ('light', 'b')]:
            if name in self.modes:
                pl.plot(downSample(self.modes[name]), color, label=name)
        pl.legend()
        self._postPlot(title, filename)

    def plotLightProjs(self, title=None, filename=None, clf=True):
        self._prePlot(clf=clf)
        pc = self.projs['light'] * abs(self.cal) * self.norms
        ns = (1.-self.projs['light']**2)**0.5 * abs(self.cal) * self.norms
        # Make a box around the light detectors
        #n, dn
        n, dn = self.lightLimits['n']
        p, dp = self.lightLimits['p']
        
        m = (self.cal!=0)*self.masks['zm']*self.masks['path']
        ll = self.assignments['light']
        for c,pop,lab in [
            ('r', ll, 'light'),
            ('k', ~ll, 'bad') ]:
            if np.sum(pop*m) > 0:
                pl.scatter(pc[pop*m], ns[pop*m], color=c, label=lab, s=3)
        dp *= self.rules['light_proj_sigma']
        dn *= self.rules['light_noise_sigma']
        pl.plot([p+dp,p+dp],[0,n+dn],'k')
        pl.plot([p+dp,p-dp],[n+dn,n+dn],'k')
        pl.plot([p-dp,p-dp],[n+dn,0],'k')
        
#        corners = [(1,1),(1,-1),(-1,-1),(-1,1),(1,1)]
#        for (x0,y0),(x1,y1) in zip(corners[:-1],corners[1:]):
#            plot([p+x0*dp,p+x1*dp],[n+y0*dn,n+y1*dn],color='k')
        pl.legend()
        # Restrict limits
        (x0,x1), (y0,y1) = pl.xlim(), pl.ylim()
        pl.xlim(0, min(x1, p+dp*3))
        pl.ylim(0, min(y1, n+dn*3))
        self._postPlot(title, filename, 'Light mode power', 'Noise power')
                                                         
    def plotDarkProjs(self, title=None, filename=None, clf=True):
        self._prePlot(clf=clf)
        m = self.masks['zm']*self.masks['path']
        # Dark projection is main signal 
        d = self.projs['dark']
        n = (1. - self.projs['light1']**2 - d**2)**0.5
        ll = self.assignments['light']
        for c,pop,lab in [
            ('b', (self.cal==0), 'no_cal'),
            ('g', self.darkCand, 'dark_cand'),
            ('k', (self.cal!=0), 'yes_cal'),
            ('r', ll, 'live') ]:
            if np.sum(pop*m) > 0:
                pl.scatter(d[pop*m], n[pop*m], color=c, label=lab, s=3)
        pl.legend()
        self._postPlot(title, filename, 'Dark CM projection',
                       'Non-mode projection')
        
    def plotAssignments(self, title=None, filename=None, clf=True):
        """
        Plot a matrix showing which detectors were identified as light, dark,
        neither, or not-present.
        """
        shape, rows, cols = self._guessDims()
        n = self.assignments['light'].shape[0]
        self._prePlot(clf)
        role = np.zeros(shape)
        d, l = self.assignments['dark'], self.assignments['light']
        role[rows, cols] = 1
        role[rows[d], cols[d]] = 2
        role[rows[l], cols[l]] = 3
        pl.imshow(role, cmap=pl.cm.gray, vmin=0., vmax=3., interpolation='nearest')
        self._postPlot(title, filename)

    def plotNoise(self, title=None, filename=None, clf=True):
        self._prePlot(clf)
        shape, rows, cols = self._guessDims()
        m = self.noise != 0
        n = self.noise[m]
        noise_max = min(n.max(), np.median(n)*4)
        nm = np.zeros(shape, 'float')
        nm[rows[m], cols[m]] = n
        pl.imshow(nm, interpolation='nearest',
                     vmin=0., vmax=noise_max)
        pl.colorbar()
        self._postPlot(title, filename)

    def standard_plots(self, prefix):
        self.plotModes(filename='%scm.png' % prefix)
        self.plotDarkProjs(filename='%sprojsd.png' % prefix)
        self.plotLightProjs(filename='%sprojsl.png' % prefix)
        self.plotAssignments(filename='%sassign.png' % prefix)
        self.plotNoise(filename='%snoise.png' % prefix)


    # depot support

    def writeToPath(self, path):
        ofile = path + 'statmode.pik'
        output = {}
        for k in self.depot_keys:
            output[k] = getattr(self, k)
        pickle.dump(output, open(ofile, 'w'))

    @classmethod
    def readFromPath(cls, path, *args):
        ifile = path + 'statmode.pik'
        input = pickle.load(open(ifile))
        o = cls(None)
        for k in o.depot_keys:
            setattr(o, k, input[k])
        return o
            
