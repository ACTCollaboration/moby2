from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from past.builtins import basestring

from moby2.libactpol import freq_space_waterfall
from moby2.libactpol import time_space_waterfall
from moby2.scripting import products
import numpy 
try:
    import pylab
except:
    pass
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import cm, patches
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
import moby2
moby2.pointing.set_bulletin_A()
ArrayData = moby2.tod.array_data.ArrayData
np = numpy

def plot_with_cuts(tod, det, interactive=True):
    if interactive: pylab.ion()
    else: pylab.ioff()
    sel = numpy.ones(tod.nsamps, dtype=bool)
    for c in tod.cuts.cuts[det]:
        sel[c[0]:c[1]] = False
    t = tod.ctime - tod.ctime[0]
    pylab.plot(t[sel], tod.data[det][sel], "k.", markersize=1)
    pylab.plot(t[~sel], tod.data[det][~sel], "r.", markersize=1)
    return pylab.gca()

class freqSpaceWaterfall( object ):
    """
    Functions to generate the waterfall data (power spectrum) and 
    shelve it to a depot
    """

    def __init__( self, tod, nfreq = 400, fmin = 0.1, fmax = 200, logx = True):
        """
        @brief Generates the waterfall matrix from the precalculated power spectrums
        @param tod    TOD to analyze.
        @param nfreq  Number of frequency points in waterfall plot.
        @param fmin   Minimum frequency to consider in waterfall
        @param fmax   Maximum frequency to consider in waterfall
        @param logx   Whether to use a logarithmic or linear scale in the x axis
        """
        self.ndet, self.ndata = tod.data.shape
        self.data = tod.data
        self.sampleTime = (tod.ctime[-1]-tod.ctime[0])/(tod.nsamps-1)
        self.name = tod.info.basename
        self.array = tod.info.array
        self.resultsDict = {}
        self.rows = tod.info.array_data["row"]
        self.Nrows = np.unique(self.rows).size
        self.cols = tod.info.array_data["col"]
        self.Ncols = np.unique(self.cols).size
        self.dets = tod.info.det_uid
        self.arrayQual = None
        self.arrayQual2 = None
        self.qual = None
        self.keys = []
        self.sort = []
        for i in range(len(self.rows)):
            self.sort.append((self.rows[i], self.cols[i]))
        self.sort = numpy.array(self.sort, dtype = [('rows', int),('cols', int)])
        self.nfreq = nfreq
        self.logx = logx
        self.matfreqs = numpy.zeros(self.nfreq)
        if (logx): 
            self.matfreqs = numpy.logspace(numpy.log10(fmin), numpy.log10(fmax), nfreq)
        else:
            self.matfreqs = numpy.linspace(fmin, fmax, nfreq)
        self.mat = freq_space_waterfall(self.data, self.matfreqs, self.sampleTime)


    def plot( self, selection = None, vmin = None, vmax = None, title = None, filename = None, 
              rowDominance = False, units = 'DAC', separators = True, show = True, logy = True,
              ratio = 1.2, size = [10.0, 10.5], dpi = None, linTickSep = None, hideYLabel = False,
              forceAll = False, **kargs):
        """
        @brief  plot function to visualize the waterfall plot.
        @param selection    bool array with selection of detectors to include in plot.
        @param vmin         Minimum value in scale of waterfall plot.
        @param vmax         Maximum value in scale of waterfall plot.
        params logy         Whether to use logarithmic scale in the Y axis.
        @param title        Title in plot.
        @param filename     Filename where to save the plot.
        @param rowDominance Sort rows in waterfall by detector rows or columns.
        @param units        Units of data in TOD (DAC or uK) for label.
        @param separators   Whether to draw a line separating detector rows or columns.
        @param show         Whether to show or not the plot on screen.
        @param ratio        Aspect ratio of the plot
        @param size         Size of the plot
        @param dpi          Resolution of the plot
        @param linTickSep   Separation of the frequency ticks for a linear space plot.
        @param hideYLabel   Do not display Y label.
        @param forceAll     Force all detectors to appear. Unselected detectors appear in black.
        """
        if show: pylab.ion()
        else: pylab.ioff()
        if self.mat is None: self.analyze()
        if selection is None: 
            sel = numpy.ones(len(self.dets), dtype = 'bool')
        elif forceAll:
            sel = numpy.ones(len(self.dets), dtype = 'bool')
            fsel = selection[self.dets]
        else: sel = selection[self.dets]

        if linTickSep is None:
            linTickSep = (self.matfreqs[-1] - self.matfreqs[0])/5
            p = numpy.power(10,numpy.floor(numpy.log10(linTickSep)))
            linTickSep = numpy.floor(linTickSep/p)*p

        if rowDominance:
            order = numpy.argsort(self.sort[sel], order = ['rows', 'cols'])
            tmp = self.rows[sel][order]
            ylabel = "Row"
        else:
            order = numpy.argsort(self.sort[sel], order = ['cols', 'rows'])
            tmp = self.cols[sel][order]
            ylabel = "Column"

        # Find frequency axis which has a resolution given by self.nfreq
        # Produce a linear scale over a logaritmic scale.
        if self.logx:
            f = numpy.log10(self.matfreqs)
            step = (f[-1]-f[0])/(float(self.nfreq-1))
            ini = numpy.floor(f[0])
            end = numpy.floor(f[-1])
            xt = (numpy.arange(ini, end+1) - f[0]) / step
            xtl = numpy.array(numpy.power(10., 
                     numpy.arange(ini, end+1, dtype = "int")), dtype = 'str')

        
        if logy: mat = numpy.log10(self.mat[sel][order]+1e-20)
        else: mat = self.mat[sel][order]
        if forceAll: fsel = fsel[order]
        if vmin is None or vmax is None: 
            if forceAll: fmat = numpy.sort(self.mat[fsel].flatten())
            else: fmat = numpy.sort(mat.flatten())
            # fmat = fmat[fmat > 0]
            ntot = len(fmat)
        if vmin is None: vmin = fmat[int(ntot*0.02)]
        if vmax is None: vmax = fmat[int(ntot*0.98)]
        sep = [0]
        z = numpy.zeros(self.nfreq)
        j = 0; yt = []; ytl = []
        yt.append(0)
        ytl.append('')
        while j < len(tmp):
            i = tmp[j]
            ini = j
            while tmp[j] == i: 
                j += 1
                if j == len(tmp): break
            assert ini != j 
            yt.append((j+ini)//2)
            ytl.append(str(i))
            sep.append(j)
            if separators:
                if j < len(tmp):
                    mat = numpy.vstack([mat[:j],z,mat[j:]])
                    tmp = numpy.hstack([tmp[:j],-1,tmp[j:]])
                    if forceAll:
                        fsel = numpy.hstack([fsel[:j], True, fsel[j:]])
            j += 1
        yt.append(j)
        ytl.append('')
            
        if forceAll: 
            mask = numpy.ones(mat.shape, dtype = bool)
            mask[fsel] = False
            mmat = numpy.ma.array(mat, mask = mask)
            m = pylab.matshow(mmat, **kargs)
        else:
            m = pylab.matshow(mat, **kargs)
        b = pylab.colorbar(shrink=0.8)
        if logy:
            b.set_label("log10("+units+"$^2$/Hz)")
        else:
            b.set_label(units+"$^2$/Hz")
        if not(hideYLabel): pylab.ylabel(ylabel)
        pylab.xlabel("Frequency [Hz]")
        if self.logx:
            m.axes.set_xticks(xt)
            m.axes.set_xticklabels(xtl)
        else:
            ini = self.matfreqs[0]-numpy.mod(self.matfreqs[0], linTickSep)+linTickSep
            end = self.matfreqs[-1]-numpy.mod(self.matfreqs[-1], linTickSep)
            f = numpy.linspace(ini, end, (end-ini)/linTickSep + 1)
            step = (self.matfreqs[-1] - self.matfreqs[0]) / self.nfreq
            xt = (f-self.matfreqs[0])/step
            xtl = numpy.array(f, dtype = str)
            m.axes.set_xticks(xt)
            m.axes.set_xticklabels(xtl)
        rat = self.nfreq / len(tmp) * ratio
        m.axes.xaxis.set_ticks_position("bottom")
        m.axes.set_yticks(yt)
        if hideYLabel: m.axes.set_yticklabels([])
        else: m.axes.set_yticklabels(ytl)
        if separators:
            for pos in sep:
                pylab.axhline(y=pos, color='black', linewidth=1)
        m.axes.set_aspect(rat)
        m.figure.set_size_inches(size[0], size[1], forward = True)
        m.set_clim(vmin = vmin, vmax = vmax)
        if "cmap" in kargs: cmap = kargs["cmap"]
        else: cmap = pylab.cm.RdYlBu_r
        cmap.set_bad([0.3, 0.3, 0.3],1.)
        m.set_cmap(cmap)
        if title is None: 
            title = "Watefall TOD %s %s " % \
                (self.name.split('.')[0], self.array)
            if rowDominance: title += "(Row Dominated)"
            else: title += "(Column Dominated)"
        pylab.title(title)
        if filename is not None: pylab.savefig(filename, dpi = dpi)
        if show: pylab.show()
        else: pylab.close()
                    


    def cumQualplot( self, filename = None, forceNew = False, nbins = 50, show = True,
                     selection = None,
                     title = None, f0 = 10., f1 = 200., Nmin = 30):
        """
        """
        self.analyze(fmin = f0, fmax = f1, logx = False)
        self.qual = numpy.zeros(len(self.mat))
        for i in range(len(self.mat)):
            p = self.mat[i][(self.matfreqs > f0)*(self.matfreqs < f1)]
            p -= p.mean()
            self.qual[i] = numpy.sqrt(p.std()/2.0/self.sampleTime)
        if selection is not None:
            self.qual = self.qual[selection]
        order = numpy.argsort(self.qual)
        self.qualDets = self.dets[order]
        self.qual = self.qual[order]
        self.qualDets = self.qualDets[self.qual != 0.0]
        self.qual = self.qual[self.qual != 0.0]
        cum = numpy.flipud(numpy.arange(len(self.qual))[Nmin:])
        q = self.qual[:-Nmin]
        pylab.plot(q, cum)
        pylab.xlabel('Noise Quality')
        pylab.ylabel('Number of Detectors')
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if show: pylab.show()
        else: pylab.close()
        return q, cum


    def plotArray( self, vmin = None, vmax = None, title = None, filename = None,
                   selection = None,
                   units = 'DAC', f0 = 10., f1 = 200., forceNew = False, show = True):
        """
        @brief Plot the quality of the power spectrum across the array in a 2D plot.
        The quality is defined as the variance of the power spectrum between 2 
        specified frequencies (default 10 and 200 Hz).
        @param vmin         Minimum value in scale of the plot.
        @param vmax         Maximum value in scale of the plot.
        @param title        Title in plot.
        @param filename     Filename where to save the plot.
        @param units        Units of data in TOD (DAC or uK) for label.
        @param f0           Minimum frequency in quality calculation.
        @param f1           Maximum frequency in quality calculation.
        @param forceNew     Whether to recalculate the quality or not.
        @param show         Whether to show or not the plot on screen.
        """
        if show: pylab.ion()
        else: pylab.ioff()
        if forceNew or self.mat is None:
            self.analyze(fmin = f0, fmax = f1, logx = False)
        if forceNew or self.qual is None:
            self.qual = numpy.zeros(len(self.mat))
            for i in range(len(self.mat)):
                p = self.mat[i][(self.matfreqs > f0)*(self.matfreqs < f1)]
                p -= p.mean()
                self.qual[i] = numpy.sqrt(p.std()/2.0)
        if selection is None:
            self.arrayQual = self.qual.reshape([self.Nrows,self.Ncol])
        else:
            q = self.qual.copy()
            q[~selection] = 0.0
            self.arrayQual = q.reshape([self.Nrows,self.Ncols])
        if vmin is None or vmax is None:
            vals = self.arrayQual.flatten()
            vals = numpy.sort(vals[vals != 0.0])
        if vmin is None: vmin = vals[int(len(vals)*0.02)]
        if vmax is None: vmax = vals[int(len(vals)*0.98)]
        m = pylab.matshow(self.arrayQual.transpose())
        m.set_clim(vmin = vmin, vmax = vmax)
        m.axes.xaxis.set_ticks_position("bottom")
        b = pylab.colorbar(shrink=0.8)
        b.formatter.set_powerlimits([-2,-2])
        #b.set_label('%s*rtsec'%units)
        b.draw_all()
        pylab.xlabel("rows")
        pylab.ylabel("cols")
        if title is None: 
            title = "Noise Quality TOD %s %s " % \
            (self.name.split('.')[0], self.array)
        pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if show: pylab.show()
        else: pylab.close()
        


class timeSpaceWaterfall( object ):
    """
    @brief Class object intended to visualize a TOD in time space as a waterfall plot
    """
    def __init__( self, tod, ntime = 1000, tmin = None, tmax = None):
        """
        @brief Initialization function for the timeSpaceWaterfall class object
        @param tod  TOD for which to produce the waterfall plot
        @param ntime  Number of time points in waterfall plot.
        @param tmin   Minimum time to consider in waterfall
        @param tmax   Maximum time to consider in waterfall
        """
        DT = (tod.ctime[-1]-tod.ctime[0])
        self.sampleTime = DT/(tod.nsamps-1)
        if tmin is None or tmin < 0.0: tmin = 0.0
        if tmax is None or tmax > DT: tmax = DT
        self.times = numpy.linspace(tmin,tmax, ntime)
        self.ntime = ntime
        self.ndet, self.ndata = tod.data.shape
        self.resultsDict = {}
        self.rows = tod.info.array_data["row"]
        self.cols = tod.info.array_data["col"]
        self.dets = tod.info.det_uid
        self.mat = time_space_waterfall(tod.data, self.times, self.sampleTime)
        self.sort = []
        for i in range(len(self.rows)):
            self.sort.append((self.rows[i], self.cols[i]))
        self.sort = numpy.array(self.sort, dtype = [('rows', int),('cols', int)])


    def plot( self, selection = None, vmin = None, vmax = None, level = 0.95, units = 'DAC', 
              title = None, rowDominance = False, separators = True, filename = None, 
              show = True):
        """
        @brief Plot function to visualize the time space waterfall plot.
        @param selection    Bool array with selection of detectors to show in plot.
        @param vmin         Minimum value in scale of waterfall plot
        @param vmax         Maximum value in scale of waterfall plot
        @param level        Fraction of values to consider in scale range (1 => min-max).
        @param units        Units of the TOD (DAC or uK).
        @param title        Title to add to the plot.
        @param rowDominance Whether to sort the waterfall by rows or columns.
        @param filename     Name of the file where to store the plot.
        @param show         Whether to display or not the plot.
        """
        if show: pylab.ion()
        else: pylab.ioff()
        if selection is None:
            sel = numpy.ones(numpy.shape(self.mat)[0], dtype = 'bool')
        else: sel = selection
        if vmin is None or vmax is None: 
            val = numpy.sort(numpy.reshape(self.mat[sel], numpy.size(self.mat[sel])))
            N = len(val)-1
        if vmin is None: vmin = val[int(N*(1.0-level)/2.0)]
        if vmax is None: vmax = val[int(N*(level+1.0)/2.0)]
        if rowDominance:
            order = numpy.argsort(self.sort[sel], order = ['rows', 'cols'])
            tmp = self.rows[sel][order]
            ylabel = "Row"
        else:
            order = numpy.argsort(self.sort[sel], order = ['cols', 'rows'])
            tmp = self.cols[sel][order]
            ylabel = "Column"
        mat = self.mat[sel][order]

        sep = [0]
        z = numpy.zeros(numpy.shape(mat)[1])
        j = 0; yt = []; ytl = []
        yt.append(0)
        ytl.append('')
        while j < len(tmp):
            i = tmp[j]
            ini = j
            while tmp[j] == i: 
                j += 1
                if j == len(tmp): break
            assert ini != j 
            yt.append((j+ini)//2)
            ytl.append(str(i))
            sep.append(j)
            if j < len(tmp):
                mat = numpy.vstack([mat[:j],z,mat[j:]])
                tmp = numpy.hstack([tmp[:j],-1,tmp[j:]])
            j += 1
        yt.append(j)
        ytl.append('')

        m = pylab.matshow(mat)
        b = pylab.colorbar(shrink=0.8)
        b.set_label(units)
        pylab.ylabel(ylabel)
        pylab.xlabel("Time [s]")

        shape = numpy.shape(self.mat[sel])
        rat = shape[1] / shape[0]*1.2
        m.axes.set_aspect(rat)

        xt = m.axes.get_xticks(); xtl = []
        st = numpy.mean(self.times[1:]-self.times[:-1])
        ti = self.times[0]
        for x in xt:
            xtl.append("%12.3f"%(x*st+ti))
        m.axes.set_xticklabels(xtl)
        m.axes.xaxis.set_ticks_position("bottom")
        m.axes.set_yticks(yt)
        m.axes.set_yticklabels(ytl)
        if separators:
            for pos in sep:
                pylab.axhline(y=pos, color='black', linewidth=1)
        m.figure.set_size_inches(10., 10.5, forward = True)
        m.set_clim(vmin = vmin, vmax = vmax)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if not(show): pylab.clf()


class scanWaterfall( object ):
    """
    @brief a waterfall plot of roughly azimuth angle versus time, stitching together the 
           common mode in pieces of left and right going scans.
    """
    def __init__( self, tod, selection = None ):
        """
        """
        if selection is not None: sel = selection
        else:
            a = numpy.min(tod.data, axis = 1)
            b = numpy.max(tod.data, axis = 1)
            sel = ~((a == 0.0)*(b == 0.0))

        cm = numpy.mean(tod.data[sel], axis = 0)
        self.cm = cm

        
        self.sampleTime = (tod.ctime[-1]-tod.ctime[0])/(tod.nsamps-1)
        T = int(1./tod.scanFreq/self.sampleTime/2)
        az = tod.az[:2*T]
        pivot = numpy.where(az == az.min())[0][0]

        i = pivot
        dir = 1
        k = 0
        self.time = [0]
        self.mat = numpy.zeros([int(tod.nsamps-pivot)//T,T])
        while i+T < tod.nsamps:
            if dir == 1: self.mat[k] = cm[i:i+T];
            else: self.mat[k] = numpy.flipud(cm[i:i+T]);
            self.time.append((self.time[-1]+T))
            i += T
            dir *= -1
            k += 1
        self.time = numpy.array(self.time)*self.sampleTime/60
        self.az = tod.az[pivot:pivot+T]*180/numpy.pi
        self.az -= self.az.mean()

            
    def plot( self, vmin = None, vmax = None, units = 'DAC', 
              title = None, filename = None, show = True):
        """
        @brief Plot function to visualize the time space waterfall plot.
        @param vmin         Minimum value in scale of waterfall plot
        @param vmax         Maximum value in scale of waterfall plot
        @param units        Units of the TOD (DAC or uK).
        @param title        Title to add to the plot.
        @param filename     Name of the file where to store the plot.
        @param show         Whether to display or not the plot.
        """
        if vmin is None: numpy.median(numpy.min(self.mat, axis = 1))
        if vmax is None: numpy.median(numpy.max(self.mat, axis = 1))

        m = pylab.matshow(self.mat)
        b = pylab.colorbar(shrink=0.8)
        b.set_label(units)
        pylab.ylabel("Time [min]")
        pylab.xlabel("dAz [deg]")

        shape = numpy.shape(self.mat)
        rat = shape[1] / shape[0] * 1.2
        m.axes.set_aspect(rat)

        m.axes.invert_yaxis()

        daz = (self.az.max()-self.az.min())/(len(self.az)-1)
        x_max = int(self.az.max())
        x = numpy.arange(2*x_max+1)-x_max
        xt = list((numpy.arange(2*x_max+1)-x_max)/daz + len(self.az)//2)
        xtl = list(numpy.array(x, dtype = 'str'))
        m.axes.set_xticks(xt)
        m.axes.set_xticklabels(xtl)
        m.axes.xaxis.set_ticks_position("bottom")
        y = numpy.arange(self.time[-1])
        yt = list(y/self.time[1])
        ytl = list(numpy.array(y, dtype = 'str'))
        m.axes.set_yticks(yt)
        m.axes.set_yticklabels(ytl)

        m.figure.set_size_inches(10., 10.5, forward = True)
        m.set_clim(vmin = vmin, vmax = vmax)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if show: pylab.show()
        else: pylab.clf()


class quality( object ):
    """
    @brief Object to cuantify the quality of the scan harmonics in the TOD.
    """
    def __init__( self, tod, f0 = 1.0, f1 = 200 ):
        """
        """
        self.name = tod.info.basename
        self.array = tod.info.array
        d, r, c = tod.listUncut()
        self.dets = numpy.array(d) 
        self.rows = numpy.array(r) 
        self.cols = numpy.array(c) 

        sel = self.dets[(self.rows>13)*(self.rows<17)*(self.cols>13)*(self.cols<17)]
        print(len(sel))
        
        f = numpy.zeros(len(sel))
        for i in range(len(sel)):
            p, nu, w = mobyUtils.power(tod.data[sel[i]], dt = tod.sampleTime)
            f[i] = tuneScanFreq(p, nu, tod.scanFreq)
            f[i] = tuneScanFreq(p, nu, f[i], scope = 0.0001)
        self.sf = numpy.median(f)
        self.f = f

        mask = generateArmonicMask(nu, self.sf, window = 6)
        sel = (nu > f0)*(nu < f1)

        print("Start arrayQual calculation")
        self.arrayQual = numpy.zeros([tod.ncol, tod.nrow])
        for i in range(tod.ndet):
            p, nu, w = mobyUtils.power(tod.data[i], dt = tod.sampleTime)
            mean1 = p[sel*mask].mean()
            mean2 = p[sel*~mask].mean()
            self.arrayQual[tod.cols[i]][tod.rows[i]] = mean1/mean2 - 1.0

    def plotQual( self, vmin = None, vmax = None, title = None, filename = None,
                   units = 'DAC', f0 = 10., f1 = 200., forceNew = False, show = True):
        """
        @brief Plot the quality of the power spectrum across the array in a 2D plot.
        The quality is defined as the variance of the power spectrum between 2 
        specified frequencies (default 10 and 200 Hz).
        @param vmin         Minimum value in scale of the plot.
        @param vmax         Maximum value in scale of the plot.
        @param title        Title in plot.
        @param filename     Filename where to save the plot.
        @param units        Units of data in TOD (DAC or uK) for label.
        @param f0           Minimum frequency in quality calculation.
        @param f1           Maximum frequency in quality calculation.
        @param forceNew     Whether to recalculate the quality or not.
        @param show         Whether to show or not the plot on screen.
        """
        if vmin is None or vmax is None:
            vals = self.arrayQual.flatten()
            vals = numpy.sort(vals[vals != 0.0])
        if vmin is None: vmin = vals[int(len(vals)*0.02)]
        if vmax is None: vmax = vals[int(len(vals)*0.98)]
        m = pylab.matshow(self.arrayQual)
        m.set_clim(vmin = vmin, vmax = vmax)
        m.axes.xaxis.set_ticks_position("bottom")
        b = pylab.colorbar(shrink=0.8)
        b.set_label(units+" rms")
        pylab.xlabel("rows")
        pylab.ylabel("cols")
        if title is None: 
            title = "Noise Quality TOD %s %s " % \
            (self.name.split('.')[0], self.array)
        pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if show: pylab.show()
        else: pylab.clf()
        




def generateArmonicMask(freqs, scanFreq, window = 10):
    """
    @brief Generates a mask that isolates those frequencies which are near a scan armonic.
    """
    w = window//2
    df = freqs[2]-freqs[1]
    index = numpy.where(numpy.mod(freqs,scanFreq) < df)[0]
    mask = numpy.zeros(len(freqs), dtype = 'bool')
    for i in index:
        if i-w < 0: mask[:i+w] = True
        elif i+w >= len(freqs): mask[i-w:] = True
        else: mask[i-w:i+w] = True
    return mask

def tuneScanFreq(p, nu, scanFreq, scope = 0.002, nsamp = 100, plot = False):
    """
    @brief Find the scan frequency by maximizing the harmonic content of a signal
    """
    df = nu[2] - nu[1]
    freqs = (numpy.arange(nsamp, dtype = 'float')/nsamp-0.5)*scope + scanFreq
    pow = numpy.zeros(nsamp)
    for i in range(nsamp):
        index = numpy.where(numpy.mod(nu,freqs[i]) < df)[0]
        pow[i] = numpy.mean(numpy.log(p[index]))
    if plot: pylab.plot(pow), pylab.show()
    mf = freqs[pow == pow.max()]
    if numpy.ndim(mf) > 0: return mf[0]
    else: return mf
        





def array_plots( param,
                 det=None, instrument = 'actpol', array = None, season = None,
                 tod=None, darks=True,
                 pmax=None, pmin=None, outrange=True,
                 param_name='', param_units='', title='',
                 display = 'show', save_name = 'newfig.png' ):
    """Plot a parameter across the array                                       
                                                                               
    Arguments:                                                                 
    ----------                                                                 
    param: parameter to plot                                                   
                                                                               
    Optional:                                                                  
    ---------                                                                  
                                                                               
    |det: list of detectors to plot (can be a list of det_uid or a tuple (row,col))    |instrument: instrument                                                    
    |array: array name ('ar1', 'pa2', 'pa3', 'pa4')                                     
    |season: observing season ('s13', 's14', 's15'...)                                
    or                                                                         
    |tod: tod object (tod.det, tod.info.array_name, tod.info.season)           
                                                                               
    pmax/pmin: range for parameter values                                      
    outrange: if True, parameters outside [pmin,pmax] will be plotted in black 
    param_name, param_units, title: informations to add to the plot            
    display: ['show', 'save'] either show the plot, or directly save it        
    save_name: if display == 'save', name of the output file                   
    """
    if det is None and array is None and tod is None:
        print("List of (dets, season and array name) or tod object must be prov\
ided")
        return 0

    det = np.asarray(det, dtype = int)
    param = np.asarray(param, dtype = float)

    if tod is not None:
        array_data = tod.info.array_data
        instrument = tod.info.instrument
        array = tod.info.array
        season = tod.info.season
    else:
        if instrument is None:
            print("Please provide the instrument if no TOD is provided")
            return 0
        if array is None:
            print("Please provide the array if no TOD is provided")
            return 0
        if season is None:
            print("Please provide the season if no TOD is provided")
            return 0
        array_data = products.get_array_data(
                                {"instrument":instrument,
                                 "array_name":array,
                                 "season":season})

    if len(det) == 2:
        Row, Col = det
        Detid = Row*32 + Col
    else:
        Detid = det

    pos, polfamily, freq = get_position( Detid, instrument, array, season )
    x, y = pos

    if pmin == None:
        pmin = param.min()
    if pmax == None:
        pmax = param.max()

    if param.min() == param.max():
        color = ['b']*param.size
    else:
        color = get_color( param, pmax = pmax, pmin = pmin, outrange = outrange )
    if np.unique(freq).size == 1:
        patchlist = get_patches( pos, color, polfamily, freq )
    else: 
        patchlist = get_patches( pos, color, polfamily, freq, radius=0.012 )


    # Add other detectors in grey
    det_uid = array_data['det_uid']
    if np.unique(freq).size == 1: 
        det_uid = det_uid[array_data['nom_freq']==np.unique(freq)]
    if Detid.size < det_uid.size:
        _dets = set( np.arange(det_uid.size) )
        _dets.difference_update( set( Detid ) )
        _dets = np.asarray( list(_dets) )
        _pos, _polfamily, _freq = get_position( _dets, instrument, array, season )
        _x, _y = _pos
        _color = np.array([[0.,0.,0.,0.3]]).repeat( _x.size, axis = 0)
        if np.unique(freq).size==1: _patchlist = get_patches( _pos, _color, _polfamily, _freq )
        else: _patchlist = get_patches( _pos, _color, _polfamily, _freq, radius=0.012 )
    else:
        _patchlist = []
    
    #Add dark detectors as black circles
    if darks:
        dd = array_data['det_uid'][array_data['det_type']=='dark_tes']
        x_dark, y_dark = array_data['sky_x'][dd], array_data['sky_y'][dd]
        for xd, yd in zip( x_dark, y_dark ):
            patchlist.append( patches.Wedge( [xd,yd], 0.02, 0., 360.,
                                             width=0.003, fc = 'k', ec='k' ) )

    
    x_lim, y_lim = get_array_plot_lims(array, season)

    plt.ioff()
    plt.figure( figsize = (10,10) )
    ax = create_plot( patchlist+_patchlist, pmin, pmax, x_lim, y_lim, array )
    set_infos( param_name, param_units, title, ax )
    
    if display == 'show':
        plt.ion()
        plt.show()
    elif display == 'save':
        plt.savefig(save_name)
        plt.close()


def get_position(Detid, instrument, array, season):
    """Return the position in the focal plane for a list of detectors
    """
    params = {
        'instrument' : instrument,
        'array_name' : array,
        'season' : season,
        }
    arraydata = moby2.scripting.get_array_data(params)
    polfamily = arraydata['pol_family']
    x = arraydata['sky_x']
    y = arraydata['sky_y']
    col = arraydata['col']
    row = arraydata['row']
    freq = arraydata['nom_freq']
    detid = row*32+col

    idx_sort = detid.argsort()
    x = x[idx_sort]
    y = y[idx_sort]
    polfamily = polfamily[idx_sort]
    freq = freq[idx_sort]
    return (x[Detid], y[Detid]), polfamily[Detid], freq[Detid]


def get_color(param, pmax=None, pmin=None, cmap = cm.RdYlBu_r, outrange = False):
    if pmax==None:
        pmax = param.max()
    if pmin==None:
        pmin = param.min()
    color = cmap( (param-pmin)/(pmax-pmin) )
    if outrange:
        color[param>pmax,:] = (0., 0., 0., 0.)
        color[param<pmin,:] = (0., 0., 0., 0.)
    return color


def get_patches(pos, color,polfamily, freq, radius=0.015):
    patchlist = []
    if len(np.unique(freq)) == 1:
        for x, y, c, pf in zip( pos[0], pos[1], color, polfamily ):
            if pf == 'A':
                theta1, theta2 = (90, 270)
            elif pf == 'B':
                theta1, theta2 = (270, 90)
            elif pf == 'X':
                theta1, theta2 = (0,360)
            patchlist.append( patches.Wedge( [x,y], radius, theta1, theta2,
                                             fc = c, ec=c ) )

    else:
        f1, f2 = np.unique(freq)[-2:]
        for x, y, c, pf, f in zip( pos[0], pos[1], color, polfamily, freq ):
            if pf == 'A':
                if f == f1:
                    theta1, theta2 = (90, 180)
                elif f == f2:
                    theta1, theta2 = (180, 270)
            elif pf == 'B':
                if f == f1:
                    theta1, theta2 = (0, 90)
                elif f == f2:
                    theta1, theta2 = (270, 360)
            elif pf == 'X':
                theta1, theta2 = (0,360)
            patchlist.append( patches.Wedge( [x,y], radius*1.2, theta1, theta2,
                                             fc = c, ec=c ) )
    
    return patchlist


def create_plot(patchlist, pmin, pmax, x_lim, y_lim, array_name):
    if pmin != pmax:
        ax1 = plt.subplot2grid((1,10),(0,0),colspan=9,aspect='equal')
        ax2 = plt.subplot2grid((1,10),(0,9))
    else:
        ax1 = plt.axes()

    for p in patchlist:
        ax1.add_patch( p )
    ax1.set_xlim( x_lim )
    ax1.set_ylim( y_lim )
    ax1.tick_params( bottom = 'off', top='off', right='off', left='off', labelbottom='off', labelleft='off', which = 'both' )

    plt.text(0.05, 0.05, 'Detector not available',
             color = 'grey', transform = ax1.transAxes)
    plt.text(0.05, 0.025, 'Detector out of range',
             color = 'black', transform = ax1.transAxes)
    plot_wafer_names(ax1, array_name)
    # PA = {'ar1':'PA1',
    #       'ar2':'PA2',
    #       'ar3':'PA3',
    #       'pa4':'PA4',
    #       "ar3_90":"PA3",
    #       "ar3_150":"PA3",
    #       }
    plt.text(0.9, 0.05, array_name,
             fontsize='xx-large',
             ha='center', va='center', transform = ax1.transAxes)
    
    # if array_name in ['ar3', 'ar5', 'ar6']:
    #     plt.text(0.1,0.95,'150', ha='center', va='center',
    #              transform=ax1.transAxes)
    #     plt.text(0.1,0.925,'90', ha='center', va='center',
    #              transform=ax1.transAxes)
    #     ax1.add_patch(
    #         patches.Wedge((0.1,0.94),0.03,0,180,transform=ax1.transAxes,edgecolor='k',fc='none')
    #         )
    #     ax1.add_patch(
    #         patches.Wedge((0.1,0.94),0.03,180,360,transform=ax1.transAxes,edgecolor='k',fc='none')
    #         )
    # elif array_name == 'ar4':
    #     plt.text(0.1,0.95,'220', ha='center', va='center',
    #              transform=ax1.transAxes)
    #     plt.text(0.1,0.925,'150', ha='center', va='center',
    #              transform=ax1.transAxes)
    #     ax1.add_patch(
    #         patches.Wedge((0.1,0.94),0.03,0,180,transform=ax1.transAxes,ec='k',fc='none')
    #         )
    #     ax1.add_patch(
    #         patches.Wedge((0.1,0.94),0.03,180,360,transform=ax1.transAxes,ec='k',fc='none')
    #         )

    if pmin != pmax:
        norm = mpl.colors.Normalize(vmin=pmin,vmax=pmax)
        cb1 = mpl.colorbar.ColorbarBase(ax2,cmap=cm.RdYlBu_r,norm=norm,orientation='vertical')
        return ax1, ax2
    else:
        return ax1


def set_infos( param_name, units, title, ax ):
    if type(ax) == tuple:
        ax1, ax2 = ax
        ax1.set_title( title, fontsize = 15 )
        ax2.set_ylabel( '%s [%s]' %(param_name, units), rotation = 270, fontsize = 20 )
    else:
        ax.set_title( title, fontsize = 15 )




def tod3D(tod, dets=None, time_resolution=2., prange=[None,None],
          sky_coords=False, pointingpar = None,
          anim_time = 10., display='show', filename=None, **kwargs):
    """3D visualization of a TOD: animation of the 2D focal plane through the TOD

    Arguments:
    ----------
    - tod: TOD object to visualize

    Optional:
    ---------
    - dets: list of detectors to plot (default will use tod.det_uid)
    - time_resolution: Tod will be downsample to reach a time resolution
                       as close as possible to the input using a power of
                       2 resampling (every 2**n samples).
                       Default is 2s
    - prange: scale
    - sky_coords: if True, project the focal plane on the sky and show the TOD animation                  in RA, DEC (must provide pointingpar)
    - pointingpar: focal plane model, eg. {'source': 'fp_file',
                   'filename': '.../RelativeOffsets/template_ar2_150529s.txt'}
    - anim_time: total time of the animated TOD in seconds
    - display: 'show' or 'save'
    - filename: only if display=='save', by default: todname.gif
    - other args to define the animation (interval, repeat, repeat_delay...)
    """

    if dets is None:
        dets = tod.det_uid


    # Re-sampling
    tres = (tod.ctime[-1] - tod.ctime[0]) / tod.nsamps
    r = time_resolution / tres
    Nresamp = 2**np.ceil(np.log2(r))
    tod_ds = tod.copy(resample=Nresamp)
    print("Downsampled tod has %i samples, with time resolution of %.2fs" %(
        tod_ds.nsamps, (tod_ds.ctime[-1] - tod_ds.ctime[0]) / tod_ds.nsamps ))


    # Get focal plane infos
    if sky_coords:
        x, y = get_sky_coords(tod_ds, pointingpar)
        x = x[dets]
        y = y[dets]
        x *= 180. / np.pi
        y *= 180. / np.pi
        if x.max() - x.min() > 180.: x[x<0] += 360.
        if y.max() - y.min() > 180.: y[y<0] += 360.
    else:
        x = tod.info.array_data['sky_x'][dets]
        y = tod.info.array_data['sky_y'][dets]
    polfamily = tod.info.array_data['pol_family'][dets]



    # Get limits fot the plot
    if sky_coords:
        x_lim = (x.min(), x.max())
        y_lim = (y.min(), y.max())
    else:
        x_lim, y_lim = get_array_plot_lims(
            tod.info.array, tod.info.season)

    pmin, pmax = prange
    if pmin == None:
        pmin = tod.data[dets].min()
    if pmax == None:
        pmax = tod.data[dets].max()

    


    plt.ioff()
    # Create animation
    fig = plt.figure(figsize=(7,8))
    ax1 = fig.add_axes([0.1,0.2,0.8,0.7])
    if sky_coords:
        ax1.tick_params(
        bottom = 'off', top='on', right='off', left='on',
        labelbottom='off', labelleft='on', labeltop='on',
        which = 'both' )
    else:
        ax1.tick_params(
        bottom = 'off', top='off', right='off', left='off',
        labelbottom='off', labelleft='off',
        which = 'both' )
    ax2 = fig.add_axes([0.1,0.1,0.8,0.1])
    ax2.plot(tod_ds.ctime - tod_ds.ctime[0], tod_ds.data[dets].T, 'b', alpha=0.1)
    ax1.set_xlim( x_lim )
    ax1.set_ylim( y_lim )
    ax2.set_xlim((0,tod_ds.ctime[-1]-tod_ds.ctime[0]))
    ax2.set_ylim((pmin,pmax))
    ax2.set_yticks((pmin,0,pmax))
    if sky_coords:
        ax1.set_xlabel('RA (degrees)')
        ax1.set_ylabel('DEC (degrees)')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('pW')
    plt.figtext(0.5, 0.95, '%s' %tod_ds.info.basename,
                ha = 'center', va = 'center',
                fontsize='xx-large')
    if not sky_coords:
        plot_wafer_names(ax1, tod_ds.info.array)

    line = ax2.axvline(0, color='r')
    color = get_color( tod_ds.data[dets,0],
                              pmax = pmax, pmin = pmin )
    colors = [ get_color( tod_ds.data[dets,i],
                          pmax = pmax, pmin = pmin ) 
               for i in range(tod_ds.nsamps) ]
    if sky_coords:
        # pos = x[:,0], y[:,0]
        patchlists = [ PatchCollection(
                get_patches( (x[:,i],y[:,i]), colors[i], polfamily ) )
                       for i in range(tod_ds.nsamps) ]
        ax1.add_collection(patchlists[0])
    else:
        pos = x, y
        patchlist = PatchCollection( get_patches( pos, colors[0], polfamily ) )
        print(patchlist.get_facecolor())
        # p = patchlist[0]
        # ax1.add_patch(p)
        ax1.add_collection( patchlist )
    
    
    def animate(i):
        if sky_coords:
            patchlist = patchlists[i]
            color = colors[i]
            patchlist.set_facecolors(color)
            patchlist.set_edgecolors(color)
            ax1.collections.pop()
            ax1.add_collection(patchlist)
        else:
            patchlist = ax1.collections[0]
            color = colors[i]
            patchlist.set_facecolors(color)
            patchlist.set_edgecolors(color)
        line.set_xdata(tod_ds.ctime[i] - tod_ds.ctime[0])
        return line

    interval = float(anim_time) / tod_ds.nsamps * 1000
    ani = animation.FuncAnimation(fig, animate,
                                  np.arange(tod_ds.nsamps),
                                  interval=interval,
                                  **kwargs)
    
    if display == 'show':
        plt.ion()
        plt.show()
        plt.draw()
    elif display == 'save':
        if filename is None:
            filename = '%s.gif' %tod_ds.info.basename
        Writer = animation.writers['imagemagick_file']
        writer = Writer(fps=1/interval*1000)
        ani.save(filename,writer=writer)
        plt.close()


        
def get_sky_coords(tod, pointingpar):
    tod.fplane = products.get_focal_plane(
        pointingpar, det_uid=tod.det_uid, tod_info=tod.info, tod=tod)
    ra, dec = moby2.pointing.get_coords(
        tod.ctime, tod.az, tod.alt, focal_plane=tod.fplane)
    return ra, dec


def get_array_plot_lims(array, season):
    if array == 'ar1':
        if season == '2013':
            xmin, xmax = -1.3094665029667607, -0.41013064910400088
            ymin, ymax = -1.220441626351696, -0.30605809090138458
        else:
            xmin, xmax = -1.3287624329466337, -0.42355235864871088
            ymin, ymax = -1.3491856371418196, -0.43444330849940416
    elif array == 'ar2':
        xmin, xmax = -0.048737122029174233, 0.85432442834243771
        ymin, ymax = -1.3379718650336261, -0.41940152960224164
    elif array == 'ar3_90' or array == 'ar3_150' or array == "ar3":
        xmin, xmax = -0.69483493785298767, 0.18037663135240153
        ymin, ymax = -0.15016921848482107,  0.80671585130610191
    elif array == 'ar4':
        if season == '2016':
            xmin, xmax = -1.3172299710057627, -0.57135351330445694
            ymin, ymax = -1.2802741932198245, -0.41247231671467965
        else:
            xmin, xmax = -0.98219689986339409, -0.23632044216208847
            ymin, ymax = -0.82135665065746755, 0.046445225847677393
    elif array == 'ar5':
        xmin, xmax = 0.2634513877020111, 1.0299899366437952
        ymin, ymax = -0.85103798177410184, 0.03837494139098363
    elif array == 'ar6':
        xmin, xmax = -0.34516751588843669, 0.39094596571401696
        ymin, ymax = 0.3670596582454384, 1.3154743040302765
    else:
        print("Array must be from the list ['ar1', 'pa2', 'ar3', 'ar4', 'ar5', 'ar6']")

    x_lim = ( xmin - (xmax-xmin)*0.1, xmax + (xmax-xmin)*0.1 )
    y_lim = ( ymin - (ymax-ymin)*0.1, ymax + (ymax-ymin)*0.1 )
    return x_lim, y_lim


def plot_wafer_names(ax, array):
    if array == 'ar1':
        plt.text(0.05, 0.2, 'W10', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.7, 0.1, 'SH2B', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.9, 0.35, 'W08', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.75, 0.85, 'SH1A', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.1, 0.85, 'W09', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.01, 0.45, 'SH2A', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
    elif array == 'ar2':
        plt.text(0.15, 0.09, 'FH3C', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.75, 0.15, 'SH4B', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.9, 0.4, 'FH6', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.6, 0.85, 'SH3B', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.1, 0.85, 'FHC1', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.01, 0.45, 'SH4A', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large',
                 rotation=90)
    elif 'ar3' in array:
        plt.text(0.5, 0.02, 'FH3', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large',
                 ha='center')
        plt.text(0.85, 0.25, 'SH1A', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.9, 0.7, 'FH4', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.5, 0.95, 'SH1B', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large',
                 ha='center')
        plt.text(0.05, 0.75, 'FH2', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')
        plt.text(0.03, 0.25, 'SH8B', color='gray',
                 transform = ax.transAxes, fontsize = 'x-large')



