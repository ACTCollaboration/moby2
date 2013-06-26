from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
This is a module of finding and subtracting the scan-synchronous
signal from the TOD data.

Classes:
Sync
"""

import os
import pickle
import numpy

import matplotlib

pylab = None
def import_pylab():
    global pylab
    import pylab as pl
    pylab = pl

import moby2
import moby2.util.log as psLib

def get_scan_info(tod=None, az=None, ctime=None):
    """
    Analyze az info (for tod, or pass in an az and ctime).

    Returns a struct with results.
    """
    if az is None:
        az = tod.az
    if ctime is None:
        ctime = tod.ctime
    output = moby2.libactpol.analyze_scan(
        numpy.asarray(az, dtype='float64'),
        numpy.asarray(ctime, dtype='float64'))
    return moby2.util.aStruct({
            'status': output[0],
            'az_speed' : output[1],
            'scan_freq' : output[2],
            'az_max' : output[3],
            'az_min' : output[4],
            'sections' : output[5],
            })


def repair_linear_vector(ctime):
    """
    Repair a vector that should be linear in index.  This fixes
    repeated samples, but does not replace missing samples.  Does not
    edit the vector in place.
    """
    dt = numpy.diff(ctime)
    dt0 = numpy.median(dt[::100])
    bad_step = abs(dt-dt0) > dt0*.3
    ok = numpy.ones(ctime.shape, 'bool')
    ok[1:] *= ~bad_step
    ok[:-1] *= ~bad_step
    t1 = ctime - ctime[0]
    p = numpy.polyfit(ok.nonzero()[0], t1[ok], 1)
    t1[~ok] = numpy.polyval(p, (~ok).nonzero()[0])
    return ctime[0] + t1

def outlier_fill(data, ok=None):
    """
    Remove extreme outliers, spline interpolating to fill the gaps.
    """
    import scipy.stats.mstats as ms
    from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
    if ok is None:
        lo, hi = ms.mquantiles(data, (0.03, 0.97))
        lo_bound = lo - (hi-lo)*.5
        hi_bound = hi + (hi-lo)*.5
        ok = (lo_bound <= data) * (data <= hi_bound)
    ## Interpolate
    x = numpy.arange(len(data))
    filler = spline1d(x[ok], data[ok])
    data[~ok] = filler(x[~ok])
    return data

def repair_pointing(tod):
    tod.ctime = repair_linear_vector(tod.ctime)
    tod.alt = outlier_fill(tod.alt)
    tod.az = outlier_fill(tod.az)


class Sync( object ):
    """
    This object represents the synchronous pickup analysis done on a TOD.

    Basic Attributes:
    ndet    - scalar, number of detectors in object.
    row     - ndet, row of each detector in object.
    col     - ndet, column of each detector in object.

    Synchronous pickup data:
    amp     - ndet, amplitude of the sinusoidal for each detector.
    phase   - ndet, phase of the sinusoidal for each detector.
    rms     - ndet, rms after removing the sinusoidal for each detector.
    phaseAz - scalar, phase of the azimuth scan.

    Analysis data:
    mask    - ndet, mask which contains information about outliers. It can take the values:
               0    Un masked detector.
              -1    Amplitude lower than "median * thrLo".
               1    Amplitude higher than "median * thrHi".
               2    Phase at less than "thrPh" from "median - 180 deg".
               3    Extended pixel (see extend function).
    params  - dictionary, parameters and logging information about the object. Items:
              todName  - name of the TOD used in the analysis.
              rows     - rows used to initialize the object.
              cols     - columns used to initialize the object.
              order    - order of the polinomial used to remove the DC drift before the fit.
              samples  - number of samples to take from TOD to fit the polinomial.
              nPeriods - number of scan periods to consider in sinusoidal fit.
              removals - detector number, amplitude and phase of sinusoidal removal.
              outlierThr - thrHi, thrLo and thrPh used for finding outliers.
    """

    @staticmethod
    def getDefaultInitParams(array):
        """
        @brief    Returns default initialization parameters from moby/params/syncDefault.par
        """
        return moby2.scripting.products.get_sync_params({'season': '2008', 'array_name': array})

    def __init__( self, tod, params = None, all = False, noFit = False, scanFreq=None):
        """
        @brief    Initialization of the sync object. Performs a fit to the synchronous signal of the given TOD.

        @param    tod    TOD.TOD object.
        @param    params dictionary containing the parameters for initialization.
        """
        import_pylab()
        self.tod = tod
        if params is None:
            params = self.getDefaultInitParams(tod.info.array)
        if all == True:
            # Get list of all det_uid for this array.
            det_uid = tod.info.array_data.select_outer()
            # But make sure we only actually process the ones in tod.
            dets = numpy.arange(tod.data.shape[0])
        else:
            # Get list of det_uid in this TOD
            det_uid = tod.info.det_uid
            # Get index of those matching our row/col selection
            dets = tod.info.array_data.select_outer({'row': params['rows'],
                                                     'col': params['cols']},
                                                    det_uid=det_uid)
        
        params['rows'] = tod.info.array_data['row'][det_uid]
        params['cols'] = tod.info.array_data['col'][det_uid]
       
        assert len(params["rows"]) > 0

        if scanFreq is None:
            scanFreq = get_scan_info(tod).scan_freq
        self.scanFreq = scanFreq

        tod_dt = tod.ctime - tod.ctime[0]
        self.sampleTime = tod_dt[-1] / (tod_dt.shape[0]-1)

        # Init empty structure
        z = numpy.zeros(len(det_uid), 'float32')
        self.sync_fit = moby2.util.aStruct({
                'amp': z.copy(), 'phase': z.copy(), 'rms': z.copy(),
                'phaseAz': 0.})
        if noFit == False:
            az = tod.az - tod.az.mean()

            sync_fit = moby2.libactpol.get_sync_amps(tod.data, dets.astype('int32'),
                                                     az, tod_dt,
                                                     self.scanFreq, params['order'],
                                                     params['samples'], params['nPeriods'])
            self.sync_fit.amp, self.sync_fit.phase, self.sync_fit.rms, \
                self.sync_fit.phaseAz = sync_fit
        else:
            pass

        self.ndet = len(det_uid)
        self.row, self.col = params['rows'], params['cols']
        self.det_uid = det_uid
        self.amp = self.sync_fit.amp
        med = numpy.median(self.sync_fit.phase)
        self.phase = moby2.util.angles.rerange(self.sync_fit.phase - med,
                                               -numpy.pi, 2*numpy.pi) + med
        self.rms = self.sync_fit.rms
        self.phaseAz = self.sync_fit.phaseAz

        self.params = params
        self.params["removals"] = [[],[],[]] # [det#, amp, phase]
        self.params["outlierThr"] = [[],[],[]] # [thrHi, thrLo, thrPh]
        self.params["todName"] = tod.info.name
        self.mask = numpy.array([0]*self.ndet)
        self.mask[numpy.isnan(self.amp)] = 1
        self.mask[numpy.isnan(self.phase)] = 1
        self.cos = numpy.cos(tod_dt*2*numpy.pi*self.scanFreq)
        self.sin = numpy.sin(tod_dt*2*numpy.pi*self.scanFreq)
        self.badcols = []

    def extend( self, delTOD = True ):
        """
        @brief   Extends the mbSync structure to cover all pixels in detector array.
        
        @return  ssEx  Synchronous pickup object defined for all detectors (extended).
        """

        ssEx = Sync(self.tod, all = True, noFit = True)

        ssEx.mask = numpy.array([3]*ssEx.ndet)
        ssEx.phaseAz = self.phaseAz
        for i in range(self.ndet):
            ssEx.amp[(ssEx.row == self.row[i])*(ssEx.col == self.col[i])] = \
                                                                    self.amp[i]
            ssEx.phase[(ssEx.row == self.row[i])*(ssEx.col == self.col[i])] = \
                                                                  self.phase[i]
            ssEx.rms[(ssEx.row == self.row[i])*(ssEx.col == self.col[i])] = \
                                                                    self.rms[i]
            ssEx.mask[(ssEx.row == self.row[i])*(ssEx.col == self.col[i])] = \
                                                                   self.mask[i]
        amp, phase, rms = ssEx.fillBlanks()
        ssEx.amp[:] = amp[:]
        ssEx.phase[:] = phase[:]
        ssEx.rms[:] = rms[:]

        zeroSel = (numpy.min(self.tod.data, axis = 1) == 0.) * \
                  (numpy.max(self.tod.data, axis = 1) == 0.)
        ssEx.mask[zeroSel] = -1

        ssEx.params["removals"] = self.params["removals"]
        ssEx.params["outlierThr"] = self.params["outlierThr"]
        
        if delTOD: self.tod = None

        return ssEx

    ## This function not updated for moby2!
    def removeSingle( self, row, col, amp, phase, plot = False):
        """
        @brief    Remove sinusoid of specified amplitude and phase from a give pixel data.

        @param  row    Specify row of detector where to remove the sinusoid.
        @param  col    Specify column of detector where to remove the sinusoid.
        @param  amp    Amplitude of the sinusoid to remove.
        @param  phase  Phase of the sinusoid to remove.
        @param  plot   Produce a plot of the TOD before and after the removal. Defaul "False".
        """

        if self.tod.data[self.tod.dets[row][col]][0] != actUtils.NO_VALUE:
            sync = amp * (numpy.cos(phase)*self.cos + numpy.sin(phase)*self.sin)
         
            if plot == True:
                filt = mbFilter.lowPassButterworth(self.tod, fc = 0.5, order = 2)
                data = mbFilter.pyFilter(self.tod.data[self.tod.dets[row][col]], filt)
         
            self.tod.data[self.tod.dets[row][col]][:] -= sync
         
            if plot == True:
                data2 = mbFilter.pyFilter(self.tod.data[self.tod.dets[row][col]], filt)
                lines = []
                lines += pylab.plot(self.tod.dt, data, 'b')
                lines += pylab.plot(self.tod.dt, data2, 'g')
                fo = matplotlib.font_manager.FontProperties(size=10)
                pylab.legend(lines, ['Before', 'After'], prop = fo)
                pylab.show()

            self.params["removals"][0].append(self.tod.dets[row][col])
            self.params["removals"][1].append(amp)
            self.params["removals"][2].append(phase)
            self.tod.abuses += [{'name':'synchronous removal', 'type':'signal subtraction', \
                 'row':row, 'col':col, 'amp':amp, 'phase':phase}]
        else:
            pslib.trace("moby", 2, "ERROR: No TOD data for removing sync signal for ROW " \
                             + str(row) + " COL " + str(col))

    def removeAll( self, useMask = True):
        """
        @brief  Remove syncronous signal from all detector in mbSync structure.

        @param  mask  Use mask to select only unmasked detectors. Default "True".
        """
        if useMask:
            mask = (self.mask==0)+(self.mask==3)
        else:
            mask = numpy.ones(self.ndet, 'bool')

        det_idx = mask[self.tod.info.det_uid].nonzero()[0].astype('int32')
        modes = numpy.array([self.cos, self.sin],
                            dtype='float32')
        amps = numpy.array([self.amp[mask] * numpy.cos(self.phase[mask]),
                            self.amp[mask] * numpy.sin(self.phase[mask])],
                            dtype='float64')
        moby2.libactpol.remove_modes(self.tod.data, det_idx, modes, amps)
        
                

    def findOutliers( self , fracAmp = None, thrldPhase = None):
        """
        @brief   Find outliers within the analyzed detectors.

        @param   fracAmp    Only accept amplitude values such that:
                            (1-fracAmp)*mAmp < Amp < (1+fracAmp)*mAmp
                            where mAmp = median(Amps)
        @param   thrldPhase  Only accept phase values such that:
                            mPhase-thrldPhase < Phase < mPhase+thrldPhase
                            where mPhase = median(Phases)
        """
        params = self.getDefaultInitParams(self.tod.info.array)
        if fracAmp is None:
            fracAmp = params['fracAmp']
        if thrldPhase is None:
            thrldPhase = params['thrldPhase']

        self.mask = numpy.zeros(self.ndet)
        zeroSel = (self.amp == 0.0)
        self.mask[zeroSel] = -2

        m_amp = numpy.median(self.amp[self.mask == 0])
        m_phase = numpy.median(self.phase[self.mask == 0])

        self.mask[self.amp > 1e5] = 1
        for c in range(32):
            sel = (self.col == c)*~zeroSel
            if len(self.col[sel]) > 2:
                m = numpy.median(self.amp[sel])
                self.mask[sel*(self.amp > m*(1+fracAmp))] = 1
                self.mask[sel*(self.amp < m*(1-fracAmp))] = 1
            elif (len(self.col[sel]) == 2) and \
                 ((self.amp[sel][0] > self.amp[sel][1]*(1+fracAmp)) or \
                  (self.amp[sel][0] < self.amp[sel][1]*(1-fracAmp))):
                if numpy.abs(self.amp[sel][0] - m_amp) > numpy.abs(self.amp[sel][1] - m_amp):
                    self.mask[numpy.where(sel)[0][0]] = 1
                else: self.mask[numpy.where(sel)[0][1]] = 1
            sel *= (self.mask == 0)
            if len(self.col[sel]) > 2:
                m = numpy.median(self.phase[sel])
                self.mask[sel*(self.phase > m+thrldPhase)] = 2
                self.mask[sel*(self.phase < m-thrldPhase)] = 2
            elif (len(self.col[sel]) == 2) and \
                 ((self.phase[sel][0] > self.phase[sel][1]+thrldPhase) or \
                  (self.phase[sel][0] < self.phase[sel][1]-thrldPhase)):
                if numpy.abs(self.phase[sel][0] - m_phase) > numpy.abs(self.phase[sel][1] - m_phase):
                    self.mask[numpy.where(sel)[0][0]] = 2
                else: self.mask[numpy.where(sel)[0][1]] = 2
            sel *= (self.mask == 0)
            if not(numpy.any(self.mask[self.col == c] == 0)):
                self.badcols.append(c)
                psLib.trace("moby", 2, "WARNING: No amplitudes left for COL " + str(c))


    def fillBlanks( self, homogenize = False ):
        """
        @brief   Puts the median of the elements of the rest of the column to the value of the 
        pixels that were labeled as outliers.

        @param  homogenize  Uses the same value to all detectors in a colunm, using the mean of the
        unmasked values for that column.
        """

        amp = self.amp - 0.0
        phase = self.phase - 0.0
        rms = self.rms - 0.0
        if homogenize:
            for ind in range(32):
                if numpy.any(self.col[self.mask == 0] == ind):
                    amp[self.col == ind]  = \
                        numpy.median(self.amp[(self.col == ind)*(self.mask == 0)])
                    phase[self.col == ind]  = \
                        numpy.median(self.phase[(self.col == ind)*(self.mask == 0)])
                    rms[self.col == ind]  = \
                        numpy.median(self.rms[(self.col == ind)*(self.mask == 0)])
                else:
                    amp[self.col == ind]  = \
                        numpy.median(self.amp[self.mask == 0])
                    phase[self.col == ind]  = \
                        numpy.median(self.phase[self.mask == 0])
                    rms[self.col == ind]  = \
                        numpy.median(self.rms[self.mask == 0])
                    psLib.trace("moby", 3, "WARNING: Missing information for COL" + str(ind) + \
                                         ", using median of whole array.")
        else: 
            #dets, rows, cols = self.tod.listUncut()
            ## This is admittedly somewhat more work...
            dets = self.tod.cuts.get_uncut()
            cols = self.tod.info.array_data['col'][self.tod.det_uid][dets]

            for ind in range(0, max(cols)+1):
                if numpy.any(self.col[self.mask == 0] == ind):
                    amp[(self.col == ind)*(self.mask != 0)]  = \
                        numpy.median(self.amp[(self.col == ind)*(self.mask == 0)])
                    phase[(self.col == ind)*(self.mask != 0)]  = \
                        numpy.median(self.phase[(self.col == ind)*(self.mask == 0)])
                    rms[(self.col == ind)*(self.mask != 0)]  = \
                        numpy.median(self.rms[(self.col == ind)*(self.mask == 0)])
                else:
                    sel = (numpy.array(cols) == ind)
                    if numpy.any(sel):
                        clm = numpy.empty(self.tod.ndata)
                        mobyLib.transposeMedian(self.tod.ctod, list(numpy.array(dets)[sel]), clm)
                        tamp, tphase = fitSine(clm, 2*numpy.pi*self.scanFreq*self.sampleTime)
                        amp[(self.col == ind)*(self.mask != 0)] = tamp
                        phase[(self.col == ind)*(self.mask != 0)] = tphase
                        rms[(self.col == ind)*(self.mask != 0)] = 0.0
                        psLib.trace("moby", 1, "WARNING: Missing information for COL %d,"
                                " fitting from live data" % ind ) 

        return amp, phase, rms

    def combine( self, ss2, badColsOnly = True ):
        """
        @brief  Combines the results from another synchronous object with the current results.
        @param  ss2  synchronous object with data to combine.
        """
        assert self.ndet == ss2.ndet

        if badColsOnly:
            sel = numpy.zeros(len(self.amp), dtype = bool)
            for c in self.badcols:
                sel += (numpy.array(self.col) == c)
        else: 
            sel = numpy.ones(len(self.amp), dtype = bool)

        A = self.amp[sel]*numpy.cos(self.phase[sel]) + ss2.amp[sel]*numpy.cos(ss2.phase[sel])
        B = self.amp[sel]*numpy.sin(self.phase[sel]) + ss2.amp[sel]*numpy.sin(ss2.phase[sel])
        self.amp[sel] = numpy.sqrt( A**2 + B**2 )
        self.phase[sel] = numpy.arctan( B/A )
        self.phase[self.amp == 0.0] = 0.0
        self.phase[sel][A < 0] += numpy.pi


    def plot2D( self, filename = None, type = "Amp", mask = True, cbmin=None, \
                cbmax=None, show = False, fill = True ):
        """
        @brief   Plots the results in a 2D matrix that represents the array.

        @param   filename   String with name where to save the plot as a .png
        @param   type   Type of plot, either "Amp", "Phase" or "RMS". Default "Amp".
        @param   mask   Use mask to plot only unmasked detectors. Masked detectors will be set to zero.
        @param   cbmin  Minimum bound for the color range.
        @param   cbmax  Maximum bound for the color range.
        @param   show   Open an X window with the plot. Default "False".
        @param   fill   Fill in masked detectors using median of unmasked ones in the same column. Def "True".
        """

        mat = numpy.array([[0.0]*32]*33)

        if fill == True:
            amp, phase, rms = self.fillBlanks()
            phase[:] -= self.phaseAz
            mask = False
        else:
            amp = self.amp - 0.0
            phase[:] = self.phase - self.phaseAz
            rms = self.rms - 0.0

        if type == "Phase":
            for i in range(self.ndet):
                if mask == True and self.mask[i] == 0:
                    mat[self.row[i]][self.col[i]] = phase[i]*180/numpy.pi
                elif mask == False:
                    mat[self.row[i]][self.col[i]] = phase[i]*180/numpy.pi
            titulo = "Synchronous Signal Phase"
            if cbmin is None:
                cbmin = numpy.floor(numpy.min(phase*180/numpy.pi))
            if cbmax is None:
                cbmax = numpy.ceil(numpy.max(phase*180/numpy.pi))
        elif type == "RMS":        
            for i in range(self.ndet):
                if mask == True and self.mask[i] == 0:
                    mat[self.row[i]][self.col[i]] = rms[i]
                elif mask == False:
                    mat[self.row[i]][self.col[i]] = rms[i]
            titulo = "Residual Noise RMS"
            if cbmin is None:
                cbmin = numpy.floor(numpy.min(rms)/100)*100
            if cbmax is None:
                cbmax = numpy.ceil(numpy.max(rms)/100)*100
        else:        
            for i in range(self.ndet):
                if mask == True and self.mask[i] == 0:
                    mat[self.row[i]][self.col[i]] = amp[i]
                elif mask == False:
                    mat[self.row[i]][self.col[i]] = amp[i]
            titulo = "Synchronous Signal Amplitude"
            if cbmin is None:
                cbmin = numpy.floor(numpy.min(amp)/100)*100
            if cbmax is None:
                cbmax = numpy.ceil(numpy.max(amp)/100)*100

        fig_width_pt = 500.0  
        inches_per_pt = 1.0/72.27               # Convert pt to inch
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*1.0              # height in inches
        fig_size =  [fig_width,fig_height]
        params = {'backend': 'ps',
                  'axes.labelsize': 10,
                  'text.fontsize': 10, 
                  'legend.fontsize': 10,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'text.usetex': False,
                  'figure.figsize': fig_size}
        im = pylab.matshow(mat, cmap=pylab.cm.jet, vmin=cbmin, vmax=cbmax)
        pylab.colorbar()      
        pylab.xlabel('Columns')
        pylab.ylabel('Rows')
        pylab.title(titulo) 
        pylab.rcParams.update(params)
#        if colmarkers:    
#            for column in colmarkers.keys():
#                pylab.axvline(x=colmarkers[column], color='black', linewidth=1)
#                pylab.axhline(y=colmarkers[column], color='black', linewidth=1)
#        pylab.axis([0,874,0,874])
        if filename is not None:
            pylab.savefig(filename)
        if show == True or filename is None:
            pylab.show()

    def writeToPath( self, path ):
        """
        @brief Stores the synchronous pickup object in Path.
        @param  path   directory path where to store the data object.
        """
        data = { 'amp': self.amp, \
                 'phase': self.phase, \
                 'rms': self.rms, \
                 'phaseAz': self.phaseAz, \
                 'mask': self.mask, \
                 'params': self.params }

        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/sync.pickle" % path
        f = file( filename, 'w' )
        p = pickle.Pickler( f )
        p.dump( data )
        f.close()

    @staticmethod
    def read_from_path(path, tod):
        """
        @brief  Reloads a stored synchronous pickup object.

        @param  path   Path to the directory where the data is stored.
        @param  tod    TOD that corresponds to the object you want to read.
        
        @return ss     Synchronous pickup object with the data read.
        """

        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/sync.pickle" % path
        f = file( filename )
        p = pickle.Unpickler( f )
        data = p.load()
        tod_name = tod
        if not isinstance(tod_name, basestring):
            tod_name = tod.info.name
        assert tod_name == data['params']['todName'], \
            "ERROR: TOD %s and stored name %s don't match" % \
            ( tod_name, data['params']['todName'] )
        ss = Sync(tod, params = data['params'], noFit = True)

        ss.amp[:] = data['amp'][:]
        ss.phase[:] = data['phase'][:]
        ss.rms[:] = data['rms'][:]
        ss.mask[:] = data['mask'][:]
        ss.phaseAz = data['phaseAz']
        ss.params = data['params']

        return ss

    readFromPath = read_from_path


def fitSine(data, wi):
    data -= numpy.mean(data)
    x = numpy.arange(len(data))
    cos = numpy.cos(wi*x)
    sin = numpy.sin(wi*x)
    c2 = numpy.dot(cos, cos)
    s2 = numpy.dot(sin, sin)
    sc = numpy.dot(cos, sin)
    cd = numpy.dot(cos, data)
    sd = numpy.dot(sin, data)
    det = c2*s2 - sc*sc
    A = (s2*cd - sc*sd)/det
    B = (c2*sd - sc*cd)/det
    amp = numpy.sqrt(A*A + B*B)
    if amp == 0: phase = 0.0
    elif A > 0: phase = numpy.arctan(B/A)
    else: phase = numpy.arctan(B/A) + numpy.pi
    return amp, phase

