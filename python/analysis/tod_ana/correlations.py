from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from past.builtins import basestring

import os, numpy, math, scipy.linalg, time, pickle
np = numpy
import moby2
import moby2.util.log as psLib
from scipy.stats import scoreatpercentile
import moby2.tod.filters as filters
from matplotlib import pyplot as plt

pylab = None
def import_pylab():
    global pylab
    import pylab as pl
    pylab = pl

#import actUtils 
#import mobyLib 
#import pyactpol
#from moby.todUtils import TOD
#from moby.todAnalysis import mbFilter
#from moby.utilities import mobyUtils, constants
#from moby.utilities.constants import *

MIN_VICTIMS = 4
# DIM_COL = 32
# DIM_ROW = 33

class correlationObject( object ):
    def __init__(self, tod, rowDominant = False, dets = None):
        """
        @brief  Correlation class: This class is intended to calculate correlations matices and plot them, 
                together with finding correlated modes and subtracting them from the data.
        """
        self.tod = tod.copy()
        self.rowDominance = rowDominant
        self.num_index = None # number of uncut detectors
        self.separators = {}  # index where row/column groups start
        self.correlation = None
        self.correlationFull = None
        self.cov = None
        self.covDiag = None
        self.varianceFigure = None
        self.quality = None
        self.log = []
        if dets is None:
            self.dets = tod.info.array_data['det_uid']
        else:
            self.dets = numpy.array(dets)
        self.Nrows = np.unique(tod.info.array_data['row']).size
        self.Ncols = np.unique(tod.info.array_data['col']).size
        self.rows = self.dets // self.Ncols
        self.cols = self.dets % self.Ncols
        self.corrDets = None
        self.commonMode = None
        self.commonModeType = None
        self.CMpower = 0

    def setTOD(self, tod ):
        """
        @brief Set TOD.TOD associated with this correlation object
        """
        self.tod = tod
        if self.commonMode is not None:
            self.commonMode.tod = tod
            
    def obtainNoiseCoeff( self, noise ):
        """
        @brief  Obtain the coefficiets given by the term (m'N^{-1}m)^{-1} in the PCG.
        @param  noise  noise object contining the precalculater inverse noise (spline) per detector.
        """
        if self.commonMode is not None:
            self.commonMode.generateNoiseCoeff( noise )
         
    def findCorrelation( self, dets = None, nstart = 0, ndata = -1, nmin = 8000 ):
        '''
        @brief   Find the correlation matrix.

        @param   dets  List of detectors to use in correlation calculations.
        '''

        if dets is None:
            dets = self.dets
        self.corrDets = dets

        ntot = self.tod.nsamps - nstart
        if ntot < nmin:
            psLib.trace('moby', 1, 'ERROR: Data selection less than nmin = %d' % nmin)
            return
        if ndata == -1: ndata = ntot
        if ndata > ntot: 
            psLib.trace('moby', 2, 'WARNING: ndata (%d) greater than ntotal (%d). ' \
                        'Reseting ndata = ntotal' % (ndata, ntot))
            ndata = ntot

        psLib.trace('moby', 1, 'Covariance calculation started')
        tic = time.time()
        data = self.tod.data[dets,nstart:nstart+ndata].copy()
        self.cov = numpy.cov(data)
        del data
        toc = time.time()
        dtime = (toc - tic)/60.0
        psLib.trace('moby', 1, 'Covariance done %g minutes' %dtime)


        self.covDiag = numpy.sqrt(numpy.diag(self.cov))
        self.varianceFigure = numpy.median(self.covDiag) 

        self.correlationFull = numpy.zeros([self.tod.det_uid.size, self.tod.det_uid.size])
        for i in range(len(dets)):
            self.correlationFull[dets[i]][dets[i]] = self.cov[i][i]/self.covDiag[i]**2
            for j in range(i+1, len(dets)):
                self.correlationFull[dets[i]][dets[j]] = self.cov[i][j]/self.covDiag[i]/self.covDiag[j]
                self.correlationFull[dets[j]][dets[i]] = self.correlationFull[dets[i]][dets[j]]

        self.reformatMatrix(rowDominant = self.rowDominance)

        self.quality = 0.0
        n = numpy.size(self.correlation, 1)
        for i in range(0, n-1):
            for j in range(i+1, n):
                self.quality += self.correlation[i][j]**2
        self.quality = numpy.sqrt(self.quality/(n*(n-1)/2.0))
                
        self.log.append({'name':'findCorrelation'})


    def reformatMatrix(self, rowDominant = None):
        """
        @brief   Reformat the order of the elements of the correlation matrix to make it row or column
                 dominant.

        @param  rowDominant   "True" for row dominance, "False" for column dominance, and "None" to 
                               toggle between the two. Default "None".
        """

        if rowDominant is None:
            self.rowDominance = not(self.rowDominance)
        else:
            self.rowDominance = rowDominant

        # now reformat the correlation to have dominant rows or columns
        # first make the index list
        psLib.trace('moby', 1, 'Reformatting the correlation matrix')
        

        dets = numpy.array(self.corrDets)
        rows = dets // self.Ncols
        cols = dets % self.Ncols

        index_list = [] 
        self.separators = {}
        if self.rowDominance:
            for row in range (0, self.Nrows):
                sep_found = False
                for col in range (0, self.Ncols):
                    index = dets[(rows == row)*(cols == col)]
                    if len(index) > 0:
                        if not sep_found:
                            self.separators[row] = len(index_list)
                            sep_found = True
                        index_list.append(index)
        else:
            for col in range (0, self.Ncols):
                sep_found = False
                for row in range (0, self.Nrows):
                    index = dets[(rows == row)*(cols == col)]
                    if len(index) > 0:
                        if not sep_found:
                            self.separators[col] = len(index_list)
                            sep_found = True
                        index_list.append(index)

        self.num_index = len(index_list)
        psLib.trace('moby', 1, 'There are '+repr(self.num_index)+' live detectors')
        self.correlation = numpy.zeros((self.num_index, self.num_index))
        for i in range(0,self.num_index):
            for j in range(0,self.num_index):
                index_i = index_list[i]
                index_j = index_list[j]
                self.correlation[i,j] = self.correlationFull[index_i,index_j]


    def plotCovMatrix(self, filename=None, show = True, separators=True, cbmin=-0.3, cbmax=0.3,
                      noTitle = False):
        """
        @brief   Plot correlation matrix

        @param  filename  address of a file (.png) where to save the figure. Default "None".
        @param  show      if "True" it will display the image. Default "True".
        @param  cbmin     set the minimum of the plot range (in correlation units from -1 to 1). 
                           Default "-0.3".
        @param  cbmin     set the maximum of the plot range (in correlation units from -1 to 1). 
                           Default "0.3".
        """
        import_pylab()
        text_x = 0
        text_y = 1.05
        pylab.rcdefaults()
        fig_width_pt = 2000.0  # Get this from LaTeX using \showthe\columnwidth
        inches_per_pt = 1.0/72.27               # Convert pt to inch
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*1.0             # height in inches
        fig_size =  [fig_width,fig_height]
        params = {'backend': 'ps',
                  'axes.labelsize': 20,
                  'text.fontsize': 20,
                  'legend.fontsize': 20,
                  'xtick.labelsize': 16,
                  'ytick.labelsize': 16,
                  'text.usetex': False,
                  'figure.figsize': fig_size}
        if noTitle:
            params['text.fontsize'] = 25
            params['xtick.labelsize'] = 20
            params['ytick.labelsize'] = 20
            text_x = -0.07
            text_y = 1.07
        pylab.rcParams.update(params)
        pylab.matshow(self.correlation, cmap=pylab.cm.jet, vmin=cbmin, vmax=cbmax)
        pylab.colorbar(shrink=0.8)

        # if bad detectors are flagged and removed from the matrix
        # then reduce the dimension of the correlation matrix accordingly
        dim_bound = self.num_index
        ax = pylab.gca()

        if separators:
            for group in self.separators.keys():
                pylab.axvline(x=self.separators[group], color='black', linewidth=1)
                pylab.axhline(y=self.separators[group], color='black', linewidth=1)
        pylab.axis([0,dim_bound,0,dim_bound])

        # generate the text on the graph
        title = 'Correlation matrix for %s %s (' %(self.tod.info.basename,
                                                   self.tod.info.array)
        if self.rowDominance:
            title += 'row dominant)'
            group_key = 'row groups: '
        else:
            title += 'column dominant)'
            group_key = 'column groups: '
        title += ' stdev: %.0f, quality: %.3f' % (math.sqrt(self.varianceFigure), self.quality)
        if not(noTitle): pylab.title(title)

        groups = list(self.separators.keys())
        groups.sort()
        for group in groups:
            group_key += '|'+repr(group)
        group_key += '|'
        pylab.text(text_x,text_y, group_key, size = "large", \
                   horizontalalignment='left', verticalalignment='top', \
                       transform=ax.transAxes)

        if filename is not None:
            pylab.savefig(filename)
            # now trim off the big (boring) borders of the plot
            try: 
                os.system('convert -trim '+filename+' '+filename)
            except:
                psLib.trace('moby', 3, "There was a problem calling imagemagick convert")
        if show == True or filename is None:
            pylab.show()
        else: pylab.close()
        pylab.rcdefaults()


    # def radialCorrelationBinning(self, \
    #                              filename = 'radial_correlation.dat', \
    #                              block_integers = False, \
    #                              block_row = False, \
    #                              block_col = False, \
    #                              make_select = True):
    #     '''
    #     block_integers is a flag that blocks integral distances from the average
    #     these are higher than the general correlation, likely from row/col correlations
    #     '''
    #     # Not updated fro moby2...
    #     greatest = (self.Ncols-1)*(self.Ncols-1)+(self.Nrows-1)*(self.Nrows-1)

    #     if make_select:
    #         selection = numpy.zeros((self.tod.det_uid.size,self.tod.det_uid.size))

    #     bin_dict = {}
    #     dets, rows, cols = self.tod.listUncut()
    #     psLib.trace('moby', 1, 'Binning the covariance radially')
    #     ii = -1; 
    #     for i in dets:
    #       ii += 1
    #       jj = -1
    #       for j in dets:
    #           jj += 1
    #           row1 = self.tod.rows[i]
    #           row2 = self.tod.rows[j]
    #           col1 = self.tod.cols[i]
    #           col2 = self.tod.cols[j]
    #           xsq = (row1-row2)*(row1-row2)
    #           ysq = (col1-col2)*(col1-col2)
    #           try:
    #             cov_val = self.correlation[ii][jj]
    #           except:
    #             self.correlation[ii][jj]

    #           row_test = True
    #           if block_row:
    #               row_test = (self.tod.rows[i] != self.tod.rows[j])
    #           col_test = True
    #           if block_col:
    #               col_test = (self.tod.cols[i] != self.tod.cols[j])

    #           if (self.tod.dets[row1, col1] != actUtils.NO_VALUE) and \
    #              (self.tod.dets[row2, col2] != actUtils.NO_VALUE):      
    #               if cov_val != -1:
    #                   if row_test and col_test:
    #                       lenindex = int(xsq+ysq)

    #                       #length = math.sqrt(float(lenindex))
    #                       arc = math.atan2(float(col2-col1), float(row2-row1))
    #                       #if make_select and (length != int(length)):
    #                       selection[i][j] = arc

    #                       if bin_dict.has_key(lenindex):
    #                           bin_dict[lenindex].append(cov_val)
    #                       else:
    #                           bin_dict[lenindex] = [ cov_val ]

    #     f_output = open(filename,'w')
    #     for i in range(0, greatest+1):
    #         if (bin_dict.has_key(i)):
    #             length = math.sqrt(i)

    #             integer_test = True
    #             if block_integers:
    #                 integer_test = (length != int(length))

    #             if (integer_test):
    #                  median = numpy.median(bin_dict[i])
    #                  mean = numpy.mean(bin_dict[i])
    #                  stdev = numpy.std(bin_dict[i])
    #                  counts = len(bin_dict[i])
    #                  output = repr(length)
    #                  output += ' '+ repr(mean)
    #                  output += ' '+ repr(median)
    #                  output += ' '+ repr(stdev)
    #                  output += ' '+ repr(counts)
    #                  f_output.write(output+'\n')
    #     if make_select:
    #         return selection
    #     else:
    #         return None



    def radialCorrelationBinning( self, bins=100 ):
        """
        @brief Compute space correlations
        
        @param bins   Bins for the correlation function:
                      - if int, number of bins
                      - if list, edges of the bins 
        """
        # Check that correlation matrix is available
        if self.correlation is None:
            print("Aborted! findCorrelation() must be performed first")
            return 0


        # Get distances between all pairs of detectors
        x_det = self.tod.info.array_data['sky_x'][self.dets]
        y_det = self.tod.info.array_data['sky_y'][self.dets]
        dist = np.sqrt( (x_det[np.newaxis,:] - x_det[:,np.newaxis])**2
                        + (y_det[np.newaxis,:] - y_det[:,np.newaxis])**2)
        assert dist.shape == self.correlation.shape, "Distance and correlation matrices do not have the same dimensions"
        

        # Select the pairs for each bin
        if type(bins) == int:
            bin_edges = np.linspace( dist.min(), dist.max(), 100 )
        elif type(bins) == list:
            bin_edges =np.asarray(bins)
        bin_center = ( bin_edges[:-1] + bin_edges[1:] ) / 2
        bins = [(b0,b1) for b0,b1 in zip(bin_edges[:-1],bin_edges[1:])]
        sel_bins = [ (dist>=b[0])*(dist<b[1]) for b in bins]

        # Average all pairs in each bin
        binned_corr = []
        low_lim = []
        up_lim = []
        N = []
        for sel in sel_bins:
            llim, med, ulim = scoreatpercentile(self.correlation[sel], [34, 50, 66])
            binned_corr.append(med)
            low_lim.append(llim)
            up_lim.append(ulim)
            N.append(sel.sum())

        self.radial_corr = np.asarray(binned_corr)
        self.radial_bins = np.asarray(bin_center)
        self.radial_corr_errlow = np.asarray(binned_corr) - np.asarray(low_lim)
        self.radial_corr_errhigh = np.asarray(up_lim) - np.asarray(binned_corr)
        

    def plotRadialCorrelation( self, xlim=(-0.1,1), ylim=(-1.1,1.1), savefig=False ):
        """
        @brief Plot the space correlation function (radialCorrelationBinning must have been called before)
        """
        if self.radial_corr is None or self.radial_bins is None:
            print("Aborted! radialCorrelationBinning() must be called first")
            return 0

        err = np.vstack([self.radial_corr_errlow, self.radial_corr_errhigh])
        if savefig:
            plt.ioff()
        plt.figure()
        plt.errorbar(self.radial_bins, self.radial_corr, err,
                     marker='o', mfc='b', ecolor='b', ls='None')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel('Pair distance [deg]')
        plt.ylabel('Correlation')
        plt.title('%s' %self.tod.info.basename)
        plt.grid()
        if savefig:
            plt.savefig(savefig)
            plt.close('all')
        else:
            plt.show()



    def write( self, data, path ):
        """
        @brief Stores the correlation object in Path.
        """
        data['todName']         = self.tod.name
        data['correlationFull'] = self.correlationFull
        data['correlation']     = self.correlation
        data['rowDominance']    = self.rowDominance
        data['separators']      = self.separators
        data['num_index']       = self.num_index
        data['quality']         = self.quality
        data['cov']             = self.cov
        data['covDiag']         = self.covDiag
        data['varianceFigure']  = self.varianceFigure
        data['dets']            = self.dets
        data['corrDets']        = self.corrDets
        data['log']             = self.log
        data['commonModeType']  = self.commonModeType
        if self.commonMode is not None: data['commonMode'] = self.commonMode.exportModes()
        else: data['commonMode'] = None

        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/correlation.pickle" % path
        f = file( filename, 'w' )
        p = pickle.Pickler( f )
        p.dump( data )
        f.close()



    def read( self, path ):
        """
        @brief  Reloads a stored correlation object.

        @param  path   Path to the directory where the data is stored.
        @param  tod    TOD that corresponds to the object you want to read.
        
        @return data   Dictionary with the data read.
        """

        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/correlation.pickle" % path
        f = file( filename )
        p = pickle.Unpickler( f )
        data = p.load()
        assert self.tod.name == data['todName'], \
            "correlations: TOD %s and stored correlations %s don't match" % \
            ( self.tod.name, data['todName'] )

        self.todName         = data['todName']
        self.correlationFull = data['correlationFull']
        self.correlation     = data['correlation']
        self.rowDominance    = data['rowDominance']
        self.separators      = data['separators']
        self.num_index       = data['num_index']
        self.quality         = data['quality']
        self.cov             = data['cov']
        self.covDiag         = data['covDiag']
        self.varianceFigure  = data['varianceFigure']
        #self.dets            = data['dets']
        self.corrDets        = data['corrDets']
        self.log             = data['log']
        #self.commonModeType  = 'multi'
        self.commonModeType  = data['commonModeType']
        if data['commonMode'] is not None: 
            self.commonMode = modeSet(self.tod, [], self.dets)
            self.commonMode.importModes(data['commonMode'])

        return data

    ################## COMMON MODE ##############################
    
    def findCommonMode( self, dets = None, useMedian = False, power = 0, filter = None ):
        """
        @brief  Calculates the common mode of the live detectors.
    
        @param   dets       List of detectors to use in common mode.
        @param   useMedian  If true use Median to calculate the common mode,
                            use the mean otherwise.
        @param   power      Power of the polynomial in case of polynomial CM.
        @param   filter     Filter to apply to common mode.
        """
        power = int(power)
        if power < 0:
            psLib.trace('moby', 1, "WARNING: power < 0 is not allowed. Setting to 0.")
            power = 0
        if power > 3:
            psLib.trace('moby', 1, "WARNING: power > 3 is not allowed. Setting to 3.")
            power = 3
        self.CMpower = power

        if dets is None:
            dets = self.dets
        if len(dets) == 0: dets = [0]

        if power == 0:
            dets = numpy.array(dets, dtype='int32')
            ##cm = numpy.empty(self.tod.ndata)
            if useMedian:
                ##mobyLib.transposeMedian(self.tod.ctod, list(dets), cm)
                cm = moby2.libactpol.data_median_axis0(self.tod.data, dets)
            else:
                ##mobyLib.transposeMean(self.tod.ctod, list(dets), cm)
                cm = moby2.libactpol.data_mean_axis0(self.tod.data, dets)
            modes = modeSet(self.tod, cm.reshape([1,len(cm)])/numpy.linalg.norm(cm), dets)
        else:
            modes = polyCommonMode(self.tod, dets, power)
        if filter is not None: modes.filterModes(filter)
        return modes

    
    ################## MODE REMOVAL ##############################
    
    def removeMode( self, mode, dets = None, kill = True, doFit = True ):
        """
        @brief  Removes a single mode from a set of detectors. The set of detectors must include 
                more than one detector.
    
        @param  tod    TOD object from where to remove the mode.
        @param  mode   mode that needs to be removed
        @param  dets   list of detectors from which to remove the modes.
        @param  doFit  choose whether to fit the mode to the detector or not.
        """
        if dets is None:
            dets = self.dets

        norm = numpy.linalg.norm(mode)
    
        n = len(dets)
        if norm != 0:
            for i in dets:
                if doFit:
                    coeff = numpy.dot(self.tod.data[i], mode)/norm**2
                else: coeff = 1.0
                if kill:
                    self.tod.data[i] -=  mode*coeff
                else:
                    self.tod.data[i] -=  mode*coeff * (n-1)/n
            self.tod.abuses += [{'name':'removeMode', 'type':'modeRemoval', 'detectors': dets}]

    
    
####################################################################################
################ DARK AND DEAD DETECTOR DECORRELATION OBJECT #######################
####################################################################################

class darkCorrObj( correlationObject ):
    '''
    @brief  Correlation object for dark detectors. Finds and removes correlations present in 
            the dark and dead detectors.
    '''
    def __init__(self, tod, rowDominant = False, dets = None):
        """
        @brief  Initialization for correlation object for dark and dead detectors.
        @param  tod           TOD related to dark correlation object.
        @param  rowDominant   dominance of the correlation matrix for plotting.
        """

        correlationObject.__init__(self, tod, rowDominant = rowDominant, dets = dets)

        self.extractedDrift = False
        self.darkModes = None
        self.darkRowModes = None
        self.darkColModes = None
        self.darkCommonMode = None


    def setTOD(self, tod ):
        """
        @brief Set TOD.TOD associated with this correlation object
        """
        correlationObject.setTOD( self, tod )
        if self.darkModes is not None:
            self.darkModes.tod = tod
        if self.darkRowModes is not None:
            self.darkRowModes.tod = tod
        if self.darkColModes is not None:
            self.darkColModes.tod = tod

    def obtainNoiseCoeff( self, noise ):
        """
        @brief  Obtain the coefficiets given by the term (m'N^{-1}m)^{-1} in the PCG.
        @param  noise  noise object contining the precalculater inverse noise (spline) per detector.
        """
        if self.darkModes is not None:
            self.darkModes.generateNoiseCoeff( noise )
        if self.darkRowModes is not None:
            self.darkRowModes.generateNoiseCoeff( noise )
        if self.darkColModes is not None:
            self.darkColModes.generateNoiseCoeff( noise )


    def removeDarkCM( self, filter = None, noise = None, forceFit = True ):
        """
        @brief   Removes the commonMode of the dark and dead detectors from the data.
        @param   filter   Filter to apply to modes.
        """
        if self.darkCommonMode is None: 
            if filter is None: filter = filters.sine2lowPass(self.tod, fc = 0.05)
            self.findDarkCM(filter = filter)
        self.darkCommonMode.removeModeSet( forceFit = forceFit, noise = noise )


    def removeDarkCorr(self, nModes = 12, forceFit = True, filter = None, noise = None):
        """
        @brief   Removes the main correlated modes found in the dark and dead detectors from the data.
        @param   nModes   Specify the number of modes to use. Default = 12.
        @param   forceFit refit existing modeSet to detector time streams for subtraction
        @param   filter   Filter to apply to modes.
        """
        if len(self.tod.dark) < nModes:
            psLib.trace('moby', 1, "ERROR: Number of dark detectors %d less than number " \
                                   "of modes to subtract %d" % (len(self.tod.dark), nModes))
            return
        try: 
            if len(self.darkModes.modes) != nModes:
                self.findDarkCorr(nModes = nModes, filter = filter)
        except:
            self.findDarkCorr(nModes = nModes, filter = filter)
        self.darkModes.removeModeSet(forceFit = forceFit, noise = noise)


    def removeDarkRowCorr(self, nModes = 4, forceFit = True, filter = None, noise = None):
        """
        @brief   Removes the main correlated modes found in the dark and dead detectors from the data.
        @param   nModes   Specify the number of modes to use. Default = 4.
        @param   forceFit refit existing modeSet to detector time streams for subtraction
        @param   filter   Filter to apply to modes.
        """
        if self.tod.nrow < nModes:
            psLib.trace('moby', 1, "ERROR: Number of rows %d less than number " \
                                   "of modes to subtract %d" % (self.tod.nrow, nModes))
            return
        try: 
            if len(self.darkRowModes.modes) != nModes:
                self.findDarkRowCorr(nModes = nModes, filter = filter)
        except:
            self.findDarkRowCorr(nModes = nModes, filter = filter)
        self.darkRowModes.removeModeSet(forceFit = forceFit, noise = noise)


    def removeDarkColCorr(self, nModes = 4, forceFit = True, filter = None, noise = None):
        """
        @brief   Removes the main correlated modes found in the dark and dead detectors from the data.
        @param   nModes   Specify the number of modes to use. Default = 4.
        @param   forceFit refit existing modeSet to detector time streams for subtraction
        @param   filter   Filter to apply to modes.
        """
        if self.tod.ncol < nModes:
            psLib.trace('moby', 1, "ERROR: Number of columns %d less than number " \
                                   "of modes to subtract %d" % (self.tod.ncol, nModes))
            return
        try: 
            if len(self.darkColModes.modes) != nModes:
                self.findDarkColCorr(nModes = nModes, filter = filter)
        except:
            self.findDarkColCorr(nModes = nModes, filter = filter)
        self.darkColModes.removeModeSet(forceFit = forceFit, noise = noise)


    def findDarkCorr(self, nModes = 12, filter = None):
        """
        @brief   Finds the main correlated modes from in the dark and dead detectors from the data and 
                 their fit to the live detectors.
        @param  nModes   Specify the number of modes to use. Default = 12.
        @param  filter   Filter to apply to modes.
        """
        if self.extractedDrift == False:
            self.extractDarkDrift()

        self.darkModes = modeSet(self.tod, self.tod.data[self.tod.dark], self.dets)
        if filter is not None: self.darkModes.filterModes(filter)
        self.darkModes.selectMainModes(nModes = nModes)


    def findDarkRowCorr(self, nModes = 4, filter = None):
        """
        @brief  Finds the row common mode of the dark and dead detectors.
        @param  nModes   Specify the number of modes to use. Default = 4.
        @param  filter   Filter to apply to modes.
        """
        if self.extractedDrift == False:
            self.extractDarkDrift()

        self.darkRowModes = modeSet(self.tod, [], self.dets)
        for row in  range(self.tod.nrow):
            li = self.tod.dark[self.tod.rows[self.tod.dark] == row]
            if len(li) > 1:
                self.darkRowModes.append(numpy.mean(self.tod.data[li], axis = 0))
            elif len(li) == 1:
                self.darkRowModes.append(self.tod.data[li])
        if filter is not None: self.darkRowModes.filterModes(filter)
        self.darkRowModes.selectMainModes(nModes = nModes)


    def findDarkColCorr(self, nModes = 4, filter = None):
        """
        @brief  Finds the column common mode of the dark and dead detectors.
        @param  nModes   Specify the number of modes to use. Default = 4.
        @param  filter   Filter to apply to modes.
        """
        if self.extractedDrift == False:
            self.extractDarkDrift()

        self.darkColModes = modeSet(self.tod, [], self.dets)
        for col in  range(self.tod.ncol):
            li = self.tod.dark[self.tod.cols[self.tod.dark] == col]
            if len(li) > 1:
                self.darkColModes.append(numpy.mean(self.tod.data[li], axis = 0))
            elif len(li) == 1:
                self.darkColModes.append(self.tod.data[li])
        if filter is not None: self.darkColModes.filterModes(filter)
        self.darkColModes.selectMainModes(nModes = nModes)


    def findDarkCM( self, filter = None ):
        """
        @brief  Finds the common mode of the dark and dead detectors.
        @param  filter  Filter to apply to dark common mode
        """
        if len(self.tod.dark) < 1:
            psLib.trace('moby', 1, "ERROR: Number of dark detectors %d less than number " \
                        "of modes to subtract %d" % (len(self.tod.dark), nModes))
            return
        self.commonMode = self.findCommonMode(dets = self.tod.dark, useMedian = False, 
                                              filter = filter)
        self.darkCommonMode = modeSet(self.tod, self.commonMode.modes, self.dets)

    def extractDarkDrift( self ):
        """
        @brief  Finds the common mode of the dark and dead detectors.
        """
        if len(self.tod.dark) < 1:
            psLib.trace('moby', 1, "ERROR: Number of dark detectors %d less than number " \
                        "of modes to subtract %d" % (len(self.tod.dark), nModes))
            return
        if self.commonMode is None:
            filt = filters.sine2lowPass(self.tod, fc = 0.05)
            self.findDarkCM(filter = filt)
        self.commonMode.removeModeSet()
        self.extractedDrift = True


    def writeToPath( self, path ):
        """
        @brief   Saves the data in the dark correlation object to path.
        @param   path  Path where to store the data
        """
        data = {}
        if self.darkModes is not None: data['darkModes'] = self.darkModes.exportModes()
        else: data['darkModes'] = None
        if self.darkRowModes is not None: data['darkRowModes'] = self.darkRowModes.exportModes()
        else: data['darkRowModes'] = None
        if self.darkColModes is not None: data['darkColModes'] = self.darkColModes.exportModes()
        else: data['darkColModes'] = None
        self.write(data, path)


    def readFromPath( path, tod ):
        """
        @brief   Reads the stored data for a dark correlation object from path.
        @param   tod   Tod for which the correlation object was originally created.
        @param   path  Path where the data was stored.
        """
        c = darkCorrObj(tod)
        data = c.read(path)
        if data['darkModes'] is not None:
            c.darkModes = modeSet(tod, [], c.dets) 
            c.darkModes.importModes(data['darkModes'])
        if data['darkRowModes'] is not None:
            c.darkRowModes = modeSet(tod, [], c.dets) 
            c.darkRowModes.importModes(data['darkRowModes'])
        if data['darkColModes'] is not None:
            c.darkColModes = modeSet(tod, [], c.dets) 
            c.darkColModes.importModes(data['darkColModes'])

        return c

    readFromPath = staticmethod(readFromPath)

####################################################################################
################### LIVE DETECTOR DECORRELATION OBJECT #############################
####################################################################################

class liveCorrObj( correlationObject ):
    def __init__(self, tod, rowDominant = False, dets = None):
        correlationObject.__init__(self, tod, rowDominant = rowDominant, dets = dets)

        self.liveSVDModes = None
        self.liveRowModes = None
        self.liveColModes = None
        self.commonModeType = 'single'

    def setTOD(self, tod ):
        """
        @brief Set TOD.TOD associated with this correlation object
        """
        correlationObject.setTOD( self, tod )
        if self.liveSVDModes is not None:
            self.liveSVDModes.tod = tod
        if self.liveRowModes is not None:
            self.liveRowModes.tod = tod
        if self.liveColModes is not None:
            self.liveColModes.tod = tod

    def obtainNoiseCoeff( self, noise ):
        """
        @brief  Obtain the coefficiets given by the term (m'N^{-1}m)^{-1} in the PCG.
        @param  noise  noise object contining the precalculater inverse noise (spline) per detector.
        """
        correlationObject.obtainNoiseCoeff( self, noise )
        if self.liveSVDModes is not None:
            self.liveSVDModes.generateNoiseCoeff( noise )
        if self.liveRowModes is not None:
            self.liveRowModes.generateNoiseCoeff( noise )
        if self.liveColModes is not None:
            self.liveColModes.generateNoiseCoeff( noise )

    def removeArrayCorr( self, noise = None, forceFindCM = False, filter = None, \
                         useMedian = False, cmtype = 'single', params = None):
        """
        @brief  Calculates and removes the common mode of the live detectors.

        @param noise         noise object to use in removal.
        @param forceFindCM   reestimate common mode even if it's already been estimated
        @param filter        Define filter to apply to common mode before subtracting it from the data.
        @param useMedian     Use the median instead of the mean for estimating the common mode.
        @param cmtype        Type of common mode to use:
                                'single': simple common mode.
                                'multi':  multi common mode obtained from sub arrays.
                                'poly':   polynomial common mode. It substracts directly
                                          without fitting the CM to each detector.
        @param params        Parameters to use depending on the type of common mode used.
        """
        if useMedian and type == 'poly':
            psLib.trace('moby', 2, 'WARNING: Cannot use median and poly type together.\n' \
                                   'Setting to single type.')
            cmtype = 'single'

        if cmtype == 'poly' and noise is not None:
            psLib.trace('moby', 2, 'WARNING: Poly CommonMode with noise feature still not implemented.')

        if type != self.commonModeType:
            self.commonModeType = cmtype
            forceFindCM = True

        if cmtype == 'single': par = {'power': 0}
        elif cmtype == 'poly': par = {'power': 2}
        elif cmtype == 'multi': par = {'nModes': 8, 'squareSize': 8}
        else: 
            psLib.trace('moby', 2, 'ERROR: unknown Common Mode type.')
            return
        if params is not None: par.update(params)

        dets, rows, cols = self.dets, self.rows, self.cols
        if (self.commonMode is None) or forceFindCM:
            psLib.trace('moby', 3, 'Finding new common mode')
            if cmtype == 'multi':
                self.commonMode = self.findMultiCommonMode( useMedian = useMedian, filter = filter, \
                                                            **par)
            else:
                self.commonMode = self.findCommonMode( dets = dets, useMedian = useMedian, \
                                                       filter = filter, **par )

        if cmtype == 'poly':
            if par['power'] == 0: self.removeMode(self.commonMode.modes[0], doFit = False)
            else: self.commonMode.removeModes()
        else:
            self.commonMode.removeModeSet( noise = noise )

        self.log.append({'name':'ArrayCorr'})


    def findMultiCommonMode( self, nModes = 8, squareSize = 8, useMedian = False, filter = None,
                             dets = None):
        """
        @brief  Divide the array given a square size finding a common mode for each square. Then select
                a specified number of modes out of the original set of common modes.
        @param  nModes      Number of modes to extract.
        @param  squareSize  Side of the square window used to divide the array. Given in number of detectors.
        @param  useMedian   Use median instead of mean in the calculation of each mode.
        @param  filter      Filter to apply to modes.
        """
        # Not updated for moby2
        if dets is None:
            alldets, allrows, allcols = self.tod.listUncut()
        else:
            alldets = list(dets)
            allrows = list(self.tod.rows[dets])
            allcols = list(self.tod.cols[dets])

        row1 = numpy.min(allrows)
        row2 = numpy.max(allrows)
        if row2 == 32:
            row2 = 31
        nrow = row2-row1+1
        nRowBlock = int(numpy.round(nrow/squareSize))
        if nRowBlock == 0:
            nRowBlock = 1

        col1 = numpy.min(allcols)
        col2 = numpy.max(allcols)
        ncol = col2-col1+1
        nColBlock = int(numpy.round(ncol/squareSize))
        if nColBlock == 0:
            nColBlock = 1

        modes = modeSet(self.tod, [], alldets)
        alldets = numpy.array(alldets)
        while nRowBlock > 0:
            nr = nrow//nRowBlock
            r = list(range(row1,row1+nr))
            row1 += nr
            nrow -= nr
            nRowBlock -= 1
            nbc = nColBlock
            tcol1 = col1
            tncol = ncol
            while nbc > 0:
                nc = tncol//nbc
                c = list(range(tcol1,tcol1+nc))
                tcol1 += nc
                tncol -= nc
                nbc -= 1

                br, bc = TOD.blockToList(r,c)
                br, bc, dets = self.tod.listFromRowCol(br, bc)
                sdets = []
                for d in dets: 
                    if numpy.any(alldets == d): sdets.append(d)
                if len(sdets) > MIN_VICTIMS:
                    sdets = numpy.array(sdets, dtype='int32')
                    #cm = numpy.empty(self.tod.ndata)
                    if useMedian:
                        cm = moby2.libactpol.data_median_axis0(self.tod.data, sdets)
                    else:
                        cm = moby2.libactpol.data_mean_axis0(self.tod.data, sdets)
                    modes.append(cm*len(sdets))

        if filter is not None: modes.filterModes(filter)
        modes.selectMainModes(nModes = nModes, plot = False)
        return modes


    def removeGaussCM( self, sigma = 10. ):
        """
        @brief  Remove a common mode calculated for every detector by averaging the detectors
                around it weighted by a gaussian of a given sigma.
        @param sigma  Sigma of the gausian used for weighting in units of detector number.
        """
        sigma = float(sigma)
        dets, rows, cols = self.tod.listUncut()
        dets = numpy.array(dets)
        data2 = self.tod.data.copy()
        for i in dets:
            xx = self.tod.rows - self.tod.rows[i]
            yy = self.tod.cols - self.tod.cols[i]
            xx = xx[dets]
            yy = yy[dets]
            w = numpy.exp(-(xx**2+yy**2)/sigma**2)
            dets2 = dets[w>0.1]
            w = w[w>0.1]
            mode = numpy.dot(data2[dets2].transpose(), w.transpose())
            norm = numpy.linalg.norm(mode)
            if norm == 0.:
                psLib.trace('moby', 2, "WARNING, zero norm for detector %d" % i )
            else:
                self.tod.data[i] -= mode * numpy.dot(data2[i], mode)/norm**2
        del data2


    def removeSVDCorr( self, nModes = 12, filter = None, ndata = 65536, plot = False, 
                       noise = None):
        """
        @brief  Finds the main modes in the data by SVD decomposition and substracts them
        @param  nModes   number of modes to subtract
        @param  filter   filter to apply to modes.
        @param  ndata    number of data points to use in calculating the covariance.
        @param  noise    noise object to use in removal.
        """
        if self.liveSVDModes is None:
            self.liveSVDModes = self.findSVDCorr( nModes = nModes, filter = filter, \
                                             ndata = ndata, plot = plot)
        self.liveSVDModes.removeModeSet( noise = noise )

        self.log.append({'name':'LiveSVDCorr'})


    def removeColCorr( self, nModes = 4, filter = None, noise = None):
        """
        @brief  Calculates and removes the column common mode of the live detectors.
        @param  nModes   number of modes to subtract
        @param  filter   filter to apply to modes.
        @param  noise    noise object to use in removal.
        """
        if self.liveColModes is None:
            self.liveColModes = self.findColCorr( nModes = nModes, filter = filter)
        self.liveColModes.removeModeSet( noise = noise )
        self.log.append({'name':'LiveColCorr'})


    def removeRowCorr( self, nModes = 4, filter = None, noise = None):
        """
        @brief  Calculates and removes the row common mode of the live detectors.
        @param  nModes   number of modes to subtract
        @param  filter   filter to apply to modes.
        @param  noise    noise object to use in removal.
        """
        if self.liveRowModes is None:
            self.liveRowModes = self.findRowCorr( nModes = nModes, filter = filter)
        self.liveRowModes.removeModeSet( noise = noise )
        self.log.append({'name':'LiveRowCorr'})


    def findSVDCorr( self, nModes = 12, filter = None, ndata = 65536, plot = False):
        """
        @brief  Finds the main modes in the data by SVD decomposition and substracts them
        @param  nModes   number of modes to subtract
        @param  filter   filter to apply to modes.
        @param  ndata    number of data points to use in calculating the covariance.
        @return liveSVDModes  SVD correlation modeSet object
        """
        liveSVDModes = modeSet(self.tod, self.tod.data[self.dets], self.dets)
        if filter is not None: liveSVDModes.filterModes(filter)
        liveSVDModes.selectMainModes( nModes = nModes, ndata = ndata, plot = plot )
        return liveSVDModes


    def findRowCorr( self, nModes = 4, filter = None):
        """
        @brief  Calculates and removes the row common mode of the live detectors.
        @param  nModes   number of modes to subtract
        @param  filter   filter to apply to modes.
        @return liveRowModes  SVD correlation modeSet object
        """
        liveRowModes = modeSet(self.tod, [], self.dets)
        for row in range(self.Nrows):
            victims = list(self.dets[self.rows == row])
            if len(victims) > MIN_VICTIMS:
                liveRowModes.append(numpy.sum(self.tod.data[victims], axis = 0))
        if filter is not None: liveRowModes.filterModes(filter)
        liveRowModes.selectMainModes( nModes = nModes )
        return liveRowModes


    def findColCorr( self, nModes = 4, filter = None):
        """
        @brief  Calculates and removes the column common mode of the live detectors.
        @param  nModes   number of modes to subtract
        @param  filter   filter to apply to modes.
        @return liveColModes  SVD correlation modeSet object
        """
        liveColModes = modeSet(self.tod, [], self.dets)
        for col in range(self.Ncols):
            victims = list(self.dets[self.cols == col])
            if len(victims) > MIN_VICTIMS:
                liveColModes.append(numpy.sum(self.tod.data[victims], axis = 0))
        if filter is not None: liveColModes.filterModes(filter)
        liveColModes.selectMainModes( nModes = nModes )
        return liveColModes


    def writeToPath( self, path ):
        """
        @brief   Saves the data in the live correlation object to path.
        @param   path  Path where to store the data
        """
        data = {}
        if self.liveRowModes is not None: data['liveRowModes'] = self.liveRowModes.exportModes()
        else: data['liveRowModes'] = None
        if self.liveColModes is not None: data['liveColModes'] = self.liveColModes.exportModes()
        else: data['liveColModes'] = None
        if self.liveSVDModes is not None: data['liveSVDModes'] = self.liveSVDModes.exportModes()
        else: data['liveSVDModes'] = None
        self.write(data, path)

    def readFromPath( path, tod ):
        """
        @brief   Read the data from a live correlation object stored in path.
        @param   tod   Tod for which the correlation object was originally created.
        @param   path  Path where the data was stored.
        """
        c = liveCorrObj(tod)
        data = c.read(path)
        if data['liveRowModes'] is not None: 
            c.liveRowModes = modeSet(tod, [], c.dets)
            c.liveRowModes.importModes(data['liveRowModes'])
        else: c.liveRowModes = None
        if data['liveColModes'] is not None: 
            c.liveColModes = modeSet(tod, [], c.dets)
            c.liveColModes.importModes(data['liveColModes'])
        else: c.liveColModes = None
        if data['liveSVDModes'] is not None: 
            c.liveSVDModes = modeSet(tod, [], c.dets)
            c.liveSVDModes.importModes(data['liveSVDModes'])
        else: c.liveSVDModes = None

        return c

    readFromPath = staticmethod(readFromPath)

####################################################################################
######################### CORRELATION MODES OBJECT #################################
####################################################################################

class modeSet( object ):
    def __init__(self, tod, modes, dets, comment = None):
        """
        @brief Modes object
        @param tod    TOD from which the modes belong.
        @param modes  modes to add to object.
        @param dets   list of uncut detectors asociated to modes.
        """
        
        self.modes = numpy.array(modes, dtype = 'float64').copy()
        self.tod = tod
        self.dets = dets
        self.coeff = []
        self.corr = []
        self.removed = False
        self.detrended = False
        self.noiseCoeff = None
        self.comment = comment
        self.renorm = None

    def setTOD( self, tod ):
        """
        @brief Set tod in object
        """
        self.tod = tod


    def append( self, mode, orthogonalize = False):
        """
        @brief  Add a new mode to a set of modes.
        @param  mode  new mode to add (could also be a set of modes)
        """
        mode = numpy.array(mode, dtype = 'float64')
        if numpy.ndim(mode) == 1: mode = mode.reshape([1, len(mode)])
        if len(self.modes) == 0:
            self.modes = mode
        else:
            if orthogonalize:
                for m in mode:
                    temp = m - numpy.dot(numpy.dot(self.modes, m), self.modes)
                    temp /= numpy.linalg.norm(temp)
                    self.modes = numpy.vstack((self.modes, temp))
            else:
                self.modes = numpy.vstack((self.modes, mode))

    def removeModeSet( self, forceFit = True, dets = None, noise = None ):
        """
        @brief   Removes a set of modes from the specified detectors by finding the least
                 mean square fit of them to the data.
        @param   forceFit   Allowes you to use a previously calculated fit to the data.
        @param   dets       Define from which detectors remove the modeSet
        @param   noise      Apply inverse noise to the modes in the PCG context
        """
        psLib.trace('moby', 5, "Removing Mode Set. ForceFit = %s"%forceFit)
        if len(self.modes) == 0:
            psLib.trace('moby', 1, "ERROR: No modes to remove.")
            return

        if dets is not None:
            forceFit = True
            useDets = dets
        else:
            useDets = self.dets
        if noise is not None:
            if self.noiseCoeff is None: self.generateNoiseCoeff( noise )
            data = self.tod.data[useDets].copy()
            print('inverseFilter')
            noise.inverseFilterTOD()
            print('fitModes')
            self.fitModes(dets = useDets)
            self.coeff = numpy.multiply( self.noiseCoeff, self.coeff )
            self.tod.data[useDets] = data[:]
            del data
            #mobyLib.removeTODModes(self.tod.ctod, list(useDets), self.modes, self.coeff)
        else:
            if self.coeff == [] or forceFit:
                self.fitModes(dets = useDets)
            #mobyLib.removeTODModes(self.tod.ctod, list(useDets), self.modes, self.coeff)
        moby2.libactpol.remove_modes(
            self.tod.data, numpy.array(useDets, dtype='int32'),
            self.modes.astype('float32'), self.coeff.astype('float64'))
        self.tod.abuses += [{'name':'removeMode', 'type':'ModeRemoval', 'detectors': useDets}]
        self.removed = True
    
    
    def fitModes( self, all = False, dets = None, useCuts = False ):
        """
        @brief   Fit a set of orthonormal modes to a specific set of detector time-streams.
        @param   all    Fit all detectors (overrides the dets specification)
        @param   dets   Manually define the detectors to fit.
        """
        if len(self.modes) == 0:
            psLib.trace('moby', 1, "ERROR: No modes to fit.")
            return

        if all: rows, cols, useDets = self.tod.listFromRowCol()
        elif dets is not None: useDets = dets
        else: useDets = self.dets

        #self.coeff = numpy.empty([len(useDets), len(self.modes)])
        #mobyLib.dotProductTOD(self.tod.ctod, list(useDets), self.modes, self.coeff, useCuts)
        self.coeff = moby2.libactpol.data_dot_modes(
            self.tod.data, numpy.array(useDets, dtype='int32'),
            numpy.asarray(self.modes, dtype='float32'),
            moby2.util.encode_c(self.tod.cuts))
        if useCuts:
            if self.renorm is None: self.getRenorm()
            self.coeff /= self.renorm
        self.coeff[numpy.where(numpy.isinf(self.coeff))] = 0.
        self.coeff[numpy.where(numpy.isnan(self.coeff))] = 0.

    def getRenorm( self, all = False, dets = None ):
        """
        @brief Obtain m^T*m considering the partial cuts per detector. This can be used
        to renormalize the fit of the modes to the TOD when cuts are present.
        """
        if len(self.modes) == 0:
            psLib.trace('moby', 1, "ERROR: No modes to fit.")
            return

        if all: rows, cols, useDets = self.tod.listFromRowCol()
        elif dets is not None: useDets = dets
        else: useDets = self.dets

        ##self.renorm = numpy.empty([len(useDets), len(self.modes)])
        ##mobyLib.getModeRenorm(self.tod.ctod, list(useDets), self.modes, self.renorm)
        self.renorm = moby2.libactpol.data_dot_modes(
            self.tod, numpy.array(useDets, dtype='int32'),
            self.modes.astype('float32'), self.tod.cuts._encode_c())

    
    def selectMainModes( self, nModes = 12, ndata = None, plot = False ):
        """
        @brief   select a set of the main modes from a set of source modes
        @param  nModes      number of main modes to select.
        """
        N = self.modes.shape[0]
        if N < nModes: nModes = N
        if N == 0: return
        elif N == 1: self.modes[0] /= numpy.linalg.norm(self.modes[0])
        elif N > 20:
            if ndata is None or ndata > self.tod.nsamps: ndata = self.tod.nsamps
            try:
                cov = numpy.cov(self.modes[:,:ndata])
                u, w, v = numpy.linalg.svd(cov, full_matrices = 0)
            except:
                cov = numpy.cov(self.modes[:,-ndata:])
                u, w, v = numpy.linalg.svd(cov, full_matrices = 0)
            self.kernel = v[:nModes]
            for i in range(nModes): self.kernel[i] /= numpy.sqrt(w[i]*len(self.modes[0]))
            self.modes = numpy.dot(self.kernel, self.modes)
        else:
            u,w,v = numpy.linalg.svd(self.modes, full_matrices = 0)
            self.modes = v[0:nModes].copy()

        if plot:
            import_pylab()
            pylab.plot(w[:nModes], '.'), pylab.show()
    

    def generateNoiseCoeff( self, noise ):
        """
        @brief generate coefficients from the (m'N^{-1}m)^{-1} term of the PCG.
        @param noise   noise object that contains the inverse noise spline fit.
        """
        delta = 1./self.tod.dt[-1]
        if numpy.rank(self.modes) == 1:
            self.noiseCoeff = numpy.zeros(len(self.dets))
            mobyLib.inverseFilterMode(noise.noise_container, self.modes.tolist(), 
                                      list(self.dets), self.noiseCoeff, delta)
        else: 
            self.noiseCoeff = numpy.zeros([len(self.modes),len(self.dets)])
            for i in range(len(self.modes)):
                mobyLib.inverseFilterMode(noise.noise_container, self.modes[i].tolist(), 
                                          list(self.dets), self.noiseCoeff[i], delta)

    
    def cleanModes( self ):
        """
        @brief   Delete from an array of modes all the ones that are zero or l.d. from each other.
        """
        if len(self.modes) == 0:
            psLib.trace('moby', 1, "ERROR: No modes to clean.")
            return
        corr = correlate(self.modes[:1000])
        s = []
        for i in range(len(self.modes)):
            if corr[i][i] == 0.0:
                s.append(i)
            else:
                for j in range(i+1, len(self.modes)):
                    if corr[i][j] == 1.0:
                        s.append(i)
        self.modes = numpy.delete(self.modes, s, 0)



    def correlateModeSet( self, plot = True, forceFit = False ):
        """
        @brief   Find the cross correlation of a set of modes fitted to the whole array.
        @param   plot   set True or False to generate or not a plot of the array correlation.
        @return  Array with the correlation.
        """
        if len(self.modes) == 0:
            psLib.trace('moby', 1, "ERROR: No modes to correlate.")
            return

        if self.coeff == [] or forceFit \
           or numpy.shape(self.coeff)[1] < self.tod.det_uid.size:
            self.fitModes(all = True)
        
        rows, cols, dets = self.tod.listFromRowCol()
        corr = numpy.zeros([self.tod.nrow, self.tod.ncol])
        for i in range(len(dets)):
            mode = numpy.dot(self.coeff[i], self.modes)
            modeNorm = numpy.linalg.norm(mode)
            dataNorm = numpy.linalg.norm(self.tod.data[i])
            corr[rows[i]][cols[i]] = numpy.dot(mode, self.tod.data[i]) / \
                                     modeNorm / dataNorm
        if plot:
            import_pylab()
            pylab.matshow(corr, vmin = -1, vmax = 1)
            pylab.colorbar(shrink=0.8)
            pylab.show()
  
        self.corr = corr
   

    def filterModes(self, filt):
        """
        @brief  Filter modes with a given filter.
        @param  filt  Filter to apply.
        """
        if len(self.modes) == 0:
            psLib.trace('moby', 1, "ERROR: No modes to filter.")
            return

        self.__detrendModes()
        self.modes = moby2.tod.filter.apply_simple(self.modes, filt)
        if filt[0] > 0.5: self.__retrendModes()
        for i in range(len(self.modes)):
            self.modes[i] -= self.modes[i].mean()
            self.modes[i] /= numpy.linalg.norm(self.modes[i])


    def getCoeffMatrix( self, modeIndex ):
        """
        @brief Gets mode coefficient array matrix for mode "modeIndex"
        """
        co = numpy.zeros(1024)
        co[self.dets] = self.coeff.transpose()[modeIndex]
        mat = numpy.reshape(co, [32, 32])
        return mat.transpose()


    def plotCoeff( self, modeIndex, title = None, vmin = None, vmax = None,
                   filename = None, show = False):
        """
        @brief  Plot fit coefficients for a given mode.
        @param  modeIndex: index of the mode for which to plot the coefficients.
        @param  title: title of the plot
        @param  vmin, vmax: minimum and maximum values in the scale.
        @param  filename: outputs to such file.
        @param  show: pops up a plot on the screen.
        """
        import_pylab()
        mat = self.getCoeffMatrix( modeIndex )
        pylab.matshow(mat, vmin = vmin, vmax = vmax)
        pylab.title(title)
        pylab.colorbar(shrink = 0.8)
        pylab.xlabel("Rows")
        pylab.ylabel("Cols")
        if filename is not None: 
            pylab.savefig(filename)
            print("Plotting coefficients")
        if show: pylab.show()
        else: pylab.close()

    def __detrendModes( self, window = 1000 ):
        """
        @brief  Detrends the modes so that they can be filtered
        """
        if ~self.detrended:
            if numpy.ndim(self.modes) == 1:
                ndata = len(self.modes)
                self.ends = numpy.zeros(2)
                self.ends[0] = numpy.median(self.modes[:window])
                self.ends[1] = numpy.median(self.modes[-window:])
                trend = numpy.arange(ndata)*(self.ends[1]-self.ends[0])/(ndata-1)
                self.modes -= trend - self.ends.mean()
            else:
                ndata = len(self.modes[0])
                self.ends = numpy.zeros([len(self.modes), 2])
                x = numpy.arange(ndata)
                for i in range(len(self.modes)):
                    self.ends[i][0] = numpy.median(self.modes[i][:window])
                    self.ends[i][1] = numpy.median(self.modes[i][-window:])
                    trend = x*(self.ends[i][1]-self.ends[i][0])/(ndata-1)
                    self.modes[i] -= trend - self.ends[i].mean()
            self.detrended = True


    def __retrendModes( self ):
        """
        @brief  Recovers the trend of the modes after detrending them  
        """
        if self.detrended:
            if numpy.ndim(self.modes) == 1:
                ndata = len(self.modes)
                trend = numpy.arange(ndata)*(self.ends[1]-self.ends[0])/(ndata-1)
                self.modes += trend - self.ends.mean()
            else:
                ndata = len(self.modes[0])
                x = numpy.arange(ndata)
                for i in range(len(self.modes)):
                    trend = x*(self.ends[i][1]-self.ends[i][0])/(ndata-1)
                    self.modes[i] += trend - self.ends[i].mean()
            self.detrended = False

    def exportModes( self ):
        """
        @brief  Export information from modeSet for storage
        """
        data = {}
        data['modes'] = self.modes
        data['dets'] = self.dets
        data['coeff'] = self.coeff
        data['corr'] = self.corr
        data['removed'] = self.removed
        data['detrended'] = self.detrended
        data["comment"] = self.comment
        return data

    def importModes( self, data ):
        """
        @brief  Import information from modeSet for storage
        @param  data  Information to add in object
        """
        self.modes = data['modes']
        self.dets = data['dets']
        self.coeff = data['coeff']
        self.corr = data['corr']
        self.removed = data['removed']
        self.detrended = data['detrended']
        self.comment = data["comment"]

    def writeToPath( self, path ):
        """
        @brief   Saves the data in the live correlation object to path.
        @param   path  Path where to store the data
        """
        data = self.exportModes()
        data['todName'] = self.tod.name

        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/modeSet.pickle" % path
        f = file( filename, 'w' )
        p = pickle.Pickler( f )
        p.dump( data )
        f.close()

    def readFromPath( path, tod ):
        """
        @brief   Read the data from a live correlation object stored in path.
        @param   tod   Tod for which the correlation object was originally created.
        @param   path  Path where the data was stored.
        """
        dets = numpy.arange(1024)
        ms = modeSet( tod, [], dets )
        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/modeSet.pickle" % path
        f = file( filename )
        p = pickle.Unpickler( f )
        data = p.load()
        assert tod.name == data['todName'], \
            "modeSet: TOD %s and stored modeSet %s don't match" % \
            ( tod.name, data['todName'] )
        ms.importModes( data )
        return ms

    readFromPath = staticmethod(readFromPath)


####################################################################################
#################### POLINOMIAL COMMOM MODE OBJECT #################################
####################################################################################

class polyCommonMode( object ):
    def __init__( self, tod, dets, power ):
        """
        """
        self.tod = tod
        self.dets = dets
        self.power = power

        A = self.__generatePoly()
        AA = numpy.dot(A, A.transpose())
        b = numpy.dot(A, self.tod.data[dets])
        self.modes =  numpy.linalg.solve(AA,b)


    def removeModes( self, dets = None):
        """
        """
        if dets is None: dets = self.dets

        A = self.__generatePoly(dets = dets)
        self.tod.data[dets] -= numpy.dot(A.transpose(), self.modes)
        self.tod.abuses += [{'name':'removeMode', 'type':'modeRemoval', 'detectors': dets}]


    def __generatePoly( self, dets = None ):
        """
        """
        if dets is None: dets = self.dets

        xx = (self.tod.rows - self.tod.rows.mean())
        yy = (self.tod.cols - self.tod.cols.mean())
        xx = xx[dets]
        yy = yy[dets]
        n = (self.power+1)*(self.power+2)/2
        A = numpy.empty([n,len(dets)])
        k = 0
        for p in range(0,self.power+1):
            for i in range(0,p+1):
                A[k] = xx**(p-i) * yy**i
                k += 1
        return A


    def filterModes(self, filt):
        """
        @brief  Filter modes with a given filter.
        @param  filt  Filter to apply.
        """
        ndata = len(self.modes[0])
        end0 = numpy.median(self.modes[0][:window])
        end1 = numpy.median(self.modes[0][-window:])
        trend = numpy.arange(ndata)*(end1-end0)/(ndata-1)
        ends_mean = (end0+end1)/2.
        self.modes[0] -= trend - ends_mean

        self.modes = moby2.tod.filter.apply_simple(self.modes, filt)

        self.modes[0] += trend - ends_mean



################## CORRELATIONS ##############################

def correlate( A ):
    """
    @brief   Correlate the vectors in an array

    @param   A   Matrix that contain the vectors to correlate.

    @return  Correlation matrix.
    """
    mat = numpy.cov(A)
    covDiag = numpy.sqrt(numpy.diag(mat))
    for i in range(len(A)):
        mat[i][i] /= covDiag[i]*covDiag[i]
        for j in range(i+1, len(A)):
            mat[i][j] /= covDiag[i]*covDiag[j]
            mat[j][i] = mat[i][j]
    return mat


def correlateMode( tod, mode):
    """
    @brief   Find the cross correlation of a mode fitted to the whole array.
    @param   tod    tod with data to correlate.
    @param   mode   mode to correlate with tod.
    @return  vector with the correlation.
    """
    assert len(mode == tod.nsamps)

    modeNorm = numpy.linalg.norm(mode)
    corr = numpy.zeros(tod.det_uid.size)
    for i in range(tod.det_uid.size):
        dataNorm = numpy.linalg.norm(tod.data[i])
        if dataNorm != 0.0:
            corr[i] = numpy.dot(mode, tod.data[i]) / modeNorm / dataNorm
        else:
            corr[i] = 0.0
    return corr

def relativeCalibration( tod, apply = True ):
    """
    @brief  Relative calibration performed using common mode without detrending
    """
    dets = tod.cuts.get_uncut()

    data1 = tod.data.copy()
    filt = filters.lowPassButterworth(tod, fc = 0.5, order = 2)
    for d in dets:
        tod.data[d] = moby2.tod.filter.apply_simple(tod.data[d], filt)
    cm = correlationObject(tod).findCommonMode()
    cm.fitModes()
    tod.data[:] = data1[:]
    del data1
    cal = numpy.ones(tod.det_uid.size)
    cal[dets] /= cm.coeff/cm.coeff.mean()
    if apply: 
        for i in dets: tod.data[i] *= cal[i]
    return cal
