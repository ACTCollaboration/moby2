from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy
import numpy as np

import moby2
import moby2.util.log as psLib

from .pcg import PCGFilter

#from moby.mapMaking import pcg
#from moby.todAnalysis import correlations, mbFilter, noise, filters
#from moby.todUtils import TOD, cuts
#from moby.database import result
#from moby import constants, rethread


class identityFilter( object ):
    """
    @brief A class to describe our pcg time domain filters
    """

    def __init__( self, noiseTOD, filterNoise = True ):
        """
        @param noiseTOD TOD.TOD with best estimate of noise
        @param filterNoise bool if False just perform decorrelation
        """
        pass

    def setTOD( self, tod ):
        """
        @brief assign this tod as the tod to be filtered
        """
        self.tod = tod

    def applyFilter(self, data=None):
        psLib.trace("moby", 6, "Indentity Filter")

class simpleFilter( PCGFilter ):

    def __init__( self, noiseTOD, filterParams):
        """
        @param noiseTOD TOD.TOD with best estimate of noise
        @param filterParams  Dictionary with parameters to use in filter and correlations, they are:
                commonModeSmoothing: smooth commom mode correlation with boxcar of this width (seconds)
                filterNoise: True or False. Whether to apply the inverse noise filter.
                noiseHighPassFc: high pass filter data with this cut-off frequency (Hz)
        """
        self.params = {'commonModeSmoothing': 0.25,
                       'filterNoise': True,
                       'noiseHighPassFc': 1.0}
        self.params.update(filterParams)
        self.dt = noiseTOD.ctime[1] - noiseTOD.ctime[0]
        self.uncut_dets = noiseTOD.cuts.get_uncut()
        
    def setTOD( self, tod ):
        """
        @brief assign this tod as the tod to be filtered
        """
        self.tod = tod

    def removeCorrelations(self, data=None):
        psLib.trace( "moby", 3, "simpleFilter.removeCorrelations begins")
        if data is None:
            data = self.tod.data
        self.filter=True
        psLib.trace( "moby", 3, "get uncut")
        psLib.trace( "moby", 3, "remove mean")
        moby2.tod.remove_mean(data=data, dets=self.uncut_dets)
        psLib.trace( "moby", 3, "computing common mode")
        cm = moby2.libactpol.data_mean_axis0(
            data, self.uncut_dets.astype('int32'))
        psLib.trace( "moby", 3, "create filter")
        cmfilt = self.params.get('commonModeFilter', 'boxcar')
        if cmfilt == 'boxcar':
            fc = 1./self.params['commonModeSmoothing']
            nw = int(self.params['commonModeSmoothing'] / self.dt)
            # Mimic behaviour of windowFilter by using only odd window widths.
            nw = nw//2*2 + 1
            filt = moby2.tod.filters2.boxcar(width_samples=nw) \
                   .getExecFor(data=data, dt=self.dt)
        elif cmfilt == 'butterworth-low':
            fc = 1./self.params['commonModeSmoothing']
            filt = moby2.tod.filters2.lowPassButterworth(fc=fc, order=4) \
                .getExecFor(data=data, dt=self.dt)
        psLib.trace( "moby", 3, "filtering common mode")
        cm = filt.applyTime(cm)
        psLib.trace( "moby", 3, "removing common mode")
        modes = cm.reshape((1,-1)).astype('float32')
        amps = numpy.ones((1,len(self.uncut_dets)), dtype='float64')
        moby2.libactpol.remove_modes(
            data, self.uncut_dets.astype('int32'), modes, amps)
        psLib.trace( "moby", 3, "simpleFilter.removeCorrelations ends")

    def noiseFilter(self, data=None):
        psLib.trace( "moby", 3, "simpleFilter.noiseFilter begins")
        if data is None:
            data = tod.data
        #d = self.tod.cuts.get_uncut()
        moby2.tod.remove_mean(data=data)
        moby2.tod.detrend_tod(data=data)
        hpf = moby2.tod.TODFilter()
        hpf.add(moby2.tod.filters.highPassButterworth,
                {'fc':self.params['noiseHighPassFc']})
        hpf.apply(self.tod)


class simpleMaskingFilter( simpleFilter ):
    """
    Like simpleFilter, but computes common mode only after masking
    data indicated in params['CMCutsFile'].
    """
    
    def removeCorrelations(self, data=None):
        self.filter=True
        tod = self.tod
        # Save original data and cuts
        d,r,c = tod.listUncut()
        tod_copy = tod.data.copy()
        cuts_ref = tod.cuts
        tod.cuts = cuts.cuts(tod.nrow, tod.ncol)
        # Get and fill cuts
        cuts.clearCuts(tod)
        cuts.cutFromFile(tod, self.params['CMCutsFile'])
        cuts.fillCuts(tod)
        # Compute common mode
        TOD.removeMean(tod)
        cm = numpy.mean(tod.data[d], axis = 0)
        # restore unmasked data and original cuts
        tod.data[:,:] = tod_copy[:,:]
        del tod_copy
        tod.setCuts(cuts_ref)
        # remove CM
        tod.data -= cm



class jawsFilter( PCGFilter ):

    def __init__( self, noiseTOD, filterParams = None):
        """
        @param noiseTOD TOD.TOD with best estimate of noise
        @param filterParams  Dictionary with parameters to use in filter 
                filterNoise: True or False. Whether to apply the inverse noise filter.
                noiseHighPassFc: high pass filter data with this cut-off frequency (Hz)
        """
        self.params = { 'filterNoise': False, \
                'useNoiseInCorr' : False, \
                'filterParams' : { 'fc':0.5, 'df':0.1 }}
        if filterParams is not None:
            self.params['filterParams'] = filterParams

    def setTOD( self, tod ):
        """
        @brief assign this tod as the tod to be filtered
        """
        self.tod = tod

    def removeCorrelations(self):
        uncutDets,uncutRows,uncutCols = self.tod.listUncut()
        cm = numpy.zeros(self.tod.ndata)
        for det in uncutDets:   # BECAREFUL HERE IF YOU USE ARRAY SLICES IT DOUBLES THE MEMORY!
            cm += self.tod.data[det]
        cm /= len(uncutDets)
        for det in uncutDets:
            self.tod.data[det] -= cm

    def noiseFilter(self):
        TOD.removeMean( self.tod )
        TOD.detrendData(self.tod)
        filt = mbFilter.mbFilter()
        filt.add(mbFilter.sine2highPass, self.params['filterParams'])
        filt.apply(self.tod)

class deluxeFilter( PCGFilter ):
    """
    A class to describe our noise filters
    """
    def __init__( self, noiseTOD, filterParams ):
        """
        @param noiseTOD TOD.TOD containing the best estimate of the noise, will be altered!
        @param filterParams  dictionary with parameters for filters and correlations, they are:
                filterNoise: True or False. Whether to apply an inverse noise filter.
                noiseParams: Dictionary with parameters of the inverse noise filter.
                useNoiseInCorr: True or False. Whether to apply the (m'N^{-1}m)^{-1} term for mode removal.
                removeLiveCorr: True or False. Whether to remove live correlation modes.
                removeArrayCorr: True or False. Whether to remove live common mode.
                arrayCorrParams: Dictionary with parameters defining the kind of common mode to use.
                removeColCorr: True or False. Whether to remove live column correlation modes.
                colCorrParams: Dictionary with parameters for the column correlations.
                removeRowCorr: True or False. Whether to remove live row correlation modes.
                rowCorrParams: Dictionary with parameters for the row correlations.
                removeSVDCorr: True or False. Whether to remove live SVD correlation modes.
                SVDCorrParams: Dictionary with parameters for the SVD correlations.
                removeDarkCorr: True or False. Whether to remove live dark correlation modes.
                darkCorrParams: Dictionary with parameters for the dark correlations.
                loadExternal: True or False. Whether to read in the noise and correlations from a depot 
                              specified in the params
                corrDepot: if loadExternal, look in this dir (under constants.CORR_DEPOT_ROOT) for 
                           correlation objects
                corrTag: location in corrDepot for correlation objects
                noiseTag: location in the noiseDepot for noise objects, yo
        """
        self.params = filterParams
        lc = None; dc = None
        if self.params['removeLiveCorr']:
            if self.params['loadExternal']:
                wd = result.setDepot( constants.CORR_DEPOT_ROOT+self.params['corrDepot'] )
                lc = result.readTODResult( correlations.liveCorrObj, noiseTOD, self.params['corrTag'] )
                result.setDepot( wd )
            else:
                lc = correlations.liveCorrObj(noiseTOD)

        if self.params['removeDarkCorr']:
            if self.params['loadExternal']:
                wd = result.setDepot(  constants.CORR_DEPOT_ROOT+self.params['corrDepot'] )
                dc = result.readTODResult( correlations.darkCorrObj, noiseTOD, self.params['corrTag'] )
                result.setDepot( wd )
            else: 
                dc = correlations.darkCorrObj(noiseTOD)

        if not self.params['loadExternal']:
            self.findAndRemoveCorrelations( lc, dc )
        if self.params['filterNoise']:
            if self.params['loadExternal']:
                wd = result.setDepot(  constants.NOISE_DEPOT_ROOT+self.params['noiseDepot'] )
                self.noise = result.readTODResult( noise.noise, noiseTOD, self.params['noiseTag'] )
                result.setDepot( wd )
                psLib.trace( "moby", 3, "Read old noise")
            else:
                self.noise = noise.noise( noiseTOD, **self.params['noiseParams'] )
                self.noise.inverseFilterTOD( useCutDet = True )
                psLib.trace( "moby", 3, "Found new noise")
            psLib.trace( "moby", 2, "Initialized noise for %s" % noiseTOD.name)
        else:
            self.noise = None

        self.condenseModes( noiseTOD, lc, dc )
        if self.params['useNoiseInCorr'] and self.params['filterNoise']:
            psLib.trace( "moby", 2, "Obtaining noise coefficient for %s" % noiseTOD.name)
            if len(self.modes.modes) > 0:
                self.modes.generateNoiseCoeff( self.noise )

    def setTOD( self, tod ):
        self.tod = tod
        if self.params['filterNoise']:
            self.noise.setTOD(tod)
        self.modes.tod = tod
        if self.polyCommonMode is not None:
            self.polyCommonMode.tod = tod

    def noiseFilter( self ):
        TOD.removeMean(self.tod)
        TOD.detrendData(self.tod, force = True)
        self.noise.inverseFilterTOD(level = self.params['filterLevel'])

    def condenseModes( self, tod, lc, dc ):
        psLib.trace( "moby", 3, "Condensing Modes." )
        dets, rows, cols = tod.listUncut()
        self.modes = correlations.modeSet(tod, [], numpy.array(dets))
        self.polyCommonMode = None
        if self.params['removeLiveCorr']:

            if self.params['removeArrayCorr']:
                if self.params['arrayCorrParams']['cmtype'] == 'poly':
                    self.polyCommonMode = lc.commonMode
                else: self.modes.append(lc.commonMode.modes)

            if self.params['removeColCorr']:
                self.modes.append(lc.liveColModes.modes)

            if self.params['removeRowCorr']:
                self.modes.append(lc.liveRowModes.modes)

            if self.params['removeSVDCorr']:
                self.modes.append(lc.liveSVDModes.modes)

        if self.params['removeDarkCorr']:
            self.modes.append(dc.darkModes.modes, orthogonalize = True)
        

    def findAndRemoveCorrelations( self, lc, dc ):
        """
        @brief Find and Remove detector-detector correlations from the TOD. Only to initialize modes.
        """
        #Remove correlations
        if self.params['removeLiveCorr']:
            psLib.trace( "moby", 3, "Removing Live Detector Correlations." )

            if self.params['removeArrayCorr']:
                psLib.trace( "moby", 3, "Removing Live CommonMode, params %s" % \
                    str(self.params['arrayCorrParams']) )
                lc.removeArrayCorr( **self.params['arrayCorrParams'] )

            if self.params['removeColCorr']:
                psLib.trace( "moby", 3, "Removing Column Correlations, params %s" % \
                    str(self.params['colCorrParams']) )
                lc.removeColCorr( **self.params['colCorrParams'] )

            if self.params['removeRowCorr']:
                psLib.trace( "moby", 3, "Removing Row Correlations, params %s" % \
                    str(self.params['rowCorrParams']) )
                lc.removeRowCorr( **self.params['rowCorrParams'] )

            if self.params['removeSVDCorr']:
                psLib.trace( "moby", 3, "Removing Live SVD Modes, params %s" % \
                    str(self.params['SVDCorrParams']) )
                lc.removeSVDCorr( **self.params['SVDCorrParams'] )

        if self.params['removeDarkCorr']:
            psLib.trace("moby", 3, "Removing Dark Detector Correlation, params %s" % \
                        str(self.params['darkCorrParams']))
            dc.removeDarkCorr( **self.params['darkCorrParams'] )

    def removeCorrelations( self ):
        """
        @brief Remove detector-detector correlations from the TOD.
        """
        psLib.trace( "moby", 3, "Removing Correlations." )
        if self.params['useNoiseInCorr'] and self.params['filterNoise']: noise = self.noise
        else: noise = None
        if self.params['removeLiveCorr']:
            if self.params['removeArrayCorr'] and self.params['arrayCorrParams']['cmtype'] == 'poly':
                if self.params['useNoiseInCorr']:
                    psLib.trace( "moby", 2, "WARNING: PolyCommonMode not implemented with noise." )
                self.polyCommonMode.removeModes()
        if len(self.modes.modes) > 0:
            self.modes.removeModeSet( noise = noise, forceFit = True )


# Export...
PCG_filter_list = {
    'simpleFilter': simpleFilter,
    'identity': identityFilter,
}
    
