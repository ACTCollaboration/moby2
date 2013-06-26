from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
This breaks all the rules.  Mostly by being very ACT specific.
"""

import os
import numpy as np
import pickle

import moby2

from moby2.util.mce import MCERunfile
from moby2.util.log import logger
from .stuff import _SimpleDetData

class IVCalibration(_SimpleDetData):
    fields = ['cut_recommendation', 'dac_to_pW', 'fraction_Rnormal']
    runfile = None
    ar_name = None

    @classmethod
    def for_runfile(cls, runfile, instrument = None):
        self = cls()
        if isinstance(runfile, basestring):
            runfile = MCERunfile(runfile)
        self.runfile = runfile
        ar_name = runfile.Item('FRAMEACQ', 'ARRAY_ID')[0].strip().lower()
        ctime = runfile.Item('FRAMEACQ', 'CTIME')
        season = int( (ctime - 949381200.) / 31557600.0 + 2000 )
        ar_data = moby2.scripting.products.get_array_data(
                        {"season":season, "array":ar_name, "instrument": instrument})

        if not 'IV' in self.runfile.data:
            raise RuntimeError("Runfile %s does not contain IV block." % \
                runfile.filename)

        cut_rec = self.get_iv_array('cut_rec_C%i', 'int').astype('bool')
        if len(cut_rec.shape) < 2 or cut_rec.shape[1] != 32:
            logger.trace(1, "IV for weird channels, det_uid might be messed up.")
        # Get the det_uid for these rows and columns
        n_row, n_col = cut_rec.shape
        rows, cols = np.ix_(np.arange(n_row), np.arange(n_col))
        rows, cols = (rows+0*cols).ravel(), (cols+0*rows).ravel()
        self.det_uid = ar_data.select_inner(
            {'row': rows, 'col': cols})
        
        # For application to TOD you need to account for DC filter gain.
        #filt_gain = runfile.ReadoutFilter().gain()

        self.cut_recommendation = cut_rec.ravel()
        self.dac_to_pW = self.get_responsivities().ravel() #/ filt_gain
        nan_bad = np.isnan(self.dac_to_pW)
        self.dac_to_pW[nan_bad] = 0.
        n_weird = (nan_bad * ~self.cut_recommendation).sum()
        if n_weird > 0:
            logger.trace(1, "%i bad calibrations are not cut_recommended." % n_weird)
        self.fraction_Rnormal = self.get_percentage_rn().ravel()
        return self

    @classmethod
    def for_tod(cls, tod_info=None, runfile=None):
        self = cls()
        if runfile is None:
            runfile = tod_info.runfile
        if isinstance(runfile, basestring):
            runfile = MCERunfile(runfile)
        self.runfile = runfile
        self.ar_name = tod_info.array

        if not 'IV' in self.runfile.data:
            raise RuntimeError("Runfile %s does not contain IV block." % \
                runfile.filename)

        cut_rec = self.get_iv_array('cut_rec_C%i', 'int').astype('bool')
        if len(cut_rec.shape) < 2 or cut_rec.shape[1] != 32:
            logger.trace(1, "IV for weird channels, det_uid might be messed up.")
        # Get the det_uid for these rows and columns
        n_row, n_col = cut_rec.shape
        rows, cols = np.ix_(np.arange(n_row), np.arange(n_col))
        rows, cols = (rows+0*cols).ravel(), (cols+0*rows).ravel()
        self.det_uid = tod_info.array_data.select_inner(
            {'row': rows, 'col': cols})
        
        # DEPRECATED For application to TOD you need to account for DC filter gain.
        #filt_gain = runfile.ReadoutFilter().gain()
        #self.dac_to_pW = self.get_responsivities().ravel() / filt_gain
        self.dac_to_pW = self.get_responsivities().ravel()

        self.cut_recommendation = cut_rec.ravel()
        nan_bad = np.isnan(self.dac_to_pW)
        self.dac_to_pW[nan_bad] = 0.
        n_weird = (nan_bad * ~self.cut_recommendation).sum()
        if n_weird > 0:
            logger.trace(1, "%i bad calibrations are not cut_recommended." % n_weird)
        self.fraction_Rnormal = self.get_percentage_rn().ravel()
        return self

    # Quick interfaces to 2D runfile data...

    def get_iv_array(self, key, type='float'):
        # pass type as a string...
        return np.array(self.runfile.Item2d('IV', key, type=type)).transpose()

    def get_squid_array(self, key, type='float'):
        # pass type as a string...
        return np.array(self.runfile.Item2d('SQUID', key, type=type)).transpose()

    # Old moby functions.

    def get_bad_pixels(self):
        cut_array = self.get_iv_array('cut_rec_C%i', type='int').astype('bool')
        return list(zip(*cut_array.nonzero()))

    def get_good_pixels(self):
        cut_array = self.get_iv_array('cut_rec_C%i', type='int').astype('bool')
        return list(zip(*(~cut_array).nonzero()))

    def get_raw_responsivities(self):
        return self.get_iv_array('Responsivity(W/DACfb)_C%i', 'float')

    def get_responsivities(self):
        # Convert to pW per unfiltered FB
        return self.get_raw_responsivities() #* 1e12 Used to be converted to pW...

    def get_squid_vphi_p2p(self):
        return self.get_squid_array('Col%i_squid_vphi_p2p')

    def get_squid_lockrange(self):
        return self.get_squid_array('Col%i_squid_lockrange')

    def get_squid_lockslope_down(self):
        return self.get_squid_array('Col%i_squid_lockslope_down')

    def get_squid_lockslope_up(self):
        return self.get_squid_array('Col%i_squid_lockslope_up')

    def get_squid_multilock(self):
        return self.get_squid_array('Col%i_squid_multilock')

    def get_percentage_rn(self):
        return self.get_iv_array('Percentage_Rn_C%i')

    def get_biaspower(self):
        return self.get_iv_array('Bias_Power(pW)_C%i')

    def get_biasvoltage(self):
        return self.get_iv_array('Bias_Voltage(uV)_C%i')

    def get_bias(self,line):
        return self.runfile.Item('HEADER', 'RB tes bias', type='float')[line-1]

    def get_mce_mode(self):
        return self.runfile.Item('HEADER', 'RB rc1 data_mode', type='int')[0]


class ACTIVCalibration(IVCalibration):
    correctForAR12008 = False

    @classmethod
    def for_tod(cls, tod_info=None, runfile=None):
        self = cls()
        if runfile is None:
            runfile = tod_info.runfile
        if isinstance(runfile, basestring):
            runfile = MCERunfile(runfile)
        self.runfile = runfile
        self.ar_name = tod_info.array

        cut_rec = self.get_iv_array('cut_rec_C%i', 'int').astype('bool')
        if cut_rec.shape[1] != 32:
            logger.trace(1, "IV for weird channels, det_uid might be messed up.")
        # Get the det_uid for these rows and columns
        n_row, n_col = cut_rec.shape
        rows, cols = np.ix_(np.arange(n_row), np.arange(n_col))
        rows, cols = (rows+0*cols).ravel(), (cols+0*rows).ravel()
        self.det_uid = tod_info.array_data.select_inner(
            {'row': rows, 'col': cols})
        
        self.cut_recommendation = cut_rec.ravel()
        self.dac_to_pW = self.get_responsivities().ravel()
        nan_bad = np.isnan(self.dac_to_pW)
        self.dac_to_pW[nan_bad] = 0.
        n_weird = (nan_bad * ~self.cut_recommendation).sum()
        if n_weird > 0:
            logger.trace(1, "%i bad calibrations are not cut_recommended." % n_weird)
        self.fraction_Rnormal = self.get_percentage_rn().ravel()

        ## Special ACT handling
        # Season 2008 mbac145 runfile data is wrong.
        if self.ar_name == 'mbac145':
            ctime = tod_info.ctime
            self.correctForAR12008 = (ctime > 1217516400) and (ctime < 1233500400)

        return self

    def get_biaspower(self):
        bp = self.get_iv_array('Bias_Power(pW)_C%i')
        if self.correctForAR12008==True:
            return bp*self.get_shunt_R()/.0007
        else:
            return bp

    def get_biasvoltage(self):
        bv = self.get_iv_array('Bias_Voltage(uV)_C%i')
        if self.correctForAR12008==True:
            return bv*self.get_shunt_R()/.0007
        else:
            return bv

    def get_responsivities(self):
        """
        These responsivities are derived the bias power and voltage together
        with the shunt resistances from the detectors module.
        Returns in units of pW/DAC
        """
        biasVoltage = self.get_biasvoltage() * 1e-6
        biasPower = self.get_biaspower() * 1e-12
        shuntArray = self.get_shunt_R()
        if shuntArray is None:
            logger.trace(0, "no shunts for array %s" % str(self.ar_name))
            return None
        denom = biasVoltage + 1*(biasVoltage==0)
        dPdI = biasVoltage - shuntArray*biasPower/denom
        dPdI[dPdI<0] = 0.
        return dPdI * self.get_amps_per_dac() * 1e12 #last factor converts to pW

    def get_amps_per_dac(self):
        return {'mbac145': 8.3033911719131563e-13,
                'mbac215': 8.3033911719131563e-13,
                'mbac280': 8.2963801464372938e-13}.get(self.ar_name, None)
        
    def get_shunt_R(self):
        from moby2.instruments import mbac
        return mbac.get_ACT_shunt_R(self.ar_name)



class ResponsivityObject( object ):
    def __init__(self, tod, responsivity = None):
        """
        @brief   Object used to store and retrieve stored responsivities.
        """
        self.tod = tod
        self.Nrows = np.unique(self.tod.info.array_data['row']).size
        self.Ncols = np.unique(self.tod.info.array_data['col']).size
        if responsivity is not None: self.resp = np.array(responsivity)
        else: self.resp = -1*np.ones([self.Nrows,self.Ncols])

    def setResponsivity( self, resp ):
        """
        @brief Set the responsivity array
        """
        self.resp = np.array(resp)

    def writeToPath( self, path ):
        """
        @brief Stores the responsivity object in a given Path
        """
        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/responsivities.pickle" % path
        f = file( filename, 'w' )
        p = pickle.Pickler( f, 2 )
        p.dump( self.resp )
        f.close()

    @staticmethod
    def readFromPath( path, tod ):
        """
        @brief Reloads a previously stored responsivity object
        """
        rp = responsivityObject(tod)

        if path[-1] == '/':
            path = path[0:-1]
        filename = "%s/responsivities.pickle" % path
        f = file( filename, 'rb' )
        p = pickle.Unpickler( f )
        rp.resp = p.load()

        return rp
