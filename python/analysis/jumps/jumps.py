from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

class JumpCorrector:
    _default_params = {
        'sample_window': 20,
        'max_jumps': 5,
        'max_frac_phi0': .05,
        }
    @staticmethod
    def find_jumps(data, thresh, step):
        """
        Find regions, of with step, where the data jumps by more than
        thresh.  Return a list of such points, centered on the middle
        the jump, roughly.
        """
        dy = data[step:] - data[:-step]
        jump_mask = (abs(dy) > thresh)
        changes = (jump_mask[1:] != jump_mask[:-1]).nonzero()[0] + step//2
        if jump_mask[0]:
            changes = np.hstack((0, changes))
        if jump_mask[-1]:
            changes = np.hstack((changes, len(data)-1))
        jump_locs = (changes[::2] + changes[1::2]) // 2
        return jump_locs

    def __init__(self):
        self.det_uid = None
        self.phi0 = None
        self.data = None
        self.params = self._default_params.copy()

    @classmethod
    def for_tod(cls, tod, runfile=None, filter_gain=None, phi0=None, params={}):
        self = cls()
        self.det_uid = tod.det_uid.copy()
        if runfile is None:
            runfile = tod.info.runfile #moby2.util.mce.MCERunfile.from_dirfile(
                #tod.info.filename, 'Extras/')
        if phi0 is None:
            rows, cols = tod.info.array_data.get_property(
                ['row','col'], det_uid=self.det_uid)
            phi0 = np.transpose(runfile.Item2dRC(
                    'HEADER', 'RB rc%i flx_quanta%%i', type='int'))[rows,cols]
        if filter_gain is None:
            filter_gain = runfile.ReadoutFilter().gain()
        self.filter_gain = filter_gain
        self.phi0 = phi0
        self.data = tod.data
        self.params.update(params)
        return self

    def find_flux_jumps(self, data=None, filter_gain=None, phi0=None,
                        dets=None):
        jumps = []
        nsamps = self.data.shape[-1]
        if dets is None:
            dets = range(self.data.shape[0])
        width = self.params['sample_window']
        for i in dets:
            jump_data = []
            phi0, data = self.phi0[i], self.data[i]
            jump_pos = self.find_jumps(data, phi0*.5*self.filter_gain,
                                       width)
            # Determine the magnitude of each shift
            for j in jump_pos:
                dy = (data[min(nsamps-1,j+width//2)] - 
                      data[max(0,j-width//2)]) / self.filter_gain
                n_phi0 = int(round(dy / phi0))
                precision = dy / phi0 - n_phi0
                jump_data.append((j, n_phi0, precision))
            jumps.append(jump_data)
        self.jumps = jumps
        return jumps

    def find_bad_channels(self, jumps=None):
        if jumps is None:
            jumps = self.jumps
        dets = self._get_dets(None)
        self.bad_channels = []
        for d in dets:
            bad = False
            if len(jumps[d]) > self.params['max_jumps']:
                bad = True
            qual = np.array([j[2] for j in jumps[d]])
            if np.any(abs(qual) > self.params['max_frac_phi0']):
                bad = True
            self.bad_channels.append(bad)
        return self.bad_channels

    def remove_jumps(self, data=None, dets=None, jumps=None):
        if data is None:
            data = self.data
        if jumps is None:
            jumps = self.jumps
        if dets is None:
            dets = range(data.shape[0])
        else:
            dets = np.asarray(dets)
            if dets.dtype == 'bool':
                dets = dets.nonzero()[0]
        for i in dets:
            if self.bad_channels[i] or len(jumps[i]) == 0:
                continue
            n = 0
            for ji in range(len(jumps[i])):
                j, dn, p = jumps[i][ji]
                if ji == len(jumps[i])-1:
                    start, stop = j, data.shape[-1]
                else:
                    start, stop = j, jumps[i][ji+1][0]
                n += dn
                data[i,start:stop] -= self.phi0[i]*n*self.filter_gain
        return data

    def _get_dets(self, dets):
        if dets is None:
            return np.arange(self.data.shape[0])
        dets = np.asarray(dets)
        if dets.dtype == 'bool':
            return dets.nonzero()[0]

    def get_cuts(self, jumps=None, dets=None, kill_unselected=True):
        if jumps is None:
            jumps = self.jumps
        dets = self._get_dets(dets)
        nsamps = self.data.shape[-1]
        width = self.params['sample_window']
        self.cuts = moby2.tod.TODCuts(nsamps=nsamps, det_uid=self.det_uid)
        for d in dets:
            if self.bad_channels[d]:
                self.cuts.set_always_cut(d)
                continue
            regions = [(max(0,j[0]-width//2),min(j[0]+width//2,nsamps-1))
                       for j in jumps[d]]
            self.cuts.add_cuts(d, regions)
        return self.cuts

    def write(self, filename):
        output = moby2.util.MobyDict()
        output['jumps'] = self.jumps
        output['phi0'] = list(self.phi0)
        output['det_uid'] = list(self.det_uid)
        output['bad_channels'] = list(self.bad_channels)
        output.write(path)

    write_to_path = write

