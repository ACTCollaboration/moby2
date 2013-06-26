from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

DEG = np.pi/180

def get_sim_tod(params, wand=None, tod=None):
    """
    Prepare a simulated source of some kind.  See code for how to set
    up parameters.

    For forwards compatibility, call parameters by keyword.

    Returns (cuts, data).  The simulated data are only defined for the
    regions marked in cuts.
    """
    if params['type'] == 'pol_swatch':
        # Produce a simple polarized gaussian source.
        sigma = params['sigma_arcmin'] / 60 * DEG
        minus_half_sigma2 = -0.5/sigma**2
        sig_func = lambda r2: np.exp(r2*minus_half_sigma2) * params['amplitude']
        mask_R = sigma * 3
        x, y = wand.get_coords(tod.fplane)
        cut_vs, sim_dat = [], []
        assert np.all(tod.fplane.det_uid == tod.det_uid)
        for i in range(len(tod.data)):
            # phi increases CCW from W.
            phi = tod.fplane.phi[i]   # in radians.
            gamma_iau = phi - np.pi/2 # CCW from N.  As recommended.
            c2, s2 = np.cos(2*gamma_iau), np.sin(2*gamma_iau)
            mask = np.zeros(len(tod.data[0]), bool)
            dat = np.zeros(len(tod.data[0]), tod.data[0].dtype)
            for (x0,y0,I,Q,U) in params['sources']:
                r2 = (x[i]-x0*DEG)**2 + (y[i]-y0*DEG)**2
                s = r2 < mask_R**2
                sig = sig_func(r2[s])
                mask += s
                dat[s] += sig * (I + Q*c2 + U*s2)
            cut_vs.append(moby2.tod.CutsVector.from_mask(mask))
            sim_dat.append(dat[mask])
        cuts = moby2.TODCuts.for_tod(tod, assign=False)
        cuts.cuts = cut_vs
        return cuts, sim_dat
    raise ValueError("Could not decode sim params %s" % params)

def add_sim_tod(tod, cuts, sim_dat):
    """
    Add the simulation data returned by get_sim_tod into a tod, in place.
    """
    for i in range(len(tod.data)):
        tod.data[i,cuts.cuts[i].get_mask()] += sim_dat[i]
