from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np
from scipy.special import j1
from time import time

import moby2

class FitParams(dict):
    """
    Manage fitting of multi-parameter beam models by providing a dictionary-
    like class that also supports ordering of selected keys.  Intended for use
    in scipy.optimize.leastsq.
    """
    def __init__(self, keys, format_string):
        dict.__init__(self)
        # These are the fit parameters, in an order.
        self.keys = [k for k in keys]
        for k in keys:
            self[k] = 0.
        self.format_string = format_string

    def pretty(self):
        # .format does not like numpy scalars
        vals = dict([(k, float(self[k])) for k in self.keys])
        return self.format_string.format(**vals)

    def tuple(self):
        """
        Code the fit parameters into an ordered set.
        """
        return tuple([self[k] for k in self.keys])
    
    def store(self, p):
        """
        Extract fit parameters from tuple or dictionary 'p' and store.
        """
        if isinstance(p, dict):
            for k in self.keys:
                if k in p:
                    self[k] = p[k]
        elif hasattr(p, '__getitem__'):
            for i,k in enumerate(self.keys):
                self[k] = p[i]
        else:
            raise TypeError('fit_param can\'t store things of type %s'%type(p))


class BeamModel:
    keys = ['dx', 'dy', 'h', 'w', 'base']
    format_string = 'x={dx:8.5f} y={dy:8.5f} h={h:.4e} w={w:.3e}'

    def __call__(self, x, y, params):
        return self.get_spatial(x, y, params)

    def init_params(self, p=None, opts={}):
        p0 = FitParams(self.keys, self.format_string)
        if not p is None:
            p0.update(p)

        # Sensible defaults.
        if p0['w'] <= 0.:
            p0['w'] = opts.get('fwhm', 1.) * (1./60)*np.pi/180

        if p0['h'] <= 0:
            p0['h'] = opts.get('amp0', 1.)

        return p0

    def get_fit_flag_sets(self):
        return [[True,  True,  False, False, True ],
                [False, False, True,  True,  False],
                [True,  True,  True,  True,  True ]]

    @staticmethod
    def get_model(model_params):
        name = model_params['name']
        if name == 'airy':
            return AiryModel()
        if name == 'gauss':
            return GaussModel()
        raise ValueError("model name %s is not known" % name)
                            
class GaussModel(BeamModel):
    def get_spatial(self, x, y, params):
        x0, y0, a, fwhm, base = params[:5]
        sig2 = fwhm*abs(fwhm) / 8 / np.log(2)
        return base + a * np.exp(-((x-x0)**2+(y-y0)**2)/2/sig2)

class AiryModel(BeamModel):
    def get_spatial(self, x, y, params):
        x0, y0, a, fwhm, base = params[:5]
        arg = (3.23268/fwhm) * ((x-x0)**2+(y-y0)**2)**.5
        arg[arg<1e-7] = 1e-7
        return base + a * (2*j1(arg)/arg)**2

