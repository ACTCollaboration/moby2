from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline as spline1d

from . import blackbody as BB

# There will be something available.  Note, before you add pyephem,
# that historically it's had a weird value for Uranus' radius.  So
# check it.
ephem_modules = []
try:
    import aephem
    ephem_modules.append('aephem')
except ImportError:
    pass


def get_planet_solid_angle(planet, ctime):
    """
    Returns the solid angle, in steradians, of the planet
    (e.g. 'uranus') at the specified ctime(s), as seen from the center
    of the Earth, which is good enough.
    """
    if len(ephem_modules) == 0:
        raise RuntimeError("No ephem modules could be loaded.")
    name = ephem_modules[0]

    if name == 'aephem':
        if np.asarray(ctime).ndim > 0:
            return np.array([get_planet_solid_angle(planet, t) for t in ctime])
        jd = aephem.ctime_to_jd(ctime)
        obs = aephem.p_orb['earth']
        obj = aephem.p_orb[planet]
        phys = aephem.p_phys[planet]
        return aephem.disc_solid_angle(jd, obs, obj, phys)


class BrightnessModel(object):
    valid_freqs = [0, 0]
    source = None

    def get_T_brightness(self, f_GHz):
        raise RuntimeError("Method not defined.")

    def _get_spectral_radiance_minus_background(self, f_GHz):
        T = self.get_T_brightness(f_GHz)
        spec_rad = BB.blackbody(T, f_GHz)
        spec_rad -= BB.blackbody(BB.T_cmb, f_GHz)
        return spec_rad
        
    def get_RJ_minus_background(self, f_GHz):
        spec_rad = self._get_spectral_radiance_minus_background(f_GHz)
        return BB.spectrumToRJ(spec_rad, f_GHz)
    
    def get_uK_minus_background(self, f_GHz):
        spec_rad = self._get_spectral_radiance_minus_background(f_GHz)
        return spec_rad / BB.DCMB_factor(f_GHz) * 1e6

class LogLambdaModel(BrightnessModel):
    p = None
    def get_T_brightness(self, f_GHz):
        if (np.any(f_GHz < self.valid_freqs[0]) or
            np.any(f_GHz >= self.valid_freqs[1])):
            raise ValueError('Model queried for planet temperatures '\
                'outside valid range of (%f,%f) GHz' % tuple(self.valid_freqs))
        phi = np.log10(f_GHz) - 2.
        return np.polyval(self.p[::-1], phi)

class UranusHasselfield2013(LogLambdaModel):
    valid_freqs = [50, 1000]
    p = (121, -78.2, 18.2)
    source = 'uranus'

class BrightnessFunction(BrightnessModel):
    func = None
    valid_freqs = [0,0]
    def get_T_brightness(self, f_GHz):
        if (np.any(f_GHz < self.valid_freqs[0]) or
            np.any(f_GHz >= self.valid_freqs[1])):
            raise ValueError('Model queried for planet temperatures '\
                'outside valid range of (%f,%f) GHz' % tuple(self.valid_freqs))
        return self.func(f_GHz)

class ESABrightnessSpline(BrightnessFunction):
    """
    Support the loading of Planck / ESA planet brightness temperature
    models.
    """
    filename = 'something.fits'
    columns = ['wave', 'T_b']
    def __init__(self):
        cfg = moby2.user_cfg
        assert 'planet_data' in cfg and 'ESA_models' in cfg['planet_data']
        import os
        full_path = os.path.join(cfg['planet_data']['ESA_models'], self.filename)
        data = moby2.util.StructDB.from_fits_table(full_path)
        x, y = [data[c] for c in self.columns]
        self.valid_freqs = [x[0], x[-1]]
        self.func = spline1d(x, y)

class UranusESA2(ESABrightnessSpline):
    filename = 'uranus_v2.fits.gz'
    source = 'uranus'
    
class NeptuneESA3(ESABrightnessSpline):
    filename = 'neptune_v3.fits.gz'
    source = 'neptune'
    
class NeptuneESA5(ESABrightnessSpline):
    filename = 'neptune_v5.fits.gz'
    source = 'neptune'

class NeptuneHybrid(BrightnessFunction):
    """
    This is a hacky correction to the ESA5 model to match Planck HFI
    measurements.  It only does so at the 2% level or so.  This is as
    good as it gets for now.
    """
    source = 'neptune'
    def __init__(self):
        self.a = NeptuneESA5()
        self.b = LogLambdaModel()
        self.b.p = [0.94525782, 0.15269391]
        self.b.valid_freqs = self.a.valid_freqs
    def get_T_brightness(self, f_GHz):
        return self.b.get_T_brightness(f_GHz) * self.a.get_T_brightness(f_GHz)

def get_brightness_model(source_name, version=None):
    if source_name == 'uranus':
        return UranusHasselfield2013()
    if source_name == 'neptune':
        return NeptuneHybrid()
    return None

