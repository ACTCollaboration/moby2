from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline as spline1d

from . import blackbody as BB

# The nice thing aephem has, and that other Python ephemeris packages
# appear to lack, is an accounting for oblateness and sub-observer
# latitude in the solid angle computation.  For Uranus, which has a
# modest difference in equatorial and polar radius but a very odd spin
# axis, this affects the solid angle at the 1-2% level.  pyephem
# reports the projected equatorial radius, in degrees, and we need to
# compute an additional correction factor to the naive pi*r**2 disk
# solid angle computation.  In the case that aephem is not present, we
# need to provide that correction somehow.

ephem_modules = []
try:
    import aephem
    ephem_modules.append('aephem')
except ImportError:
    pass
try:
    import ephem
    ephem_modules.append('ephem')
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

    if name == 'ephem':
        act = moby2.ephem.ACTEphem()
        single = np.asarray(ctime).ndim == 0
        if single:
            ctime = [ctime]
        results = []
        for t in ctime:
            act.set_timestamp(t)
            u = act.get_object(planet.capitalize())
            r = u.radius
            omega = np.pi * u.radius**2
            results.append(omega)
        # Corrections?
        if planet in ['uranus', 'jupiter']:
            efactor = solid_angle_correction(np.asarray(ctime), planet)
            results = [om*e for om,e in zip(results, efactor)]
        else:
            raise RuntimeError("Solid angle correction factor not "
                               "available for '%s'" % planet)
        if single:
            return results[0]
        return np.array(results)


def solid_angle_correction(t, target):
    """Returns the approximate sub-observer latitude, in degrees, of
    planetary target at ctime t (which can be an array or scalar).
    The current model gets Uranus and Jupiter to <.1%; Saturn to ~.2%
    (but you'd want to account for the rings anyway).

    """
    lims = [1356998400, 1893456000]
    p = {
        'uranus': [ 3.88233016e+00,  4.81513918e+01, -2.54042514e+00, -1.20158753e+00,
                    -2.46918964e-04,  6.21477912e+00,
                    1577836800, 0.047481307441877396],
        'jupiter': [ 2.21912058e-02, -7.21533904e-02, -2.11716337e+00,  2.79381276e+00,
                     1.40789074e-03,  5.36641642e-01,
                     1577836800, 0.1435630105925818],
        'saturn': [ 1.80574507e+00, -8.90907531e+00,  3.67589879e+01, -2.85961149e+01,
                    1.55770490e-03,  1.65458760e-01,
                    1577836800, 0.22899679859188238]
    }[target]
    t0, e = p[-2:]
    m, b, A, B, omega_m, omega_b = p[:-2]
    dt = (t - t0) / (86400*365.24)
    assert np.all(lims[0] <= np.asarray(t)) and np.all(np.asarray(t) <= lims[1])
    omega = omega_m * dt + omega_b
    lat = dt*m + b + A*np.cos(omega*dt) + B*np.sin(omega*dt)
    return (1. + np.sin(lat*np.pi/180)**2 * e)**-.5


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

class JupiterSimple(BrightnessFunction):
    def get_T_brightness(self, f_GHz):
        # Rough wmap, only valid for LF array!!
        p = [ 1.07794677, 111.55898163 ]
        return np.polyval(p, f_GHz)

class JupiterSimple(BrightnessFunction):
    def get_T_brightness(self, f_GHz):
        # Rough wmap, only valid for LF array!!
        p = [ 1.07794677, 111.55898163 ]
        return np.polyval(p, f_GHz)

class JupiterPlanck2013(BrightnessFunction):
    """This is the expression from "Planck 2013 results V - LFI
    calibration".  Note the caption says 187.8K - 4.5 K cm^-1 lambda,
    but the plot has a slope of -4.5 K mm^-1 so the caption must be a
    typo.

    This is only shown from 100 down to 25 GHz.  The WMAP point at 90
    GHz is substantially above the line so we probably need something
    different outside of the LF bands.

    """
    def get_T_brightness(self, f_GHz):
        lamb = BB.c / f_GHz * (1e-9*1e3)  # in mm
        T_RJ = np.polyval([-4.5, 187.8], lamb)
        return BB.RJToBlackbody(T_RJ, f_GHz)


def get_brightness_model(source_name, version=None):
    if source_name == 'uranus':
        return UranusHasselfield2013()
    if source_name == 'neptune':
        return NeptuneHybrid()
    if source_name == 'jupiter':
        return JupiterPlanck2013()
    return None
