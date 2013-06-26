from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os, pickle

import moby2

from . import beam_models as beamModels
from . import util

# I need these.
from numpy import *
from pylab import *

SIGMA_TO_FWHM = sqrt(8.*log(2.))


class BeamObsList(list):
    def __init__(self, filename=None, obs=None, read=True, loadSource=False):
        list.__init__(self)
        self.filename = filename
        if self.filename is not None and read:
            self.read(loadSource=loadSource)
        if obs is not None:
            self += obs

    def write(self, filename=None, force=False):
        if filename is None: filename = self.filename
        if os.path.exists(filename) and not force:
            raise RuntimeError('output file %s exists, use force=True to overwrite.' % filename)
        output = [x.encode() for x in self]
        pickle.dump(output, open(filename, 'w'))

    def read(self, filename=None, loadSource=False):
        if filename is None: filename = self.filename
        encoded = pickle.load(open(filename))
        for e in encoded:
            b = BeamObs()
            b.decode(e, loadSource=loadSource)
            self.append(b)

    def select(self, mask):
        new_obs = BeamObsList()
        mask = np.asarray(mask)
        if mask.dtype == 'bool':
            if len(mask) != len(self):
                raise ValueError('expected selection mask of same length as BeamObsList')
            mask = mask.nonzero()[0]
        if mask.dtype != 'int':
            raise ValueError('expected selection mask or integer indices into BeamObsList')
        for i in mask:
            new_obs.append(self[i])
        return new_obs
        
    def sort(self, k1=None, k2=None, data=None, reverse=False):
        """
        Sort entries.  Pass either k1 and k2 (as in .get) or a data vector.
        """
        if data is None:
            data = self.get(k1, k2)
        items = [x for x in self]
        while self: self.pop()
        for _,i in sorted(zip(data, items), reverse=reverse):
            self.append(i)

    def get(self, k1, k2=None):
        """
        Create array of values from particular data in each member observation.

        e.g. to get an array of all the gaussian azimuth offsets:
            BeamObsList.get('gaussian', 'dx')
        """
        if k2 is None:
            keys = list(self[0][k1].keys())
            output = {}
            for k in list(self[0][k1].keys()):
                output[k] = self.get(k1, k)
            return output
        return array([o[k1][k2] for o in self])

    @classmethod
    def from_file(cls, filename):
        self = cls()
        self.read(filename)
        return self

def xyz_for_mask(map, mask):
    idx = mask.nonzero()
    return map.x[0, idx[1]], map.y[idx[0], 0], map.data[mask]

class BeamObs(dict):

    def __init__(self, filename=None, mapfile=None, map=None, read=True,
                 loadSource=False, savePlots=None, **kwargs):
        """
        Create a beam analysis object.

        To associate this object with a map, pass in map=... or mapfile=...

        Use filename=... to set the analysis input/output filename, with
        read=True|False.
        """
        dict.__init__(self)
        self['init'] = {
            'source': None,
            'center0': kwargs.get('center', (0.,0.)),
            'radius0': kwargs.get('radius', 0.1),
            }
        self['flags'] = {}
        self.map = None
        if map is not None:
            self.map = map
            self['source'] = 'map'
        elif mapfile is not None:
            self._loadMap(mapfile)
        self.savePlots = savePlots
        self.filename = filename
        if self.filename is not None and read:
            self.read(loadSource=loadSource)

    #### IO

    def write(self, filename=None, force=False):
        if filename is None: filename = self.filename
        if os.path.exists(filename):
            raise RuntimeError('output file %s exists, use force=True to overwrite.' % filename)
        output = self.encode()
        pickle.dump(output, open(filename, 'w'))

    def read(self, filename=None, loadSource=False):
        if filename is None: filename = self.filename
        input = pickle.load(open(filename))
        self.decode(input, loadSource=loadSource)

    def encode(self):
        output = {}
        for k in list(self.keys()):
            output[k] = self[k]
        return output

    def decode(self, input, loadSource=False):
        for k, v in zip(list(input.keys()), list(input.values())):
            self[k] = v
        if loadSource:
            self._loadMap()

    #### low-level map and mask handling

    def _loadMap(self, filename=None):
        if filename is None:
            sw = self['init']['source'].split(':')
            if sw[0] == 'file':
                filename = sw[1]
            else:
                raise RuntimeError('cannot load source map for %s' % self['init']['source'])
        self._clearMap()
        self.mapfile = filename
        self.map = moby2.mapping.fits_map.spaceMap(self.mapfile)
        self['init']['source'] = 'file:%s' % self.mapfile

    def _clearMap(self):
        for k in ['map', '_cmask', '_wmask']:
            if hasattr(self, k): delattr(self, k)
        
    def _checkMap(self):
        if not hasattr(self, 'map'):
            raise RuntimeError('BeamObs needs a map for this operation.')

    def _getWeightMask(self, regenerate=False):
        if regenerate or not hasattr(self, '_wmask'):
            self._wmask = (self.map.data != 0)
        return self._wmask.copy()

    def _getCentralMask(self, regenerate=False):
        if regenerate or not hasattr(self, '_cmask'):
            if 'basic' in self:
                x0, y0 = self['basic']['center']
            else:
                x0, y0 = self['init']['center0']
            r0 = self['init']['radius0']
            self._cmask = self.map.circleMask(r0, x0=x0, y0=y0)
        return self._cmask.copy()

    def _getMask(self, regenerate=False):
        return self._getWeightMask(regenerate=regenerate) * \
               self._getCentralMask(regenerate=regenerate)

    #### Plotting

    def _savefig(self, tag):
        if self.savePlots == 'show':
            return show()
        else:
            return savefig('%s%s.png' % (self.savePlots, tag))


    #### Retrieval...

    def reloadMap(self):
        map_src = self.get('init',{}).get('source')
        if map_src is None:
            self.map = None
        if map_src.startswith('file:'):
            filename = map_src[len('file:'):]
            self.map = moby2.mapping.fits_map.spaceMap(filename)
        else:
            raise ValueError("could not load map from init source = %s" % map_src)
        return self.map

    #### Analysis
    # the fit* routines compute parameters based on the map
    #

    def fit(self, **kwargs):
        return self.fitBasic(**kwargs) and \
               self.fitGaussian(**kwargs) and \
               self.fitNoise(**kwargs) and \
               self.fitPeak(**kwargs)
    
    def fitBasic(self, inverted=None, regenerate=False):
        """
        Perform naive analysis of map features.

        Returns True on success, and puts results in self['basic'].
        """
        if 'basic' in self and not regenerate:
            return True
        self._checkMap()
        
        m = self._getMask(regenerate=regenerate)
        self['basic'] = {}
        B = self['basic']  # alias
        B.update({
            'mask_pix': sum(m),
            'mask_frac': float(sum(m)) / sum(self._getCentralMask()),
            })

        self['flags']['sparse'] = B['mask_frac'] < 0.1
        if self['flags']['sparse']:
            # Don't even bother
            return False

        # Get basic amplitude level features, including a guess of inversion
        x, y, d = xyz_for_mask(self.map, m)
        mx, mn = d.max(), d.min()
        if inverted is None:
            inverted = median(d) > (mx + mn)/2
        B.update({
            'cal_sign': 1. - 2 * float(int(inverted)),
            'level': d.mean(),
            'amp': mx - mn,
            })

        # Pick a peak
        peak_idx = (d * B['cal_sign']).argmax()
        x0, y0 = x[peak_idx], y[peak_idx]
        B['center'] = (x0, y0)
        # Trigger regeneration of central mask
        self._getCentralMask(regenerate=True)

        # Annular mask to estimate power in peak.
        m0 = self._getWeightMask()
        r0 = self['init']['radius0']
        m1, m2 = self.map.circleMask(r0/2, x0=x0, y0=y0), \
                 self.map.circleMask(r0, x0=x0, y0=y0)
        m2 *= m0*~m1
        m1 *= m0
        # Base level from outer ring
        base = self.map.data[m2].mean()
        dx, dy = self.map.x[0,1] - self.map.x[0,0], self.map.y[1,0] - self.map.y[0,0]
        dA = dx*dy
        # Integrate inner ring to estimate flux in beam
        flux = dA*sum(self.map.data[m1] - base)*B['cal_sign']
        # Under gaussian assumption, estimate sigma and thus fwhm.
        B['sigma'] = sqrt(flux / (2.*pi*B['amp']))
        B['fwhm'] = B['sigma'] * SIGMA_TO_FWHM
        return True

    def fitGaussian(self, regenerate=False):
        """
        Simple gaussian fit; good for determining beam center and base level.

        Returns True on success, and puts results in self['gaussian'].
        """
        self._checkMap()
        if 'basic' not in self:
            self.fitBasic()
        m = self._getMask()
        mi = m.nonzero()

        # Data and initial parameters for fit.
        B = self['basic']
        x, y, z = xyz_for_mask(self.map, m)
        z = z * B['cal_sign']
        p = list(B['center']) + [B['amp'], B['sigma'], B['sigma'],
                                 B['level']*B['cal_sign'], 0.]

        # Fit first without angle, then full model
        args = (beamModels.gauss_model, x, y, z)
        for flags in [[1,1,0,0,0,0,0],
                      [1,1,1,1,1,1,0],
                      [1,1,1,1,1,1,1]]:
            p, ok = moby2.util.fitting.multi_leastsq(
                beamModels.bg_resid, p, args=args,
                free=np.array(flags, 'bool'))

        p = list(p)
        # Keep fwhm positive
        p[3], p[4] = abs(p[3]), abs(p[4])
        # Force fwhm_a > fwhm_b
        if p[3] < p[4]:
            p[3], p[4], p[6] = p[4], p[3], p[6] - 90.
        # Keep angle in +-90
        p[6] = fmod(fmod(p[6], 180.) + 180, 180.)
            
        # Residual RMS
        z0 = beamModels.gauss_model(x, y, p)
        rms = (z - z0).std()

        G = {
            'center': tuple(p[0:2]),
            'amp': p[2],
            'fwhm_a': p[3] * SIGMA_TO_FWHM,
            'fwhm_b': p[4] * SIGMA_TO_FWHM,
            'level': p[5],
            'angle': p[6],
            'dx': p[0],
            'dy': p[1],
            'p': p,
            'resid': rms,
            }
        G['solid_angle'] = 2.*pi*p[3]*p[4]*(pi/180.)**2

        self['gaussian'] = G

        if self.savePlots is not None:
            clf()
            subplot(211)
            r = ((x-p[0])**2+(y-p[1])**2)**0.5 * 60.
            scatter(r, z, marker='+')
            ylabel('Inner beam')
            subplot(212)
            z0 = beamModels.gauss_model(x, y, p)
            scatter(r, z-z0, marker='+')
            xlabel('Radius (arcmin)')
            ylabel('Fit residual')
            self._savefig('gaussfit')
        if any(isnan(p)):
            return False
        return True

    def fitPeak(self, regenerate=False):
        """
        
        To estimate peak, fit airy model only to inner-most points
        (defined as being inside the Gaussian FW ellipse).  Do not
        allow the background level to vary, however.
        
        Returns True on success, and puts results in self['peak'].
        """
        self._checkMap()
        if 'gaussian' not in self:
            self.fitGaussian()

        # Mask in the top half of the beam, according to gauss.
        B = self['gaussian']
        x, y = self.map.x - B['dx'], self.map.y - B['dy']
        xe, ye = beamModels.elliptify((x+y*0).ravel(), (y+x*0).ravel(),
                                      B['fwhm_a']/SIGMA_TO_FWHM,
                                      B['fwhm_b']/SIGMA_TO_FWHM,
                                      B['angle'])
        m = self._getMask() * (xe**2+ye**2 < 1.).reshape(self.map.data.shape)
        if sum(m) < 4:
            return False
        if sum(m) < 10:
            print('Warning: fitting peak with only %i points.' % sum(m))
        x, y, z = xyz_for_mask(self.map, m)
        z = z * self['basic']['cal_sign']

        # Start fit at gaussian result.
        fs = beamModels.airy_fwhm_scale / beamModels.gauss_fwhm_scale
        p = [B['dx'], B['dy'], B['amp'], B['fwhm_a']*fs, B['fwhm_b']*fs,
             B['level'], B['angle']]

        # Fit without varying background level or center or angle
        args = (beamModels.airy_model, x, y, z)
        for flags in [[0,0,1,1,1,0,0]]:
            p, ok = moby2.util.fitting.multi_leastsq(
                beamModels.bg_resid, p, args=args,
                free=np.array(flags, 'bool'))

        p = list(p)
        # Force fwhm_a > fwhm_b
        if p[3] < p[4]:
            p[3], p[4], p[6] = p[4], p[3], p[6] - 90.
        # Keep angle in +-90
        p[6] = fmod(fmod(p[6], 180.) + 180, 180.)

        # Fit residuals and quality
        z0 = beamModels.airy_model(x, y, p)
        quality = (z-z0).std() / self['noise']['rms']

        fs = beamModels.airy_fwhm_scale
        P = {
            'center': tuple(p[0:2]),
            'amp': p[2],
            'fwhm_a': p[3] * fs,
            'fwhm_b': p[4] * fs,
            'angle': p[6],
            'p': p,
            'quality': quality,
            }
        self['peak'] = P

        if self.savePlots is not None:
            clf()
            subplot(211)
            r = ((x-p[0])**2+(y-p[1])**2)**0.5 * 60.
            scatter(r, z, marker='+')
            ylabel('Inner beam')
            subplot(212)
            z0 = beamModels.airy_model(x, y, p)
            scatter(r, z-z0, marker='+')
            xlabel('Radius (arcmin)')
            ylabel('Fit residual')
            self._savefig('peakfit')
            clf()
            s = self.map.data.shape
            mm = ones(s, 'bool')
            x, y, z = xyz_for_mask(self.map, mm)
            r = (z - beamModels.airy_model(x, y, p)).reshape(s)
            subplot(211)
            self.map.imshow(data=r)
            self.map.zoom(.02)
            colorbar()
            subplot(212)
            self.map.imshow(data=m)
            self.map.zoom(.02)
            colorbar()
            savefig('peakfit_r.png')
            print(r[m].mean())
        if any(isnan(p)):
            return False
        return True

    def fitNoise(self, regenerate=False):
        """
        Estimate map white noise level in annulus.

        Returns True on success, and puts results in self['noise'].
        """
        self._checkMap()
        m0 = self._getWeightMask()
        r, (x0, y0) = self['init']['radius0'], self['gaussian']['center']
        m1, m2 = self.map.circleMask(r/2, x0=x0, y0=y0), \
                 self.map.circleMask(r, x0=x0, y0=y0)
        self['noise'] = {
            'rms': self.map.data[m2*~m1*m0].std()
            }
        return True

    #### Auxiliary information
    # the find* routines accumulate observation information from misc. databases

    def find(self, basename=None):
        return self.findInfo(basename=basename) and self.findWeather() and \
               self.findEphemeris()

    def findInfo(self, basename=None):
        """
        Load interesting todInfo data for this observation.

        Sets self['info'] and returns True on success.
        """
        from moby2.instruments import actpol
        tdb = actpol.TODDatabase()
        if 'info' in self:
            return True
        if basename is None:
            basename = util.extract_basename(self['init']['source'])
        rec = tdb.get_record(basename)
        self['info'] = {}
        for k in rec.fields: # fi._keys:
            self['info'][k] = getattr(rec, k)
        return True
        
    def findWeather(self):
        """
        Load APEX weather data for this observation.
        
        Sets self['weather'] and returns True on success.
        """
        if 'weather' in self:
            return True
        from moby.weather.APEX_weather import get_weather
        if 'info' not in self:
            self.findInfo
        self['weather'] = {}
        keys = ['radiometer', 'temperature', 'pressure',
                'humidity', 'dewpoint', 'windspeed', 'winddirection']
        t0, t1 = self['info']['ctime_start'], self['info']['ctime_end']
        for k in keys:
            self['weather'][k] = get_weather(t0, t1, k)
        self['weather']['ok'] = not any(isnan(array(list(self['weather'].values()))))
        return True

    def findEphemeris(self, sourceName=None):
        if sourceName is None:
            if 'info' not in self:
                self.findInfo()
            sourceName = self['info']['obs_detail']
        sourceName = sourceName.capitalize()
        # Use these old ephemerides; evaluate at center time of observation
        from calibration import calibration_support as cs
        try:
            e = cs.planetEphemeris(target=sourceName)
        except IOError:
            print('Could not find ephemeris data for %s' % sourceName)
            return False
        t = 0.5*(self['info']['ctime_start'] + self['info']['ctime_end'])
        E = dict([ (k, float(e.get(k, t))) for k in list(e.interpData.keys()) ])
        # Store the planet radii as well
        p = cs.ptsrcPhysicalInformation(target=sourceName)
        r, f = p.data['equatorial_radius'], p.data['flattening']
        E['flattening'] = f
        E['equatorial_radius'] = r
        E['polar_radius'] = r*(1-f)
        # Compute corrections due to oblateness aspect
        lat = E['subelat'] * pi/180.      # geodetic sub earth latitude
        ae = (1-f)                        #oblateness aspect
        # This probably has a real name...
        E['second_radius'] = r*ae / sqrt(sin(lat)**2 + ae**2*cos(lat**2))
        # Solid angle, steradians
        E['solid_angle'] = pi*r*E['second_radius'] / E['earthdist']**2

        self['eph'] = E
        return True


if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser(usage='%prog [options] filename[, filename ...]')
    o.add_option('--list',action='store_true',
                 help='Load files into a list of BeamObsList objects.')
    o.add_option('--update',
                 help='Recompute the specified stage for the loaded archive.')
    opts, args = o.parse_args()

    if opts.list:
        bols = [BeamObsList(f) for f in args]
        print('variable "bols" has list of BeamObsLists')
        if opts.update is not None:
            print('Recomputing "%s"' % opts.update)
            for i, bol in enumerate(bols):
                print(args[i])
                for b in bol:
                    if opts.update in b:
                        b.pop(opts.update)
                    b.find()
                bol.write(args[i] + 'x')
    else:
        for filename in args:
            b = BeamObs(mapfile=filename)
            b.fit()
            #b.find()
            print(filename.split('/')[-1])
            for d,k,s in [('gaussian','amp', 1.),
                          ('gaussian','fwhm_a', 60.),
                          ('gaussian','fwhm_b', 60.), 
                          ('gaussian','angle', 1.),
                          ('noise','rms', 1.),
                          ('basic','mask_frac', 1.)]:
                print('%-10s %10.6f' % (k+':', s*b[d][k]))
            print()

