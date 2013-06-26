from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Implementation of a global pointing model.

The global pointing model accounts for difference between true horizon
coordinates and the naive horizon coordinates one might compute from
the telescope azimuth and altitude encoders.

The focal plane coordinates are defined as xi, eta with xi to the
right (parallel to az) and eta up (parallel to alt).  These are mapped
into a 3d right-handed cartesian system with (x,y,z) = (eta, xi, 1 -
eps) where eps is probably pretty small.

Horizon coordinates have x up, y east, z north.  This simplifies the
transformation from focal plane to horizon.
"""

import numpy as np

from . import quaternions as quat

DEG = np.pi/180
ARCMIN = DEG/60
ARCSEC = DEG/60

DEFAULT_ALT_SAG_PIVOT = 45. * DEG

_pretty_units = {
    '_default': ('arcmin', ARCMIN),
    'alt_sag_alt': ('deg', DEG),
    'alt_sag': ('arcmin per deg', ARCMIN/DEG),
    'alt_cant_sag': ('arcmin per deg', ARCMIN/DEG),
    'alt_rot': ('arcmin per deg', ARCMIN/DEG),
    }

_units_of = lambda k: _pretty_units.get(k, _pretty_units['_default'])


class GlobalModel(dict):
    """
    We work with quaternion rotations, like libactpol.

    Most routines accept (az,alt) coordinates in radians, which can be
    scalars or vectors of equal length.
    """

    # Version -- should be incremented if the meaning of any of the
    # keys changes.  This can be used to prevent stored parameters
    # from being used in an invalid way.
    model_version = 2

    # Keys this model understands
    model_keys = [
        'az0',
        'alt0',
        'tilt_cosaz',
        'tilt_sinaz',
        'alt_cant',
        'alt_sag',
        'alt_sag_alt',
        'alt_cant_sag',
        'alt_rot',
        'ar_dx',
        'ar_xi',
        'ar_eta',
        ]

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        if not 'alt_sag_alt' in self:
            self['alt_sag_alt'] = DEFAULT_ALT_SAG_PIVOT
        for k in self.model_keys:
            if not k in self:
                self[k] = 0.

    def to_string(self, pretty_units=True, nonzero_only=False):
        output = ''
        scode = '%%%is' % max(list(map(len, self.model_keys)))
        for k in self.model_keys:
            if nonzero_only and self[k] == 0:
                continue
            o = scode % k + '  '
            if pretty_units:
                us, uu = _units_of(k)
                o += '%11.3f %s' % (self[k]/uu, us)
            else:
                o += '%11.5e' % self[k]
            output += o +'\n'
        return output

    def __repr__(self):
        return '%s %s' % (self.__class__.__name__, dict.__repr__(self))

    def get_q(self, az, alt):
        """
        Return the quaternion rotation that takes libactpol focal
        plane coordinates to the corrected (az_corr, alt_corr), given
        the encoder (az,alt).
        """
        tilt_cosaz, tilt_sinaz, alt_cant, az0, alt0, \
            alt_sag, alt_sag_alt, alt_cant_sag, \
            ar_dx, alt_rot = \
            [self.get(key, 0.) for key in 
             ['tilt_cosaz', 'tilt_sinaz', 'alt_cant', 'az0', 'alt0',
              'alt_sag', 'alt_sag_alt', 'alt_cant_sag',
              'ar_dx', 'alt_rot']]
        # Apply sag to alt and cant in alt.
        alt_sagged = alt + (alt-alt0-alt_sag_alt) * alt_sag - alt0
        alt_cant_sagged = alt_cant + (alt-alt0-alt_sag_alt) * alt_cant_sag
        # Convert base tilt into rotation and tilt
        base_az = np.arctan2(tilt_sinaz, tilt_cosaz)
        base_tilt = (tilt_sinaz**2 + tilt_cosaz**2)**.5
        # Linear coupling of altitude into azimuth
        az_twist = alt_rot * (alt-alt0)
        # And then it's just a matter of...
        return quat.rots(((1, -base_az),
                          (2, base_tilt),
                          (1, base_az),
                          (1, -(az-az0-az_twist)),
                          (3, alt_cant_sagged),
                          (2, alt_sagged),
                          (3, -alt_cant_sagged),
                          (1, -ar_dx),
                          (1, -self['ar_xi']),
                          (2, self['ar_eta']),
                          ))

    def get_naive_q(self, az, alt):
        """
        Return the quaternion rotation that takes z^ to the
        uncorrected (az,alt).
        """
        return quat.rots(((1, -az),
                          (2, alt),
                          (1, -self['ar_xi']),
                          (2, self['ar_eta'])))
        
    def get_diff_q(self, az, alt):
        """
        Get the difference rotation between the full model computation
        and the naive rotation.  This can be applied as a correction
        to a vector that is already in naive horizon coordinates.
        """
        return quat.mul(quat.conj(self.get_naive_q(az,alt)),
                        self.get_q(az, alt))

    def get_focal_plane_offset(self, az, alt):
        """
        Obtain the focal plane shifts (dx, dy) induced by this model
        at the specified encoder position (az,alt).  dx and dy are
        taken roughly parallel to +az and +alt respectively.
        """
        diffq = self.get_diff_q(az, alt)
        z = np.array([0,0,0,1.])
        n = quat.rotate(diffq, z)
        return n[...,2], n[...,1]

    # Replace get_focal_plane_offset...
    def get_pointing_shift(self, az, alt):
        return [-x for x in self.get_source_shift(az, alt)]

    def get_source_shift(self, az, alt):
        return self.get_focal_plane_offset(az, alt)

    def get_indexed_reference(self, keys):
        """
        Get a child object that can get / set parameters in the parent
        object according to integer indices rather than keyword.  This
        is especially useful in fitting routines that like to pass in
        tuples of parameters.

        keys should be a list of names of model parameters that will
        be associated with numerical indices [0,1,...] in the child
        object.

        E.g.:

           >>> indexed = my_global_model.get_indexed_reference(['az0', 'alt0'])
           >>> print indexed.values()
           (0., 0.)
           >>> my_global_model['az0'] = .001
           >>> print indexed.values()
           (0., 0.001)
           >>> indexed.update((0.002, 0.003))
           >>> print my_global_model['alt0']
           0.002
        """
        class IndexedGlobalModelReference:
            def __init__(self, parent, keys):
                self.parent = parent
                self.keys = keys
                for k in keys:
                    if not k in parent:
                        print('Warning, %s is not a parent parameter.' % k)
            def update(self, values):
                if len(values) != len(self.keys):
                    raise ValueError("provided %i of %i expected values" % \
                        (len(values), len(self.keys)))
                for k, v in zip(self.keys, values):
                    self.parent[k] = v
            def values(self):
                return tuple([self.parent[k] for k in self.keys])
            def get_focal_plane_offset(self, az, alt):
                return self.parent.get_focal_plane_offset(az, alt)
        return IndexedGlobalModelReference(self, keys)

    @classmethod
    def naive(cls):
        return cls()

