from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np

import moby2
import moby2.libactpol as libactpol
from moby2.util import encode_c, angles

from . import base as pointing_base
from . import coords
from .focal_plane import FocalPlane

import moby2.util.log as logger


class ArrayWand:
    """
    Encapsulation of pointing data necessary for ArrayWand pointing
    approximation.
    """
    
    proj_type = 'ra_dec'
    polarized = False
    fcenter = None
    nsamps = None

    ra = None
    dec = None
    cos_dec = None
    sin_dec = None
    cos_gamma = None
    sin_gamma = None
    cos_2gamma = None
    sin_2gamma = None
    cos_2hwp = None
    sin_2hwp = None

    x = None
    y = None

    def __init__(self, ctime, az, alt, fplane,
                 proj_type='tangent', polarized=False,
                 ref_coord=None, weather=None,
                 hwp_angle=None,
                 ):

        self.proj_type = proj_type
        self.polarized = polarized

        # All wands need the rotation angle
        fields = ['cos_gamma', 'sin_gamma']

        # ... and some need projection-specific coordinates
        if self.proj_type == 'tangent':
            fields.extend(['y', 'z'])
        else:
            fields.extend(['ra', 'dec', 'cos_dec', 'sin_dec'])

        # ... and some are polarized
        if self.polarized:
            fields.extend(['cos_2gamma', 'sin_2gamma'])

        # Rotation to the reference coordinate... can in principle be an array
        if ref_coord is None:
            ref_coord = (0., 0., 0.)
        elif len(ref_coord) == 2:
            ref_coord = tuple(ref_coord) + (0., )
        Q = coords.qrot_eq_to_eqo(*ref_coord)

        # Split focal plane about its center, point the center.
        self.fcenter, _ = fplane.split_mean()
        wand_data = pointing_base.get_coords(
            ctime, az, alt, fields=fields,
            focal_plane=self.fcenter,
            final_q=Q, weather=weather)

        if hwp_angle is not None:
            assert (hwp_angle is not True and hwp_angle is not False)
            self.cos_2hwp = np.cos(2*np.array(hwp_angle, np.float64))
            self.sin_2hwp = np.sin(2*np.array(hwp_angle, np.float64))

        # Correct RA to keep it in a single branch.
        for data, field in zip(wand_data, fields):
            if not field in ['ra']:
                continue
            data[:] = angles.rerange(data, data[0][0] - np.pi, 2*np.pi)

        # Store nsamps and the vectors
        self.nsamps = len(wand_data[0][0])
        for f,w in zip(fields, wand_data):
            setattr(self, f, w[0])

    def _encode_c(self):
        # Must match conventions in C interface
        if self.proj_type in ['ra_dec', 'ra_sindec']:
            return (self.proj_type, int(self.polarized),
                    int(self.cos_2hwp is not None),
                    self.cos_gamma, self.sin_gamma,
                    self.cos_2gamma, self.sin_2gamma,
                    self.ra, self.dec, self.cos_dec, self.sin_dec,
                    self.cos_2hwp, self.sin_2hwp)
        elif self.proj_type == 'tangent':
            return (self.proj_type, int(self.polarized),
                    int(self.cos_2hwp is not None),
                    self.cos_gamma, self.sin_gamma,
                    self.cos_2gamma, self.sin_2gamma,
                    self.y, self.z,
                    self.cos_2hwp, self.sin_2hwp)

    def get_coords(self, focal_plane, get_coords=True, get_rot=False,
                   get_pol=False):
        """
        Get coordinates at positions in focal_plane which should be a
        FocalPlane object.

        Returns tuple of pointing vectors; a subset of:

            x, y, cos_gamma, sin_gamma, cos_2gamma, sin_2gamma

        If get_coords==False, x and y are excluded.

        If get_rot==False, cos_gamma and sin_gamma are excluded.

        If get_pol==False, cos_2gamma and sin_2gamma are excluded.
        """
        # FIXME: this should be a rotation
        fp = FocalPlane(focal_plane.x - self.fcenter.x[0],
                        focal_plane.y - self.fcenter.y[0],
                        focal_plane.phi)
        return libactpol.wand_get_coords(encode_c(self),
                                         encode_c(fp),
                                         int(get_coords),
                                         int(get_rot),
                                         int(get_pol))

    get_pointing = get_coords

    @classmethod
    def for_tod(cls, tod, coords='tangent', ref_coord=(0., 0., 0.),
                polarized=False, hwp=None):
        """
        Like the default constructor, but boresight pointing is taken
        from tod.
        """
        # Split focal plane about its center, point the center.
        fcenter, _ = tod.fplane.split_mean()
        assert abs(fcenter.phi[0] - np.pi/2) < 1e-7
        if hwp is None and hasattr(tod, 'hwp_angle'):
            hwp = tod.hwp_angle
        self = cls(tod.ctime, tod.az, tod.alt, fcenter, 
                   proj_type=coords, polarized=polarized, ref_coord=ref_coord,
                   hwp_angle=hwp)
        return self

    @classmethod
    def for_tod_source_coords(cls, tod, ref_coord=(0., 0.), polarized=False,
                              warning_distance=10./60, scan_coords=False,
                              hwp=None
                              ):
        """
        Get ArrayWand for tod, in tangent plane coordinates centered
        on the provided ref_coord=(ra,dec), and possibly including a
        rotation that makes the X axis parallel to the scan direction.

        A warning is issued if the TOD does not scan within
        warning_distance (degrees) of the ref_coord.
        """
        if isinstance(ref_coord, basestring):
            ra_src, dec_src = moby2.ephem.get_source_coords(
                ref_coord, tod.ctime.mean())
            gamma_in = 0.
        elif len(ref_coord) == 3:
            ra_src, dec_src, gamma_in = ref_coord
        elif len(ref_coord) == 2:
            ra_src, dec_src = ref_coord
            gamma_in = 0.
        else:
            raise ValueError("ref_coord must consist of 2 or 3 values, or a string to decode.")
        # Preliminary computation, no rotation
        wand0 = cls.for_tod(tod, coords='tangent', polarized=polarized,
                            ref_coord=(ra_src, dec_src, gamma_in),
                            hwp=hwp)
        if not scan_coords:
            return wand0
        # Get index of sample that is nearest the ref_coord
        distance = (wand0.y**2+wand0.z**2)**.5
        ic0 = np.argmin(distance)
        assert tod.ctime[ic0] != 0  # To be dealt with elsewhere...
        if distance[ic0] * (180/np.pi) > warning_distance:
            logger.trace('moby', 0, 'Source at [%f, %f] is not in %s' %
                         (ra_src, dec_src, tod.info.name))
        gamma_src = np.arctan2(wand0.sin_gamma[ic0], wand0.cos_gamma[ic0])
        # Recompute with this rotation
        return cls.for_tod(
            tod, coords='tangent', polarized=polarized,
            ref_coord=(ra_src, dec_src, gamma_in+gamma_src),
            hwp=hwp)


class WandProjector:
    def __init__(self, wand, proj, fplane):
        self.wand = wand
        self.proj = proj
        # FIXME: this should be a rotation
        fp = FocalPlane(fplane.x - wand.fcenter.x[0],
                        fplane.y - wand.fcenter.y[0],
                        fplane.phi)
        self.fplane = fp

    def deproject_from_map(self, map_vect, output=None, cuts=None):
        if output is None:
            output = np.zeros((len(self.fplane.x), self.wand.nsamps),
                              dtype='float32')
        pointing_base.project_map_to_data(map_vect,
                                          cuts,
                                          self.wand,
                                          self.proj,
                                          self.fplane,
                                          data_out=output)
        return output

    def project_to_map(self, tod_vect,
                       output=None, weights_output=None,
                       cuts=None):
        if output is None:
            output = self.proj.zeros()
        if weights_output is True:
            weights_output = self.proj.zeros(dtype='int32')
        pointing_base.project_data_to_map(tod_vect,
                                          cuts,
                                          self.wand,
                                          self.proj,
                                          self.fplane,
                                          map_out=output,
                                          weights_out=weights_output)
        return output, weights_output

    # Interface for pcg

    def projectMapToTOD(self, tod, map_in, cuts=True):
        cuts = {True: tod.cuts, False: None}.get(cuts, cuts)
        return self.deproject_from_map(map_in.data, output=tod.data, cuts=cuts)

    def projectTODToMap(self, map_out, tod, cuts=True):
        cuts = {True: tod.cuts, False: None}.get(cuts, cuts)
        return self.project_to_map(tod.data,
                                   output=map_out.data,
                                   weights_output=map_out.weight,
                                   cuts=cuts)

class PolWandProjector:
    def __init__(self, wand, proj, fplane):
        self.wand = wand
        self.proj = proj
        # FIXME: this should be a rotation
        fp = FocalPlane(fplane.x - wand.fcenter.x[0],
                        fplane.y - wand.fcenter.y[0],
                        fplane.phi)
        self.fplane = fp

    def deproject_from_map(self, maps_iqu, output=None, cuts=None):
        if output is None:
            output = np.zeros((len(self.fplane.x), self.wand.nsamps),
                              dtype='float32')
        pointing_base.project_IQU_to_data(maps_iqu,
                                          cuts,
                                          self.wand,
                                          self.proj,
                                          self.fplane,
                                          data_out=output)
        return output

    def project_to_map(self, tod_vect,
                       iqu_output=None, weights_output=None,
                       cuts=None):
        if iqu_output is None:
            iqu_output = [self.proj.zeros() for x in 'iqu']
        if weights_output is True:
            weights_output = self.proj.zeros(dtype='int32')
        pointing_base.project_data_to_IQU(tod_vect,
                                          cuts,
                                          self.wand,
                                          self.proj,
                                          self.fplane,
                                          iqu_output=iqu_output,
                                          weights_out=weights_output)
        return iqu_output, weights_output

    # Interface for moby-style pcg not implemented!
