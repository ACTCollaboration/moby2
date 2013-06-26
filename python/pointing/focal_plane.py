from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np

import moby2
from . import base

class FocalPlane:

    """
    The focal plane coordinate system is defined relative to a nominal
    boresight aligned with z^.  x^ corresponds roughly to +azimuth,
    and y^ corresponds to +altitude.  Yes, it's sort of left-handed.

    The angle "phi" describes the polarization angle and is defined,
    roughly, as the angle between the polarization axis and x^,
    increasing towards y^.

    The "mask" is a boolean array indicating which channels have valid
    offsets.  "det_uid" gives the detector id associated with
    each offset.
    """

    x = None
    y = None
    phi = None
    mask = None
    det_uid = None

    def __init__(self, x=None, y=None, phi=None, mask=None, det_uid=None):
        if x is not None:
            x = self.to_offset_vector(x)
            if phi is None:
                phi = 0*x
            if mask is None:
                mask = np.ones(x.shape, 'bool')
        self.x = self.to_offset_vector(x)
        self.y = self.to_offset_vector(y)
        self.phi = self.to_offset_vector(phi)
        self.mask = np.array(mask, dtype='bool')
        if det_uid is not None:
            self.det_uid = np.array(det_uid, dtype='int')

    def get_coords(self, ctime, az, alt, fields=None, final_q=None, weather=None):
        """
        Compute high-precision detector coordinates at specified ctime
        and boresight az and alt.  All arguments are passed directly
        to pointing.get_coords.
        """
        return base.get_coords(ctime, az, alt, focal_plane=self, fields=fields,
                               final_q=final_q, weather=weather)

    @staticmethod
    def to_offset_vector(x):
        return np.array(x, dtype='float64')

    def _encode_c(self, dets=None):
        x, y, phi = [self.to_offset_vector(_x)
                     for _x in [self.x, self.y, self.phi]]
        if dets is not None:
            x, y, phi = x[dets], y[dets], phi[dets]
        # Must match conventions in C interface
        return (x, y, phi)

    def subset(self, idx=None, det_uid=None, superset_ok=True):
        """
        Get a FocalPlane object that is a subset of this one.

        By default, all the valid pointings are selected.  Otherwise,
        indices are taken from idx, or matched to det_uid.
        """
        if idx is None:
            if det_uid is None:
                idx = np.ones(len(self.x), 'bool')
                mask = np.ones(idx.shape, 'bool')
            else:
                ldu = list(self.det_uid)
                idx = []
                for i in det_uid:
                    if i in ldu:
                        idx.append((True, ldu.index(i)))
                    else:
                        idx.append((False, 0))
                mask, idx = list(map(np.array, list(zip(*idx))))
                #idx = [ldu.index(i) for i in det_uid]
        else:
            idx = np.asarray(idx)
            mask = None
        out = self.__class__()
        for k in ['x', 'y', 'phi', 'mask', 'det_uid']:
            v = getattr(self, k)
            if v is not None:
                v = np.asarray(v)[idx]
                if mask is not None:
                    v[~mask] = 0
            setattr(out, k, v)
        if det_uid is not None:
            out.det_uid = np.array(det_uid)
        return out

    def copy(self):
        """
        Returns a FocalPlane that is a copy of self.
        """
        return self.subset()

    def split_mean(self):
        """
        Returns a pair of FocalPlane objects, the first of which gives
        the approximate position of the center of the array, and the
        second of which gives the positions of the detectors relative
        to that center.
        """
        mask = self.mask
        if mask is None:
            mask = np.ones(len(self.x), 'bool')
        x0, y0 = self.x[mask].mean(), self.y[mask].mean()
        # This implementation is an approximation... needs checking.
        center = self.__class__([x0], [y0], [np.pi/2], [True])
        children = self.__class__(np.asarray(self.x)-x0,
                                  np.asarray(self.y)-y0,
                                  np.asarray(self.phi).copy(),
                                  np.asarray(mask).copy())
        return center, children

    @classmethod
    def from_super_po_file(cls, filename):
        """
        Load pointing offsets from ACT-style multiOffset file.
        """
        if isinstance(filename, basestring):
            filename = open(filename)
        xya = np.zeros((3,33,32), 'float')
        on = False
        while True:
            line = filename.readline()
            if line == '': break
            w = line.split()
            if not on:
                on = (len(w) > 0) and (w[0] == 'begin_offsets')
                continue
            if w[0] == 'end_offsets':
                break
            r, c, x, y = [cast(_x) for cast,_x in zip([int,int,float,float],w)]
            xya[:,int(w[0]),int(w[1])] = [float(w[2]), float(w[3]), 1.]
        self = cls()
        self.x = xya[0].ravel()
        self.y = xya[1].ravel()
        self.phi = self.x * 0.
        self.mask = xya[2].ravel().astype('int').astype('bool')
        self.det_uid = np.arange(33*32)
        return self

    @classmethod
    def from_ACT_po_file(cls, filename):
        """
        Load pointing offsets from ACT-style pointingOffset file.
        """
        if isinstance(filename, basestring):
            filename = open(filename)
        header = filename.readline().split()
        nr, nc = int(header[2]), int(header[5])
        xy = np.zeros((2,nr,nc), 'float')

        # burn 2
        filename.readline()
        filename.readline()
        for r in range(nr):
            filename.readline()
            for c in range(nc):
                xy[:,r,c] = list(map(float, filename.readline().split()))
        self = cls()
        xy.shape = (2,-1)
        self.mask = abs(xy[0]) < 100
        xy[:,~self.mask] = 0.
        self.y, self.x = xy
        self.phi = self.x * 0.
        self.det_uid = np.arange(33*32)
        return self

