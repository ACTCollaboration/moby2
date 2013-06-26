from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
from moby2 import libactpol
import numpy as np

class GridPixelization:
    """
    Defines a mapping between a cartesian plane and a finite grid of
    pixels.  This information is passed to C-level projection
    routines.

    The approach is similar to the FITS conversion system.  The pixel
    index i along an axis is related to the cartesian coordinate x by

       i0 = pix + (x-x0)/dx

    Integral pixel indices are obtained by rounding, and the rounded
    values are invalid if they do not fall in [0,n0-1].

    We do this for two axes.  Axis 0 is the "fast" dimension for associated
    2-d array data.

    Note that we don't support rotations here.  Rotations can usually
    be handled in the main pointing code, at 0 cost.

    For CEA projections, dx0 and dy0 must be computed to take into
    account the fiducial latitude; dx0 and dy0 are not direct analogs
    to the FITS CDELT keyword.
    """

    n0 = 1
    n1 = 1
    pix0 = 0
    pix1 = 0
    x0 = 0.
    x1 = 0.
    dx0 = 1.
    dy0 = 1.
    
    def _encode_c(self):
        # Must match conventions in C interface
        return (self.n0, self.n1, self.pix0, self.pix1,
                self.x0, self.x1, self.dx0, self.dx1)

    def get_pix(self, x0, x1):
        return (x0 - self.x0) / self.dx0 + self.pix0, \
            (x1 - self.x1) / self.dx1 + self.pix1

    def zeros(self, dtype='float32'):
        return np.zeros((self.n1, self.n0), dtype=dtype)

    @classmethod
    def forFitsMap(cls, map0=None, pixn=None):
        if map0.fheader['CTYPE1'].strip() == 'X':
            pass
        elif map0.fheader['CTYPE1'].strip() == 'RA---TAN':
            pass
        else:
            raise ValueError('only implemented for linear and TAN maps (got'\
                ' %s)' % map0.fheader['CTYPE1'])
        DEG = np.pi/180
        self = cls()
        if pixn is None:
            pixn = map0.pixn
        g = pixn.grid
        self.dx0, self.dx1 = -g['CDELT1']*DEG, g['CDELT2']*DEG
        g = pixn.mapgrid
        self.pix0, self.pix1 = [int(g['CRPIX'+i]) - 1 for i in '12']
        self.n0, self.n1 = [g['NAXIS'+i] for i in '12']
        g = pixn.projection
        self.x0, self.x1 = [g['CRVAL'+i]*DEG for i in '12']
        return self
