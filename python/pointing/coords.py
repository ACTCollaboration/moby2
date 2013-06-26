from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np

from . import quaternions as quat

Q_GAL = None

def eq_to_q(ra, dec):
    """
    Get quaternion representation of point(s) at (ra,dec).  Then you
    can rotate it ('em).  Inverse of q_to_eq.
    """
    z = np.sin(dec)
    c = np.cos(dec)
    x, y = np.cos(ra)*c, np.sin(ra)*c
    del c
    out = np.empty(x.shape + (4,))
    out[...,0] = 0.
    out[...,1] = x
    out[...,2] = y
    out[...,3] = z
    return out

def q_to_eq(q):
    """
    Get (ra, dec) coordinates corresponding to quaternion(s) in q.
    Inverse of eq_to_q.
    """
    x, y, z = q[...,1], q[...,2], q[...,3]
    return np.arctan2(y,x), np.arcsin(z)

def qrot_fp_to_eq(ra, dec, gamma=0.):
    """
    Returns the quaternion rotation that takes vectors in the "focal
    plane" (here defined as vectors near z^, with x^ "up" and y^
    "right").

    The angle gamma produces an initial rotation, clockwise in the
    (x,y) plane, prior to the rotation from z^ to the point
    corresponding to (ra, dec).
    """
    return quat.rots([(3, ra), (2, np.pi/2 - dec), (3, np.pi + gamma)])

def qrot_eq_to_fp(ra, dec, gamma=0.):
    """
    Inverse of qrot_fp_to_eq.
    """
    return quat.conj(qrot_fp_to_eq(ra, dec, gamma))

def qrot_eqo_to_eq(ra, dec, gamma=0.):
    """
    Returns the quaternion rotation that takes a vector pointing to
    the equatorial origin (1,0,0) to the equatorial point (ra, dec).

    The angle gamma produces an initial rotation, clockwise in the
    (y,z) plane, prior to the rotation from x^ to the point
    corresponding to (ra, dec).
    """
    return quat.rots([(3, ra), (2, -dec), (1, gamma)])

def qrot_eq_to_eqo(ra, dec, gamma=0.):
    """
    Inverse of qrot_eqo_to_eq.
    """
    return quat.conj(qrot_eqo_to_eq(ra, dec, gamma))

def qrot_J2000_to_gal():
    """
    Returns the quaternion rotation that takes a vector in J2000
    equatorial coordinates to a vector in Galactic coordinates.
    """
    global Q_GAL
    if Q_GAL is None:
        DEG = np.pi/180.
        ra, dec, lon_pos = 192.859508*DEG, 27.128336*DEG, 122.932*DEG
        Q_GAL = quat.rots([(3, lon_pos-np.pi), (2, -np.pi/2+dec), (3, -ra)])
    return Q_GAL

def J2000_to_gal(ra, dec):
    """
    Returns galactic (l, b) associated with J2000 (ra, dec).
    """
    Q = qrot_J2000_to_gal()
    return q_to_eq(quat.rotate(Q, eq_to_q(ra, dec)))

def galactic_to_J2000(l, b):
    """
    Returns J2000 (ra, dec) associated with galactic (l, b).
    """
    Q = quat.conj(qrot_J2000_to_gal())
    return q_to_eq(quat.rotate(Q, eq_to_q(l, b)))
    
