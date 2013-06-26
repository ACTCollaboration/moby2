from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np
from numpy import pi, abs, fmod, modf

# Handling of branch cuts and angle conversions

_ANGLE_PERIOD = 360.
_ANGLE_BRANCH = -180.

def rerange(x, lo=_ANGLE_BRANCH, period=_ANGLE_PERIOD):
    """
    Remove multiples of period to bring x into interval [lo, lo+period).
    """
    return fmod(fmod(x - lo, period)+period,period) + lo

def is_within(x, lo, hi, period=_ANGLE_PERIOD):
    """
    Returns true if some branch of x satisfies lo <= x < hi.
    
    x can be an array.
    """
    # Note the rerange numpifies the angle sufficiently that we can
    # use np.array-style boolean operations (the '*' operator for
    # AND).
    x = rerange(x, lo, period=period)
    return (lo <= x)*(x < hi)


def to_deg(x):
    return 180.*x/pi

def to_rad(x):
    return pi*x/180


def from_sexagesimal(x, hours=False):
    """
    Given string of the form "dd:mm:ss.sss" or "hh:mm:ss.sss" (when
    hours=True), convert to decimal degrees.
    """
    w = x.split(':')
    w[0] = w[0].strip()
    s = 1.
    if w[0][0]=='-':
        s=-1.
        w[0] = w[0][1:]
    y = 0.
    c = 1.
    for yy in w:
        y += float(yy) * c
        c /= 60.
    if hours: y*= 15.
    return y*s

def to_sexagesimal(x, hours=False, hms=False):
    """
    Given x in decimal degrees, convert to sexagesimal degrees or hours.
    """
    s = ''
    if x < 0:
        s = '-'
        x *= -1
    if hours or hms:
        x = x / 15.
    di, i = modf(x)
    dm, m = modf(abs(di)*60.)
    if hms:
        return s+'%ih%02im%05.2fs' % (int(i), int(m), dm*60)
    return s+'%i:%02i:%05.2f' % (int(i), int(m), dm*60)

if __name__ == '__main__':
    # Test conversion
    z_set = [0.0166, 4.543123, 0.12, -12., -78.32, -0.0166]
    for z_f in z_set:
        z_s = to_sexagesimal(z_f)
        z_ff = from_sexagesimal(z_s)
        print('%12.6f %14s %12.6f' % (z_f, z_s, z_ff))
    
