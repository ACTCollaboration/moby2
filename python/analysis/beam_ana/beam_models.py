from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
from numpy import *
from scipy.special import j1


# Coordinates

def elliptify(x, y, wx, wy, angle):
    """
    Transform x,y into coordinates on the ellipse specified by
    semi-axis lengths wx,wy oriented at angle to the x-axis.

    This will, for example map the ellipse itself to the unit
    circle.
    """
    g = array([[cos(angle), sin(angle)], [-sin(angle), cos(angle)]])
    v = vstack((x,y))
    om = array([1./wx, 1./wy])
    arg = dot(g, v) * om.reshape(-1, 1)
    return arg[0], arg[1]


# 2D gaussian with baseline

def gauss_model(x, y, params):
    x0, y0, h, wx, wy, base, angle = params
    angle *= pi/180.
    X = 1.
    if abs(angle) > 4*360:
        X = 100.
    g = array([[cos(angle), sin(angle)], [-sin(angle), cos(angle)]])
    xx, yy = elliptify(x-x0, y-y0, wx, wy, angle)
    arg = xx**2 + yy**2
    return X * (h * exp(-arg / 2) + base)

# 2D airy function with baseline

def airy_model(x, y, params):
    x0, y0, h, wx, wy, base, angle = params
    angle *= pi/180.
    if abs(angle) > 4*360:
        X = 100.
    else:
        X = 1.
    g = array([[cos(angle), sin(angle)], [-sin(angle), cos(angle)]])
    v = vstack((x-x0,y-y0))
    om = array([1./wx**2, 1./wy**2])
    arg = sqrt(dot(om, dot(g, v)**2))
    arg[(abs(arg)<1e-9).nonzero()] = 1e-9
    return X * (h * (2.*j1(arg)/arg)**2 + base)


# generic residualizer for beam models

def bg_resid(*args):
    model, x, y, z = args[1:]
    return z - model(x, y, args[0])


# Models and their parameters

gauss_fwhm_scale =  2.*sqrt(2.*log(2))
airy_fwhm_scale = 2.*1.61634

models = {
    'gauss': {
    'model_function': gauss_model,
    'params': ['x0', 'y0', 'h', 'wx', 'wy', 'base', 'angle'],
    'defaults': [0., 0., 0.01, 0.01, 0.01, 0., 0.],
    'scales': [1., 1., 1., gauss_fwhm_scale, gauss_fwhm_scale, 1., 1.],
    'formats': ['%10.4f', '%10.4f', '%10.6f', '%10.6f', '%10.6f', '%8.4f', '%10.3f'],
    'fit_sequence': [[1,1,1,1,1,1,0], [1,1,1,1,1,1,1]],
    },
    'airy': {
    'model_function': airy_model,
    'params': ['x0', 'y0', 'h', 'wx', 'wy', 'base', 'angle'],
    'defaults': [0., 0., 0.01, 0.01, 0.01, 0., 0.],
    'scales': [1., 1., 1., airy_fwhm_scale, airy_fwhm_scale, 1., 1.],
    'formats': ['%10.4f', '%10.4f', '%10.6f', '%10.6f', '%10.6f', '%8.4f', '%10.3f'],
    'fit_sequence': [[1,1,1,1,1,1,0], [1,1,1,1,1,1,1]],
    }
    }

