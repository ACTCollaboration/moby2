from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np

import moby2.libactpol as libactpol
from moby2.util import encode_c, get_user_config
from . import quaternions as quat

def get_coords(ctime, az, alt,
               fields=None,
               focal_plane=None,
               final_q=None,
               weather=None):
    """
    Get "exact" J2000 pointing for the telescope.

    ctime, az, alt are vectors of unix time-stamp, boresight azimuth
    and altitude (radians).  These should have dtype=float64 and be
    the same length.

    fields is a list of strings that indicate what coordinate vectors
    to return.  Note that no matter how many you pass in, they are all
    getting computed.  Allowed fields are:

       'ra', 'dec', 'gamma'  -- ra, dec, and polarization angle
       'cos_dec', 'sin_dec', 'cos_gamma', 'sin_gamma','cos_2gamma',
         'sin_2gamma'        -- pre-computed trig
       'x', 'y', 'z'         -- cartesian coordinates of (ra, dec)

    focal_plane is either None (in which case the boresight
    coordinates are converted) or a FocalPlane object (in which case
    the coordinates of each detector are computed).

    The final_q is an optional quaternion that is applied as a final
    rotation to the J2000 coordinates.  This can be used to select an
    origin for the output coordinate system (the native spherical
    coordinates).  Use "get_tracking_quaternion" and friends to
    generate one.
    """
    ctime, az, alt = [np.asarray(_t, dtype='float64') for _t in [ctime, az, alt]]
    return libactpol.get_coords(ctime, az, alt,
                                fields,
                                encode_c(focal_plane),
                                final_q,
                                encode_c(weather))


def get_coords_inverse(ctime, az, alt, ra, dec,
                       final_q=None,
                       weather=None):
    return libactpol.get_coords_inverse(ctime, az, alt, ra, dec,
                                        final_q, encode_c(weather))


def set_bulletin_A(mjd=None, dut1=None, dx=None, dy=None, params=None):
    """
    Update the IERS bulletin A parameters in libactpol.  mjd, dut1,
    dx, dy must be arrays of equal length.  mjd should be integral and
    sequential.
    """
    if all([thing is None for thing in [dx, dy, dut1, mjd, params]]):
        params = get_user_config().get('bulletin_A_settings', {})
    if params is not None:
        data = np.loadtxt(params['filename'], unpack=1)
        data = [encode_c(data[i].copy()) for i in params['columns']]
        mjd = data[0].astype('int32')
        dx, dy, dut1 = data[1:]
    try:
        libactpol.set_bulletinA(int(mjd[0]), int(mjd[-1]), dut1, dx, dy)
    except:
        raise RuntimeError("libactpol rejected bulletin A data.")
    return mjd, dut1, dx, dy


## Where do these belong?

def project_data_to_map(data, cuts, wand, proj, fplane,
                        map_out=None, weights_out=None,
                        dets=None):
    if map_out is None:
        map_out = proj.zeros(dtype='float32')
    if weights_out is None:
        weights_out = proj.zeros(dtype='int32')
    if cuts is not None:
        cuts = encode_c(cuts)
    libactpol.wand_project_data_to_map(
        data,
        encode_c(cuts),
        encode_c(wand),
        encode_c(proj),
        encode_c(fplane),
        map_out, None, None,
        weights_out)
    return map_out, weights_out

def project_map_to_data(map_in, cuts, wand, proj, fplane,
                        data_out=None,
                        dets=None):
    if data_out is None:
        data_out = np.zeros((fplane.x.shape[0], wand.ra.shape[0]), 'float32')
    libactpol.wand_project_map_to_data(
        map_in,
        None, #Q
        None, #U
        encode_c(cuts),
        encode_c(wand),
        encode_c(proj),
        encode_c(fplane),
        data_out)
    return data_out

def project_data_to_IQU(data, cuts, wand, proj, fplane,
                        iqu_output=None, weights_out=None,
                        dets=None):
    if iqu_output is None:
        iqu_output = [proj.zeros(dtype='float32') for x in 'iqu']
    if weights_out is None:
        weights_out = proj.zeros(dtype='int32')
    if cuts is not None:
        cuts = encode_c(cuts)
    libactpol.wand_project_data_to_map(
        data,
        encode_c(cuts),
        encode_c(wand),
        encode_c(proj),
        encode_c(fplane),
        iqu_output[0], iqu_output[1], iqu_output[2],
        weights_out)
    return iqu_output, weights_out

def project_IQU_to_data(maps_iqu, cuts, wand, proj, fplane,
                        data_out=None,
                        dets=None):
    if data_out is None:
        data_out = np.zeros((fplane.x.shape[0], wand.ra.shape[0]), 'float32')
    libactpol.wand_project_map_to_data(
        maps_iqu[0], maps_iqu[1], maps_iqu[2],
        encode_c(cuts),
        encode_c(wand),
        encode_c(proj),
        encode_c(fplane),
        data_out)
    return data_out
