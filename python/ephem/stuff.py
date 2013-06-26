from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np
import astropy
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Longitude, Latitude, Angle

import moby2
import moby2.util.angles as angles
import moby2.util.log as psLib

from .actEphem import ACTEphem

calibration_sources = [
    ('Jupiter','planet'),
    ('Mars', 'planet'),
    ('Venus', 'planet'),
    ('Saturn', 'planet'),
    ('Uranus', 'planet'),
    ('Neptune', 'planet'),
    ('J2000 83.63 22.01', 'source'),      #Tau A
    ('J2000 201.36 -43.01', 'source'),    #Cen A
    ('J2000 1.55 -6.39', 'source'),       #QSO B0003-066
    ('J2000 187.69 12.40', 'source'),     #JPB2009
    ('J2000 187.27 2.05', 'source'),      #WMAP 170
    ('J2000 164.62 1.57', 'source'),      #4C 01.28
    ('J2000 237.36 2.62', 'source'),      #QSO B1546+0246
    ('J2000 133.69 20.11', 'source'),     #QSO J0854+2006'
    ]

DEG = np.pi / 180

def get_sources_in_patch(ctime=None, ra_lims=None, dec_lims=None, tod=None,
                         source_list=None, add_source_list=None):
    """
    Find sources in a region of RA and dec.  The usual planets will be
    checked, unless a list of particular objects is passed in
    source_list.  The RA and dec. limits will be determined from the
    tod argument unless they are passed in explicitly (as (min,max)
    tuples).  Since these things tend not to move very fast, we do the
    computation at a single ctime (which defaults to
    tod.ctime.mean()).

    Returns a list of matches, with each match a tuple (source_name, ra, dec).
    """
    if tod is not None:
        mask = tod.pointing_mask
    if ra_lims is None or dec_lims is None:
        wand_eq = moby2.pointing.ArrayWand.for_tod(tod, coords='ra_dec')
        # Get ra, dec limits covered by this TOD.  Just use the array center, buffered by ~.5 deg.
        cos_dec0 = np.cos(wand_eq.dec.mean())
        if ra_lims is None:
            ra_lims = (wand_eq.ra[mask].min() - 0.5*DEG/cos_dec0,
                       wand_eq.ra[mask].max() + 0.5*DEG/cos_dec0)
        if dec_lims is None:
            dec_lims = (wand_eq.dec[mask].min() - 0.5*DEG, 
                        wand_eq.dec[mask].max() + 0.5*DEG)
    if ctime is None:
        ctime = tod.ctime[mask].mean()
    if source_list is None:
        source_list = calibration_sources
    if add_source_list is not None:
        source_list.extend(add_source_list)
    
    ephem = ACTEphem()
    ephem.set_ctime(ctime)

    matches = []

    for source_name, source_type in source_list:
        s_ra, s_dec = get_source_coords(source_name, ctime)
        psLib.trace('moby', 3, 'Source = %s,  RA = %f, Dec = %f' %
                    (source_name, s_ra/DEG, s_dec/DEG))
        if angles.is_within(s_ra, ra_lims[0], ra_lims[1], 2*np.pi) and \
                angles.is_within(s_dec, dec_lims[0], dec_lims[1], 2*np.pi):
            matches.append((source_name, float(s_ra), float(s_dec),
                            source_type))

    return matches

def get_sources_in_tod(tod, source_list=None, site=None, pointing_shift=None):
    """
    Find sources in a TOD. The usual planets will be                       
    checked, unless a list of particular objects is passed in      
    source_list.
    Since these things tend not to move very fast, we do the                       
    computation at a single ctime (which defaults to                                    
    tod.ctime.mean()).
    site argument can be passed as an astropy EarthLocation object,
    otherwise ACT loction will be used

    Returns a list of matches, with each match a tuple (source_name, ra, dec). 
    """
    if site is None:
        # Use ACT site
        site = EarthLocation(lon=-1.183116812024908 * u.rad,
                             lat=-0.40070141631911815 * u.rad,
                             height=5188 * u.m)

    if isinstance(source_list, basestring):
        sources = astropy.io.ascii.read(
            source_list,
            names=['name', 'ra', 'dec'])

    # if pointing_shift is not None:
    #     print np.deg2rad(pointing_shift[0]), np.deg2rad(pointing_shift[1] )
        # tod.az -= np.deg2rad(pointing_shift[0])
        # tod.alt -= np.deg2rad(pointing_shift[1])

    time = Time(tod.ctime.mean(), format='unix')
    min_az = Longitude(tod.az.min(), unit=u.rad)
    max_az = Longitude(tod.az.max(), unit=u.rad)
    throw = Angle(tod.az.max() - tod.az.min(), unit=u.rad)
    if np.median(tod.alt) > np.pi/2:
        print("TOD elevation is higher than 90 deg")
        return []
    mean_alt = Latitude(np.median(tod.alt), unit=u.rad)
    min_coord = SkyCoord(min_az, mean_alt, frame='altaz', location=site)
    max_coord = SkyCoord(max_az, mean_alt, frame='altaz', location=site)
     
    source_coords = SkyCoord(
        ra=sources['ra'], dec=sources['dec'], unit=u.deg,
        location=site,obstime=time)

    altaz = source_coords.transform_to('altaz')
    sel_alt = np.abs( mean_alt - altaz.alt ) < 5*u.deg
    sel_az = np.logical_and(
        altaz.az - min_coord.az > -2.*u.deg,
        max_coord.az - altaz.az > -2.*u.deg)
    sel = np.logical_and(sel_alt,sel_az)

    matches = []
    for source in sources[sel]:
        matches.append(
            (source['name'],
             np.deg2rad(source['ra']),
             np.deg2rad(source['dec']),
             'source'))

    return matches


def get_source_coords(source_desc, ctimes=None):
    """
    Get equatorial coordinates of a source.  source_desc should be one
    of:

    - a planet name, e.g. "Saturn"
    
    - a string of the form "J2000 <ra> <dec>", with <ra> and <dec>
      representing equatorial coordinates in degrees.

    - a tuple of the form ("J2000", ra, dec), where ra and dec are floats
      with units of radians.
    
    - instead of J2000, I guess you can use "HOR" too, in which case
      the arguments are azimuth and altitude.
    
    For planetary sources, or horizon coordinates, ctimes must be
    provided.
    """
    destring = False
    if isinstance(source_desc, basestring):
        # Convert string to tuple.
        words = source_desc.split()
        source_desc = (words[0],) + tuple(map(angles.to_rad, list(map(float, words[1:]))))
        destring = True
    else:
        return source_desc[1], source_desc[2]
        # Store original shape, then convert to 1-d array.
    orig_shape = np.asarray(ctimes).shape
    ctimes = np.array(ctimes).ravel()
    az, alt = np.zeros(ctimes.shape), np.zeros(ctimes.shape)

    if source_desc[0] == 'J2000':
        return source_desc[1], source_desc[2]
    elif source_desc[0] == 'HOR':
        az[:], alt[:] = angles.to_rad(source_desc[1]), angles.to_rad(source_desc[2])
    else:
        # For high precision work, we can only trust the az/alt
        # returned by actEphem.  This is because the equatorial
        # coordinates it knows about either include refraction
        # correction, or are geocentric.  We need topo-centric
        # coordinates (because sometimes Mars is nearby) and we want
        # to do our own atmospheric correction.

        # Set pressure to 0 to suppress refraction correction
        planet = source_desc[0].capitalize()
        ae = moby2.ephem.ACTEphem()
        ae.site.pressure = 0.
        for i in range(ctimes.shape[0]):
            ae.set_ctime(ctimes[i], fixDUT1=True)
            obj = ae.get_object(planet)
            az[i], alt[i] = obj.az, obj.alt
        
    ra, dec = moby2.pointing.get_coords(ctimes, az, alt, weather=(0.,0.,0.,0.))
    ra.shape, dec.shape = orig_shape, orig_shape
    return ra, dec
