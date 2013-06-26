.. -*- mode: rst; mode: auto-fill -*-

========
Pointing
========

Coordinate Systems
==================

``libactpol`` (which is an independent work upon which ``moby2`` depends)
provides astrometric conversions from detector coordinates to various
celestial systems, including telescope horizon coordinates, J2000
equatorial coordinates, and galactic coordinates.  ``moby2`` makes these
computations available at the python level.

For astrometry, ``moby2`` uses radians as the native units for angles
and the unix time stamp ("ctime") for time vectors.

Detector positions
==================

The focal plane coordinates used by ``moby2`` differ slightly from the
conventions used in libactpol.  In both cases the boresight is taken
to lie along z^.  In ``moby2``, the x^ direction is such that a detector
with positive x coordinate will project to a point on the sky that is
to the right of the boresight.  Detectors with positive y^ will lie
above the boresight, when projected onto the sky.

Detector positions are encapsulated in :doc:`FocalPlane objects
<focal_plane>`.


High-precision Pointing
=======================

``libactpol`` performs high precision computations that take the focal
plane position and orientation of detectors to J2000 coordinates.

To access these transformations in ``moby2``, use the pointing module.  In
``moby2`` the telescope boresight coordinates are described by azimuth and
altitude angles, and a unix timestamp ("ctime").  The J2000 RA and dec
can be obtained like this:

.. code-block:: python

  import moby2
  import numpy as np
  
  ctime = np.array([1234560000.])
  az = np.array([10.]) * np.pi/180
  alt = np.array([50.]) * np.pi/180
  ra, dec = moby2.pointing.get_coords(ctime, az, alt)

If you have detectors in a focal plane, you can obtain the RA, dec of
each detector at each time.  First initialize a FocalPlane object
somehow, then pass the FocalPlane to get_coords.

.. code-block:: python

  fp = moby2.pointing.FocalPlane(x=[0., 0.003, 0.], y=[0., 0., 0.006])
  ra, dec = moby2.pointing.get_coords(ctime, az, alt, focal_plane=fp)
  
Alternately, you can call get_coords on the FocalPlane object
directly:

.. code-block:: python

  fp = moby2.pointing.FocalPlane(x=[0., 0.003, 0.], y=[0., 0., 0.006])
  ra, dec = fp.get_coords(ctime, az, alt)


Accelerated Pointing
====================

Because the exact pointing is expensive to compute, ``moby2`` implements
an approximation that can determine the pointing of individual
detectors very rapidly if the pointing of the array center is already
known.  The approximation is called :doc:`ArrayWand <array_wand>`.

Note that the approximation degrades for detectors that are far (> 1
degree) from the focal plane center, or in parts of spherical
coordinate systems that are highly curved (i.e. near the poles).
Users may want to check the fidelity for their particular application.
For wands generated in tangent plane coordinates, or at J2000
latitudes visible from ACT, the accuracy should be better than 1".

To create an ArrayWand, use one of the factory classmethods.  For
planet mapping, for example, we often want tangent plane coordinates
centered on an ephemeris source position:

.. code-block:: python

  # Load a TOD
  tod = moby2.TOD.from_dirfile(filename)
  # Set the default focal plane
  tod.fplane = moby2.scripting.products.get_focal_plane(
       {'source': 'template'}, tod_info=tod.info)
  # Get a wand centered on a particular RA, dec:
  saturn_ra_dec = (.13, 0.01)
  wand = moby2.pointing.ArrayWand.for_tod_source_coords(tod,
       ref_coord=saturn_ra_dec)

Then, to get the detector pointing, pass the focal plane (which does
not need to match the one used initially) to wand.get_coords:
   
.. code-block:: python

  x, y = wand.get_coords(tod.fplane)


Pointing in Mapping Applications
================================

To project back and forth between a map and the time domain, create
the appropriate projection object using the following incantation.  It
is assumed that you already have a map, a focal plane and a pointing
wand in hand:

.. code-block:: python

  # Load tod...
  # Get a focal_plane...
  # Get wand...
  # Create map...
  grid = moby2.pointing.GridPixelization.forFitsMap(map)
  proj = moby2.pointing.WandProjector(wand, grid, focal_plane)

With this projector, you can project from any time-domain data with
the same size and dtype as tod.data into any map with the same size
and dtype as map.data:

.. code-block:: python

  # Project TOD into map
  proj.project_to_map(tod.data, map.data, map.weight)
  # Project map to TOD
  proj.deproject_from_map(map.data, tod.data)

In both cases, the data are added into the target array, without
first zeroing it.

To cut certain samples or channels from the projection, pass a TODCuts
object to the ``cuts`` keyword, for example ``cuts=tod.cuts``.


Alternative spherical coordinate systems
========================================

If one wants to use an output coordinate system that is not ICRF /
J2000, an arbitrary final rotation may be defined to be applied to the
libactpol computations.  This is specified using a quaternion, because
libactpol uses quaternion rotations internally.

This is useful, for example, if galactic coordinates or tangent plane
coordinates are desired.

To get a quaternion that rotates the coordinate system so that a
particular (RA, dec) will be mapped to the celestial north poll, for
example, use:

.. code-block:: python

  # Where is the telescope pointing?
  ra, dec = moby2.pointing.get_coords(ctime, az, alt)
  # Get rotation that re-centers on that point.
  gamma = 0.
  final_q = moby2.pointing.get_tracking_quaternions(ra, dec, gamma)
  # Check it...
  ra1, dec1 = moby2.pointing.get_coords(ctime, az, alt, final_q=final_q)
  

