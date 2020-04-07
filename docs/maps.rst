.. -*- mode: rst ; mode: auto-fill -*-

==========================
Maps - fitsMap and friends
==========================

``fitsMap`` is a bit of a mess, but has some powerful features and the
TOD-to-map projection routines are set up to work with it.  The module
can handle common ACT WCS, and can read and write FITS files as well
as create new, blank maps.

There are two main kinds of ``fitsMap``: ``spaceMap`` for real space
maps and ``freqMap`` for spatial frequency domain maps.  Some of the
useful interfaces these expose are described below.

Because FITS uses degrees, so do we.

Basic operation
===============

The first argument to the ``spaceMap`` for ``freqMap`` constructors is
assumed to be a filename.  If a file with the same root but extension
'.weight.fits' is present, it is loaded as the map weight.  The
weight=... keyword can be set to the filename of an alternative
weights map, or may be set to False to suppress tracking of map
weight.

.. code-block:: python
  
  import moby2
  map0 = moby2.mapping.fits_map.spaceMap('source.fits')

Inspect data:

.. code-block:: python

  print map0.data.shape
  (2299, 4841)

Alter data, save map:

.. code-block:: python
  
  map0.data[:] = 0.
  map0.write('source_zero.fits')

The ``fitsMap.imshow`` method renders the map on the current pylab
axes.  It accepts several useful options; see source.

.. code-block:: python

  import moby2
  import pylab as pl
  
  map0 = moby2.mapping.fits_map.spaceMap('source.fits')
  map0.imshow(units='arcmin')
  map0.zoom(10.)
  pl.savefig('map0.png')

To create a new tangent plane map:

.. code-block:: python

  x_lims = (-.5,.5)
  y_lims = (-.5,.5)
  pitches = (30./3600, 30./3600)
  map1 = moby2.mapping.fits_map.simpleMap(x_lims, y_lims, pitches)


Coordinates and pixelization
============================

WCS and pixelization information (basically the FITS header) are
encoded in the ``pixn`` attribute.

.. code-block:: python

  print map.pixn.bounds()
  (array([-2.42760417,  2.61510417]), array([-1.52864583,  0.86614583]))

Note that bounds returns an array of X positions followed by an array
of Y positions.

Converting between pixel index and WCS coordinate is achieved using
methods of the form ``<sys1>_to_<sys2>``.  For example:

.. code-block:: python

  print map0.pixn.sky_to_IJ(0.,0.)
  (2330.0, 1467.0)

The token "IJ" is used to indicate the pixel coordinate system,
(column, row).  Note that I is the column index, and J is the row
index.  But fitsMap.data and fitsMap.weight are 2-d arrays indexed by
[row,column].  In this case, for example, the coordinate (0., 0.) is
associated with map0.data[1467,2330].

For simple projections such as TAN and CEA, the spaceMap.x and
spaceMap.y vectors may be a simpler way to quickly compute map
coordinates at each pixel.  The ``spaceMap.radii(center=(0.,0.))``
function returns the distances of each pixel to the indicated
coordinate.

Copies and extracted maps
=========================

A map can be copied by calling its copy method, e.g. ``map1 =
map0.copy()``.

The extract method can be used to cut out a small section of a map,
and return a new map with the coordinate system appropriately updated.

.. code-block:: python

  map1 = map0.extract((-.2, 2), (-.1,.1), coords='sky')
  print map1.pixn.bounds()
  (array([-0.20052083,  1.99947917]), array([-0.10052083,  0.10052083]))


Fourier space
=============

You can get a fourier space representation of the map by evaluating:

.. code-block:: python

  fmap0 = map0.fftMap()

The fmap0.data will be complex valued.  To plot it, be sure to pass in
the absolute value or something:

.. code-block:: python

  fmap0.imshow(data=abs(fmap0.data), units='ell')
  fmap0.zoom(10000)
  pl.savefig('beam_ell_space.png')

Note that fmap0.x, fmap0.y, fmap0.radii() work with freqMaps.  But
they probably return k vector coordinates that are conjugate to
degrees.  To convert to multipole ell, multiply by 180/pi.

You can filter a map and transform it back to real space.  For
example, to filter out signal below ell of 1000:

.. code-block:: python

  map0 = moby2.mapping.fits_map.spaceMap('source.fits')
  fmap = map0.fftMap()
  ell = fmap0.radii() * 180 / np.pi
  fmap.data *= (ell > 1000)
  map1 = fmap.ifftMap()


Sub-pixel operations
====================

Using Fourier techniques, we can...

Shift an image by a non-integral number of pixels.  Note this will
cycle pixels from one side of the map into the other side -- so be
careful.

.. code-block:: python

  # Load a map
  map0 = moby2.mapping.fits_map.spaceMap('source.fits')
  # Find the position of the peak or something...
  x0, y0 = -0.05, +0.02
  # Get the shifted map.
  map1 = moby2.mapping.fits_map.shiftImage(map0, -x0, -y0)

