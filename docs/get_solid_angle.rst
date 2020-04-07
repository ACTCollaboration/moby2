.. -*- mode: rst ; mode: auto-fill -*-

===================================================
``get_solid_angle`` - Compute basic beam parameters
===================================================

This script lives in analysis.beam_ana.  Its purpose is to analyze
planet maps and produce some basic beam statistics, such as FWHM and
solid angle.  This code should give fairly high fidelity results,
suitable for precision calibration work but somewhat less
sophisticated than a full beam analysis.

Analyzing a set of maps
=======================

Get a recent config file from
http://abs-inthe.princeton.edu/~mhasse/distribution/beams/

To get the latest version (for which the following information will
hopefully apply), execute this in your working directory:

.. code-block:: bash
  
 wget http://abs-inthe.princeton.edu/~mhasse/distribution/beams/get_solid_angle.zip
 unzip get_solid_angle.zip

Copy the input file from the unzipped folder.  Edit the file to point
to your beam maps.  If you want the code to make frequency
band-dependent decisions on the maps, the filename should contain the
frequency code (e.g. f090, f150, etc.).

.. code-block:: bash
  
  moby2 get_solid_angle get_solid_angle.in

This produces a few plots for each beam map, and an output file
`table.fits` with the analysis results for each map.  The FITS file
can be loaded through StructDB and reduced further, e.g. to get an
average solid angle and associated error estimate.

Configuration options
=====================

The main groups in the configuration file are:

* ``moby_options``
* ``source_maps``
* ``analysis``
* ``output``


Configuration: ``source_maps``
------------------------------

``'source': <map_list_args>``

  This is passed to get_map_list.  The most probably thing here is
  something along the lines of ``('glob', './maps/map*.fits')``.


``'basename_extractor': <extractor>``

  Specify the extractor to use for creating a map name from the map
  filenames.  If you're globbing for files, and the path/filename
  contain both the map basename and the frequency, then use
  'standard_freq'.  If you're working with combined files where the
  path/filename doesn't contain a TOD name, you might just want to use
  'filename'.


Configuration: ``analysis``
------------------------------

``'mask_radii_arcmin': <tuple>``

  This is a tuple of 3 radii, in arcminutes, which are used to rough
  out the data.  The first number is the radius within which the rough
  solid angle is computed, after subtracting the baseline estimate
  from data between the first and second radius.  The third radius
  number is used to limit the data considered in all subsequent
  analyses and plotting.


``'wing_fit_radii_arcmin': <2-tuple>``

  The 1/theta^3 wing is fit to all data between the two radii
  specified here.  This is one you might want to specify by frequency
  range, e.g.: ``{'f090': (2.0, 4.0), 'f150': (1.5, 3.0), 'f220': (3.,
  5.)}``.

``'wing_fit_ellipticity': <boolean>``

  When the wing is fit, you can do it with all the data points and a
  model that includes ellipticity, or you can bin up the points and
  disregard ellipticity.  For low signal-to-noise situations you
  probably want to set this to False.

``'force_frequency': <freq_string>``

  In the case that the frequency code (e.g. f090) cannot be determined
  from the filename, you can set it here.

``'recenter': <boolean>``

  If True, map coords are adjusted to remove the beam centroid.
  Defaults to True.
