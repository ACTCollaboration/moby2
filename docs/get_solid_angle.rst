.. -*- mode: rst ; mode: auto-fill -*-

===================================================
``get_solid_angle`` - Compute basic beam parameters
===================================================

This script lives in ``analysis.beam_ana``.  Its purpose is to analyze
planet maps and produce some basic beam statistics, such as FWHM and
solid angle.  This code should give fairly high fidelity results,
suitable for precision calibration work but somewhat less
sophisticated than a full beam analysis.

Analyzing a set of maps
=======================

An example configuration file can be found below.  To use it as is,
create an ascii file called map_list.txt with the path to each map you
want to analyze.  If the band code (e.g. f150) and TOD name are not
obvious from the filename, then you might need to use a 2-column
format that gives the basename and frequency separately.  See detailed
configuration options below.

To run the script:

.. code-block:: bash
  
  moby2 get_solid_angle get_solid_angle.in

This produces a few plots for each beam map, and an output file
`table.fits` with the analysis results for each map.  The FITS file
can be loaded through StructDB and reduced further, e.g. to get an
average solid angle and associated error estimate.

Some basic summarization can be triggered through the output.summaries
configuration settings; see below.


Configuration options
=====================

For a starting point, see this example: 
:download:`get_solid_angle.in <params/get_solid_angle.in>`.

The main groups in the configuration file are:

* ``moby_options``
* ``source_maps``
* ``analysis``
* ``output``

The configuration file is parsed by the ``solid_angle``
:meth:`~moby2.analysis.beam_ana.util.driver` method.


``source_maps`` block
---------------------

``'source': <map_list_args>``

  This is a tuple that is used as the ``params`` for
  :meth:`~moby2.analysis.beam_ana.util.get_map_list`.  An example is
  ``('file', 'map_list.txt')``.


``'basename_extractor': <extractor>``

  Specify the extractor to use for creating a map name from the map
  filenames.  If you're globbing for files, and the path/filename
  contain both the map basename and the frequency, then use
  ``'standard_freq'``.  If you're working with combined files where the
  path/filename doesn't contain a TOD name, you might just want to use
  'filename'.


``analysis`` block
------------------

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


``output`` block
----------------

``'prefix'': <path-like>``

  Sets the prefix for output (note you need a trailing ``/`` if it's a
  directory).  A bunch of plots, a log of the fit output, and
  ``table.fits`` will be written there.

``'summaries': <list of dicts>``

  A list of instructions for summary computations.  Each dict is used
  as parameters for the
  :meth:`~moby2.analysis.beam_ana.solid_angle.summary_op` function.
  Probably the only think you need to worry about is the key ``'select'``,
  the value of which will be passed to
  :meth:`~moby2.analysis.beam_ana.util.get_obs_mask` to restrict the
  set of maps included in each average.


Reference
=========

Selected methods from beam_ana.solid_angle
''''''''''''''''''''''''''''''''''''''''''

.. autofunction:: moby2.analysis.beam_ana.solid_angle.driver

.. autofunction:: moby2.analysis.beam_ana.solid_angle.summary_op

Selected methods from beam_ana.util
'''''''''''''''''''''''''''''''''''

.. autofunction:: moby2.analysis.beam_ana.util.get_map_list

.. autofunction:: moby2.analysis.beam_ana.util.get_obs_mask
