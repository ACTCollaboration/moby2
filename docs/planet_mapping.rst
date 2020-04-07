.. -*- mode: rst; mode: auto-fill -*-

=======================================
``planet_mapping`` - Mapping of planets
=======================================

The scripts ``bin/planetCuts`` and ``bin/planetMap`` can be used to
process planet observation time streams to produce maps.

When running these scripts: the template configuration files should be
copied into a working directory and edited for the purpose at hand.
Below, we first describe the two scripts and how they are invoked.
After that there is a description of important configuration
parameters.

.. Link test: :ref:`get-focal-plane`.

planetCuts
==========

``planetCuts`` processes one or more merged TOD and for each one
writes some important results to disk.  These include:

* Results to be consumed by planetMap (or other interested agents):

  * Time-stream cuts (TODCuts).
  * Time-stream calibration (Calibration).
  * Common mode set, with couplings to the calibrated data
    (CommonMode).

* Summary and debugging plots and statistics to be inspected by the
  analyzer:

  * Plots of the common mode, dark modes, surviving detectors, etc.
  * Statistics on how many detectors survived successive cuts.

The program is invoked as

.. code-block:: bash

  planetCuts [-o <output_prefix> ] [-v <verbosity>] [-i] \
    planetCuts.in [tod_filename, ...]

Some useful options:

* -i: causes interactive debugging, which means that python exceptions
  cause a python debugger session to be initiated in the context of
  the fault (so local variables can be inspected).
* -v: set the verbosity level.  Higher numbers are more verbose; 4
   gives really quite a lot of detail.
* -o: override the output prefix from the configuration file.

The ``tod_filename`` can often just be the tod_name
(e.g. 1234567890.1234568000.ar2) that you're interested in.  Note that
planetCuts won't necessarily read the TODs from the command line --
that's a setting you can manipulate.

The next sections are about the planetCuts configuration file.

Configuration file overview
---------------------------

The top-level entries in the planetCuts configuration file are:

* moby2_options - common block specifying default verbosity and
  logging behaviors.
* output -- where to write outputs, and what tag names to give to
  various results.
* source_tods -- how to find TODs, restrictions on what detectors or
  samples should be loaded.
* pointing -- how to load the pointing information for a TOD.
* time_constants -- how to load the time constants.
* calibration -- what calibration operations should be performed.
* planet_cuts -- parameters for the various cuts.
* common_mode -- parameters for computing the common mode and for the
  calibration / cuts associated with common_mode correlation.

Input data and meta-data
------------------------

To specify the TOD list, see :ref:`get-tod-list`.


Tuning the cut options and parameters
-------------------------------------

Control of output data locations and formats
--------------------------------------------

By default, the cut results are written to a depot in ``./depot/``,
under generic tags "position_cuts" and "planet_cuts".  The tags and
depot indicated in the ``output`` block of the parameters file should
agree with those in the planetMap configuration file.

The ``planetCuts`` script is typically run on a single merged TOD,
though multiple TODs can be specified on the command line.  The cut
specifications are given in a configuration file like :download:`planetCuts.in
</_moby2/params/planetCuts.in>`.


planetMap
=========

``planetMap`` reprocesses the TOD and creates a maximum likelihood
map.  The map is written in FITS format, along with associated
weights.  Presently it works best when making maps in "source_scan"
coordinates.  These are tangent plane coordinates centered on the
ephemeris position of the planet, rotated so that the X direction is
parallel to the azimuth scan.

Use a configuration file like :download:`planetMap.in
</_moby2/params/planetMap.in>`.  TODs should be passed in on the
command line, or through a tod_list file.  Multiple TODs will be
mapped simultaneously into a single map... so you probably just want
to pass a single TOD.

The planet mapper consumes the output of ``planetCuts``; make sure the
depot and tag settings in the ``tod_cuts`` and ``position_cuts``
blocks point to the right place.

Invocation:

.. code-block:: bash

  planetMap [-o <output_prefix> ] [-v <verbosity>] [-i] \
    planetMap.in [tod_filename, ...]


Configuration file parameters
=============================

The ``planetCuts.in`` and ``planetMap.in`` parameter files provide
important inputs, such as detector positions, cuts file locations,
etc.  The syntax is meant to be quite flexible, to provide the user
with the ability to jam in their own data whenever it might be
convenient.

A description of certain configuration blocks is outlined below.

Relative detector positions
---------------------------

To use an FPFitFile for relative detector offsets, the ``pointing``
block should look like

.. code-block:: python

  pointing = {
      'source': 'fp_file',
      'filename': 'path/to/fp_output.txt',
  }


Once moby2 has been configured with templates, the pointing block can
simply read:

.. code-block:: python

  pointing = {
      'source': 'template',
  }

To load the ``det_uid``, ``x`` and ``y`` positions from an ascii file,
provide the filename and the indices of the columns containing the
required data:

.. code-block:: python

  pointing = {
      'source': 'columns_file',
      'filename': 'path/to/columns_data.txt',
      'columns': [0,1,3],
  }

When loading from a columns_file, the ``x`` and ``y`` coordinates
should be provided in radians.


Detector time constants
-----------------------

To specify detector time constants, use the same basic syntax as for
the relative detector positions.  But put it in the ``time_constants``
block.  E.g., to load time constants from an FPFitFile:

.. code-block:: python

  time_constants = {
      'source': 'fp_file',
      'filename': 'path/to/fp_output.txt',
  }

When providing time constants in a columns_file, provide the time
constant in seconds.

Calibration
-----------

The ``calibration`` block is currently ignored by ``planetCuts``, but
is essential to ``planetMap``.  It describes how the mapper should
convert from readout DAC units to physical units.  For planet maps
we've historically just gone to detector pW, possibly with a flat
field applied.  Calibration consists of a series of operations, which
are specified in a list.  For example, to take IV calibration from the
TOD's runfile and then apply the template flat field:

.. code-block:: python

  calibration = {
      'cal_steps': [
          # Use runfile IV responsivity to start
          {'name': 'IV analysis',
           'type': 'iv', 
           'source': 'data', }, 
          # Apply a flatfield
          {'name': 'Flat field',
           'type': 'flatfield',
           'source': 'template',
           },
          ]
      }

No block is mandatory; the calibration factor will default to 1 if no
steps are given.  Other useful steps:

.. code-block:: python
  
  # Calibration, or recalibration, in a calgc-style / ACTDict file:
  {'name': 'my calgc',
   'type': 'cal_dict',
   'filename': 'my_calgc.dict',
   }
  
  # Calibration factors in an asciifile (provide columns for det_uid
  # and cal factor):
  {'name': 'my recal',
   'type': 'columns_file',
   'filename': 'my_recal.txt',
   'columns': [0,1],
   }

Note that inserting arbitrary calibration factors can be used to cope
with unexpected calibration weirdness, such as correcting the sign of
the detector response (so we don't get negative signal planets).
Applying a cal factor of 0 at this stage will cause detectors to be
ignored in the light/dark mode analysis, effectively cutting them from
the map.
   
Cuts
----

``planet_map`` loads separate cuts for the planet position and for
general sample and detector masking.  These can be loaded from a depot
(the default) or can be specified as files.  Currently moby2 defaults
to the ACT cuts format.

Loading cuts from the depot is achived with a configuration block
like:

.. code-block:: python
  
  tod_cuts = {
      'source': 'depot',
      'depot': { 'path': './depot',
                 'act_depot': False },
      'tag': 'planet_cuts',
  }

If you want instead to load the cuts directly from some file:

.. code-block:: python
  
  tod_cuts = {
      'source': 'file',
      'filename': 'path/to/cuts.txt',
  }

Output Maps
-----------

I'm pretty sure we can only write a single map right now, but there's
infrastructure ready to go for multiple maps.  The map parameters are
specified as list entries in the `maps` configuration block.  For
example:

.. code-block:: python
  
  maps = [
      ('source', {
              'coords': 'source_scan',
              'pitch': (3.75/3600),
              'center': (0., 0.),
              'size': (2., 1.5),
              # If source coords could be ambiguous, provide them here
              'source_name': None,
              }
       )]

The parameter ``pitch`` is the pixel pitch, in degrees.  It defaults
to 3.75 arcseconds.

The ``center`` and ``size`` parameters should either both be provided,
or both commented out / set to None.  When not provided, the map will
be made large enough to contain all the data.

The ``coords`` parameter can (at present) be set to:

* ``source``: tangent plane at source position, with X anti-parallel to
  J2000 RA.
* ``source_scan``: like ``source``, but rotated so X is parallel to
  the scan direction

If you want the map to be centered on something other than the
program's best guess, specify ``source_name``, using one of the following
forms:

* ``'source_name': 'Saturn'``
* ``'source_name': ('J2000', 171.23, 5.19)`` -- with RA and dec in
  degrees.

