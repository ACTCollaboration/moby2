.. -*- mode: rst ; mode: auto-fill -*-

=======================================
``fp_fit`` - Focal plane fitting module
=======================================

``fp_fit`` is a driver program for the fp_fit analysis module that can
be used to determine the relative detector offsets and optical time
constants for an array of detectors, by fitting a beam model to time
stream data of a planet observation.

In the single TOD fitting mode, ``fp_fit`` updates detector parameters
by stepping through a sequence of operations specified by the user.

The operations supported by ``fp_fit`` are:

* ``load_tod`` - load detector data.  The ``fit_coarse`` and
  ``fit_model`` steps require detector data.
* ``fit_coarse`` - obtain rough detector positions by projecting each
  detector time-stream into a map using the TOD boresight pointing.
* ``fit_model`` - fit a beam model to the detector time-stream data to
  obtain detector positions and time constants.
* ``load_fits`` - load detector positions from a file (the output from
  ``fit_coarse`` may have been saved to a file previously, for
  example).
* ``match_to_template`` - load a pointing template of some kind, and
  try to model the current detector offsets by shifting and possibly
  distorting the template.
* ``apply_cuts`` - reject position fits if fit parameters, such as the
  apparent source amplitude, do not fall within a certain range.
* ``average_fits`` - combine multiple fits and obtain average detector
  positions and time constants.

To analyze a single planet observation, ``fp_fit`` is typically run
several times, with different configuration files, to advance the fit.


Fitting a single planet observation
===================================

If a planet (or other point source) is bright enough to show up with
high significance in each detector of a single TOD, the ``fp_fit``
program can be used to determine the detector pointing offset and time
constant.

The default configuration file ``fp_fit.in`` causes ``fp_fit`` to do a
full, minimal-priors analysis of a single TOD whose filename is passed
in on the command line:

.. code-block:: bash
  
  fp_fit my_fp_fit.in path/to/tod

The results of the fit are written to files in ``./fp_output/``, in
the ASCII format used by the ``FPFitFile`` class.

Note that ``fp_fit`` accepts some standard ``moby2`` scripting
options, such as ``-v <verbosity>``, and ``-o <output_prefix>``.


Useful ``fp_fit`` configuration parameters
------------------------------------------

``'source_tods': <instructions>``

  This can be used to provide a file with a TOD list, instead of
  passing TODs on the command line.  

``'forks': <n_threads>``

  Time-consuming operations can be accelerated through trivial
  parallelization by setting ``forks`` to the number of threads.  The
  program does not currently use MPI or openMP; it simply forks child
  processes and manages IPC using the file system.

``'prefetch': <True|False>``

  On a system with sufficient RAM, ``prefetch=True`` can be used to
  load a TOD in the background while another TOD is being analyzed.


Using an existing pointing template to guide fitting
----------------------------------------------------

To increase the yield of time stream fitting, a template can be used
to provide a good starting point for the position fitter.

This mode is configured by changing the ``steps`` setting in the
``fp_fit`` section of the config file.

.. code-block:: python
  
  ...
  ## Default steps:
  #'steps': [ 'load_tod', 'fit_coarse', 'apply_cuts', 'fit_model',
  #           'apply_cuts'], 
  # Match to template:
  'steps': [ 'load_tod', 'fit_coarse', 'apply_cuts',
             'match_to_template', 'fit_model', 'apply_cuts'],
  ...

The ``match_to_template`` step uses parameters described in the
``template_opts`` section of the configuration block:
 
.. code-block:: python
  
  ...
  'template_opts': {
    'source': 'file',
    'filename': 'my_template.txt',
    'shift': True,
    'rotation': True,
    'scale': True,
    'shear': False,
  }
  ...

In addition to the template filename, the parameters ``shift``,
``rotation``, ``scale``, and ``shear`` specify which transformations
are permitted when fitting the template to the coarse detector
positions.

When fitting with a realistic template (i.e. one based on previous
observations rather than a mock-up of the focal plane), only the
``shift`` should be needed to get a good fit.

The ``summary_output`` parameter can be used to specify a filename
where the resulting model fit parameters (i.e., the shifts, rotation
angles, etc.) should be written.

In this mode, you can speed things up by configuring the
``coarse_fit`` step to ignore most of the detectors:

.. code-block:: python
  
  'coarse_opts': {
  ...
    # Fit one tenth of the detectors.
    'decimation': 10


Combining several pointing fits to get a pointing template
==========================================================

Set ``source_tods`` to read the list of relevant TODs from a file.
(It's probably enough for the list to contain basenames, since the
TODs will not get loaded):

.. code-block:: python

  'source_tods' = {
    ...
    # List of filenames / basenames as ('tod_list', filename, column)
    'source': ('tod_list', 'my_tod_list.txt', 0),


Set the ``steps`` parameter to load existing fits, and to perform averaging:

.. code-block:: python

  'fp_fit': {
     ...
     'steps': ['load_fits', 'average_fits'],

The ``average_opts`` parameters may also need tweaking.  The averaging
includes outlier rejection, so ideal parameters may depend on how many
TODs you are combining.

The ``average_fits`` step also produces plots of the detector
positions and histograms of the positions and time constants.


Option description for individual operations
============================================

This section is not complete; we focus on opaque / tweakable /
dangerous parameters.

``load_tod``
------------

* ``detrend: <detrend?>`` - boolean; causes removal of the trend line
  from each detector.
* ``scale: <rescale_factor>`` - float; number by which to rescale the
  data (brings plot units and planet amplitudes into some
  comprehensible range).
* ``poly_order: <order>`` - integer; order of polynomial to remove.
  Comment out to suppress polynomial subtraction.
* ``poly_decimation: <n>`` - integer; fit polynomial to every nth
  sample of TOD only; acceleration.  100 is fine.
* ``guess_sign: <guess?>`` - boolean; flip the sign of a time stream
  if its mean lies closer to its maximum than to its minimum.
* ``common_mode: {'source': <source_spec>, 'fit': <fit?>}`` - For
  common mode removal.

  * ``<source_spec>``: Only ``('file', <filename>)`` is currently
    supported.  The file must be an ascii file with a single column
    and the same number of samples as the TOD.
  * ``<fit?>`` - when set to True, the common mode will be fit to the
    TOD prior to subtraction.  False isn't that useful unless you have
    a calibrated and flat-fielded TOD, which we probably don't at this
    stage.


``fit_coarse``
--------------

* ``resolution: <delta>`` - sets the map resolution, in arcminutes.
* ``decimation: <n>`` - perform coarse fit on only every nth detector;
  used on high signal TODs in combination with ``match_to_template``.
* ``filter_band_pass: (<f_center>, <f_sigma>)`` - causes a gaussian
  band pass filter to be applied to the TOD (for the ``fit_coarse``
  mapping step only).  A nice way to kill atmosphere and white noise
  to pick out the planetary signal.

