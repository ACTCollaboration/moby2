.. -*- mode: rst ; mode: auto-fill -*-

Calibration Objects
===================

Overview
--------

The main difficulty in calibrating a TOD is not the multiplication, but the
actual retrieval and construction of the correct calibration numbers to use.
The swiss-knife function for assembling calibration vectors is
``scripting.get_calibration()``.  The goal of this function is that all
calibration can be specified in a calibration parameter block, suitable for
``MobyDict``-style configuration files.  Example:

.. code-block:: python

  # Import moby2 and get the local filebase
  import moby2
  fb = moby2.scripting.get_filebase()
  # Load a particular planet TOD from season 2.
  filename = fb.get_full_path('1412925756.1412925775.ar2')
  tod = moby2.scripting.get_tod({'filename': filename})
  # Get a calibration object, requesting the basic IV cal stored with the data.
  cal = moby2.scripting.get_calibration({'type': 'iv', 'source': 'data'}, tod=tod)
  # Apply the calibration.
  cal_mask, cal_val = cal.get_property('cal', det_uid=tod.det_uid)
  tod.data *= cal_val[:,None]

This example is quite generic aside from the call to ``get_calibration``,
where a particular calibration step has been requested:

.. code-block:: python

  {'type': 'iv',
   'source': 'data'}


To do other kinds of TOD calibration, see the recipes below.


Recipes
-------

Here we provide examples of calibration parameter blocks that may be passed as
the first argument to ``get_calibration``.

The code currently accepts a few other weird `type=...` options.  Note that
most parameters that accept a `filename=...` argument will decode the filename
relative to a depot if you specify a depot through the `'depot'` key.

**Multiple calibration steps:** To apply several calibrations in sequence,
  pass a list instead of a dict as the first argument of ``get_calibration``.
  Each item of the list must also be a valid calibration parameter block.

.. code-block:: python

  [{'type': 'iv',
    'source': 'data'},
   {'type': 'constant',
    'value': 1e7}]

**IV calibration:** Usually you will want `source='data'`, but if you have an
  MCE IV analysis file (or other runfile) that you want to use you can put
  that filename into attribute `filename` and set `source='runfile'`.

.. code-block:: python

  {'type': 'iv',
   'source': 'data'}

**Constant value:** For example, to scale from pW to approximate uK units.

.. code-block:: python

  {'type': 'constant',
   'value': '1e7'}

**ArrayData parameter:** This exists primarily for the purposes of
  automatically applying the ``optical_sign`` that is recorded in the
  ArrayData.  But in principle you can look up any value stored in ArrayData.

.. code-block:: python

  {'type': 'array_data_parameter',
   'parameter': 'optical_sign'}

**ASCII data:** Must be indexed by det_uid.  Provide the column index for the
  det_uid and calibration number columns.

.. code-block:: python

  {'type': 'columns_file',
   'columns': {'det_uid': 0, 'cal': 1},
   'filename': 'my_flatfield.txt'}

**ABSCal ASCII data:** For looking up the TOD basename in a list and applying
  a single recalibration factor for all detectors.

.. code-block:: python

  {'type': 'abscal_columns_file',
   'columns': {'det_uid': 0, 'cal': 1},
   'filename': 'my_per_tod_calibration.txt'}

**Referenced calibration:** This is basically designed to load a bias_step
  calibration, given a file that associates time ranges of acquisition to a
  particular calibration result.  See ``scripting.get_referenced_calibration``
  for more details, or seek advice.

.. code-block:: python

  {'type': 'referenced_cal',
   'assignments': {'depot': {'label': 'actpol_shared'},
                   'structure': './BiasStepTimes/intervals_{season}_{array}_{tag}.txt',
                   'tag': '150131',
                   'columns': [0,1,4]},
   'library': {'depot': {'label': 'actpol_egrace1'},
               'structure': 'calibration/bias_step/{tag}.cal',
               'type': 'cal_dict',
               }
   }


**Removing the readout filter gain:** Divide out the DC gain of the MCE's
  readout filter.  Note that this is not necessary in combination with
  `type='iv'`.  But it will often appear in combination with bias_step-derived
  calibartion results.

.. code-block:: python

  {'type': 'remove_readout_filter_gain'}

