.. -*- mode: rst ; mode: auto-fill -*-

Database access
===============

The ``tod_info`` Table
----------------------

The ``tod_info`` table stores lots of useful information about the merged
TODs.  To get some tod_info, run:

.. code-block:: python

  from moby2.instruments import actpol
  db = actpol.TODDatabase()

  db.get_record('1378762587.1378762626.ar1')

To see a summary of all the information in a TODRecord, use the
``to_string`` method:

.. code-block:: shell

  print tod_record.to_string()

  mce_data_acq_id : 6864L
  array           : 'AR1'
  mean_alt        : 50.0411
  rms_alt         : 0.00013
  min_az          : 276.919
  max_az          : 281.514
  scan_speed      : 1.32151
  pwv             : None
  min_dec         : -13.0941
  max_dec         : -10.3552
  min_RA_south    : 204.329
  min_RA_north    : 205.45
  max_RA_south    : 206.397
  max_RA_north    : 207.518
  ctime_start     : 1378762588L
  ctime_end       : 1378763082L
  obs_type        : 'planet'
  obs_detail      : 'venus'
  obs_drift       : 'setting'
  malfunction     : None
  datafile_id     : 20597L
  basename        : '1378762587.1378762626.ar1'
  ctime           : 1378762587L
  merged_ctime    : 1378762626L

Each field listed is available as an attribute in the TODRecord
object.  All fields except for datafile_id, basename, ctime, and
merged_ctime can be used for selecting TODs using keyword arguments to
the ``select_tods`` function as described in the next section.


Selecting TODs
--------------

The ``select_tods`` method returns a TODRecordList of TODRecords.

.. code-block:: python

  # Get a list of TODRecord instances
  tod_list = db.select_tods(obs_detail='venus')

  # Print basename of last item
  print tod_list[-1].basename  

The TODRecordList is like a normal list, but an additional
``get_fields`` method returns an array of values for a desired
attribute (or list of attributes):

.. code-block:: python

  ctime = tod_list.get_fields('ctime_start')
  min_dec, max_dec = tod_list.get_fields(['min_dec', 'max_dec'])
  print 'Venus mean declination: ', (min_dec + max_dec).mean() / 2

Keyword arguments passed to the ``select_tods`` method are interpreted
as expressions to match in the tod_info table.  In addition to simple
matching to a single value:

* If a **list** is passed, then records will be return that match any
  one of the values listed.
* If a **tuple** is passed, then it should have a length two and will
  be interpreted as describing a range of values to match.

For ctime fields, you can pass in a date string like '2013-09-21
12:00:00' instead of an integer or float and it will be converted
automatically.

Examples:

.. code-block:: python

  # Get all TODs for deep fields 1 or 2
  tod_list = db.select_tods(obs_detail=['deep1','deep2'])

  # Get all planet TODs with 1234560000 <= ctime_start < 1234570000
  tod_list = db.select_tods(obs_type='planet',
                            ctime_start=(1234560000, 1234570000))

  # Get all stare TODs from Sep 15 UTC:
  tod_list = db.select_tods(obs_type='planet',
                            ctime_start=('2013-09-15 00:00:00',
			                 '2013-09-16 00:00:00'))


Getting filenames
-----------------

Datafile locations are stored in the ACTpol database.

.. code-block:: python

  from moby2.instruments import actpol
  filebase = actpol.Filebase()

  # Find the filename associated with my Venus observation
  filename = filebase.filename_from_name('1378762587.1378762626.ar1', single=True)
  print filename

The ``actpol.Filebase`` configuration can be passed in to the
constructor, or provided in the .moby2 file with a block like:

.. code-block:: shell

  # Config for actpol.Filebase
  actpol_filebase = {
      'manifest_config': '/ACTpol/etc/manifest.conf',
      'nodes': [('merlin','raid')
                ('mce1','raid',None,'prefix5')],
  }
  
  # Use actpol.Filebase as the default filebase
  filebase = {
      'type': 'actpol_manifest',
  }

The second block is used by scripts to decide how to instantiate a
default filebase.


Finding Calibration and Housekeeping Data
-----------------------------------------

SQUID tuning, IV curve and bias step data from the MCE are not merged
into dirfiles.  These acquisitions are recorded in the
``mce_data_acq`` table of the database and the data are stored in the
raw data archive.

To find calibration data, use ``select_acqs``.  The important things
to tell it are the data ``suffix`` and ``array``.  You can also filter
on ``ctime``, as one would when using ``select_tods``.

Popular values for suffix:

* ``iv``: IV curves.  The IV analysis (extension .out) is included
  with the data.
* ``bc1_step``, ``bc2_step``, ``bc3_step``: Bias steps.
* ``sqtune``: Analysis block from SQUID setup.
* ``RCs_ssa``, ``RCs_sq2servo``, ``RCs_sq1servo``, ``RCs_sq1ramp``,
  ``RCs_sq1rampc``: SQUID tuning raw data.

Example:

.. code-block:: shell

  from moby2.instruments import actpol
  database = actpol.TODDatabase()
  filebase = actpol.Filebase()
  
  ctimes = ('2013-09-18 00:00:00', '2013-09-18 08:00:00')
  ar = 'ar1'
  recs = database.select_acqs(suffix='iv', ctime=ctimes,
                              array=ar)
  
  print filebase.filename_from_name(recs[0].mce_name, array=ar)


Auto-doc
--------

.. autoclass:: moby2.instruments.actpol.TODRecord
   :members:

.. autoclass:: moby2.instruments.actpol.TODDatabase
   :members:

.. autoclass:: moby2.instruments.actpol.Filebase
   :members:

