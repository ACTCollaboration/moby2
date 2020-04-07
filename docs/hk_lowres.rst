.. -*- mode: rst ; mode: auto-fill -*-

Low resolution housekeeping data
================================

Overview
--------

For studies on long (>1 minute) time scales, the housekeeping data recorded by
AMCP are stored at reduced resolution in text format and can be loaded by
moby2.  The interface is based on the one used to load :doc:`APEX weather data
<weather>`.

The data that this module makes available must be compiled from the
full-resolution dataset.  This occurs automatically on raid for the current
season.  For past seasons the most up-to-date archive is probably on feynman.


Recipes
-------

See :doc:`APEX weather data<weather>` for more examples of how this sort of thing
works.  An example, though, is:

.. code-block:: python

  import moby2
  chan_el = moby2.aux_data.hk.HKChannel('Enc_El_Deg')
  el, dt = chan_el.get_nearest(1417046400)
  print (el, dt)
  (59.990000000000002, 33.0)

This means that the nearest data to ctime 1417046400 was 33 seconds later, and
at that time the Enc_El_Deg field was reading 59.99 degrees.

Note that the "dt" result effectively carries the validity of the data.  The
user must always check the time difference between requested and returned
data, to confirm it is within acceptable tolerances.  AMCP is not always
running, and even when it is it might have been set to not acquire the data
you are asking for.  So check dt to make sure it's within the tolerance of
your application.

A more complex example: Look up a bunch of values at once, and screen them for
validity using a boolean mask based on dt.

.. code-block:: python

  # Initialize the looker-up.
  import moby2
  import numpy as np
  chan_el = moby2.aux_data.hk.HKChannel('Enc_El_Deg')
  
  # Create a vector of ctimes covering October 2014, in 10 minute spacing.
  t0 = moby2.util.ctime.from_string('2014-10-01 0:0:0')
  t1 = moby2.util.ctime.from_string('2014-11-01 0:0:0')
  ctimes = np.arange(t0, t1, 10*60)
  
  # Look-up the elevation and observation time offset at each time.  Since 
  # we pass in a vector of ctimes, the returned values el and dt are also
  # vectors.
  el, dt = chan_el.get_nearest(ctimes)
  
  # Keep only the data that are less than 20 minutes out of date.
  mask = (abs(dt) < 20*60)
  el, dt = el[mask], dt[mask]



Configuration
-------------

For sensible default behavior, your .moby2 file should contain a block
with the path to the local APEX data archive:

**On galahad/raid:**

.. code-block:: shell

  APEX_weather = {
      'targets_directory': '/dataz/aux_data/hk_lowres/',
      'ctime_range': [1370000000, 2000000000],
  }

**On feynman:**

.. code-block:: shell

  HK_lowres = {
      'targets_directory': '/mnt/act3/users/mhasse/depots/mhasse0/hk_lowres/',
      'ctime_range': [1370000000, 2000000000],
  }

**On tigercpu/della:**

.. code-block:: shell

  HK_lowres = {
      'targets_directory':
      '/projects/ACT/mhasse/depots/actpol_shared/aux_data/hk_lowres',
      'ctime_range': [1370000000, 2000000000],
  }


Auto-doc
--------

.. autoclass:: moby2.aux_data.hk.HKChannel
   :members:

