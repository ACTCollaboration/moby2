.. -*- mode: rst ; mode: auto-fill -*-

APEX weather data
=================

Recipes
-------

The simplest way to get the PWV at a particular time is:

.. code-block:: python

  import moby2
  apex_pwv = moby2.aux_data.apex.WeatherChannel('radiometer')
  ctime = 1378850000  # Sep 10 2013 at 21:53:20 UTC
  print apex_pwv.get_nearest(ctime)
  (2.6075330000000001, -20.0)

In this example, the PWV was about 2.6 mm at a a time 20 seconds
before the requested ctime.

If you've established that there are good weather readings in some
time range that you are interested in, you can get the average APEX
readings over that time range:

.. code-block:: python

  # A five minute interval
  ctimes = 1378850000, 1378850000 + 5*60
  # Get (n_readings, mean PWV, stdev PWV):
  print apex_pwv.get_average(ctimes[0], ctimes[1])
  (5, 2.4839551999999996, 0.11184779471299382)


Configuration
-------------

For sensible default behavior, your .moby2 file should contain a block
with the path to the local APEX data archive:

.. code-block:: shell

  APEX_weather = {
      'targets_directory': '/usr/local/share/actpol/apex_weather/targets/',
      'ctime_range': [1370000000, 1420000000],
  }


Auto-doc
--------

.. autoclass:: moby2.aux_data.apex.WeatherChannel
   :members:

