.. -*- mode: rst ; mode: auto-fill -*-

TODInfo objects
===============

.. py:class:: TODInfo

The point of the TODInfo object is to provide some basic information
about a TOD, so that decisions can be made about how to handle it.
The TODInfo class provides centralized methods for determining
properties like the observing season, array name, etc. from the
dirfile filename and runfile and stuff.

Basic attributes:

- ``instrument``: string identifying telescope; 'actpol'
- ``array_name``: string identifying the array, e.g. 'ar1'
- ``season``: string identifying the epoch of observation; for
  ACT/ACTpol we use 'season2013' and so on.
- ``ctime``: unix timestamp corresponding to the approximate start
  time of the TOD.
- ``name``: short string identifying the TOD.  This will often just be
  the ``basename`` of the file from which the TOD was loaded,
  e.g. 1234560000.12345678043.ar1.  But in cases where the TOD wasn't
  loaded from a file, or whatever, it can be set according to some
  other formula.
- ``array_data``: an :doc:`ArrayData object <array_data>`, which
  contains numerous static features of the array (such as MCE row and
  column assignments, bias line arrangement, hex names).

For TODs loaded from a file, TODInfo may also provide:

- ``filename``: the filename from which the TOD was loaded, if any.
- ``basename``: typically this is the last bit of the filename.
- ``sample_index``: the offset of the present data relative to the
  start of the file from which it was loaded.
- ``downsample_level``: the factor by which the data have been
  resampled.
- ``runfile``: MCERunfile object, with the MCE configuration.
