.. -*- mode: rst ; mode: auto-fill -*-

Working with time-ordered data (TOD)
====================================

The TOD class is the basic container for detector and pointing data.
Populating a TOD with data can be done a number of different ways, so
loading TODs is normally done through the moby2.scripting.get_tod
function, which is able to decode complicated parameters and handle
the peculiarities of different experiments.

Quick start
-----------

To load all of the detector and pointing data for a dirfile:

.. code-block:: python

  import moby2
  filename = '1234560000.1234569999.ar1'
  tod = moby2.scripting.get_tod({'filename': filename})

Look at the data and pointing info:

.. code-block:: python

  >>> print tod.data.shape
  (1056, 256000)
  >>> print tod.data[0]         # First detector data
  [1.34e3, 1.33e3, ...]
  >>> print tod.az
  [...]
  >>> print tod.alt
  [...]

To load only certain samples, use the ``start`` and ``end`` arguments.
These have the same semantics as array indices [start:end], including
the handling of negative arguments:

.. code-block:: python

  # Load samples 200 to 1200 samples
  tod = moby2.TOD.from_dirfile({'filename': filename,
                                'start': 200, 'end': 1200})

  # Load entire TOD except for the last 400 samples:
  tod = moby2.TOD.from_dirfile({'filename', 'end': -400})

To select particular detectors, use ``det_uid`` or ``rows`` and
``cols`` with ``read_block``:

.. code-block:: python

  # Load detector channel 123
  tod = moby2.TOD.from_dirfile({'filename': filename, 'det_uid': [123]})

  # Load multiplexing row 2
  tod = moby2.TOD.from_dirfile({'filename': filename, 'rows': [2]})

  # Load only (row 2, col 3) and (row 8, col 0) by passing block_read=False
  tod = moby2.TOD.from_dirfile({'filename': filename,
                                'rows': [2,8], 'cols': [3,0],
                                'block_read': False})

Full class documentation
------------------------

.. module:: moby2


.. autoclass:: TOD
   :members:
   :undoc-members:
