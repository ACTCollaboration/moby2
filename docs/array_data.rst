.. -*- mode: rst ; mode: auto-fill -*-

ArrayData objects
=================

Overview
--------

For basic TOD loading and manipulation, moby2 requires basic information about
an experiment's detector arrays.  These are provided through a simple FITS
table, stored in the shared repository, for each instrument, season, and array.

The ArrayData database is intended to store idealized information rather than
final analysis products.  For example, ArrayData contains detector focal plane
positions, but these are appropriate for plotting and not for map-making.

When a TOD is loaded, the ArrayData object is stored in member `TOD.info.array_data`:

.. code-block:: python

  >>> tod.info.array_data['det_uid']
  array([   0,    1,    2, ..., 1053, 1054, 1055])

But if you wanted, you could instead load ArrayData for a particular
season/array like this:

.. code-block:: python

  >>> array_data = moby2.scripting.get_array_data(
  ... {'season': 2014, 'array_name': 'ar2'})

The ArrayData object currently inherits from dict and so you can grab whatever
vector you want by name:

.. code-block:: python

  >>> print array_data['det_uid']
  [   0    1    2 ..., 1053 1054 1055]
  >>> print array_data['sky_x']
  [ 0.42267853  0.44379987  0.39369126 ...,  0.          0.          0.        ]

To grab multiple vectors and extract only a preferred set of `det_uid`, use
the `get_property()` call:

.. code-block:: python

  >>> array_data.get_property(['sky_x','sky_y'], det_uid=[0,1])
  [array([ 0.42267853,  0.44379987]), array([-1.30226105, -1.17007225])]

ArrayData also supports some database-like functionality.  This allows one to
select detectors matching some criteria.  For example, the det_uid of the
row=3 detectors:

.. code-block:: python

  >>> array_data.select_outer({'row': [3]})
  array([ 96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
         109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
         122, 123, 124, 125, 126, 127])


In order to find the FITS source files, the ``.moby2`` file must contain a
pointer to the ``actpol_shared`` depot.


Data
----

Here we list some of the data contained in the ArrayData files.

**det_uid**- The detector UID.

**col**, **row** - The multiplexing column and row.

**wafer** - The wafer name.

**det_type** - The detector type (e.g., 'tes', 'dark_tes', 'dark_squid', 'johnson', 'no_conn').

**nom_freq** - The nominal detector frequency band.  (This a low-precision
  number that serves as a label for separating the bands in a multi-chroic
  array.  I.e., we don't want more than 1 or 2 sigfigs here.  We currently
  use 150. and 90. to describe the ACTPol bands.)

**sky_x**, **sky_y**, **sky_angle** - Position of the detector in the focal
  plane (degrees), and polarization angle (degrees CCW from east as projected
  onto the sky).  **These are idealized.  Do not use for mapping.** Co-located
  detector pairs should have identical sky_x and sky_y values.  Positions are
  only valid for det_type='tes'.

**optical_sign** - Either +1 or -1 depending on how the optical signal is
  coupled into the readout.

**pol_family** - An arbitrary division of detectors into categories 'A' and
  'B' (or 'X' for undefined).  The two detectors in a co-located pair will
  have different values for this parameter.

**bias_name**, **bias_line** - The MCE bias line assigned to the detector.

