.. -*- mode: rst ; mode: auto-fill -*-

==========
FocalPlane
==========

FocalPlane objects can be constructed by hand, by passing in the
positions and polarization angles of the detectors.

Factory function
================

The ``moby2.scripting.products`` module provides a function for
getting best FocalPlane in many cases:

.. autofunction:: moby2.scripting.products.get_focal_plane

Typical use cases are:

1. You have a TOD and want the system default template:

.. code-block:: python

  params = {'source': 'template'}
  fplane = moby2.scripting.products.get_focal_plane(params, tod_info=tod.info)
  
2. You want to load an FPFitFile as a FocalPlane:

.. code-block:: python

  params = {'source': 'fp_file', 'filename': 'my_fp_fit.txt'}
  fplane = moby2.scripting.products.get_focal_plane(params)



Documented methods
==================


.. autoclass:: moby2.pointing.FocalPlane
   :members:

