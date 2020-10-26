.. -*- mode: rst ; mode: auto-fill -*-

=======================================================
``refit_planet_cal`` - refit absolute calibration model
=======================================================

This script lives in ``analysis.det_sens``.  Its purpose is to
re-analyze planet peak height results and produce a full absolute
calibration (AbsCal) model for multiple seasons, arrays, and bands.

This code consumes the output from one or more fit_planet_cal runs,
allowing the user to combine results and apply them to different time
ranges.

Running the code
================

After setting up the configuration file (see below), run:

.. code-block:: bash
  
  moby2 refit_planet_cal fit_planet_cal.in

This will call the driver program in
:meth:`~moby2.analysis.det_sens.refit_planet_cal.main`.

Running that will produce plots and a calibration model file, of the
form ``models_<tag>_cand.fits``.  When the model is acceptable, rename
to ``models_<tag>.fits`` and run::

  moby2 refit_planet_cal fit_planet_cal.in --write-abscal

This will cause the AbsCal files (hdf5 and txt) to be written out.

Configuration options
=====================

Configuration options are described in the example config file, 
