.. -*- mode: rst ; mode: auto-fill -*-

===================================================
``fit_planet_cal`` - Fit absolute calibration model
===================================================

This script lives in ``analysis.det_sens``.  Its purpose is to analyze
planet peak height results and produce a calibration model for a
particular data set (e.g. a single season of a single array and band).
This code operates on a list of planet peak heights, along with a
solid angle estimate, and does not process planet maps directly.

Running the fitter
==================

After setting up the configuration file, run:

.. code-block:: bash
  
  moby2 fit_planet_cal fit_planet_cal.in

This will call the driver program in
:meth:`~moby2.analysis.det_sens.fit_planet_cal.main`.

The input file will need modification to specify appropriate data
selection and cuts.


Configuration options
=====================

You will need to provide some configuration options!
