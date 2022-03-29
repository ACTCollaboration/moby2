.. -*- mode: rst ; mode: auto-fill -*-

===================================================
``fit_planet_cal`` - Fit absolute calibration model
===================================================

.. py:module:: moby2.analysis.det_sens.fit_planet_cal

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
:meth:`main`.

The input file will need modification to specify appropriate data
selection and cuts.


Configuration options
=====================

The configuration file for fit_planet_cal is fairly dense.  An example
can be found here: :download:`fit_planet_cal.in
<params/fit_planet_cal.in>`.

You might want to refer to the docs/source for :meth:`load_amplitudes`

Auto-doc
========

.. autofunction:: main

.. autofunction:: load_amplitudes
                  
