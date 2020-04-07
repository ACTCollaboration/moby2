.. -*- mode: rst ; mode: auto-fill -*-

=======================================================
``det_sens`` - Tools for detector and array sensitivity
=======================================================

The ``det_sens`` module lives in python/analysis/det_sens.  This
module includes supporting code for:

* TOD power spectrum computation
* Blackbody and RJ spectral radiance calculations
* Planetary brightness models and relevant ephemeris info


Getting the Array Sensitivity from Planet Observations
======================================================

It is assumed that you have fp_fit results for one or more
observations of some planet for which a reasonable brightness model
exists (that probably means Uranus or Neptune).

Configuration file sets are available here:
http://abs-inthe.princeton.edu/~mhasse/distribution/array_sens/

To get the latest version (for which the following information will
hopefully apply), execute this in your working directory:

.. code-block:: bash
  
 wget http://abs-inthe.princeton.edu/~mhasse/distribution/array_sens/array_sens.zip


For each of the TODs of interest, run the following 3 scripts.  It's
assumed that $B has been set to the tod_id:

.. code-block:: bash
  
  moby2 get_fpfit_planetcal get_fpfit_planetcal.cfg $B
  
  moby2 get_tod_spec get_tod_spec.cfg $B 
  
  moby2 get_cal_noise get_cal_noise.cfg $B

If that goes well, you can combine the results into an array sensitivity:

.. code-block:: bash

  moby2 get_array_sens get_cal_noise.cfg

See details below for configuration options, debugging, etc.


Script configuration details
============================

get_fpfit_planetcal
-------------------

get_tod_spec
-------------------

get_cal_noise
-------------------

get_array_sens
-------------------


Other useful sub-modules
========================

blackbody
---------

This module contains functions for converting between Blackbody, RJ,
and CMB spectra.  To see the docstring, run:

.. code-block:: python

  >>> from moby2.analysis import det_sens
  >>> help(det_sens.blackbody)

Here's an example, for getting the conversion from RJ to CMB uK:

.. code-block:: python

  >>> from moby2.analysis import det_sens
  >>> BB = det_sens.blackbody
  >>> f = 150.
  >>> BB.spectrumToDCMB(BB.rayleigh(1., f), f)
  1.7351210431239954

