.. -*- mode: rst ; mode: auto-fill -*-

.. py:module:: moby2.analysis.det_sens

=======================================================
``det_sens`` - Tools for detector and array sensitivity
=======================================================

The ``det_sens`` module lives in python/analysis/det_sens.  This
module includes supporting code for:

* TOD power spectrum computation
* Blackbody and RJ spectral radiance calculations (see `blackbody`_).
* Planetary brightness models and relevant ephemeris info


Getting the Array Sensitivity from Planet Observations
======================================================

This section describes using observations of Uranus (or possibly some
other planet, if a reasonable brightness model exists) to measure
array sensitivities.  Going in to this you will need to already have:

- a list of observations, perhaps in a file called list.txt
- fp_fit results for those observations
- glitch cuts and planet position cuts for those observations
- approximate solid angles and band centers for each array and band.

Some template configuration files are found here:
::download:`get_solid_angle.in <params/get_solid_angle.in>`

There are 4 or so steps to follow.


Convert fp_fit to absolute calibration results
----------------------------------------------

The fp_fit results store a peak height, in DAC units.  This can be
used to determine the DAC to uK calibration, using assumptions about
the band center and the beam solid angle.  Before running, update
those values in the configuration file (along with the location of the
fp_fit results).  Then run:

.. code-block:: bash
  
  moby2 get_fpfit_planetcal get_fpfit_planetcal.cfg
  

Get calibrated TOD power spectra
---------------------------------

Now on each TOD of interest, run the code that measures the spectrum.
This is the most computationally expensive step, though it only takes
a few core seconds per TOD.

.. code-block:: bash

  moby2 get_tod_spec get_tod_spec.cfg

Extract calibrated noise levels
-------------------------------

This code will analyze the spectra and extract the white noise level
and store it for the next step.

.. code-block:: bash

  moby2 get_cal_noise get_cal_noise.cfg


If that goes well, you can combine the results into an array sensitivity:

.. code-block:: bash

  moby2 get_array_sens get_cal_noise.cfg


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
and CMB spectra.

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

Some key functions are auto-documented here.  Throughout this module,
the argument ``T`` is a temperature in Kelvin, and ``f`` is a
frequency in GHz.  The units of radiance are non-standard; see
docstrings.  The symbol ``DCMB`` (and also ``dT``) is used to denote
CMB anisotropy temperature units, where a temperature :math:`\delta T` at
frequency :math:`f` corresponds to spectral radiance

.. math::

   S_\nu = \left.\frac{\partial B_\nu(T)}{\partial T}
   \right|_{T=T_{\textrm CMB}} \delta T



.. autofunction:: moby2.analysis.det_sens.blackbody.blackbody_x
.. autofunction:: moby2.analysis.det_sens.blackbody.rayleigh
.. autofunction:: moby2.analysis.det_sens.blackbody.blackbody
.. autofunction:: moby2.analysis.det_sens.blackbody.spectrumToBrightness
.. autofunction:: moby2.analysis.det_sens.blackbody.spectrumToRJ
.. autofunction:: moby2.analysis.det_sens.blackbody.blackbodyToRJ
.. autofunction:: moby2.analysis.det_sens.blackbody.RJToBlackbody
.. autofunction:: moby2.analysis.det_sens.blackbody.DCMB_factor
.. autofunction:: moby2.analysis.det_sens.blackbody.DCMB
.. autofunction:: moby2.analysis.det_sens.blackbody.spectrumToDCMB
.. autofunction:: moby2.analysis.det_sens.blackbody.RJToDCMB
