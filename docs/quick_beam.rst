.. -*- mode: rst ; mode: auto-fill -*-

==========================================
``quick_beam_fit`` - gaussian model fitter
==========================================

This program performs an elliptical gaussian fit to the signal in
planet map(s).  The output can be analyzed to assess focus quality.

Analyze maps
============

Use a parameters file like :download:`quick_beam.in
</_moby2/params/quick_beam.in>`.  The list of maps to analyze can be
passed in on the command line, or specified in a file.

.. code-block:: bash

  quick_beam_fit quick_beam.in


Output looks like:

.. code-block:: bash
  
  #             peak_x  peak_y   peak_h  fwhm_a  fwhm_b   angle  baseline    noise
  source1.fits  0.0001  0.0001   0.2730  0.0234  0.0213  117.38  8.161e-04 5.492e-04
  source2.fits -0.0013 -0.0027   0.2700  0.0234  0.0215  116.85  8.004e-04 5.425e-04

The position of the peak in the map is given, in degrees, by
``peak_x`` and ``peak_y``.  The major and minor FWHM of the gaussian
model are ``fwhm_a`` and ``fwhm_b``, in degrees.  The orientation of
the major axis, measured in degrees from X towards Y, is ``angle``.
The ``baseline`` and ``noise`` indicate the approximate mean map level
and white noise level, respectively.

The code also creates plots with annotations.  The plotting can be
configured to some extent, to change the zoom factor or set linear /
log color scale.


Focus assessment
----------------

If evaluating beam quality in a planetfest situation, the product of
the major and minor FWHM gives a good proxy for beam solid angle
(especially if you scale it by pi/2).  In the plots, this number is
labeled "Omega_G".  The peak amplitude is another useful indicator of
beam quality.

