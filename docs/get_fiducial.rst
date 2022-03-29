.. -*- mode: rst ; mode: auto-fill -*-

==============================================
``get_fiducial`` - Flat field and fiducial set
==============================================

This script lives in ``analysis.fp_fit``.  Its purpose is to analyze
peak heights from fp_fit operations, calibrated to pW, and determine
the average relative optical couplings (proportional to uK / pW) for
the detectors in an array.

The code also allows the user to identify a subset of the detectors as
"fiducial", meaning that their relative optical calibration seems to
be stable enough to act as an anchor for absolute calibration
purposes.


Running an analysis
===================

First pass
----------

Use a configuration file like the example, here:
::download:`get_fiducial.in <params/get_fiducial.in>`.

Before running the code, update at least these settings:

``fpfit_reduced``
    The pattern for the path to reduced fp_fit results for each
    observation; for example:
    ``'../200423/fp_output1/reduced/{tod_id}_fp.fits'``

``tod_list``
    The path to an ascii file listing the observations you want to
    include (with tod_id in the first column).

We will update the cut parameters after looking at the preliminary
output.  So go ahead and run:

.. code-block:: bash

  mkdir depot
  moby2 get_fiducial get_fiducial.in

On the first pass, the code will write calibration results (copied
from the IV) into depot/.  This calibration will be re-used on
subsequent runs.  The flat field numbers will be relative to this
baseline calibration.


Iteration: set the rms_cut
--------------------------

To iterate, first see if you want to adjust the "rms_cut" setting.
Look at the "10_rms_hist" plots, and adjust rms_cut in the parameter
file to exclude things above the cut line.

Iteration: set the fiducial_selection
-------------------------------------

We mark detectors as fiducial if they lie within some bounded region
in the plane of average flatfield correction ("amp") and scatter (or
fractional rms, "frms").  The plots "{tag}_20_{fcode}.png" show each
detector's position in this plane, and dashed lines that cut out a
rectangle in this space.  The boundaries of the rectangle are set in
the parameter file under "fiducial_selection", which is a dictionary
index by fcode.  A typical entry looks like this::

  {'fcode': 'f039',
   'amp': (0.001, 0.200),
   'frms': (0.001, 0.10),
  }

Outputs
-------

Once you have converged on a reasonable bounding box, look at the
textual output from the program and assess whether the number of
fiducial / flatfielded detectors are acceptable.


Configuration options
=====================

All configuration options are members of a block called ``get_fid``.
In most cases these are well-documented in the example configuration
file provided above.
