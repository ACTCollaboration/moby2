.. -*- mode: rst; mode: auto-fill -*-

==========================================
``pathologies`` - TOD pathologies analysis
==========================================

The ``pathologies`` code analyzes TODs and produces statistics
indicative of TOD quality.  It classifies detectors as live, dark, or
dead, and produces TODCuts objects identifying bad data segments.

As of this writing, ``pathologies`` is not an analysis module.  But it
eventually will be.  Use the binary ``process_cuts`` to analyze one or
more TODs.

Getting cuts for a TOD
======================

This program has not been completely updated for moby2.  But it
works.  Run it like this:

.. code-block:: bash

  process_cuts <tod_filename> <params_file>

You will need two parameter files:

1. One like :download:`process_cuts.in
</_moby2/params/process_cuts.in>`, which gives instructions to
``process_cuts``.  The path to this file must be provided to the
``process_cuts`` script on the command line.

2. One like :download:`cut_params.in </_moby2/params/cut_params.in>`,
which sets some thresholds for various cuts.  The path to this file
must be provided in the ``process_cuts.in`` parameter called
"cutParams".

You will probably want to make sure that the depot and output
locations (for cuts_file, e.g.) exist prior to running the program.
