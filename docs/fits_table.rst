.. -*- mode: rst ; mode: auto-fill -*-

======================================
FITS data table support: MobyFitsTable
======================================

Quick start
===========

Loading FITS table data
-----------------------

To load a FITS table into a ``StructDB`` (which inherits from the
numpy structured ndarray):

.. code-block:: python

  >>> import moby2
  >>> ardata = moby2.util.StructDB.from_fits_table(filename)
  >>> print ard[0]
  (0, 0, 0, 150.0, -0.9482451509415123, -1.287722644556525, 'W10', 'LN_BIAS_08', 12, 'dark_squid', 1, 'tesdatar00c00')
  >>> print ard['det_uid']
  [   0    1    2 ..., 1053 1054 1055]

(This example was generated using the ``ArrayData`` file for ACTPol
PA1, season 2.)  Pass keyword index=[1,2,3...] to .from_fits_table to
select a different HDU.

Writing FITS table data
-----------------------

If your data are already in a ``StructDB``, you can just call
``StructDB.to_fits_table(filename)``.  But they probably aren't.  So:

.. code-block:: python

  >>> import moby2
  >>> data = [('row', [0,1,2,3,0,1,2,3]),
  ...         ('col', [0,0,0,0,1,1,1,1])]
  >>> header_info = [
  ...     ('instrum', 'ACTPol'),
  ...     ('season', '2014'),
  ...     ('array', 'PA1'),
  ...     ]
  >>> formats = {'row': '%2i', 'col': '%2i'}
  >>> ftw = moby2.util.moby_fits.MobyFitsTableWriter(
  ...    data=data, header=header_info, formats=formats)
  >>> ftw.write('my_table.fits')

The MobyFitsTable code will add FITS header entries to store the
printf-style format codes (e.g., '%2i').

Convert FITS table to ASCII
---------------------------

Use the ``moby2_fitsdump`` script to dump a FITS table to the
terminal.  This will respect the MobyFitsTable printf codes, if
present.

.. code-block:: bash

  moby_fitsdump my_table.fits
  # row col
   0  0
   1  0
   2  0
  ...  

