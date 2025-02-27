from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
FITS-based replacement for columnar ascii data.  We want to use FITS
because it provides a standardized encapsulation of field names and
data types.

The advantage of ASCII is portability, which includes portability to
the eyes of a human that just wants to browse or assess the data at a
glance.

So the goal here is to provide a base class that can write columnar
data to FITS files in a way that includes enough metadata for accurate
ASCII dumping in the absence of custom handling by the child class.
"""

import os, sys
try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits       # requires pyfits or astropy

import numpy as np

# List from smallest to largest so that when we loop through them and
# ask "is this good enough?", we will match the smallest type.

numpy_to_fits = [
    # Start with unsigned integers.  The original data type.
    (np.uint8,  'B', None),
    (np.uint16, 'J', np.dtype('int32')),
    (np.uint32, 'K', np.dtype('int64')),
    # Now signed integers.
    (np.int8,   'I', np.int16),
    (np.int16,  'I', None),
    (np.int32,  'J', None),
    (np.int64,  'K', None),
    # And floats
    (np.float32, 'E', None),
    (np.float64, 'D', None),
    # Hack for bools; store as byte.
    (np.bool_,    'B', None),
    (bool,    'B', None),
]

_numpy_to_fits = {
    'float64': 'D',
    'float32': 'E',
    'int64': 'K',
    'int32': 'J',
    'int16': 'I',
#    'int8': 'A', #int8 not fully supported... 'A' loads as a
                  # character string and 'S' (TSBYTE) does not work
    'uint8': 'B',
    'bool': 'L',
    }

numpy_to_format = {
    'float64': '%15.8e',
    'float32': '%11.4e',
    'int64': '%11i',
    'int32': '%11i',
    'int16': '%6i',
    'int8': '%3i',
    'uint8': '%3i',
    'bool': '%1i',
    }

def get_bintable_hdu(columns):
    return fits.BinTableHDU.from_columns(columns)

class MobyFitsTable:
    """
    MobyFitsTable provides routines for reading and writing simple
    columnar table data to FITS files.  Classes that encapsulate
    simple columnar data can subclass MobyFitsTable to acquire the
    read_fits, write_fits, write_ascii routines.

    The subclass should maintain a list of field names, in
    self.fields.  To customize formatting, self.formats should be a
    dict with printf style format codes for each field (or just the
    tricky fields).

    In addition to the FITS standard table entries, the class reads
    and writes a number of additional keys (such as printf-style
    format codes) that increase portability.
    """

    _mobyfits_container = 'attr'

    def _mobyfits_data(self, field):
        if self._mobyfits_container == 'attr':
            return getattr(self, field)
        elif self._mobyfits_container == 'dict':
            return self[field]
        raise ValueError

    def _mobyfits_fitscode(self, field):
        if field in getattr(self, 'fits_formats', {}):
            return self.fits_formats[field]
        return get_fitscode(data=self._mobyfits_data(field))

    def _mobyfits_fmtcode(self, field):
        fmt = getattr(self, 'formats', {}).get(field)
        if fmt is not None:
            return fmt
        return get_fmtcode(data=self._mobyfits_data(field))

    def write_fits(self, filename=None, fields=None, formats=None, fits_formats=None,
                   clobber=True):
        """
        Write the data to a FITS file, including printf formatting
        information for easy ascii dumping.

        'fields' is a list of fields to write; this will default to
        self.fields.  'fits_formats' and 'formats' also default to
        self.<that>, but if these formats are missing they will be
        guessed based on the data.

        This function returns the resulting astropy.io.fits BinTableHDU, so
        pass filename=None if the HDU is what you really want.
        """
        if fields is None:
            fields = self.fields
        if formats is None:
            formats = {}
        else:
            formats = formats.copy()
        if fits_formats is None:
            fits_formats = {}
        else:
            fits_formats = fits_formats.copy()
        # Get stored or default fits_formats and ASCII formats
        for name in fields:
            if not name in fits_formats:
                fits_formats[name] = self._mobyfits_fitscode(name)
            if not name in formats:
                formats[name] = self._mobyfits_fmtcode(name)
        # Create fits columns
        fits_columns = []
        for name in fields:
            fits_columns.append(fits.Column(name=name,
                                            format=fits_formats[name],
                                            array=self._mobyfits_data(name)))
        tbhdu = get_bintable_hdu(fits_columns)
        # Add printf formatting strings
        for i, name in enumerate(fields):
            tbhdu.header.update('MFFMT%i' % i, formats[name])

        if filename is not None:
            if clobber and os.path.exists(filename):
                os.remove(filename)
            tbhdu.writeto(filename)
        return tbhdu

    @classmethod
    def from_fits(cls, filename, quiet=False):
        """
        Instantiate an object of this class, and load the associated
        columnar data from a FITS file.
        """
        self = cls()
        hdulist = fits.open(filename)
        hdu = hdulist[1]
        fields = []
        formats = {}
        # Do some lower level matching of column names to get their
        # FITS id numbers.
        n_col = hdu.header['TFIELDS']
        fits_idx_map = {}
        for i in range(1,n_col+1):
            fits_idx_map[hdu.header['TTYPE%i' % i]] = i
        # Now check for Moby format codes.  Yes, we are too lazy to
        # translate the FITS standard formats.
        for cidx, column in enumerate(hdu.columns):
            fields.append(column.name)
            data = hdu.data.field(column.name)
            if self._mobyfits_container == 'attr':
                setattr(self, column.name, data)
            if self._mobyfits_container == 'dict':
                self[column.name] = data
            key = 'MDISP%i' % fits_idx_map[column.name]
            if key in hdu.header:
                formats[column.name] = hdu.header[key]
            else:
                if not quiet:
                    print("Key %s not found in header, setting default "\
                        "format for %s." % (key, column.name))
                formats[column.name] = get_fmtcode(data)
        self.fields = fields
        self.formats = formats
        return self
    
    def write_ascii(self, filename=None, fields=None, formats=None,
                    header=True, clobber=True):
        """
        Write a MobyFitsTable in columnar ASCII.  filename should be
        a string or an open file handle.
        """
        if filename is None:
            filename = sys.stdout
        if fields is None:
            fields = self.fields
        if formats is None:
            formats = {}
        else:
            formats = formats.copy()
        n_recs = len(self._mobyfits_data(fields[0]))
        
        formats = [formats.get(f, self._mobyfits_fmtcode(f)) for f in fields]
        data_refs = [self._mobyfits_data(f) for f in fields]
        # Init the file
        if isinstance(filename, basestring):
            if os.path.exists(filename):
                if clobber:
                    os.path.remove(filename)
                else:
                    raise RuntimeError("output file %s exists; remove or "\
                        "pass clobber=True")
            fout = open(filename, 'w')
        else:
            fout = filename
        # Write data
        for row in range(n_recs):
            tokens = [fmt % d[row] for fmt,d in zip(formats, data_refs)]
            if header:
                # Pad, when possible, to keep header aligned with data
                header, didx = '#', 0
                for f, t in zip(fields, tokens):
                    didx += len(t)
                    header += ' ' * max(1,didx-len(header)-len(f)) + f
                    didx += 1
                fout.write(header + '\n')
                # Mark header as written
                header = False
            fout.write(' '.join(tokens) + '\n')
        del fout


def get_fitscode(data=None, dtype=None):
    """
    Get the FITS format code appropriate for data.dtype (or dtype if provided directly).
    """
    if dtype is None:
        dtype = data.dtype
    fits_type = None
    #for k,v in numpy_to_fits.items():
    for (template_type, code, container_dtype) in numpy_to_fits:
        if dtype.type is template_type:
            fits_type = code
            break
    if fits_type is None:
        # String?
        if dtype.kind == 'S':
            fits_type = '%iA' % dtype.itemsize
    if fits_type is None:
        raise ValueError("Could not determine FITS type code for "\
            "numpy.dtype=%s" % str(dtype))
    return fits_type

def get_fmtcode(data):
    """
    Get a printf-style format string for the numpy array 'data'.
    Float types are always returned in exponential notation with 8
    digit precision; formats for integer and string types are widened
    or narrowed to fit the data.
    """
    if data.dtype.kind in ['i', 'u']:
        # Integer-like
        dig_count = int(np.any(data<0))
        dig_count += int(np.ceil(np.log10(abs(data[np.isfinite(data)].max()+.1))))
        return '%' + str(dig_count) + 'i'
    elif data.dtype.kind in ['b']:
        return '%1i'
    elif data.dtype.kind in ['S', 'U']:
        # One might be tempted to use itemsize... for S that works but
        # for U it's 4 times too big because apparently utf-32 is used
        # internally (?).
        # n = data.dtype.itemsize
        n = max(map(len, data))
        return '%' + str(n) + 's'
    elif data.dtype.kind == 'f':
        # Floats
        return '%15.8e'
    raise ValueError("Could not create default printf format for "\
        "numpy.dtype=%s" % str(data.dtype))

def _retyped_for_write(data):
    """
    If we're in python3, convert the ndarray data to the
    python2-compatible data type for writing to FITS.  This probably
    just means that arrays of unicode strings are converted to arrays
    of the "bytes" type; numpy "S" rather than "U".

    Returns the unmodified array, or else a new array or view that has
    conversion applied.
    """
    if sys.version_info.major == 2:
        return data
    if data.dtype.kind == 'U':
        return data.astype('S')
    return data

def _retyped_for_read(data):
    """
    If we're in python3, convert the ndarray data from the
    python2-compatible data type encoded in FITS.  This probably just
    means that arrays of "bytes" type are converted to unicode
    strings; numpy "U" rather than "S".

    Returns the unmodified array, or else a new array or view the
    required conversion applied.
    """
    if sys.version_info.major == 2:
        return data
    if data.dtype.kind == 'S':
        return data.astype('U')
    return data

class MobyFitsTableWriter:
    def __init__(self, data=None, formats={}, key_order=None,
                 header=None):
        self.columns = []
        self.formats = formats.copy()
        self.header = header
        if self.header is None:
            self.header = []
        if data is None:
            return

        # First form a list of (name, vector) pairs.
        if isinstance(data, np.ndarray):
            vectors = [(name, data[name]) for name in data.dtype.names]
        elif isinstance(data, list) or isinstance(data, tuple):
            vectors = [(name, np.asarray(d)) for name,d in data]
        elif isinstance(data, dict):
            if key_order is None:
                key_order = list(data.keys())
            vectors = [(name, np.asarray(data[name])) for name in key_order]
        else:
            raise ValueError("Unknown data format passed to MobyFitsTableWriter.")

        # Convert vectors to fits Column objects; be sure to apply the
        # python2/3 compatibility transformation.
        for name, d in vectors:
            d = _retyped_for_write(d)
            self.columns.append(fits.Column(
                name=name,
                format=get_fitscode(d),
                array=d))

    def get_hdu(self):
        tbhdu = get_bintable_hdu(self.columns)
        # Add printf format codes.
        for i in range(1,tbhdu.header['TFIELDS']+1):
            key = tbhdu.header['TTYPE%i' % i]
            if key in self.formats:
                tbhdu.header['MDISP%i' % i] = self.formats[key]
        for h in self.header:
            if len(h) == 2:
                k, v = h
                tbhdu.header[k] = v
            elif len(h) == 3:
                k, v, c = h
                tbhdu.header[k] = (v, c)
            else:
                raise ValueError("Header data should be tuples of "\
                    "(card,value) or (card,value,comment).")
        return tbhdu
    def write(self, filename, clobber=False):
        tbhdu = self.get_hdu()
        if filename is not None:
            if clobber and os.path.exists(filename):
                os.remove(filename)
            tbhdu.writeto(filename)
        return tbhdu

class MobyFitsTableReader:
    def __init__(self, filename, index=1):
        hdulist = fits.open(filename)
        self.hdu = hdulist[index]
    def get_header(self):
        return self.hdu.header
    def get_formats(self, keys=None):
        if keys is None:
            keys = self.hdu.data.dtype.names
        formats = {}
        for i in range(1,self.hdu.header['TFIELDS']+1):
            key = self.hdu.header['TTYPE%i' % i]
            card = 'MDISP%i' % i
            if key in keys and card in self.hdu.header:
                formats[key] = self.hdu.header[card]
        return formats
    def get_array(self, keys=None):
        return self.hdu.data
    def get_data(self, keys=None):
        if keys is None:
            keys = self.hdu.data.dtype.names
        output = []
        for k in keys:
            output.append((k, _retyped_for_read(self.hdu.data[k])))
        return output
