from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

import os
import glob
import numpy as np
from moby2.libactpol import get_nearest

class AuxChannel:
    def get_nearest(self, ctime):
        """
        Finds reading closest to specified ctime.  Returns
        (channel_value, dtime), where dtime is the (ctime_of_reading -
        ctime).

        If an array (or list) of ctimes is passed in, arrays of the
        channel_value and dtime are returned.
        """
        if np.asarray(ctime).ndim == 0:
            ii = np.argmin(abs(self.data[0] - ctime))
            return self.data[1][ii], self.data[0][ii] - ctime
        vals, dtimes = np.transpose( [get_nearest(self.data, t) for t in ctime] )
        return vals, dtimes

    def get_average(self, ctime_start, ctime_end):
        """
        Average channel data from times between ctime_start and
        ctime_end.  Returns (n, mean, std), where n is the number of
        readings within the given time range.
        """
        s = (ctime_start <= self.data[0]) * (self.data[0] < ctime_end)
        vals = self.data[1][s]
        if len(vals) == 0:
            return 0, 0., 0.
        return len(vals), vals.mean(), vals.std()

    def get_interp(self, spline=True):
        """
        Return an interpolator.  Care should be taken to avoid regions
        where the data are sparse, and especially regions outside of
        the dataset time range.
        """
        import scipy.interpolate as interp
        if spline:
            return interp.InterpolatedUnivariateSpline(self.data[0], self.data[1])
        return interp.interp1d(self.data[0], self.data[1])

    def _load_ascii_weather(self, prefix, ctime_range=None):
        """Load data from moby standard weather archives."""
        if ctime_range is None:
            ctime_range = (0, 1e10)

        self.sources = []
        self.ctime_range = ctime_range
        epoch_start, epoch_stop = ctime_range
        self.data = np.zeros((2,0))

        t_low = epoch_start   # this will prevent us from accidental overlap
        for t0, t1, filename in get_weather_files(prefix):
            if t1 < t_low or t0 > epoch_stop:
                self.sources.append((filename, 0))
                continue
            d = read_simple_data(filename, epoch_start, epoch_stop)
            self.sources.append((filename, len(d)))
            if len(d):
                self.data = np.hstack((self.data, d))
            t_low = t1


def seek_ctime(ifile, ctime):
    """
    Binary search for ctime in first column of ascii file.  Stops with
    file pointer at the start of the first row having t >= ctime.

    ifile is an open file object, or a filename.
    
    Returns the file.
    """
    if isinstance(ifile, basestring):
        iofile = open(ifile)
    step = 10000       # bytes
    bounds = [0, None] # upper and lower offset boundaries
    # First search upwards
    for counter in range(50):
        if bounds[1] is not None:
            step = (bounds[1] - bounds[0])/2
        else:
            step *= 4
        offset = bounds[0] + step
        ifile.seek(offset)
        new_bounds = [b for b in bounds]
        if offset != 0:  # Eat a \n ?
            offset += len(ifile.readline())
        line = ifile.readline()  # Read a real line
        # Check it
        if len(line) == 0:
            ifile.seek(0, 2)
            new_bounds[1] = ifile.tell()
        else:
            t = float(line.split()[0]) 
            if t >= ctime:
                new_bounds[1] = offset
            else:
                new_bounds[0] = offset
        if new_bounds == bounds:
            break
        bounds = new_bounds
    else:
        print('seek_ctime died in oscillation.  Using last lower bound.')
    ifile.seek(new_bounds[0])
    # Finish with a linear search
    while 1:
        offset = ifile.tell()
        line = ifile.readline()
        if len(line) == 0 or float(line.split()[0]) >= ctime:
            break
    ifile.seek(offset)
    return ifile
    

def read_simple_data(filename, epoch_start, epoch_stop):
    fin = open(filename, 'r')
    seek_ctime(fin, epoch_start)
    rows = []
    for line in fin:
        entry = line.split()
        date = float(entry[0])
        if ((date >= epoch_start) and (date < epoch_stop)):
            rows.append((date, float(entry[1])))
        elif (date >= epoch_stop):
            break
    return np.transpose(rows)


def get_weather_files(prefix):
    """
    Find weather data files in target_dir.  Filenames must have format

         <prefix>_<start>_<end>.dat

    Where start and end are integer ctimes and prefix can have path
    components in it.  Returns a list of tuples
         (filename, start, end)

    """
    files = glob.glob('%s_*.dat' % prefix)
    data = []
    for f in files:
        try:
            t0, t1 = f.split('/')[-1].rstrip('.dat').split('_')[-2:]
            data.append((int(t0), int(t1), f))
        except:
            raise RuntimeError(
                'get_weather_files: error parsing filename "%s"' % f)
    return sorted(data)
