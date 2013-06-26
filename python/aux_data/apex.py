from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

import numpy as np
import os, glob

from . import tools

class WeatherChannel(tools.AuxChannel):
    """
    Class for loading and reduction of APEX weather data.
    """
    def __init__(self, channel_name, ctime_range=None, source_dir=None):
        """
        channel_name should be one of:

            radiometer
            pressure
            temperature
            dewpoint
            humidity
            windspeed
            winddirection
        """
        cfg = moby2.user_cfg.get('APEX_weather', {})
        if source_dir is None:
            source_dir = cfg.get('targets_directory')
        if ctime_range is None:
            ctime_range = cfg.get('ctime_range')
        if ctime_range is None:
            ctime_range = (0, 1e10)

        self.sources = []
        self.ctime_range = ctime_range
        epoch_start, epoch_stop = ctime_range

        self.data = np.zeros((2,0))
        t_low = epoch_start   # this will prevent us from accidental overlap
        for t0, t1, filename in get_weather_files(source_dir, channel_name):
            if t1 < t_low or t0 > epoch_stop:
                self.sources.append((filename, 0))
                continue
            fin = open(filename, 'r')
            if t0 < t_low:
                tools.seek_ctime(fin, epoch_start)

            if min(epoch_stop, t1) - max(t_low, t0) > 86400*30:
                # This call really is only faster for large numbers of
                # rows...  there is major overhead to the call for some
                # reason, and since it needs to read to the end of the
                # file, it does way more IO than necessary in all but
                # the most extreme cases.
                data = np.loadtxt(fin, unpack=1)
                self.data = np.hstack((self.data, np.array(data)))
            else:
                data = [[],[]]
                for line in fin:
                    entry = line.split()
                    date = float(entry[0])
                    if ((date >= epoch_start) and (date < epoch_stop)):
                        data[0].append(date)
                        data[1].append(float(entry[1]))
                    if (date > epoch_stop):
                        break
                self.data = np.hstack((self.data, np.array(data)))
            self.sources.append((filename, len(data[0])))
            t_low = t1

if 0:
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
        vals, dtimes = np.transpose( list(map(self.get_nearest, ctime)) )
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


class Radiometer(WeatherChannel):
    """
    Load APEX 'radiometer' (PWV at zenith in mm) data.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'radiometer', ctime_range=ctime_range,
                                source_dir=source_dir)

class Pressure(WeatherChannel):
    """
    Load APEX 'pressure', in mBar.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'pressure', ctime_range=ctime_range,
                                source_dir=source_dir)

class Temperature(WeatherChannel):
    """
    Load APEX 'temperature', in C.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'temperature', ctime_range=ctime_range,
                                source_dir=source_dir)

class Dewpoint(WeatherChannel):
    """
    Load APEX 'temperature', in C.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'dewpoint', ctime_range=ctime_range,
                                source_dir=source_dir)

class Humidity(WeatherChannel):
    """
    Load APEX 'dewpoint', in C.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'humidity', ctime_range=ctime_range,
                                source_dir=source_dir)

class WindSpeed(WeatherChannel):
    """
    Load APEX 'windspeed', in m/s.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'windspeed', ctime_range=ctime_range,
                                source_dir=source_dir)

class WindDirection(WeatherChannel):
    """
    Load APEX 'winddirection', in degrees east of North.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'winddirection', ctime_range=ctime_range,
                                source_dir=source_dir)



"""
Support routines for WeatherChannel.
"""

def get_weather_files(target_dir, channel):
    """
    Find weather data files in target_dir.  Filenames must have format

         APEX_<channel>_<start>_<end>.dat

    Where start and end are ctimes.  Returns a list of tuples
         (filename, start, end)
    """
    if not os.path.exists(target_dir):
        raise RuntimeError('get_weather_files: directory "%s" does not exist.' % target_dir)
            
    files = glob.glob('%s/APEX_%s_*.dat' % (target_dir, channel))
    data = []
    for f in files:
        try:
            chan, t0, t1 = f.split('/')[-1].rstrip('.dat').split('_')[1:4]
            data.append((int(t0), int(t1), f))
        except:
            raise RuntimeError('get_weather_files: error parsing filename "%s"' % f)
    return sorted(data)


