from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

import numpy as np
import os

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
        self._load_ascii_weather(os.path.join(source_dir, 'APEX_%s' % channel_name))


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
