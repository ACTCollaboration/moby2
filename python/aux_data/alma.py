from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

import numpy as np
import os, glob

from . import tools

class WeatherChannel(tools.AuxChannel):
    """
    Base class for loading a FITS table of (time,weather) information.
    """
    def __init__(self, channel_name, ctime_range=None, source_dir=None):
        """
        channel_name should be one of:

            pwv
        """
        cfg = moby2.user_cfg.get('ALMA_weather', {})
        if source_dir is None:
            source_dir = cfg.get('targets_directory')
        if source_dir is None:
            source_dir = os.path.join(moby2.user_cfg.get('aux_data', '/'), 'alma_weather')
        if ctime_range is None:
            ctime_range = cfg.get('ctime_range')
        if ctime_range is None:
            ctime_range = (0, 1e10)

        self.ctime_range = ctime_range
        epoch_start, epoch_stop = ctime_range

        columns = [[], []]
        for filename in sorted(glob.glob('%s/*%s*.fits' % (source_dir, channel_name))):
            db = moby2.util.StructDB.from_fits_table(filename)
            t, y = db['ctime'], db[channel_name]
            mask = (epoch_start <= t) * (t < epoch_stop)
            [c.append(_c) for c,_c in zip(columns, [t, y])]
        columns = list(map(np.hstack, columns))
        self.data = np.array(columns)

class Radiometer(WeatherChannel):
    """
    Load ALMA 'radiometer' (PWV at zenith in mm) data.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        WeatherChannel.__init__(self, 'pwv', ctime_range=ctime_range,
                                source_dir=source_dir)
