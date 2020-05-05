from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

import os

from . import tools


class ACTRadiometer(tools.AuxChannel):
    """
    Class providing access to ACT radiometer archives.

    You need to have aux_data defined in .moby2.
    """
    def __init__(self, ctime_range=None, source_dir=None):
        aux_data = moby2.user_cfg.get('aux_data')
        if source_dir is None:
            assert aux_data is not None  # Pass in source_dir or set .moby2:aux_data
            source_dir = os.path.join(aux_data, 'act_radiometer')
            assert os.path.exists(source_dir)
        self._load_ascii_weather(os.path.join(source_dir, 'ACT_radiometer'))

    def apex_recal(self, p):
        """
        Apply a recalibration to ACT radiometer measurements (p) that make
        them more consistent with APEX readings.
        """
        a, m1, m2 = [ 1.40004431,  1.15260892,  0.2157575 ]
        dx = p - 1.4
        return (a + m1*dx + m2*abs(dx))

    def get_nearest_apex_recal(self, t):
        # Call super class get_nearest so that user can set
        # self.get_nearest = self.get_nearest_apex_recal if that's
        # what they need.
        p, dt = moby2.aux_data.tools.AuxChannel.get_nearest(self, t)
        return self.apex_recal(p), dt
