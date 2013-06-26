from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
#
# Parse the config file when the module loads, and take some decisive
# actions before too many other modules get loaded.
#
from . import util
user_cfg = util.get_user_config()
#
# 1) Alter sys.path to force perferred versions of some modules?
#
path_inserts = user_cfg.get_deep(('environment', 'python_path_inserts'), [])
if len(path_inserts) > 0:
    import sys
    for x in path_inserts[::-1]:
        sys.path.insert(0, x)
#
# 2) Set the matplotlib backend first, before anyone has a chance to
#    mess it up.
#
if user_cfg.get_deep(('plotting_opts', 'force_non_interactive')):
    from .util import noninteractive_plots

# Import C routines.
from . import libactpol

# Import the rest of the alphabet.
from . import aux_data
from . import detectors
from . import ephem
from . import instruments
from . import pointing
from . import tod
from . import mapping
from . import util
from . import scripting
from .scripting import products

# Objects that we want at root-level:
from .util import DirfileManager
from .detectors import Calibration
from .tod import TOD
from .tod import TODCuts
from .tod import TODInfo

package_init_filename = __file__

