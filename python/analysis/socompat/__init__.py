# This analysis module support the use of moby2 routines to populate
# SO data structures.  We do not need to support Python2 in this
# sub-module.

# Basic requirements.
import sotodlib

# Local.
from .tod import *
from .db import *
from .metadata import *

import os
if not os.environ.get('MOBY2_NO_SO'):
    # Unless MOBY2_NO_SO is set in the environment, register all our
    # metadata loader functions.  This alters data structures inside
    # sotodlib.
    register_loaders()
