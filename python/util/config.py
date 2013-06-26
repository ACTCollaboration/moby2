from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os

from .moby_dict import MobyDict

class MobyConfig(MobyDict):
    filename = None
    def __init__(self, filename=None):
        super(MobyDict, self).__init__()
        if filename is None:
            # Look for .moby2... but not that hard.
            filename = os.getenv('HOME', '')
            filename = os.path.join(filename, '.moby2')
            if not os.path.exists(filename):
                filename = None
        if filename is not None:
            self.filename = filename
            self.update_from_file(filename)

user_config = None

def set_user_config(filename=None):
    global user_config
    user_config = MobyConfig(filename=filename)

def get_user_config():
    global user_config
    if user_config is None:
        set_user_config()
    return user_config
    
