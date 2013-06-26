from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
import moby2

def package_data_filename(path):
    root = os.path.split(moby2.package_init_filename)[0]
    return os.path.join(root, 'data', path)

