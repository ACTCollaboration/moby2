from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import sys

def main(args):
    cat = moby2.scripting.get_obs_catalog()
    cat.to_column_file(sys.stdout)
