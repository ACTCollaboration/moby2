from __future__ import print_function

import os
BASE = os.getenv('MOBY2_PREFIX')
if BASE == None or BASE == '':
    print()
    print('Set MOBY2_PREFIX in environment or in moby2_config.py')
    print()
    raise RuntimeError

include_dirs = [BASE + '/include']
library_dirs = [BASE + '/lib', '/usr/lib64/atlas/']
libraries = ['actpol']

def get_bash_lines():
    print('export PATH=%s/bin:${PATH}' % BASE)
    print('export PYTHONPATH=$PYTHONPATH:%s/lib64/python2.7/site-packages:'\
        '%s/lib/python2.7/site-packages' % (BASE, BASE))
    print('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s/lib64:${BASE}/lib' % \
        BASE)


