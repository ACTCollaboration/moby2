from __future__ import print_function

from distutils.core import setup, Extension

VERSION = '0.1'

import os
import glob
import sys
try:
    import actpol_config
except ImportError:
    print("Failed to load actpol_config.")
    print("The search path is:", sys.path)
    sys.exit(10)

import numpy
includes = actpol_config.include_dirs + [numpy.get_include()]
libdirs = actpol_config.library_dirs + []
library = ['fftw3f', 'gslcblas', 'gsl', 'lapack'] + actpol_config.libraries

sources = glob.glob('src/*.c')
headers = glob.glob('src/*.h')

compile_args = ['-std=c99','-fopenmp', '-Wno-strict-aliasing']
extra_link_args = ['-fopenmp']

## Everything in bin is a script
scripts = [x for x in glob.glob('bin/*') if x[-1] != '~']


# Set C extension to live at in libactpol package.

module1 = Extension('moby2/libactpol',
                    sources=sources,
                    depends=headers,
                    extra_compile_args=compile_args,
                    extra_link_args=extra_link_args,
                    include_dirs=includes,
                    library_dirs=libdirs,
                    libraries=library)

setup (name = 'moby2',
       version = VERSION,
       description = 'ACTpol support',
       package_dir = {'moby2': 'python'},
       ext_modules = [module1],
       scripts = scripts,
       packages = ['moby2',
                   'moby2.analysis',
                   'moby2.analysis.det_sens',
                   'moby2.analysis.fp_fit',
                   'moby2.analysis.beam_ana',
                   'moby2.analysis.jumps',
                   'moby2.analysis.tod_ana',
                   'moby2.analysis.hwp',
                   'moby2.aux_data',
                   'moby2.detectors',
                   'moby2.ephem',
                   'moby2.instruments',
                   'moby2.instruments.actpol',
                   'moby2.instruments.mbac',
                   'moby2.mapping',
                   'moby2.pointing',
                   'moby2.scripting',
                   'moby2.tod',
                   'moby2.util',
                   ])
