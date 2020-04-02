from __future__ import print_function

from distutils.core import setup, Extension

VERSION = '0.1'

import glob
import sys

import numpy
includes = [numpy.get_include()]
libdirs = []
library = ['fftw3f', 'gslcblas', 'gsl', 'lapack', 'actpol']

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

analysis_modules =  [
    'det_sens',
    'fp_fit',
    'beam_ana',
    'jumps',
    'tod_ana',
    'hwp',
]
if sys.version_info.major >= 3:
    analysis_modules.append('socompat')

setup (name = 'moby2',
       version = VERSION,
       description = 'ACTpol support',
       package_dir = {'moby2': 'python'},
       ext_modules = [module1],
       scripts = scripts,
       packages = (
           ['moby2',
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
            'moby2.analysis'] +
           ['moby2.analysis.%s' % p for p in analysis_modules]
       )
   )
