import sys, os, shutil

for f in sys.argv[1:]:
    d, fn = os.path.split(f)
    tf = os.path.join(d, '__patch_' + fn)
    fout = open(tf, 'w')
    preambled = False
    for line in open(f):
        if not preambled:
            if not line.startswith('#!'):
                preambled = True
                fout.write('from __future__ import print_function\n')
                fout.write('from __future__ import absolute_import\n')
                fout.write('from past.builtins import basestring\n')
        fout.write(line)
    fout.close()
    shutil.move(tf, f)

