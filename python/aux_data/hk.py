from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

import numpy as np
import os, glob

from . import tools

class HKChannel(tools.AuxChannel):
    """
    Class for loading and looking up low-resolution HK archival data.
    """
    def __init__(self, channel_name, ctime_range=None, source_dir=None):
        """
        channel_name is probably one of:
          Enc_Az_Astro_Deg
          Enc_El_Deg
        """
        cfg = moby2.util.config.get_user_config().get('HK_lowres', {})
        if source_dir is None:
            source_dir = cfg.get('targets_directory')
        if source_dir is None:
            source_dir = os.path.join(moby2.user_cfg.get('aux_data', '/'), 'hk_lowres')
        if ctime_range is None:
            ctime_range = cfg.get('ctime_range')
        if ctime_range is None:
            ctime_range = (0, 1e10)

        self.sources = []
        self.ctime_range = ctime_range
        epoch_start, epoch_stop = ctime_range

        self.data = np.zeros((2,0))
        t_low = epoch_start   # this will prevent us from accidental overlap
        for t0, t1, filename in get_hk_files(source_dir, channel_name):
            if t1 < t_low or t0 > epoch_stop:
                self.sources.append((filename, 0))
                continue
            fin = open(filename, 'r')
            if t0 < t_low:
                tools.seek_ctime(fin, epoch_start)

            if min(epoch_stop, t1) - max(t_low, t0) > 86400*30:
                # This call really is only faster for large numbers of
                # rows...  there is major overhead to the call for some
                # reason, and since it needs to read to the end of the
                # file, it does way more IO than necessary in all but
                # the most extreme cases.
                data = np.loadtxt(fin, unpack=1)
                self.data = np.hstack((self.data, np.array(data)))
            else:
                data = [[],[]]
                for line in fin:
                    entry = line.split()
                    date = float(entry[0])
                    if ((date >= epoch_start) and (date < epoch_stop)):
                        data[0].append(date)
                        data[1].append(float(entry[1]))
                    if (date > epoch_stop):
                        break
                self.data = np.hstack((self.data, np.array(data)))
            self.sources.append((filename, len(data[0])))
            t_low = t1


def get_hk_files(target_dir, channel):
    """
    Find HK data files in target_dir.  Filenames must have format

         <field>_<start>_<end>.dat

    Where start and end are ctimes.  Returns a list of tuples
         (filename, start, end)
    """
    if not os.path.exists(target_dir):
        raise RuntimeError('get_hk_files: directory "%s" does not exist.' % target_dir)
            
    files = glob.glob('%s/%s_*.dat' % (target_dir, channel))
    data = []
    for f in files:
        try:
            tokens = os.path.split(f)[1].rstrip('.dat').split('_')
            chan, t0, t1 = '_'.join(tokens[:-2]), tokens[-2], tokens[-1]
            data.append((int(t0), int(t1), f))
        except:
            raise RuntimeError('get_hk_files: error parsing filename "%s"' % f)
    return sorted(data)


def get_hk_fields(ctime_range=None, source_dir=None):
    cfg = moby2.util.config.get_user_config().get('HK_lowres', {})
    if source_dir is None:
        source_dir = cfg.get('targets_directory')
    if ctime_range is None:
        ctime_range = cfg.get('ctime_range')
    if ctime_range is None:
        ctime_range = (0, 1e10)

    files = glob.glob('%s/*.dat' % (source_dir))
    fields = []
    for f in files:
        try:
            tokens = os.path.split(f)[1].rstrip('.dat').split('_')
            field, t0, t1 = '_'.join(tokens[:-2]), tokens[-2], tokens[-1]
            if field in fields:
                continue
            t0, t1 = max(int(t0), ctime_range[0]), min(int(t1), ctime_range[1])
            if (t0 <= t1):
                fields.append(field)
        except:
            raise RuntimeError('get_hk_files: error parsing filename "%s"' % f)
    return fields


