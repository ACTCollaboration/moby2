from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import re
import glob
import numpy as np

#
# Helper functions for extracting patterns from filenames.
#

def extract_basename(text):
    """
    Find a TOD name in text (ctime.ctime.ar1) and return it.
    """
    m = re.search('[0-9]+\.[0-9]+\.ar[123456789]', text)
    return m.group(0)

def extract_basename_freq(text):
    return extract_basename(text) + '_' + extract_freq(text)

def extract_basename_splits(text):
    """
    Take a moby2 map name, with splits extension, and return a
    combined basename of the form ctime.ctime.ar1_split.
    """
    m = re.search('([0-9]+\.[0-9]+\.ar[123456789]).*split([0-9]*).fits', text)
    return '_'.join(m.groups())

def extract_basename_filename(text):
    return text.split('/')[-1]

class BasenameExtractor:
    @classmethod
    def from_params(cls, params):
        self = cls()
        bex = None
        if isinstance(params, basestring):
            self._ex_name = params
            if params == 'standard':
                bex = extract_basename
            elif params == 'standard_freq':
                bex = extract_basename_freq
            elif params == 'splits':
                bex = extract_basename_splits
            elif params == 'filename':
                bex = extract_basename_filename
        if bex is None:
            raise ValueError("Could not determine extractor from parameters: %s" % repr(params))
        self._ex = bex
        return self

    def extract(self, items):
        if isinstance(items, basestring):
            try:
                return self._ex(items)
            except Exception as e:
                print()
                print('** Failed to run extractor "%s" on string "%s"' % (
                    self._ex_name, items))
                print()
                raise e
        return list(map(self.extract, items))

def extract_freq(text):
    """
    Find a frequency band identifier in filename text.  Look for the form f###.
    """
    m = re.search('f[0-9]+', text)
    return m.group(0)

    
def get_map_list(params, args=[], data={}):
    """
    Get a list of maps and their associated names, somehow.  This is
    intended for aid in scripting.  The argument "params" is a dict
    describing where to get the maps.

    Returns a tuple of two lists, (basenames, filenames).

    The most important entry in the dict is 'source', which should
    look like one of:

    'source': ('file', 'map_list.txt')
     - to read filenames from column 0 of map_list.txt

    'source': ('file', 'map_list.txt', 5)
     - to read filenames from column 5 of map_list.txt

    'source': ('file', 'map_list.txt', (0, 1))
     - to read filenames from column 1 and basenames from column 0 of map_list.txt

    'source': ('glob', 'pattern*.fits')

    'source': ('list', ['file1.fits', 'file2.fits'])

    'source': ('list', [('map1', 'path/to/file1.fits'),
                        ('map2', 'path/to/file2.fits')])

    'source': ('command_line', )

    If the source you specify does not provide basenames somehow, the
    basenames will default to be the same as the filenames.  If you
    want to extract some standard basename from the filenames
    automatically, set params['basename_extractor']:

    'basename_extractor': 'standard'
    -- to pull out ctime.ctime.array basename.

    'basename_extractor': 'splits'
    -- to pull out ctime.ctime.array_split00 sub-map name.

    'basename_extractor': 'filename'
    -- to simply remove the leading/path/elements/.

    'basename_extractor': 'standard_freq'
    -- find a basename and frequency code and connect them with _.

    Finally, you can sort the whole thing by basename on its way out.
    Just set

    'sort': True
    """

    # Get a pre-list and generate map names?
    # ...

    # Get the map list.
    msource = params.get('source')
    if msource[0] == 'glob':
        pattern = msource[1]
        map_files = glob.glob(pattern.format(**data))
        basenames = map_files
    elif msource[0] == 'command_line':
        map_files = args
        basenames = map_files
    elif msource[0] == 'file':
        filename = msource[1].format(**data)
        columns = msource[2]
        if not hasattr(columns, '__len__'):
            columns = [columns]
        if len(columns) == 1:
            columns = columns + columns  # duplicate for basename
        data = []
        for line in open(filename):
            w = line.split()
            if len(w) == 0 or w[0][0] == '#': continue
            data.append([w[c] for c in columns])
        basenames, map_files = list(zip(*data))
    elif msource[0] == 'list':
        items = []
        for item in params.get('map_list', []):
            if isinstance(item, basestring):
                item = (item,item)
            items.append((item[0].format(**data), item[1].format(**data)))
        basenames, map_files = list(zip(*data))
    else:
        raise ValueError('Not set up to handle source_maps source=%s' % msource[0])

    # Apply a basename extractor?
    bex = params.get('basename_extractor')
    if bex is not None:
        basenames = BasenameExtractor.from_params(bex).extract(basenames)

    # Sort by something?
    sort_code = params.get('sort')
    if sort_code is True:
        basenames, map_files = list(zip(*sorted(zip(basenames, map_files))))
    
    return basenames, map_files


def get_obs_mask(obs_db, rules):
    """
    Refine obs_db, which is a StructDB (or really any structured numpy
    array).

    rules is a list of simple selection rules to apply to the vectors
    in obs_db.  The results of the rule applications are and-ed
    together to get a mask.  Rules must be of one of the following forms::

        (key, 'in', [item1,item2])
        (key, 'range', (lo, hi))
        (key, 'range_mod', (lo, hi, modulus))

    """
    mask = np.ones(len(obs_db), bool)
    counts = []
    for rule in rules:
        x = obs_db[rule[0]]
        if rule[1] == 'in':
            mask *= [_x in rule[2] for _x in x]
        elif rule[1] == 'range':
            lo, hi = rule[2]
            mask *= (lo <= x) * (x < hi)
        elif rule[1] == 'range_mod':
            lo, hi, modulus = rule[2]
            x = (x-lo) % modulus
            mask *= (0 <= x) * (x < (hi-lo))
        counts.append(mask.sum())
    return mask, counts
