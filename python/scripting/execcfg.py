from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2

"""
Get configuration variables from a user-specified python script
executed in a restricted, but also enhanced, environment.  Perform
string substitution to expand things like depot paths and tod.info
variables.
"""

EXEC_HEADER = """
def namespace(k):
    g = globals()
    g.update(g[k])
def include(_filename):
    g = globals()
    _filename = _filename.format(**g)
    exec(open(_filename).read(), g)

"""

EXEC_FOOTER = """
del include, namespace
"""

def load_vars_exec(filename, data={}, do_substitutions=True):
    text = EXEC_HEADER + open(filename).read() + EXEC_FOOTER
    exec(text, data)
    del data['__builtins__']
    if do_substitutions:
        data = recursive_substitution(data)
    return data

def recursive_substitution(data, base_data=None):
    if base_data is None:
        # Keep calling until it stops changing.
        assert(isinstance(data, dict))
        while True:
            changes, data = recursive_substitution(data, data)
            if not changes:
                break
        return data

    if isinstance(data, basestring):
        new_data = data.format(**base_data)
        return (new_data != data), new_data
    if isinstance(data, list):
        changes, data = list(zip(*[recursive_substitution(v, base_data) for v in data]))
        changes = any(changes)
        return changes, list(data)
    if isinstance(data, tuple):
        changes, data = list(zip(*[recursive_substitution(v, base_data) for v in data]))
        changes = any(changes)
        return changes, tuple(data)
    if isinstance(data, dict):
        changes = False
        for k in data:
            c, v = recursive_substitution(data[k], base_data)
            changes = changes or c
            data[k] = v
        return changes, data
    return False, data

class InputChooser:
    def __init__(self, cfg=None):
        if cfg is None:
            cfg = moby2.user_cfg.get('execcfg', {})
        self.depots = dict(cfg.get('depots', {}))
        self.special_blocks = cfg.get('special_processors', {})
        self.tod_info_block = cfg.get('special_processors', {}).get('tod_info')

    def get_special(self, special, params):
        if not special in self.special_blocks:
            raise RuntimeError(
                'No handler for "%s" in .moby2:execcfg.' % special)
        data = {'depots': self.depots}
        data.update(params)
        return load_vars_exec(self.special_blocks[special].format(**data), data)
        
    def get_config(self, filename, tod_id=None, data=None):
        if data is None:
            _data_in = {}
        else:
            _data_in = data
        data = {}
        data.update(_data_in)
        data['depots'] = self.depots
        if tod_id is not None:
            if isinstance(tod_id, basestring):
                tod_id = load_vars_exec(self.tod_info_block.format(**data),
                                        {'id': tod_id})['tod_info']
            data['tod_info'] = tod_id
        # Do substitutions on the filename, including tags, within
        # namespace('depots').
        fn_kw = data.copy()
        fn_kw.update(data['depots'])
        filename = filename.format(**fn_kw)
        # Return only a single variable?
        k = None
        if ':' in filename:
            filename, k = filename.split(':', 1)
        output = load_vars_exec(filename, data)
        if k is not None:
            output = output[k]
        return output
