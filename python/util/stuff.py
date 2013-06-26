from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
For quick construction of structs.
"""

class aStruct:
    def __init__(self, items=None):
        self._keys = []
        if items is not None:
            self._update(items)
    def _update(self, items=None):
        if isinstance(items, dict):
            items = list(items.items())
        for k, v in items:
            setattr(self, k, v)
            if not k in self._keys:
                self._keys.append(k)
    def __repr__(self):
        return '<%s@%#x: %s>' % (
            self.__class__, id(self), ', '.join(['%s=%s' %(k,repr(getattr(self, k)))
                                                 for k in self._keys]))
        
def encode_c(x, insist=False):
    if x is None:
        return None
    if hasattr(x, '_encode_c'):
        return x._encode_c()
    if insist:
        raise RuntimeError("object %s does not support _encode_c" % repr(x))
    return x



if __name__ == '__main__':
    qs = aStruct({'n_fish': 5, 'name': 'monger'})
    qs._update([('name', 'jim')])
    qs.n_fish += 1
    print(qs.name, qs.n_fish)

