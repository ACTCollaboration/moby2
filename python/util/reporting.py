from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
from .moby_dict import MobyDict

class ProgressDict(MobyDict):

    def __init__(self, filename=None, aggressive=True, dummy=False, read=False):
        super(ProgressDict, self).__init__(self)
        self.filename = filename
        self.aggressive = aggressive
        self.dummy = dummy
        self.key_order = []
        if read:
            self.dummy = True # Prevent clobbering
            self.update_from_file(self.filename)

    def write_to_file(self, filename=None):
        if filename is None:
            filename = self.filename
        if self.dummy:
            return
        return MobyDict.write_to_file(self, filename)

    def __setitem__(self, key, value):
        MobyDict.__setitem__(self, key, value)
        if not key in self.key_order:
            self.key_order.append(key)
        if self.aggressive:
            self.write_to_file(self.filename)

    def exitWithError(self, message, code=1):
        self['exit_code'] = code
        self['exit_message'] = message
        self.write_to_file(self.filename)
        sys.exit(code)

