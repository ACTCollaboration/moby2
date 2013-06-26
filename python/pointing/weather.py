from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
class ACTpolWeather:
    temperature = 0.    # C
    pressure = 550.     # mbar
    humidity = 0.2      #
    lapse_rate = 0.0065 # K/m
    def __init__(self, **kwargs):
        for k,v in list(kwargs.items()):
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise ValueError("initializer %s not known" % k)

    def _encode_c(self):
        return (self.temperature, self.pressure, self.humidity, self.lapse_rate)

    @classmethod
    def no_atmosphere(cls):
        return cls(pressure=0.)
    
