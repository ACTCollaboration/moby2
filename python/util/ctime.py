from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import time
import calendar

#
# Date/Time Utilities
#

#DEFAULT_FORMAT = '%b %d %H:%M:%S %Y UT'
DEFAULT_FORMAT = '%Y-%m-%d %H:%M:%S'

def to_string(t, format=DEFAULT_FORMAT):
    """
    @brief Translate a ctime (int) to a formatted string.
    @param ctime int: ctime to translate
    @return str: formatted time string
    """
    return time.strftime(format, time.gmtime(t))

def from_string(s, format=DEFAULT_FORMAT):
    """
    @brief Translate a formatted time string to a ctime 
    @param asctime str: formatted time string (see mobyManifest.timeFormat)
    @return int: ctime
    """
    return calendar.timegm(time.strptime(s, format))

if __name__ == '__main__':
    t0 = 1234560000
    s = to_string(t0)
    t1 = from_string(s)
    print(t0)
    print(s)
    print(t1)
