from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

import numpy
import math, time, datetime, calendar
import ephem, pytz

from moby2.util import constants

def _fixed( ra, dec ):
    f = ephem.FixedBody()
    f._ra = ra
    f._dec = dec
    return f

class ACTEphem:

    _objects = {
        'Sun'       : ephem.Sun(),
        'Moon'      : ephem.Moon(),
        'Mercury'   : ephem.Mercury(),
        'Venus'     : ephem.Venus(),
        'Mars'      : ephem.Mars(),
        'Jupiter'   : ephem.Jupiter(),
        'Saturn'    : ephem.Saturn(),
        'Uranus'    : ephem.Uranus(),
        'Neptune'   : ephem.Neptune(),

        'Abell 85'      : _fixed( '00:41:50', '-09:18:07' ),
        'Abell 478'     : _fixed( '04:13:26', '+10:27:58' ),
        'Abell 3404'    : _fixed( '06:45:29', '-54:13:08' ),
        '1ES0657'       : _fixed( '06:58:31', '-55:56:49' ),
        'Abell 754'     : _fixed( '09:08:50', '-09:38:12' ),
        'Abell 1689'    : _fixed( '13:11:30', '-01:20:07' ),
        'Abell 3558'    : _fixed( '13:27:58', '-31:29:17' ),
        'Abell 3571'    : _fixed( '13:47:28', '-32:51:14' ),
        'Abell 2163'    : _fixed( '16:15:49', '-06:09:00' ),
        'Abell 3667'    : _fixed( '20:12:31', '-56:49:55' ),
        'WMAP J0540 5416' : _fixed( '05:40:32', '-54:19:15'),
        'Tau_a'         : _fixed( '05:34:32', '+22:00:52' ),
        'Cen_a'         : _fixed( '13:25:27', '-43:01:09' ),
    }

    def listObjects():
        """
        @brief static method to find sources with positions known to actEphem
        @return List: object strings 
        """
        return list(ACTEphem._objects.keys())

    listObjects = staticmethod(listObjects)

    def __init__( self ):
        site = ephem.Observer()
        site.lat =  str(constants.site['latitude_deg'])
        site.long = str(constants.site['longitude_deg'])
        site.elev =     constants.site['elevation_m']

        # Use the pyephem optical refraction model, with a 1% correction for mm-wave
        # refraction.  This gives results correct to better than 1" provided P_water is
        # between 0.5 and 1.5 mbar (typically, that means 6-25% relative humidity at 273 Kelvin).
        # Plan is to use constant typical temp/pressure for now (April 2008).
        #
        # Note that pressure at ACT is 0.990 of that at APEX, given the elevations (5188, 5107
        # meters, respectively) and the 8780 m density scale height of the US Standard
        # Atmosphere between 4 and 10 km altitude.
        site.temp = -2.0 # Celsius, a typical number
        site.pressure = 548.0 # mbar.
        OPTICAL_MM_FUDGE_FACTOR = 1.01
        site.pressure *= OPTICAL_MM_FUDGE_FACTOR
        # Roughly refraction in optical at 548 mb is same as millimeter refraction at 553 mb
        # Expect 548 mb at ACT site as typical value

        site.horizon = '32'
        self.site = site

    def set_timestamp( self, timestamp, fixDUT1=False ):
        if fixDUT1:
            dut1 = dUT1Getter()(timestamp)
        else:
            dut1 = 0.
        d = datetime.datetime.utcfromtimestamp( timestamp + dut1 )
        t = d.year, d.month, d.day, d.hour, d.minute, d.second+d.microsecond*1e-6
        self.site.date = ephem.date( t )

    set_ctime = set_timestamp

    def sidereal_time( self ):
        return self.site.sidereal_time()

    def altaz_to_radec( self, alt, az ):
        return self.site.radec_of( az*math.pi/180., alt*math.pi/180. )

    def radec_to_altaz( self, ra, dec ):
        o = _fixed( ra*math.pi/180., dec*math.pi/180. )
        o.compute( self.site )
        return o.alt, o.az

    def get_all_object_names( self ):
        return list(self._objects.keys())

    def get_object( self, name ):
        p = ACTEphem._objects[name]
        p.compute( self.site )
        return p

    def get_alt_crossing_times( self, name, alt, start, stop ):
        assert start < stop

        self.site.horizon = alt*math.pi/180.

        # add an extra day before and after
        year, month, day = start.triple()
        date = ephem.date( (year,month,int(day)) ) - 1.
        ndays = int(math.ceil(stop-start)) + 2

        rise_time = []
        set_time = []

        for i in range(ndays):
            self.site.date = date + i
            p = self.get_object( name )
            if p.rise_time is not None:
                rise_time.append( p.rise_time )
            if p.set_time is not None:
                set_time.append( p.set_time )

        if len(rise_time) == 0 and len(set_time) == 0:
            return None

        def insert_missing_risets( name, riset ):
            for i in range(len(riset)-1):
                if riset[i+1]-riset[i] > 1.8:
                    avg = ephem.date( 0.5*(riset[i+1]+riset[i]) )
                    print("found missing %s time between %s, %s" % (name, riset[i], riset[i+1]))
                    print(" - inserting ", avg)
                    riset.insert( i+1, avg )

        insert_missing_risets( 'rise', rise_time )
        insert_missing_risets( 'set', set_time )

        assert abs(len(rise_time)-len(set_time)) < 2

        result = []
        for rise in rise_time:
            if start < rise and rise < stop:
                self.site.date = rise
                p = self.get_object( name )
                result.append( (rise,0,p.az) )
        for set in set_time:
            if start < set and set < stop:
                self.site.date = set
                p = self.get_object( name )
                result.append( (set,1,p.az) )

        if len(result) > 0:
            result.sort()
            return result

        return None

    def object_hours_per_night( self, name, year, month, day ):
        date = ephem.date( (year,month,day) )

        rise_time = []
        set_time = []

        self.site.date = date
        for i in range(2):
            p = self.get_object( name )
            if p.neverup:
                return 0.
            if p.circumpolar:
                return 12.
            rise_time.append( p.rise_time )
            set_time.append( p.set_time )
            self.site.date += 1.

        if set_time[0] is None or set_time[1] is None:
            print(rise_time[0], set_time[0])
            print(rise_time[1], set_time[1])

        if rise_time[1]-rise_time[0] > 1.8:
            print(name, 'bad rise time', rise_time[1], float(rise_time[1]), date)
            rise_time[1] -= 1
        if set_time[1]-set_time[0] > 1.8:
            print(name, 'bad set time', set_time[1], float(set_time[1]), date)
            set_time[1] -= 1

        uptimes = []
        if rise_time[0] < set_time[0]:
            uptimes.append( (rise_time[0],set_time[0]) )
            uptimes.append( (rise_time[1],set_time[1]) )
        else:
            uptimes.append( (date,set_time[0]) )
            uptimes.append( (rise_time[0],set_time[1]) )
            uptimes.append( (rise_time[1],date+2) )

        # night is defined as 10pm-10am CST
        night_start = date + (22.+4.)*ephem.hour # 4==CST
        night_stop = night_start + 12.*ephem.hour
        t = 0.
        for rise,set in uptimes:
            if not rise < set:
                #print rise_times.keys()
                #print set_times.keys()
                #print date
                #print uptimes
                #print float(rise), float(set), set-rise
                return -1.
            a = max(night_start,rise)
            b = min(night_stop,set)
            t += max(0.,24.*(b-a))
        if t == 12.:
            print(list(rise_times.keys()))
            print(list(set_times.keys()))
        return t

#def get_planet_from_name( name ):
#    p = {
#        'sun'       : ephem.Sun(),
#        'mars'      : ephem.Mars(),
#        'jupiter'   : ephem.Jupiter(),
#        'saturn'    : ephem.Saturn(),
#        'neptune'   : ephem.Neptune(),
#    }
#    return p.get( name )
#
#def ctime_to_ephemdate( ctime ):
#    s = time.gmtime( ctime )
#    t = s.tm_year, s.tm_mon, s.tm_mday, s.tm_hour, s.tm_min, s.tm_sec
#    return ephem.date( t )
#
#def planet_altaz_at_ctime( name, ctime ):
#    planet = get_planet_from_name( name )
#    obs = ACTObserver()
#    obs.date = ctime_to_ephemdate( ctime )
#    planet.compute( obs )
#    return planet.alt, planet.az
#
#def hits_sun( time, alt, az ):
#
#    obs = ACTObserver()
#    sun = ephem.Sun()
#
#    d2r = math.pi/180.
#    fail = False
#    min_sep = 180.
#    for i in range(len(time)):
#        obs.date = ctime_to_ephemdate( time[i] )
#        sun.compute( obs )
#        # sep is in degrees
#        sep = ephem.separation( (sun.az,sun.alt), (az[i]*d2r,alt[i]*d2r) )
#        #min_sep = min( min_sep, sep )
#        if sep < 10.*d2r:
#            print sun.ra, sun.dec, alt[i], az[i]
#            fail = True
#            break
#    #if min_sep < 10.*math.pi/180.:
#    #    fail = True
#    #print min_sep
#
#    return fail

def ephemdate_to_datetime( d ):
    t = d.tuple()
    fsec,isec = math.modf(t[5])
    return datetime.datetime( t[0], t[1], t[2], t[3], t[4], \
        int(isec), int(fsec*1e6), tzinfo=pytz.utc )

def UTC_to_CLT( d ):
    u = ephemdate_to_datetime( d )
    c = u.astimezone( pytz.timezone('America/Santiago') )
    return c.strftime( "%a %b %d %H:%M:%S %Z %Y" )

def ephemdate_to_ctime( d ):
    dt = ephemdate_to_datetime( d )
    return calendar.timegm( dt.timetuple() )


# Local sunrise and sunset

def sun_times(date):
    """
    Return the sunrise and sunset times at the site on the specified
    date.  Times are just horizon crossing, they do not account for
    mountains being in the way or wahtever.

    If date is a float or int, it is interpreted as a ctime and two
    ctimes are returned.  If date is any other format, it is handed
    straight to pyephem and two ephem.Dates are returned.

    Conveniently, sunset at our site is always within a few hours
    _before_ 0:00 UTC, and sunrise is always well after midnight.  So
    there will be no ambiguity if you always ask about 0:00 <date>.
    """
    
    obs = ephem.Observer()
    obs.lat = '-22.9585'
    obs.long= '-67.7876'
    sun = ephem.Sun()
    ctime = hasattr(date, '__float__')
    if ctime:
        d = datetime.datetime.utcfromtimestamp(date)
        t = d.year, d.month, d.day, d.hour, d.minute, d.second+d.microsecond*1e-6
        date = ephem.date(t)
    obs.date = date
    rise_set = [obs.next_rising(sun), obs.next_setting(sun)]
    if ctime:
        return [calendar.timegm(ephemdate_to_datetime(d).timetuple()) for d in rise_set]
    return rise_set



if __name__ == '__main__':

    #x = ACTEphem()

    #x.set_timestamp( time.time() )
    #print x.date
    #print datetime.datetime.utcnow()

#    for name in x._planets:
#        print name, x.planet_hours_per_night( name, 2007, 12, 1 )

#    elevation = 32.
#    start = ephem.date((2007,3,8,22+4))
#    stop = ephem.date((2007,3,9,10+4))
#    print x.get_alt_crossing_times( 'Saturn', elevation, start, stop )

    x = ACTEphem()
    elevation = 60.
    start = ephem.date((2007,7,13,18))
    stop = ephem.date((2007,7,14,18))

    rule = "-"*72
    print(rule)
    print("Planet crossings at elevation=%f" % elevation, "between")
    print(UTC_to_CLT(start), "and", UTC_to_CLT(stop))
    print(rule)
    riset_str = ["rises","sets"]
    for planetname in ["Mars","Jupiter","Neptune","Uranus"]:
        xings = x.get_alt_crossing_times( planetname, elevation, start, stop )
        print("%s:" % planetname)
        if xings is None:
            print("None")
        else:
            for xing in xings:
                print(ephemdate_to_ctime(xing[0]), UTC_to_CLT(xing[0]), \
                    riset_str[xing[1]], "at azimuth", xing[2])
    print(rule)

def mjd_to_ctime(mjd):
    return (mjd - 40587.)*86400.

def ctime_to_mjd(ctime):
    return float(ctime)/86400. + 40587.

class dUT1Getter:
    def __init__(self, filename=None):
        if filename is None:
            import moby2
            user_config = moby2.util.get_user_config()
            filename = user_config.get('bulletin_A_settings')['filename']
        self.filename = filename
        data = numpy.loadtxt(filename, unpack=1)
        self.mjd0 = data[0][0]
        assert(numpy.all(numpy.diff(data[0]) == 1))
        self.xa, self.yz, self.dUT1 = data[1:4]
    def __call__(self, times):
        idx = numpy.floor(ctime_to_mjd(times) - self.mjd0).astype('int')
        if numpy.any(idx >= len(self.dUT1)):
            raise RuntimeError(
                "User requested dUT1 for MJD beyond end of dUT1 source file (%s)." %
                self.filename)
        return self.dUT1[idx]
        
