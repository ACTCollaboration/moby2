from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Create connection to MySQL database.  Adapted from ACTManifest, so
that the same manifest.conf -style config file can be used.
"""

try:
    import MySQLdb
    from MySQLdb.constants import CR
    import socket
except:
    MySQLdb = None
    #print('Could not import MySQLdb... some DB functions not available.')

import os, re

def tunnel_active( host, port ):
    sd = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sd.connect((host, port))
    except socket.error:
        return False
    sd.close()
    return True

def create_tunnel_to_mysql_server(host, port, mysqlPort,
                                  user=None, identity=None):
    command = "/usr/bin/ssh"
    if user is not None:
        command += " -l %s" % user
    if identity is not None:
        command += " -i %s" % identity
    command += " -x -L %d:localhost:%d %s -N" % (port,mysqlPort,host)
    if not tunnel_active('localhost', port):
        #print "tunnel inactive -- starting"
        command = "/usr/bin/ssh"
        if user is not None:
            command += " -l %s" % user
        if identity is not None:
            command += " -i %s" % identity
        command += " -x -L %d:localhost:%d %s -N" % (port,mysqlPort,host)
        pipe = os.popen(command,'r')
        pipe.close()

def mysql_bool( s ):
    if s == 'Y':
        return True
    elif s == 'N':
        return False
    return None

def bool2mysql( b ):
    if b: return 'Y'
    return 'N'

def version():
    """Return current cvs version"""
    version = '$Name:  $'                 # version from cvs
    mat = re.search(r"^\$[N]ame:\s*(\S+)\s*\$$", version)
    if mat:
        version = mat.groups()[0]
    else:
        version = "(NOCVS)"
    return version

from .moby_dict import MobyDict

class _Config(dict):
    def __init__(self, filename=None):
        line_no = 0
        if filename is not None:
            line_no += 1
            for line in open(filename):
                line = line.strip()
                if len(line)== 0 or line[0] == '#':
                    continue
                if not '=' in line:
                    raise RuntimeError("Could not parse line %i of %s" % (line_no, filename))
                ei = line.index('=')
                key, value = line[:ei].strip(), line[ei+1:].strip()
                self[key] = value
    def get(self, key, type=None):
        val = dict.get(self, key)
        if type is not None:
            val = type(val)
        return val

    def exists(self, key):
        return key in self

class ManifestConnector(object):

    def __init__( self, user=None, db=None, server=None,
                  configFile=None, tunnel=None, configValues={},
                  tablePrefix=None):
        """Setup the connection to the specified db as \"user\".  Needed parameters
        are read from configFile using actConfig; the dictionary configValues may
        be used to augment or replace this set of parameters"""

        if not configFile and not configValues:
            configFile = '/ACTpol/etc/manifest.conf'  # Just a guess.
                        
        config = _Config(configFile)
        config.update(configValues)

        if db is None:
            if config.exists("db"):
                db = config.get( "db" )
            else:
                db = "act_manifest"

        if not server:
            if config.exists("server"):
                server = config.get("server")
            else:
                server = "localhost"
        if not user:
            user = config.get("user")
        passwd = config.get("passwd:%s" % user)
        mysqlPort = config.get("port", int)

        if not tunnel:
            if config.exists("tunnel"):
                tunnel = config.get("tunnel")
            else:
                tunnel = "no"

        if tablePrefix is None:
            if config.exists('table_prefix'):
                tablePrefix = config.get('table_prefix')
            else:
                tablePrefix = ''

        if tunnel == "yes":
            connect_host = '127.0.0.1'
            port = 9000
            if config.exists('tunnel_user'):
                tunnel_user = config.get('tunnel_user')
            else:
                tunnel_user = None
            if config.exists('tunnel_identity'):
                tunnel_identity = config.get('tunnel_identity')
            else:
                tunnel_identity = None
            create_tunnel_to_mysql_server( server, port, mysqlPort,\
                user=tunnel_user, identity=tunnel_identity )
        else:
            connect_host = server
            port = mysqlPort

        self.db = None
        self.table_prefix = tablePrefix
        self.connect_kwargs = {
            'host'      : connect_host,
            'port'      : port,
            'user'      : user,
            'passwd'    : passwd,
            'db'        : db,
        }
        self.connect()
    
    def connect( self ):
        self.close()

        try:
            self.db = MySQLdb.connect( **self.connect_kwargs )
        except MySQLdb.OperationalError as e:
            # raise unless it's a connection error
            errno = e.args[0]
            if errno != CR.CONNECTION_ERROR and \
               errno != CR.CONN_HOST_ERROR:
                raise
            else:
                print('(mysqldb connection error, %s)' % e.args[1])
    def connected( self ):
        try:
            self.db.ping()
            return True
        except:
            return False

    def close( self ):
        if self.connected():
            self.db.commit()
            self.db.close()
        self.db = None


if __name__ == '__main__':
    # Test
    m = ManifestConnector(configFile='/u/mhasse/manifest_density/2013/manifest.conf')

