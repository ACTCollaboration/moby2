from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
ACTpol-specific database classes.
"""

from moby2.util.database import TODRecordList, convert_ctimes
from moby2.util.filebase import MobyFilebase

import os

class TODRecord:
    fields = None
    def __init__(self):
        self.fields = []

    @classmethod
    def from_database(cls, tdb, mce_id=None):
        """
        Load TOD information from ACTpol manifest over database
        connection tdb, which must be an actpol.TODDatabase instance.
        mce_id should be the mce_data_acq_id, or basename.
        """
        self = cls()
        c = tdb.db.cursor()

        datafile_id = None
        basename = None

        #Convert basename to mce_data_acq_id?
        if isinstance(mce_id, basestring):
            basename = mce_id.strip('/').split('/')[-1]
            # Does this match a datafile?
            q = 'select id from %s where producer="merlin" and '\
                'name="%s" or name="%s.zip"' % \
                (tdb._table('datafiles'), basename, basename)
            if c.execute(q) > 0:
                datafile_id = c.fetchone()[0]
            # To guarantee we can find this in tod_info, get the
            # mce_data_acq_id.
            basename_tokens = basename.split('.')
            if len(basename_tokens) != 3:
                raise ValueError("Could not decode basename in '%s'" % mce_id)
            ctime, ar = int(basename_tokens[0]), basename_tokens[2]
            q = 'select id from %s where '\
                'ctime=%i and array="%s"' % \
                (tdb._table('mce_data_acq'), ctime, ar)
            if c.execute(q) <= 0:
                raise ValueError("mce_data_acq entry not found for "\
                    "ctime=%i, array=%s" % (ctime, ar))
            mce_id = c.fetchone()[0]
        else:
            if c.execute('select ctime,array from %s where '
                         'id=%i' % (tdb._table('mce_data_acq'), mce_id)) == 0:
                raise valueError("mce_data_acq entry not found for "\
                    "id=%i" % mce_id)
            ctime, ar = c.fetchone()

        # Start with tod_info fields
        q = 'select * from %s where mce_data_acq_id=%i' % \
            (tdb._table('tod_info'), mce_id)
        n = c.execute(q)
        if n == 0:
            return None
        for i, v in enumerate(c.fetchall()[-1]):
            k = c.description[i][0]
            self.fields.append(k)
            setattr(self, k, v)
        # And mce_data_acq fields
        for k, v in [('ctime', ctime), ('array', ar)]:
            if not k in self.fields:
                self.fields.append(k)
            setattr(self, k, v)

        # We want to set basename and datafile_id.  At this point, either:
        #
        # - user passed in a basename, and maybe we also found a
        #   datafile_id for it.  Set .basename and .datafile_id to
        #   those values, regarless of tod_info contents.
        #
        # - user passed in an mce_data_acq.id, so we should take
        #   tod_info.datafile_id and load its basename from datafiles.
        if basename is None:
            datafile_id = self.datafile_id
            if datafile_id is not None and c.execute(
                'select name from %s where id=%i' %
                (tdb._table('datafiles'), datafile_id)) > 0:
                basename = c.fetchone()[0]

        if basename is not None and basename.endswith('.zip'):
            basename = basename[:-4]

        self.basename = basename
        self.fields.append('basename')
        self.datafile_id = datafile_id
        return self

    @classmethod
    def from_cursor(cls, cursor):
        self = cls()
        row = cursor.fetchone()
        for i, v in enumerate(row):
            k = cursor.description[i][0]
            self.fields.append(k)
            setattr(self, k, v)
        return self

    def to_string(self):
        """
        Return textual summary of data in this record.
        """
        n = max([len(k) for k in self.fields])
        fmt = '%%%is : %%s' % -n
        text = ''
        for k in self.fields:
            text += fmt % (k, repr(getattr(self, k))) + '\n'
        return text

class AcquisitionRecord(TODRecord):

    @classmethod
    def from_database(cls, tdb, mce_id, ar=None):
        """
        Load ACQ information from ACTpol manifest over database
        connection db, which must be an actpol.TODDatabase instance.
        mce_id should be the mce_data_acq_id.  Alternately mce_id may be the
        MCE filename (<ctime>_<suffix> format), but this requires the
        array name to also be passed in.
        """
        c = tdb.db.cursor()

        #Convert basename to mce_data_acq_id?
        if isinstance(mce_id, basestring):
            basename_tokens = mce_id.strip('/').split('/')[-1].split('_')
            if len(basename_tokens) != 2:
                raise ValueError("Could not decode MCE filename in '%s'" % mce_id)
            ctime, suffix = int(basename_tokens[0]), basename_tokens[1]
            q = 'select *, concat(ctime,"_",suffix) as mce_name '\
                'from %s where '\
                'ctime=%i and suffix=%s and array="%s"' % \
                (tdb._table('mce_data_acq'), ctime, suffix, ar)
        else:
            q = 'select *, concat(ctime,"_",suffix) as mce_name '\
                'from %s where '\
                'id=%i' % (tdb._table('mce_data_acq'), mce_id)
        n = c.execute(q)
        if n == 0:
            return cls()
        return cls.from_cursor(c)


class HKRecord(TODRecord):

    @classmethod
    def from_database(cls, tdb, merlin_id, ar=None):
        """
        Load HK archive information from ACTpol manifest over database
        connection db, which must be an actpol.TODDatabase instance.
        merlin_id should be the merlin.id.  Alternately merlin_id may be the
        basename.
        """
        c = tdb.db.cursor()
        tables = {'merlin': tdb._table('merlin'),
                  'datafiles': tdb._table('datafiles')}

        q = 'select {merlin}.*, {datafiles}.name as basename '\
            'from {merlin} join {datafiles} on '\
            '{merlin}.datafile_id = {datafiles}.id '
        #Convert basename to merlin.id?
        if isinstance(merlin_id, basestring):
            basename = merlin_id.strip('/').split('/')[-1]
            q += 'where {datafiles}.name="{basename}"'
        else:
            q += 'where {merlin}.id={id_}'
        n = c.execute(q.format(basename=basename,id_=merlin_id,**tables))
        if n == 0:
            return cls()
        return cls.from_cursor(c)


def get_actpol_db(config_file=None, season=None):
    import moby2
    if config_file is None:
        db_config = moby2.user_cfg.get('actpol_manifest', {})
        config_file = db_config.get('config_file')
    if config_file is None:
        config_file = '/ACTpol/etc/manifest.conf'
        if not os.path.exists(config_file):
            raise RuntimeError("Location of manifest.conf not known.  "\
                "Provide it in .moby2.")
    manifest = moby2.util.mysql.ManifestConnector(configFile=config_file)
    return manifest.db, manifest.table_prefix

class TODDatabase:
    ctime_fields = ['ctime', 'ctime_start', 'ctime_end', 'merged_ctime']

    table_prefix = ''

    def __init__(self, config_file=None, season=None):
        self.db, self.table_prefix = get_actpol_db(config_file, season)

    def _table(self, table):
        return self.table_prefix + table

    def get_record(self, acq_id):
        """
        Return a TODRecord for the acq_id specified, which must be either
        the mce_data_acq_id or a basename/filename.
        """
        return TODRecord.from_database(self, acq_id)

    def select_tods(self, fields=None, clauses=[], id_only=False,
                    order=None,
                    **kwargs):
        """
        Select records from ACTpol TOD database.  Returns either a
        list of records, or a list of the mce_data_acq_ids if
        id_only=True.

        If you are only interested in certain fields from the tod_info
        table, pass fields=[...].  This bypasses the individual record
        queries that get the basename, etc.  It's way faster.  In this
        case the raw DB row data are returned.

        To select records using tod_info fields, pass in strings in
        the "clauses" list, or use kwargs to specify values or ranges
        for particular tod_info fields.

        E.g.:
            db = actpol.TODDatabase()

            # Match some values exactly:

            tods = db.select_tods(array='AR1', obs_detail='saturn')

            # Match a range of values -- pass a tuple:

            tods = db.select_tods(ctime_start=(1234560000,1234570000))

            # Match any from a list of values -- pass a list:

            tods = db.select_tods(obs_detail=['venus', 'saturn'])
            
        """
        clauses = [c for c in clauses] + self._assemble_clauses(kwargs)
        if fields is None:
            q = 'select mce_data_acq_id from %s ' % \
                self._table('tod_info')
        else:
            q = 'select %s from %s ' % (','.join(fields),
                                        self._table('tod_info'))
        if len(clauses) > 0:
            q += 'where ' + ' and '.join(clauses)
        if order is not None:
            if isinstance(order, basestring):
                order = [order]
            q += ' order by ' + ','.join(order)
        c = self.db.cursor()
        try:
            c.execute(q)
        except:
            print('Exception executing query: ', q)
            print()
            raise
        rows = c.fetchall()
        if fields is not None:
            return rows
        if id_only :
            # Return list of ids
            return [r[0] for r in rows]
        # Look up details
        recs = TODRecordList()
        for id_, in rows:
            recs.append(TODRecord.from_database(self, id_))
        return recs

    def select_acqs(self, fields=None, clauses=[], id_only=False,
                    order=None,
                    **kwargs):
        """
        Find acquisitions in the mce_data_acq table.  This handles
        kw_args in the same way as select_tods.

        Returns a TODRecordList of AcquisitionRecord objects.  If
        that's too slow, get raw DB rows by passing id_only=True or
        passing a list of desired fields in fields.

        Example::

           tdb = TODDatabase()
           recs = tdb.select_acqs(array='ar1', suffix='iv',
                                  ctime=(1379703312, 1379789712))
           prints recs[0].to_string()

        """
        clauses = [c for c in clauses] + self._assemble_clauses(kwargs)
        if fields is None:
            q = 'select id from %s ' % self._table('mce_data_acq')
        else:
            q = 'select %s from %s ' % (','.join(fields),
                                        self._table('mce_data_acq'))
        if len(clauses) > 0:
            q += 'where ' + ' and '.join(clauses)
        if order is not None:
            if isinstance(order, basestring):
                order = [order]
            q += ' order by ' + ','.join(order)
        c = self.db.cursor()
        try:
            c.execute(q)
        except:
            print('Exception executing query: ', q)
            print()
            raise
        rows = c.fetchall()
        if fields is not None:
            return rows
        if id_only:
            # Return list of ids
            return [r[0] for r in rows]
        # Look up details
        recs = TODRecordList()
        for id_, in rows:
            recs.append(AcquisitionRecord.from_database(self, id_))
        return recs

    def select_hk(self, clauses=[], order=None, **kwargs):
        """
        Find HK dirfiles in the merlin table.  This handles
        kw_args in the same way as select_tods.

        Example:

           tdb = TODDatabase()
           recs = tdb.select_hk(ctime=(1379703312, 1379789712))
           prints recs[0].to_string()
        """
        tables = {'merlin': self._table('merlin'),
                  'datafiles': self._table('datafiles')}
        q = 'select {merlin}.id, ctime, ctime_end, merged_ctime, '\
            'serialnum_start, serialnum_end, flags, datafile_id, '\
            '{datafiles}.name as basename '\
            'from {merlin} join {datafiles} on datafile_id={datafiles}.id '.\
            format(**tables)
        clauses = [c for c in clauses] + self._assemble_clauses(kwargs)
        clauses.append('(archive_type="hk" and succeeded="Y")')
        q += 'where ' + ' and '.join(clauses)
        if order is not None:
            if isinstance(order, basestring):
                order = [order]
            q += ' order by ' + ','.join(order)
        c = self.db.cursor()
        try:
            n_rows = c.execute(q)
        except:
            print('Exception executing query: ', q)
            print()
            raise
        # Look up details
        recs = TODRecordList()
        for i in range(n_rows):
            recs.append(HKRecord.from_cursor(c))
        return recs

    def _assemble_clauses(self, kwargs):
        def db_test(k, v):
            if v is None:
                return 'isnull(%s)' % k
            return '%s=%s' % (k, repr(v))
        clauses = []
        for k, v in list(kwargs.items()):
            if k in self.ctime_fields:
                # Convert string dates to unix timestamps
                v = convert_ctimes(v)
            if isinstance(v, list):
                clauses.append('(' + 
                               ' or '.join([db_test(k, x) for x in v]) +
                               ')')
            elif isinstance(v, tuple):
                clauses.append('(%s <= %s) and (%s < %s)' % \
                                   (repr(v[0]), k, k, repr(v[1])))
            else:
                clauses.append(db_test(k, v))
        return clauses

class Filebase(MobyFilebase):
    """
    Provide filenames for merged TODs in simple colossus storage.
    
    Requires that a connection be made to the ACTpol manifest
    database.
    """
    table_prefix = ''

    def __init__(self, manifest_config=None, nodes=None, db=None):
        """
        Settings are loaded, by default, from moby2.user_cfg.  Overrides:
        
        db is a _mysql.connection to use instead of loading
        actManifest.
                           
        manifest_config is the filename of ACTpol manifest config file
        to use.
        
        nodes is a list of storage_nodes to check for files.  Provide
        each node as a tuple, with entries
          (producer, node_name, root, structure)

        The "root" and "structure" are optional.  "root" specifies an
        alternative mount point and defaults to the database value.
        "structure" is used to describe the sub-directory structure of
        the storage node.  This be:

           'flat' / None: No structure, all files in <root>/<producer>

           'prefix5': Files in <root>/<producer>/<prefix>, where
                prefix is the first 5 characters of the file name.

        E.g. ('merlin', 'raid')
             ('mce1', 'mce1', None, 'prefix5')s
        """
        # Get defaults from moby config
        import moby2
        params = moby2.user_cfg.get('actpol_filebase', {})
        if nodes is None:
            nodes = params.get('nodes', [])
        if manifest_config is None:
            manifest_config = params.get('manifest_config', None)
        self.db = db
        if db is None and manifest_config is not None:
            self.db, self.table_prefix = get_actpol_db(manifest_config)
        self.nodes = {}
        self.node_clause = {}
        if nodes is not None:
            for node in nodes:
                if isinstance(node, basestring):
                    self._add_node(node)
                else:
                    self._add_node(*node)
            
    def _table(self, table):
        return self.table_prefix + table

    def _add_node(self, producer, name, root=None, structure=None):
        c = self.db.cursor()
        if c.execute('select id, root from %s where name="%s"' %
                     (self._table('storage_nodes'), name)) == 0:
            raise ValueError('Storage node "%s" not found' % name)
        id_, sroot = c.fetchone()
        if root is None:
            root = sroot
        node_cfg = self.nodes.get(producer, {})
        node_cfg[id_] = (name, root, structure)
        self.node_clause[producer] = ' or '.join(['node_id=%i' % nid 
                                                  for nid in list(node_cfg.keys())])
        self.nodes[producer] = node_cfg

    def _get_filename(self, producer, node_id, name):
        _, root, structure = self.nodes[producer][node_id]
        if structure is None:
            return os.path.join(root, producer, name)
        if structure == 'first5':
            return os.path.join(root, producer, name[:5], name)
        raise ValueError("Unknown storage node structure '%s'" % structure)

    def _get_file_copies(self, disk_name, id_, producer, single=False):
        if id_ is None:
            if single: return None
            return []
        c = self.db.cursor()
        # Find file copies on nodes we know about
        n = c.execute('select node_id from %s where '\
                          'datafile_id=%i and '\
                          'has_file="Y" and (%s)' % \
                          (self._table('file_copies'),
                           id_, self.node_clause[producer]))
        # Construct paths
        filenames = []
        for node_id, in c.fetchall():
            filenames.append(self._get_filename(producer, node_id, disk_name))
        if single:
            if len(filenames) == 0:
                return None
            return filenames[0]
        return filenames

    def id_from_name(self, name, producer='merlin', get_disk_name=False):
        c = self.db.cursor()
        # Find entry in datafiles table
        n = c.execute('select id, name from %s where ' 
                      'producer="%s" and ' 
                      'name="%s" or name="%s.zip"' % 
                      (self._table('datafiles'),
                       producer, name, name))
        if n == 0:
            return None
        if get_disk_name:
            return c.fetchone()
        else:
            return c.fetchone()[0]
        

    def filename_from_name(self, name, single=False, array=None):
        """
        Find files on present storage nodes and return a list of full
        paths.  This works for TODs, MCE data, and merged
        housekeeping.  "name" is the either:
        
        * TOD basename, in the form <ctime>.<ctime>.<ar>
        * MCE filename, <ctime>_suffix[extension]
        * Housekeeping basename, <ctime>.<ctime>.hk
        
        For MCE filenames, the "array" must also be specified.

        If single=True then instead of a list, one of the filenames is
        returned, or None is returned if there were no matches.
        """
        if array is None:
            producer = 'merlin'
        else:
            producer = {
                'ar1': 'mce1',
                'ar2': 'mce2',
                'ar3': 'mce3',
                }.get(array.lower(), array)
        datafile_id, disk_name = self.id_from_name(name, producer=producer,
                                                   get_disk_name=True)
        return self._get_file_copies(disk_name, datafile_id, producer, single=single)

    def filename_from_id(self, id_, single=False):
        """
        Same behavior as filename_from_name, except that instead of a
        basename the user should pass in the datafile_id.
        """
        c = self.db.cursor()
        # Find entry in datafiles table
        n = c.execute('select name,producer from %s where id=%i' % 
                      (self._table('datafiles'), id_))
        if n == 0:
            return None
        datafile_name, producer = c.fetchone()
        return self._get_file_copies(datafile_name, id_, producer, single=single)
        
