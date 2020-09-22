"""This submodule provides functions to convert key moby2 data
structures into databases compatible with the sotodlib Context system.
It also helps to index ACTPol metadata archives (such as directories
full of cuts and calibration files) for use with the Context system.

"""

import moby2
from sotodlib.core import metadata

import os
import re

__all__ = ['make_detdb', 'make_obsdb', 'make_obsfiledb', 'make_cuts_db', 'make_cal_db']


TOD_ID_PAT =  '[0-9]{10}\.[0-9]{10}\.ar.'


def make_detdb():
    """
    Convert recent ACTPol ArrayData to an metadata.DetDb.
    """
    db = metadata.DetDb()
    db.create_table('base', [
        "`name` varchar(8)",
        "`array_name` varchar(8)",
        "`old_det_id` integer",
        "`band` varchar(8)",
        "`optical_sign` float",
        "`optics_tube` integer",
        "`det_type` varchar(8)",
        "`wafer` varchar(8)",
        "`pol_family` varchar(8)",
    ])
    db.create_table('mce', [
        "`row` integer",
        "`col` integer",
        "`bias_line` integer",
        "`bias_name` varchar(16)",
    ])
    db.create_table('sky', [
        "`x` float",
        "`y` float",
        "`angle` float",
    ])
    db.create_table('array', [
        "`x` float",
        "`y` float",
        "`angle` float",
    ])

    # The only reason to avoid pa1/2/3 is the weird wiring change in
    # pa1 from s13-s14.  We just need to decide whether to have a
    # fixed det_uid change location (that's what the ArrayData do) or
    # to have a fixed (new) det_uid hold location but change MCE
    # address.

    targets = [
        ('pa4', 's17', {'array_name': 'pa4', 'optics_tube': 1}),
        ('pa5', 's17', {'array_name': 'pa5', 'optics_tube': 2}),
        ('pa6', 's17', {'array_name': 'pa6', 'optics_tube': 3}),
        ('pa7', 's20', {'array_name': 'pa7', 'optics_tube': 3}),
    ]

    for pa, scode, info in targets:
        adata = moby2.scripting.products.get_array_data_new(
            {'array': pa, 'season': scode})
        base_map = [
            (k, k) for k in ['optical_sign', 'det_type', 'pol_family']
        ]
        other_maps = [
            ('mce', [
                (k, k) for k in ['row', 'col', 'bias_line', 'bias_name']]),
            ('sky', [
                (f'sky_{k}', k) for k in ['x', 'y', 'angle']]),
            ('array', [
                (f'array_{k}', k) for k in ['x', 'y', 'angle']]),
        ]

        for u in adata['det_uid']:
            # The `tolist` here is used to convert from numpy scalar
            # to native Python types -- otherwise sqlite gets
            # confused.
            new_uid = '%s_%04i' % (pa, u)
            # base table:
            data = {dest: adata[src][u].tolist() for src, dest in base_map}
            data['old_det_id'] = u.tolist()
            data['band'] = 'f%03i' % float(adata['nom_freq'][u])
            data.update(info)
            db.add_props('base', new_uid, **data,
                         commit=False)
            for table_name, table_map in other_maps:
                data = {dest: adata[src][u].tolist()
                        for src, dest in table_map}
                db.add_props(table_name, new_uid, **data,
                             commit=False)
    db.conn.commit()
    return db


def make_obsdb(cat=None):
    """Convert an ObsCatalog to an ObsDb.  If a catalog is not provided,
    the default is loaded.

    """
    if cat is None:
        cat = moby2.scripting.get_obs_catalog()
    obsdb = metadata.ObsDb()
    obsdb.add_obs_columns([
        'timestamp float',
        'pa string',
        'obs_type string',
        'obs_detail string',
    ])
    for i in cat['tod_name'].argsort():
        tags = []
        row = cat[i]
        obs_id = row['tod_id']
        data = {
            'timestamp': float(row['ctime'])
        }
        data.update({k: row[k] for k in ['pa', 'obs_type', 'obs_detail']})
        if row['obs_type'] == 'planet':
            tags.extend(['planet', row['obs_detail']])
        obsdb.update_obs(obs_id, data=data, tags=tags, commit=False)

    obsdb.conn.commit()
    return obsdb


def make_obsfiledb(cat=None, filebase=None, detdb=None, db_in=None,
                   ignore_duplicates=True):
    """Make an ObsFileDb using the moby2 filebase and an ObsCatalog.  If a
    catalog is not provided, the default is loaded.

    """
    if cat is None:
        cat = moby2.scripting.get_obs_catalog()
    if detdb is None:
        detdb = make_detdb()
    ignore = []
    if db_in is None:
        obsfiledb = metadata.ObsFileDb()
        pas = sorted(list(set(cat['pa'])))
        for pa in pas:
            names = detdb.dets(props={'array_name': pa})['name']
            obsfiledb.add_detset(pa, names)
    else:
        obsfiledb = db_in
        if ignore_duplicates:
            ignore = db_in.get_obs()
    if filebase is None:
        filebase = moby2.scripting.get_filebase()
    for row in cat:
        obs_id = row['tod_id']
        if obs_id in ignore:
            continue
        f = filebase.get_full_path(obs_id)
        if f is None:
            continue
        detset = row['pa']
        obsfiledb.add_obsfile(f, obs_id, detset, None)
    return obsfiledb


def _cuts_and_cal_helper(root_dir, loader, restrictions, db_in, re_suffix):
    scheme = metadata.ManifestScheme()\
             .add_exact_match('obs:obs_id')\
             .add_data_field('loader')
    # Additional restrictions...
    for k in restrictions:
        scheme.add_data_field(k)
    if db_in is None:
        db = metadata.ManifestDb(scheme=scheme)
        ignore = []
    else:
        db = db_in
        ignore = [r[0] for r in db.conn.execute('select distinct `obs:obs_id` from map')]
    product_re = re.compile('(%s)%s' % (TOD_ID_PAT, re_suffix))
    entry = dict(restrictions)
    entry['loader'] = loader
    for root, dirs, files in os.walk(root_dir):
        for f in files:
            m = product_re.fullmatch(f)
            if m is None:
                continue
            entry['obs:obs_id'] = m.group(1)
            if entry['obs:obs_id'] in ignore:
                continue
            db.add_entry(entry, filename=os.path.join(root, f))
    return db


def make_cuts_db(root_dir, loader=None, restrictions={},
                 db_in=None):
    """Scan root_dir for cuts results and add them to ManifestDb.

    """
    if loader is None:
        loader = 'actpol_cuts'
    return _cuts_and_cal_helper(root_dir, loader, restrictions, db_in, '\.cuts')


def make_cal_db(root_dir, loader=None, restrictions={},
                 db_in=None):
    """Scan root_dir for cal results and add them to ManifestDb.

    """
    if loader is None:
        loader = 'actpol_cal'
    return _cuts_and_cal_helper(root_dir, loader, restrictions, db_in, '\.cal')
