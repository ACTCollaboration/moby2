"""This submodule provides functions to convert key moby2 data
structures into databases compatible with the sotodlib Context system.
It also helps to index ACTPol metadata archives (such as directories
full of cuts and calibration files) for use with the Context system.

"""

import moby2
from moby2.analysis import socompat
from sotodlib.core import metadata
from sotodlib import io

import argparse
import os
import re
import time
import yaml
import h5py

from tqdm import tqdm
import numpy as np

__all__ = ['make_detdb', 'make_obsdb', 'make_obsfiledb', 'make_cuts_db', 'make_cal_db',
           'make_pointofs_hdf', 'process_cuts_release', 'write_context']


TOD_ID_PAT =  '[0-9]{10}\.[0-9]{10}\.ar.'


def make_detdb():
    """
    Convert recent ACTPol ArrayData to an metadata.DetDb.
    """
    db = metadata.DetDb()
    db.create_table('base', [
        "`readout_id` varchar(16)",
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
            (k, k) for k in ['optical_sign', 'det_type', 'wafer', 'pol_family']
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
            data['readout_id'] = new_uid
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
        'scode string',
        'duration float',
        'obs_type string',
        'obs_detail string',
    ])
    for i in cat['tod_name'].argsort():
        tags = []
        row = cat[i]
        obs_id = row['tod_name']
        data = {
            'timestamp': float(row['ctime']),
            'duration': float(row['duration']),
        }
        data.update({k: row[k] for k in ['pa', 'obs_type', 'obs_detail', 'scode']})
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
    for row in tqdm(cat):
        obs_id = row['tod_name']
        if obs_id in ignore:
            continue
        f = filebase.get_full_path(obs_id)
        if f is None:
            continue
        detset = row['pa']
        obsfiledb.add_obsfile(f, obs_id, detset, None)
    return obsfiledb


def make_pointofs_hdf(offsets_file, hdf_out, dataset='pointofs',
                      obs_list=None, hdf_relout=None):
    # Read offsets ...
    data = moby2.util.StructDB.from_column_file(
        offsets_file, [('obs:obs_id', 0),
                       ('dx', 5),
                       ('dy', 6),
                       ('gamma', 6)])
    data['gamma'] = 0.

    if obs_list is not None:
        # Restrict ...
        idx = data.select_inner({'obs:obs_id': obs_list})
        assert(np.all(idx > 0))
        data = data[idx]

    if hdf_relout is None:
        hdf_relout = hdf_out
    with h5py.File(hdf_out, 'a') as h:
        rs = metadata.ResultSet(data.dtype.names)
        rs.rows.extend(list(data))
        io.metadata.write_dataset(rs, h, dataset)

    scheme = metadata.ManifestScheme()\
             .add_range_match('obs:timestamp')\
             .add_data_field('loader')\
             .add_data_field('dataset')
    man = metadata.ManifestDb(scheme=scheme)
    man.add_entry({'obs:timestamp': (0,2e9),
                   'dataset': dataset,
                   'loader': 'actpol_pointofs'},
                  hdf_relout,
                  commit=False)
    man.conn.commit()
    return man

def make_abscal_hdf(offsets_file, hdf_out, dataset='abscal'):
    if hdf_out is None:
        # Try to open offsets_file as hdf ...
        hdf_out = offsets_file
        offsets_file = None
    # Read offsets ...
    data = moby2.util.StructDB.from_column_file(
        offsets_file, [('obs:obs_id', 0),
                       ('dx', 5),
                       ('dy', 6),
                       ('gamma', 6)])
    data['gamma'] = 0.

    with h5py.File(hdf_out, 'a') as h:
        rs = metadata.ResultSet(data.dtype.names)
        rs.rows.extend(list(data))
        io.metadata.write_dataset(rs, h, dataset)

    scheme = metadata.ManifestScheme()\
             .add_range_match('obs:timestamp')\
             .add_data_field('loader')\
             .add_data_field('dataset')
    man = metadata.ManifestDb(scheme=scheme)
    man.add_entry({'obs:timestamp': (0,2e9),
                   'dataset': dataset,
                   'loader': 'actpol_pointofs'},
                  hdf_out,
                  commit=False)
    man.conn.commit()
    return man


def _cuts_and_cal_helper(root_dir, loader, restrictions, db_in, re_suffix,
                         source_prefix):
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
            print(f, m)
            if m is None:
                continue
            entry['obs:obs_id'] = m.group(1)
            if entry['obs:obs_id'] in ignore:
                continue
            db.add_entry(entry, filename=os.path.join(source_prefix, root, f))
    return db


def make_cuts_db(root_dir, loader=None, restrictions={},
                 db_in=None, source_prefix=''):
    """Scan root_dir for cuts results and add them to ManifestDb.

    """
    if loader is None:
        loader = 'actpol_cuts'
    return _cuts_and_cal_helper(root_dir, loader, restrictions, db_in, '\.cuts',
                                source_prefix=source_prefix)


def make_cal_db(root_dir, loader=None, restrictions={},
                 db_in=None, source_prefix=''):
    """Scan root_dir for cal results and add them to ManifestDb.

    """
    if loader is None:
        loader = 'actpol_cal'
    return _cuts_and_cal_helper(root_dir, loader, restrictions, db_in, '\.cal',
                                source_prefix=source_prefix)

def process_cuts_release(release_filename, temp_dir='temp/',
                         output_dir='./',
                         output_pattern='metadata_{category}.sqlite'):
    """
    Process a release file from cutslib.
    """
    if isinstance(release_filename, dict):
        cutsc = release_filename
    else:
        cutsc = yaml.safe_load(open(release_filename).read())

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    db_files = {k: os.path.join(output_dir, output_pattern).format(category=k)
                for k in ['cal', 'cuts', 'pcuts']}

    cuts_map = {}
    for k in cutsc['tags'].keys():
        #pa4_f150_s17_c11
        pa, fcode, scode, x = k.split('_', 3)
        key = f'{pa}_{scode}'
        cuts_map[key] = cuts_map.get(key, [])
        cuts_map[key].append((fcode, k))

    # Make temporary dbs -- this efficiently walks each tag tree.
    print('Making temporary dbs for each item...')
    for key, cuts in cuts_map.items():
        for fcode, k in cuts:
            print(f'  {fcode} {k}')
            for t in ['tag_out', 'tag_planet']:
                print(f'    {t}')
                temp_db_file = 'temp/_%s_%s.sqlite' % (k, t)
                if os.path.exists(temp_db_file):
                    print(f'        skipping because {temp_db_file} exists')
                    continue
                _tag = cutsc['tags'][k][t]
                base = '{depot}/TODCuts/{tag}/'.format(
                    depot=cutsc['depot'], tag=_tag)
                db = socompat.make_cuts_db(base)
                db.to_file(temp_db_file)
            for t in ['tag_cal']:
                print(f'    {t}')
                temp_db_file = 'temp/_%s_%s.sqlite' % (k, t)
                if os.path.exists(temp_db_file):
                    print(f'        skipping because {temp_db_file} exists')
                    continue
                _tag = cutsc['tags'][k][t]
                base = '{depot}/Calibration/{tag}/'.format(
                    depot=cutsc['depot'], tag=_tag)
                db = socompat.make_cal_db(base)
                db.to_file(temp_db_file)

    # Join them together.
    scheme = metadata.ManifestScheme() \
                     .add_exact_match('obs:obs_id') \
                     .add_data_field('dets:band') \
                     .add_data_field('loader')
    scutsdb = metadata.ManifestDb(scheme=scheme)
    pcutsdb = metadata.ManifestDb(scheme=scheme)
    caldb = metadata.ManifestDb(scheme=scheme)

    print()
    print('Joining temp dbs together')
    for key, cuts in cuts_map.items():
        for fcode, k in cuts:
            print(f'  {fcode} {k}')
            for t, db, loader in [('tag_out', scutsdb, 'actpol_cuts'),
                                  ('tag_planet', pcutsdb, 'actpol_cuts'),
                                  ('tag_cal', caldb, 'actpol_cal')]:
                db_in = metadata.ManifestDb('temp/_%s_%s.sqlite' % (k, t))
                c_in = db_in.conn.execute('select `obs:obs_id`,name from map join files on map.file_id=files.id')
                c = db.conn.cursor()
                for row in tqdm(c_in):
                    obs_id, filename = row
                    _ = c.execute('insert into files (name) values (?)', (filename,))
                    fid = c.lastrowid
                    _ = c.execute('insert into map (`obs:obs_id`, `dets:band`, `loader`, `file_id`) '
                              'values (?,?,?,?)', (obs_id, fcode, loader, fid))
                    # Cursors are expensive.
    #                db.add_entry({'obs:obs_id': obs_id,
    #                              'dets:band': fcode,
    #                              'loader': loader},
    #                             filename=filename, commit=False)

    scutsdb.to_file(db_files['cuts'])
    pcutsdb.to_file(db_files['pcuts'])
    caldb.to_file(db_files['cal'])
    return db_files

def write_context(filename, metadata_root='./'):
    context_info = {
        'tags': {
            'metadata_root': metadata_root,
            },
        'imports': [
            'sotodlib.io.metadata',
            'moby2.analysis.socompat',
        ],
        'obsdb': '{metadata_root}/obsdb.sqlite',
        'detdb': '{metadata_root}/detdb.sqlite',
        'obsfiledb': '{metadata_root}/obsfiledb.sqlite',
        'obs_colon_tags': ['band', 'name'],
        'obs_loader_type': 'actpol_moby2',
        'metadata': []
    }
    md = yaml.safe_load("""
metadata:
  - db: '{metadata_root}/focalplane.sqlite'
    name: 'focal_plane'
  - db: '{metadata_root}/pointofs.sqlite'
    name: 'boresight_offset'
  - db: '{metadata_root}/metadata_abscal.sqlite'
    name: 'abscal&cal'
  - db: '{metadata_root}/metadata_cal.sqlite'
    name: 'relcal&cal'
  - db: '{metadata_root}/metadata_cuts.sqlite'
    name: 'glitch_flags&flags'
  - db: '{metadata_root}/metadata_pcuts.sqlite'
    name: 'source_flags&flags'
  - db: '{metadata_root}/timeconst.sqlite'
    name: timeconst&

""")
    context_info['metadata'].extend(md['metadata'])
    yaml.dump(context_info, open(filename, 'w'), sort_keys=False)
