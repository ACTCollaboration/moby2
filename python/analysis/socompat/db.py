# Convert certain ACTPol metadata to simonsobs database formats.
#

import moby2
import sotoddb

__all__ = ['make_detdb', 'make_obsdb']


def make_detdb():
    """
    Convert recent ACTPol ArrayData to an sotoddb.DetDB.
    """
    db = sotoddb.DetDB()
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
#        ('pa7', 's20', {'pa': 'pa7', 'optics_tube': 3}),
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


def make_obsdb():
    db = sotoddb.ObsDB()
    cat = moby2.scripting.get_obs_catalog()
    c = db.conn.cursor()
    s = cat['tod_name'] != '1410143132.1018'
    for i in s.nonzero()[0]:
        c.execute('insert into obs (obs_id, timestamp) values (?,?)',
                  (cat[i]['tod_name'],
                   float(cat[i]['ctime'])))
    db.conn.commit()
    return db
