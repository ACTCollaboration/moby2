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

def load_tod_list(filename, fmt=None):
    if fmt is None:
        ext_map = {'.txt': 'txt',
                   '.sqlite': 'obsdb',
                   '.gz': 'obsdb',
                   '.fits': 'structdb'
               }
        fmt = ext_map[os.path.splitext(filename)[1]]
    if fmt == 'txt':
        tod_id = moby2.util.ascii.read_columns(filename)[0]
    elif fmt == 'obsdb':
        db = metadata.ObsDb.from_file(filename)
        tod_id = db.get()['obs_id']
    elif fmt == 'structdb':
        db = moby2.util.StructDB.from_fits_table(filename)
        tod_id = db['tod_id']
    return tod_id

def get_obs_catalog(catalog):
    if catalog:
        return moby2.util.StructDB.from_fits_table(catalog)
    return moby2.scripting.get_obs_catalog()

def relativify_paths(srcdir, odir):
    """Returns (srcdir_from_here, srcdir_relpath).  If srcdir is and odir
    are both relative, then srcdir_from_here is simply srcdir but
    srcdir_relpath is a prefix that can be used to reach srcdir
    starting from the odir.  Otherwise, srcdir_from_here is the
    absolute path to srcdir and srcdir_relpath is ''.

    """
    rel = ''
    if srcdir[0] == '/':
        pass
    elif odir[0] == '/':
        srcdir = os.path.abspath(srcdir)
    else:
        rel = os.path.relpath('./', odir)
    return (srcdir, rel)


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--tod-list', '-i')
    parser.add_argument('--output-dir', '-o', default='./')
    parser.add_argument('--force', action='store_true')

    sp = parser.add_subparsers(dest='module', metavar='module', help=
                               "Module to activate \"%(prog)s {module} -h\" for more help")
    
    p = sp.add_parser('obsdb', help="Write out an obsdb obsfiledb.")
    p.add_argument('--catalog', help="Alternative obs_catalog (FITS) to use.")

    p = sp.add_parser('obsfiledb', help="Write out an obsfiledb.")
    p.add_argument('--catalog', help="Alternative obs_catalog (FITS) to use.")

    p = sp.add_parser('detdb', help="Write out a detdb.")

    p = sp.add_parser('scan-hdf', help="Generically scan HDF5 files(s) looking for "
                      "obs_ids.")
    p.add_argument('sources', action='append', default=[])
    p.add_argument('--db-file', default='generic.sqlite')
    p.add_argument('--update', action='store_true', help=
                   "If output file exists, update it.")

    p = sp.add_parser('pointofs', help="Write out a per-TOD pointofs archive.")
    p.add_argument('infile')
    p.add_argument('h5file', help="Path of HDF5 output (relative to output dir).")
    p.add_argument('--dataset', default='pointofs')
    p.add_argument('--db-file', default='pointofs.sqlite')

    p = sp.add_parser('abscal', help="Write out a per-TOD pointofs archive.")
    p.add_argument('h5file')
    p.add_argument('--dataset', default='abscal')
    p.add_argument('--db-file', default='abscal.sqlite')

    p = sp.add_parser('timeconst', help="Write out a timeconst loader archive.")
    p.add_argument('scan', nargs='+', help='Top level dir to scan.')
    p.add_argument('--db-file', default='timeconst.sqlite')

    p = sp.add_parser('focalplane', help="Write out a per-TOD pointofs archive.")
    p.add_argument('spec_file', help="Filename to focal plane definitions.")
    p.add_argument('h5file', help="Path of HDF5 output (relative to output dir)." )
    p.add_argument('--db-file', default='focalplane.sqlite')

    p = sp.add_parser('cuts_release', help="Process a cuts release into cuts+cal archives.")
    p.add_argument('release_file')

    p = sp.add_parser('cuts_dir', help=
                      "Create an index for a single cuts depot result.")
    p.add_argument('src_dir', help="Path to cuts dir (relativity will be preserved).")
    p.add_argument('--db-file', default='cuts.sqlite')
    p.add_argument('--subset', nargs=2, action='append', default=[], help=
                   "State a detector subset to which this applies; e.g. "
                   "--subset dets:band f150")
    p.add_argument('--update', action='store_true', help=
                   "If output file exists, update it.")

    p = sp.add_parser('cal_dir', help=
                      "Create an index for a single cal depot result.")
    p.add_argument('src_dir', help="Path to cuts dir (relativity will be preserved).")
    p.add_argument('--db-file', default='cuts.sqlite')
    p.add_argument('--subset', nargs=2, action='append', default=[], help=
                   "State a detector subset to which this applies; e.g. "
                   "--subset dets:band f150")
    p.add_argument('--update', action='store_true', help=
                   "If output file exists, update it.")


    p = sp.add_parser('context', help="Write a context.yaml file.")
    p.add_argument('context_file')

    return parser


def _checkfile(filename, args, parser=None, updatable=False):
    full_path = os.path.join(args.output_dir, filename)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    err_str = None
    if os.path.exists(full_path):
        if updatable:
            if not args.update:
                err_str = f'Found existing {full_path}, pass --update to modify it.'

        elif not args.force:
            err_str = f'Found existing {full_path}, pass --force to regenerate.'

    if err_str:
        if parser is not None:
            return parser.error(err_str)
        return print(err_str)

    print(f'Creating {full_path}...')
    return full_path


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = get_parser()
    args = parser.parse_args(args=args)

    if args.module == None:
        parser.error('Select a submodule.')

    elif args.module == 'obsdb':
        fn1 = _checkfile('obsdb.sqlite', args, parser=parser)
        cat = get_obs_catalog(args.catalog)
        if args.tod_list:
            tods = load_tod_list(args.tod_list)
            print(f'Restricting TOD list to {len(tods)} items from {args.tod_list}...')
            idx = cat.select_inner({'tod_id': tods})
            n_bad = (idx<0).sum()
            if n_bad:
                print(' -- warning, did not match %i of %i tods from input list.' % (
                    n_bad, len(idx)))
                idx = idx[idx>=0]
            cat = cat[idx]
        socompat.make_obsdb(cat=cat).to_file(fn1)

    elif args.module == 'obsfiledb':
        fn1 = _checkfile('obsfiledb.sqlite', args)
        cat = get_obs_catalog(args.catalog)
        if args.tod_list:
            tods = load_tod_list(args.tod_list)
            print(f'Restricting TOD list to {len(tods)} items from {args.tod_list}...')
            idx = cat.select_inner({'tod_id': tods})
            n_bad = (idx<0).sum()
            if n_bad:
                print(' -- warning, did not match %i of %i tods from input list.' % (
                    n_bad, len(idx)))
                idx = idx[idx>0]
            cat = cat[idx]
        socompat.make_obsfiledb(cat=cat).to_file(fn1)
        
    elif args.module == 'detdb':
        fn1 = _checkfile('detdb.sqlite', args)
        socompat.make_detdb().to_file(fn1)

    elif args.module == 'scan-hdf':
        fn1 = _checkfile(args.db_file, args, parser=parser, updatable=True)

        if args.tod_list:
            tod_list = load_tod_list(args.tod_list)
        else:
            tod_list = None

        if os.path.exists(fn1):
            db = metadata.ManifestDb.from_file(fn1)
        else:
            scheme = metadata.ManifestScheme()\
                             .add_data_field('dataset')\
                             .add_exact_match('obs:obs_id')
            db = metadata.ManifestDb(scheme=scheme)

        for source_file in args.sources:
            print(f'Scanning {source_file}...')
            with h5py.File(source_file, 'r') as h:
                n_added = 0
                for k in h.keys():
                    if tod_list is None or k in tod_list:
                        db.add_entry({'dataset': k, 'obs:obs_id': k},
                                     source_file, replace=True)
                        n_added += 1
                print(f' ... found {n_added} entries to keep')
        db.to_file(fn1)

    elif args.module == 'pointofs':
        if args.tod_list:
            tods = load_tod_list(args.tod_list)
        else:
            tods = None
        # Clean this up ... what if hdf already exists, elsewhere, etc.
        fn1 = args.infile
        fn2 = os.path.join(args.output_dir, args.h5file)
        fn3 = _checkfile(args.db_file,
                         args, parser=parser)

        proddb = socompat.make_pointofs_hdf(fn1, fn2, dataset=args.dataset,
                                            obs_list=tods, hdf_relout=args.h5file)
        proddb.to_file(fn3)

    elif args.module == 'abscal':
        fn1 = _checkfile(args.db_file, args, parser=parser)
        proddb = socompat.metadata.get_abscal_proddb(args.h5file, dataset=args.dataset)
        proddb.to_file(fn1)

    elif args.module == 'timeconst':
        fn1 = _checkfile(args.db_file, args, parser=parser)
        if args.tod_list:
            tod_list = load_tod_list(args.tod_list)
        else:
            tod_list = None
        # Scan some directories.
        if len(args.scan) == 0:
            parser.error('No directories specified (use --scan)')

        scheme = metadata.ManifestScheme()\
                 .add_exact_match('obs:obs_id')\
                 .add_data_field('loader')\
                 .add_data_field('pa_hint')
        db = metadata.ManifestDb(scheme=scheme)
        entry = {
            'loader': 'actpol_timeconst'
        }
        TOD_ID_PAT =  '[0-9]{10}\.[0-9]{10}\.ar.'
        product_re = re.compile('(%s)\.tau' % (TOD_ID_PAT, ))
        for root_dir in args.scan:
            print(f'Working on {root_dir} ...')
            for root, dirs, files in os.walk(root_dir):
                if len(files):
                    print(f'  looking at {len(files)} in {root}')
                for f in files:
                    m = product_re.fullmatch(f)
                    if m is None:
                        continue
                    entry['obs:obs_id'] = m.group(1)
                    entry['pa_hint'] = 'pa' + m.group(1)[-1]
                    if tod_list is None or entry['obs:obs_id'] in tod_list:
                        db.add_entry(entry, filename=os.path.join(root, f))
        db.to_file(fn1)

    elif args.module == 'focalplane':
        fn1 = _checkfile(args.db_file, args, parser=parser)
        # For a bit of generality, request a map from array and
        # time_range to offset and polarization files.
        spec = yaml.safe_load(open(args.spec_file, 'rb'))
        # Prepare the output database...
        scheme = metadata.ManifestScheme()\
                         .add_data_field('dataset')\
                         .add_range_match('obs:timestamp')\
                         .add_exact_match('obs:pa')
        db = metadata.ManifestDb(scheme=scheme)
        # Write results to hdf5
        hdf_out = os.path.join(args.output_dir, args.h5file)
        with h5py.File(hdf_out, 'a') as h:
            for row in spec['table']:
                pa, t0, t1, pos_file, pol_file = row
                dset = f'{pa}_{t0}_{t1}'
                aman = socompat.metadata.load_detoffsets_file(
                    os.path.join(spec['prefix'], pos_file),
                    os.path.join(spec['prefix'], pol_file), pa=pa)
                # Convert to ResultSet and write out.
                rs = metadata.ResultSet(keys=['dets:readout_id', 'xi', 'eta', 'gamma'])
                for i, d in enumerate(aman.dets.vals):
                    rs.rows.append([d, aman['xi'][i], aman['eta'][i], aman['gamma'][i]])
                io.metadata.write_dataset(rs, h, dset, overwrite=args.force)
                db.add_entry({'dataset': dset, 'obs:pa': pa, 'obs:timestamp': (t0, t1)},
                             args.h5file)
        db.to_file(fn1)

    elif args.module == 'cuts_release':
        socompat.process_cuts_release(args.release_file, output_dir=args.output_dir)

    elif args.module == 'cuts_dir':
        src_dir, src_prefix = relativify_paths(args.src_dir, args.output_dir)
        if src_prefix != '':
            print(f'output_dir and src_dir are both relative, so target files '
                  f'will be prefixed with {src_prefix}')

        fn1 = _checkfile(args.db_file, args, parser=parser, updatable=True)
        db = None
        if os.path.exists(fn1):
            db = metadata.ManifestDb.from_file(fn1)
        subset = dict(args.subset)
        db = socompat.make_cuts_db(src_dir, db_in=db, source_prefix=src_prefix,
                                   restrictions=subset)
        db.to_file(fn1)

    elif args.module == 'cal_dir':
        src_dir, src_prefix = relativify_paths(args.src_dir, args.output_dir)
        if src_prefix != '':
            print(f'output_dir and src_dir are both relative, so target files '
                  f'will be prefixed with {src_prefix}')
        fn1 = _checkfile(args.db_file, args, parser=parser, updatable=True)
        db = None
        if os.path.exists(fn1):
            db = metadata.ManifestDb.from_file(fn1)
        subset = dict(args.subset)
        db = socompat.make_cal_db(src_dir, db_in=db, source_prefix=src_prefix,
                                  restrictions=subset)
        db.to_file(fn1)

    elif args.module == 'context':
        fn1 = _checkfile('context.yaml', args)
        socompat.write_context(fn1)

    else:
        parser.error(f'Module "{args.module}" not implemented.')
