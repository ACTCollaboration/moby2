from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Functions to help scripts find information needed for analyzing TODs,
mapping, etc.

The model here is that these routines handle all the details of
decoding parameters to arrive at a conclusion about what files to
load; the calling script should need to do only minimal inspection of
the parameters before making the call.

To support different instruments / eras, the "instrument" option is
used to indicate a paradigm for finding and reading data.  When not
provided, instrument will be sought in the user's config file.
"""

import os, glob
import numpy as np
import warnings

import moby2
import moby2.mapping.fits_map as fits_map
import moby2.instruments as instrument
from moby2.util import mce

DEFAULT_INSTRUMENT = 'unknown'
shared_libraries = {"actpol":"actpol_shared", "mbac":"act_shared"}

from . import execcfg
input_chooser = execcfg.InputChooser()
exec_chooser = input_chooser.get_config

def translate_array_name(array, instrument):
    if instrument == "mbac":
        if "mbac" in array: return array
        else: return {"ar1":"mbac145","ar2":"mbac215","ar3":"mbac280"}[array.lower()]
    elif instrument == "actpol":
        if "actpol" in array: return {"actpol1":"ar1",
                                      "actpol2":"ar2",
                                      "actpol3":"ar3",
                                      "actpol4":"ar4",}[array]
        if 'pa' in array:
            return {"pa1":"ar1",
                    "pa2":"ar2",
                    "pa3":"ar3",
                    "pa4":"ar4",}[array]
        else: return array.lower()
    else:
        raise ValueError("Unknown instrument for array translation")

def trace(level, msg):
    moby2.util.log.trace('moby', level, msg)

def warn_deprecated(message):
    warnings.warn(message, Warning, stacklevel=2)

def _get_instrument(params={}, tod_info=None, instrument=None):
    """
    Return the instrument identification string, which may have been
    passed by the user through the 'instrument' argument or in the
    params dict, or be stored in the TODInfo, but might just need to
    be looked up in the user config file.
    """
    if instrument is not None:
        return instrument
    if params is not None and 'instrument' in params:
        return params['instrument']
    if tod_info is not None and tod_info.instrument is not None:
        return tod_info.instrument
    return moby2.util.get_user_config().get('default_instrument',
                                            DEFAULT_INSTRUMENT)

def _get_season_array_name(params = {}, tod_info = None, instrument = None):
    """
    Returns (season, array_name), based on the params dict, TODInfo
    tod_info, and instrument string.

    params['season'] and params['array_name'] are used if they are not
    None.  Otherwise, tod_info is inspected.
    """
    instrument = _get_instrument(params, tod_info, instrument)

    # Season:
    season = params.get('season')
    if season is None and tod_info is not None:
        season = tod_info.season

    # Array name:
    array_name = params.get('array_name')
    if array_name is None and tod_info is not None:
        array_name = tod_info.array

    return season, translate_array_name(array_name, instrument)

def get_tod_info_new(params):
    if not 'instrument' in params:
        params['instrument'] = _get_instrument(params)
    return input_chooser.get_special('tod_info', {'id': params['tod_id']})['tod_info']

def get_tod_info(params, tod=None):
    """
    Get TODInfo object for arguments.

    If tod is None, a tod_info is generated based entirely on the
    params dict, which will probably need to contain 'basename' and
    possibly 'season', 'array', ...

    If tod is not None, it should be a string (filename), a TODInfo
    object, or a TOD object with a valid info member.
    """
    if os.path.isfile(params.get("filename")):
        filename = params.get("filename")
    else:
        fb = get_filebase()
        filename = fb.filename_from_name(os.path.basename(params.get("filename")),
                                         single=True)
    if tod is None:
        self = moby2.TODInfo()
        self.filename = filename
        # Note .runfile can fail if user is passing in a basename
        # instead of a filename... just issue a warning in this case.
        try:
            self.runfile = mce.MCERunfile.from_dirfile(filename,
                                                       prefix='Extras/')
        except:
            trace(4, 'Did not load runfile for %s' % filename)
            self.runfile = None
            
        if filename.endswith('.zip'):
            filename = filename[:-4]
        self.basename = filename.strip('/').split('/')[-1]
        self.name = self.basename
        self.tod_id = self.basename  #new!
        self.instrument = _get_instrument(params)
        try:
            # Assumes ctime.ctime.array or ctime.ctime-and-it's-2007.
            pieces = self.basename.split('.')
            self.ctime = int(pieces[0])
            self.day_str, self.time_str = moby2.util.ctime.to_string(
                self.ctime).split()
            # Rough but good enough; uses ~Feb 1 as a pivot.
            self.season = '%i' % int((self.ctime - 949381200.)/31557600.0
                                           + 2000)
            if len(pieces) > 2:
                ar = pieces[2]
            else:
                ar = 'ar1'
            self.array = translate_array_name(ar, self.instrument)
        except:
            self.ctime = None
            self.day_str, self.time_str = None, None
            self.array = 'unknown'
            self.instrument = 'unknown'
            self.season = 'unknown'
            raise
 
        self.array_data = get_array_data({}, tod_info=self)
        return self
    
    if isinstance(tod, moby2.TOD) and getattr(tod, 'info', None) is not None:
        return tod.info

    if isinstance(tod, moby2.tod.TODInfo):
        return tod
        
    raise RuntimeError("Could not convert %s to a TODInfo." % tod)

def get_array_data_new(params):
    params = params.copy()
    if 'tod_id' in params:
        tod_info = input_chooser.get_special('tod_info', {'id': params['tod_id']})
        tod_info.update(params)
        params = tod_info
    if not 'instrument' in params:
        params['instrument'] = _get_instrument(params)
    ad_params = input_chooser.get_special('array_data', params)
    return moby2.tod.ArrayData.from_fits_table(ad_params['filename']).to_struct_db()

def get_array_data(params, tod_info=None):
    """
    Get array data for arrays.  params should be a dict
    """
    instrument = _get_instrument(params, tod_info=tod_info)
    try:
        label = shared_libraries[instrument]
    except:
        raise _UnhandledInstrument(instrument)
    depot_params = params.get('depot', {'label': label}) 
    depot = moby2.scripting.get_depot(depot_params)
    season, array = _get_season_array_name(params, tod_info, instrument)
    filename = os.path.join(depot.depot_path, 'ArrayData/%s/%s/default.fits'%(season,array))
    return moby2.tod.ArrayData.from_fits_table(filename)

def get_tod(*params, **kw_params):
    """Get a TOD.  Parameters are backend specific... but there is only
    one backend.  See moby2.instruments.actpol.get_tod.

    Historically the call had to be made like this:

       get_tod({'filename': 'tod.zip', 'det_uid': [1,2,3]})

    But now you can pass any number of such dictionaries:

       get_tod({'filename': filename}, {'det_uid': [1,2,3]})

    Or pass a mix of dicts and keyword args (the latter will
    override):

       get_tod({'filename': filename}, det_uid=[1,2,3])

    And in fact the first argument can be a string, which will be
    jammed in as 'filename':

       get_tod(filename, det_uid=[1,2,3])

    Returns the TOD object if all goes well.
    """
    all_params = {}
    for ip, p in enumerate(params):
        if ip == 0 and isinstance(p, str):
            p = {'filename': p}
        all_params.update(p)
    all_params.update(kw_params)

    instrument = _get_instrument(all_params)
    if instrument in ["actpol","mbac"]:
        from moby2.instruments.actpol import get_tod as get_act_tod
        return get_act_tod(**all_params)
    else:
        raise _UnhandledInstrument(instrument)
    
class _UnhandledInstrument(Exception):
    def __init__(self, instrument):
        self.value = instrument
    def __str__(self):
        return 'Instrument "%s" is not handled by this method.' % \
            (self.value)


def get_depot(params, default_path=None):
    if params is None:
        if default_path is None:
            return None
        return moby2.util.depot.Depot(default_path)
    if isinstance(params, basestring):
        params = {'path': params}
    if params.get('act_depot', False):
        return moby2.util.Depot.for_ACT(params.get('path', None))
    if params.get('label') is not None and params.get('path') is None:
        # Look up the depot by label in the moby config
        depot_list = moby2.util.get_user_config().get('depots', [])
        if not params['label'] in depot_list:
            raise ValueError("Unknown depot label '%s'; update your"\
                " .moby2 file?" % params['label'])
        params = depot_list[params['label']]
    return moby2.util.depot.Depot(params.get('path', None))


def get_product_filename(params, cls=None, **kwargs):
    """
    Perform standard decoding of params, a dict containing indications
    of where to find a file.  The cases to support are:

    1. Absolute path
    params = { 'filename': 'path/to/file.txt' }

    2. File relative to depot:
    params = { 'depot': '/path/to/depot',
               'filename': 'location/within/depot.txt' }

    3. Tagged result in depot (uses depot_structure from class/object cls):
    params = { 'depot' : '/path/tp/depot',
               'tag': 'results_2000' }
    or
    params = { 'depot' : {'label':'my_label'},
               'tag': 'results_2000' }
       where depot label is defined in .moby2


    4. Tagged result in depot with custom depot structure:
    params = { 'depot' : '/path/tp/depot',
               'structure': 'MyClass/{season}/{tag}.txt',
               'tag': 'results_2000' }

    The kwargs are all passed to depot.get_full_path.  cls (the class
    for the target object) is likely to be needed if 'tag' is used and
    a 'structure' override is not provided.
    """
    structure, tag = params.get('structure'), params.get('tag')
    if structure is None and tag is None:
        structure, tag = '{tag}', params.get('filename')

    depot = params.get('depot')
    if depot is None:
        filename = params.get('filename')
    else:
        depot = get_depot(params.get('depot'))
        filename = depot.get_full_path(cls=cls, tag=tag, structure=structure,
                                       **kwargs)
    return filename


# map_from_spec should not be in this module.
def map_from_spec(params, tod_limits=None,
                  dtype='float32', wtype='int32'):
    """
    Get map based on params dict and possibly on the pointing limits
    tod_limits = ((xmin, xmax),(ymin,umax))
    """
    if params.get('center') is not None:
        x0, y0 = params['center']
        dx, dy = params['size']
        x0, x1 = x0-dx/2, x0+dx/2
        y0, y1 = y0-dy/2, y0+dy/2
    else:
        if 'size' in params:
            print('Warning: map "size" is ignored if "center" is not specified.')
        x0, x1 = min(tod_limits[0]), max(tod_limits[0])
        y0, y1 = min(tod_limits[1]), max(tod_limits[1])
    delt = params.get('pitch', 3.75/3600)
    
    pixn = fits_map.fitsPixelization(
        projection=fits_map.LinearProjection(),
        grid=fits_map.fitsGrid(cdelt=(delt,delt)))
    pixn.mapgrid = pixn.grid.getMapGrid((x0,x1),(y0,y1))
    map0 = fits_map.spaceMap(pixn=pixn, dtype=dtype, wtype=wtype)
    return map0


def get_tod_list(params, cmdline_args=None):
    """
    Obtain a list of TOD filenames, as specifed in the dict "params".
    """
    # Get a list of TOD filenames or basenames or other instructions
    src, src_args = params['source'][0], params['source'][1:]
    if src == 'command_line':
        tod_list = cmdline_args
    elif src == 'tod_list':
        # src_args = (filename, column_index)
        tod_list = []
        for line in open(src_args[0]):
            words = line.split()
            if len(words) == 0 or words[0][0] == '#':
                continue
            tod_list.append(words[src_args[1]])
    elif src == 'glob':
        # User provides a pattern to glob.
        tod_list = glob.glob(src_args[0])
    elif src == 'list':
        # User provides a list inline.
        tod_list = src_args[0]
    else:
        raise ValueError('Unknown tod_list spec "%s"' % src)

    # Convert list into filenames, if it's not already.
    if params.get('use_database', False):
        warn_deprecated('instead of scripting.get_tod(use_database=True), '\
            'use scripting.get_tod(use_filebase=True)')
        if 'use_filebase' not in params:
            params['use_filebase'] = params['use_database']

    if params.get('use_filebase', False):
        mfb = get_filebase()
        for i,t in enumerate(tod_list):
            if not os.path.exists(t):
                t = mfb.filename_from_name(t, single=True)
                if t is not None:
                    tod_list[i] = t

    return tod_list

def get_filebase(params=None):
    """
    Decode "params" dictionary and return a MobyFilebase object that
    can turn TOD references into filenames.
    """
    if params is None:
        params = moby2.user_cfg.get('filebase', {})

    fb_type = params.get('type')
    fb_par = params.get('params', {})
    if fb_type is None:
        fb = moby2.util.filebase.NullFilebase()
    elif fb_type == 'actpol_manifest':
        from moby2.instruments import actpol
        fb = actpol.Filebase(**fb_par)
    elif fb_type == 'filesystem':
        fb = moby2.util.filebase.FilesystemFilebase(**fb_par)
    elif fb_type == 'multi':
        fbs = [get_filebase(p) for p in params.get('components')]
        fb = moby2.util.filebase.MultiFilebase(fbs)
    else:
        raise ValueError("Unknown filebase type %s" % fb_type)
    return fb

def get_tod_id(tod_id=None, tod_info=None, tod=None, filename=None):
    if tod_id is not None:
        return tod_id
    if tod_info is not None:
        if isinstance(tod_info, basestring):
            return tod_info
        return tod_info.name
    if isinstance(tod, basestring):
        return get_tod_id(filename=tod)
    if hasattr(tod, 'info'):
        return tod.info.tod_id
    if filename is not None:
        basename = os.path.split(filename)[1]
        for term in ['.zip', '.fits', '.hdf']:
            if basename.endswith(term):
                return basename[:-len(term)]
        return basename
    return None

def get_detector_offsets(params, tod_info=None, det_uid=None, tod=None):
    """
    Decode "params" dictionary and return a FocalPlane object.
    Depending on params, tod_info may be inspected.  If det_uid is not
    None, then the FocalPlane will be reduced to match det_uid before
    being returned.

    The filename from which to load the offsets is specified by
    params['filename'] and related parameters, as decoded by
    get_product_filename.

    File format is specified by params['format'] (though
    params['source'] will work if 'format' is not specified), which
    defaults to 'ascii'.  Formats:
    
      'ascii'     -- load from ascii file.  A column map must be
                     provided in params['columns'], of the form
                     [('basename', 0),('x', 1),('y',2),('mask',3)].
                     The 'mask' column is optional, and will be used
                     to either exclude or include results depending on
                     whether the value in that column matches the
                     value of params['match_bad'], or
                     params['match_good'], respectively.

      'fp_file'   -- load from FPFitFile (output format of fp_fit).

      'po_file'   -- load from ACT-style pointing offsets file.  Provide
                     the filename in params['filename'].
    """
    if '_execcfg' in params:
        tod_id = get_tod_id(tod=tod, tod_info=tod_info)
        params = exec_chooser(params['_execcfg'], tod_id=tod_id)

    format = params.get('format', params.get('source', 'ascii'))
    filename = get_product_filename(params)

    if format == 'po_file':
        fplane = moby2.pointing.FocalPlane.from_ACT_po_file(filename)

    elif format == 'super_po_file':
        fplane = moby2.pointing.FocalPlane.from_super_po_file(filename)

    elif format == 'fp_file':
        import moby2.analysis.fp_fit as fp_fit
        fits = fp_fit.FPFitFile.from_columns_file(filename)
        fplane = moby2.pointing.FocalPlane(
            x=fits.x0, y=fits.y0, mask=fits.ok, det_uid=fits.det_uid)

    elif is_ascii_columns(format):
        col_defs = params['columns']
        cdb = moby2.util.StructDB.from_column_file(
            filename, col_defs)

        if 'mask' in cdb.dtype.names:
            mask = np.ones(len(cdb), bool)
            if 'match_good' in params:
                mask *= (cdb['mask'] == params['match_good'])
            if 'match_bad' in params:
                mask *= (cdb['mask'] != params['match_bad'])
            cdb = cdb[mask]
        fplane = moby2.pointing.FocalPlane(
            x=cdb['x'], y=cdb['y'], det_uid=cdb['det_uid'])

    elif format == 'h5':
        dset_name = params.get('dataset', 'focalplane')
        cdb = moby2.util.StructDB.from_hdf(filename, dataset=dset_name)
        # subselect for this array / season...

        fplane = moby2.pointing.FocalPlane(
            x=cdb['x'], y=cdb['y'], det_uid=cdb['det_uid'])
        
    else:
        raise ValueError('Unknown focal plane source "%s"' % source)

    # Match the detectors to the TOD?
    if det_uid is not None:
        fplane = fplane.subset(det_uid=det_uid)

    return fplane.subset(det_uid=det_uid)

def get_focal_plane(params, tod_info=None, det_uid=None, tod=None):
    """
    Decode "params" dictionary and return a FocalPlane object.
    Optionally add polarization angle information and pointing
    offsets.

    params is a dict with blocks describing different actions::

      params = {
        'detector_offsets': {...},
        'polarization_angles': {...},
        'pointing_shifts': {...},
      }

    Parameters in detector_offsets are passed to get_detector_offsets;
    parameters in pol_source (if not None) are passed to
    get_polarization angles.  Parameters in shift_generator are passed
    to get_pointing_offset.
    
    E.g.::
    
      params = {
        'detector_offsets': {
           'format': 'fp_file',
           'filename': 'template_ar1_130821s.txt' },
      }

    Depending on params, tod_info and stuff may be inspected.

    """
    if '_execcfg' in params:
        tod_id = get_tod_id(tod=tod, tod_info=tod_info)
        params = exec_chooser(params['_execcfg'], tod_id=tod_id)

    if 'detector_offsets' in params:
        # New style focal plane description...
        do_pars = params['detector_offsets']
    else:
        # One big jumble
        do_pars = params 
    # This should be a FocalPlane object.
    fplane = get_detector_offsets(do_pars, tod_info=tod_info,
                                  det_uid=det_uid, tod=tod)
        
    # Load polarization as well?
    pol_pars = params.get('polarization_angles')
    if 'pol_source' in params:
        warn_deprecated('pointing description block should use "polarization_angles" '
                        'instead of "pol_source"')
        if pol_pars is None:
            pol_pars = params.get('pol_source')
    if not pol_pars is None:
        pol = get_polarization_angles(pol_pars, tod=tod, tod_info=tod_info)
        pol_ok, pol_ang = pol.get_property('pol_angle', det_uid=fplane.det_uid)
        fplane.mask *= pol_ok
        fplane.phi = pol_ang

    # Apply a global pointing model offset thing?
    shift_pars = params.get('pointing_shifts')
    if 'shift_generator' in params:
        warn_deprecated('pointing description block should use "pointing_shifts" '
                        'instead of "shift_generator"')
        if shift_pars is None:
            shift_pars = params.get('shift_generator')
    if shift_pars is not None:
        dx_dy = get_pointing_offset(
            shift_pars, tod=tod, tod_info=tod_info, source_offset=False)
        if dx_dy is None:
            raise RuntimeError("Could not determine source offset!")
        fplane.x += dx_dy[0] * np.pi/180
        fplane.y += dx_dy[1] * np.pi/180

    # Match the detectors to the TOD?
    if det_uid is not None:
        fplane = fplane.subset(det_uid=det_uid)

    return fplane.subset(det_uid=det_uid)


def get_time_constants(params, tod_info=None):
    """
    Decode "params" dictionary and return a TimeConstants object.
    Depending on params, tod_info may be inspected.

    Very different things are done depending on the value of
    params['source']:
    
      'template'  -- load season/array default, taking season/array
                     from tod_info if they are not given in params.

      'file'      -- load from params['filename']
    """
    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'],
                              tod_id=get_tod_id(tod_info=tod_info))

    cls = moby2.detectors.TimeConstants
    filename = get_product_filename(params, cls=cls, tod_info=tod_info)

    source = params.get('source', params.get('format'))

    if source in [None, 'file', 'depot']:
        return cls.read_from_path(filename)

    elif is_ascii_columns(source):
        return moby2.detectors.TimeConstants.from_columns_file(
            filename, fields=['det_uid','tau'], columns=params['columns'])

    elif source == 'fp_file':
        import moby2.analysis.fp_fit as fp_fit
        fits = fp_fit.FPFitFile.from_columns_file(filename)
        tc = moby2.detectors.TimeConstants()
        tc.det_uid = fits.det_uid
        tc.tau = fits.tau
        return tc

    elif source == 'act_file':
        f3dB = params.get("f3dB", False) 
        return moby2.detectors.TimeConstants.from_ACT_file(filename, f3dB)

    elif source == 'constant':
        tc = cls()
        tc.det_uid = np.arange(params['n_det'])
        tc.tau = params['value'] + tc.det_uid*0.
        return tc

    raise ValueError('unknown time constants source "%s"'% params.get('source'))


def get_polarization_angles(params, tod=None, tod_info=None):
    if '_execcfg' in params:
        tod_id = get_tod_id(tod=tod, tod_info=tod_info)
        params = exec_chooser(params['_execcfg'], tod_id=tod_id)

    cls = moby2.detectors.PolarizationAngles
    filename = get_product_filename(params, cls)
    polang = cls.read_from_path(
        filename,
        format=params.get('format'),
        columns=params.get('columns'),
        units=params.get('units'),
        convention=params.get('convention'),
        bad_pol_value=params.get('fail_value'))
    if params.get('mirror', False):
        polang.pol_angle[:] = np.pi/2 - polang.pol_angle[:]
    return polang
        

def get_hwp_angles(params, tod=None):
    """
    Load HWP position angles.  If tod is passed in, the
    angles will be trimmed to match the samples loaded in TOD.
    """
    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'],
                              tod_id=get_tod_id(tod=tod))

    trimmed = False
    source = params.get('source', None)
    if source is None:
        return None
    elif source == 'binary':
        filename = get_product_filename(params, tod=tod)
        dtype = params.get('dtype', 'float32')
        angles = np.fromfile(filename, dtype=dtype)
    elif source == 'tod':
        # Yuck.
        from moby2.instruments import actpol
        angles = actpol.get_hwp_angles(params, tod=tod)
        trimmed = True
    elif source == 'dirfile':
        field_name = params.get('field', 'Hwp_Angle')
        filename = get_product_filename(params, tod=tod)
        df = moby2.util.DirfileManager(filename)
        info = df.get_frame_count(field_name)
        angles = df.load_channel(field_name, 0, info.n_samples)
    else:
        raise RuntimeError("HWP angle source='%s' unknown" % source)
    # Trim to match TOD?
    if not (trimmed or params.get('no_trimming') or (tod is None)):
        offset, count = tod.info.sample_index, tod.nsamps
        angles = angles[offset:offset+count]
    # Final conversions?
    angles = angles * params.get('rescale_to_degrees', 1.)
    return angles


def get_pointing_offset(params, tod=None, tod_info=None,
                        ctime=None, az=None, alt=None,
                        source_offset=False, model_cache=False):
    """
    Get focal plane shifts d_xi, d_eta, that should be added to the
    detector relative offsets prior to mapping.  That statement sets
    the sign convention.  Pass source_offset=True to flip the sign.

    E.g.:
        {'source': 'file',
         'filename': ...,
         'columns': [0,2,3],
         'rescale_degrees': 60.,
         }
    or
        {'source': 'model',
         'filename': ...,
         }
    You can also pass in a list of such blocks, to sum shifts.

    Note that passing a tod_info object as tod is often good enough.
    Yes, that interface is stupid.
    """
    # Allow combination of multiple sources (e.g. a model plus a
    # per-file shift) if the params come in as a list of param dicts.
    if isinstance(params, list):
        x, y = 0., 0.
        for p in params:
            dx, dy = get_pointing_offset(p, tod=tod, tod_info=tod_info,
                                         ctime=ctime, az=az,
                                         alt=alt, source_offset=source_offset,
                                         model_cache=model_cache)
            x, y = x + dx, y + dy
        return x, y

    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'],
                              tod_id=get_tod_id(tod=tod, tod_info=tod_info))
        return get_pointing_offset(params, tod=tod, tod_info=tod_info,
                                   ctime=ctime, az=az, alt=alt,
                                   source_offset=source_offset,
                                   model_cache=model_cache)

    filename = get_product_filename(params)
    if params.get('source') == 'file':
        # Convert tod to a basename?
        key = get_tod_id(tod=tod, tod_info=tod_info)
        offset_data = moby2.util.ascii.read_columns(
            filename, columns=params.get('columns',[0,1,2]))
        try:
            idx = list(offset_data[0]).index(key)
        except ValueError:
            return None
        offset = offset_data[1][idx], offset_data[2][idx]
        offset = [x*params.get('rescale_degrees',1.) for x in offset]
    elif params.get('source') == 'model':
        gm = None
        if isinstance(model_cache, dict):
            gm = model_cache.get(filename)
        elif model_cache is False:
            pass
        else:
            raise RuntimeError("model_cache should be either False or a dict.")
        if gm is None:
            ppar = moby2.util.MobyDict.from_file(filename) 
            model_version = ppar.get('model_version', 1)
            gm = moby2.pointing.GlobalModel()
            if gm.model_version != model_version:
                raise RuntimeError("Pointing shift model version mismatch.")
            gm.update(dict(list(zip(ppar['params'], ppar['values']))))
            if model_cache is not False:
                model_cache[filename] = gm
        if az is None:
            assert(tod is not None)   # Need a TOD, or else pass in az,alt.
            az = tod.az.mean()
        if alt is None:
            assert(tod is not None)   # Need a TOD, or else pass in az,alt.
            alt = tod.alt.mean()
        dx, dy = gm.get_pointing_shift(az, alt)
        offset = dx * 180/np.pi, dy * 180/np.pi
    elif params.get('source') == 'constant':
        offset = params['value']
    else:
        raise ValueError("could not decode source offset instructions: %s" %\
            str(params))
    if source_offset and offset is not None:
        offset = -offset[0], -offset[1]
    return offset

def get_weighted_modes(params, tod=None, match_det_uid=False):
    """
    Load time-ordered data modes with their per detector weights.
    """
    fmt = params['format']
    if fmt == 'moby_modes':
        modes_file = products.get_product_filename(params, tod=tod)
        weights = moby2.util.StructDB.from_fits_table(modes_file, index=1)
        modes = moby2.util.StructDB.from_fits_table(modes_file, index=2)
        det_uid = None
        if match_det_uid is True:
            det_uid = tod.det_uid
        elif hasattr(match_det_uid, '__getitem__'):
            det_uid = match_det_uid
        if det_uid is None:
            return modes, weights, None
        widx = weights.select_inner({'det_uid': tod.det_uid})
        return modes, weights[widx], (widx>=0)
    else:
        raise ValueError("I do not understand modes_format='%s'" % modes_format)

def get_cuts(params, tod=None, raw=True):
    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'], tod_id=get_tod_id(tod=tod))

    filename = get_product_filename(params, cls=moby2.TODCuts, tod=tod)
    format, source = params.get('format'), params.get('source')
    if format is None:
        format = source
    if format is None and any([filename.endswith(x) for x in ['.h5', '.hdf', '.hdf5']]):
        format = 'flags'

    if format in ['file', 'depot', None]:
        cuts = moby2.TODCuts.read_from_path(filename)
    elif format == 'flags':
        flags = get_tod_flags(params, tod=tod)
        cuts = flags.get_cuts(params['flag_name'])
    elif format == 'act_file':
        cuts = moby2.TODCuts.from_act_cuts_file(filename)
    elif format == 'actpol_file':
        cuts = moby2.TODCuts.from_actpol_cuts_file(filename)
    else:
        raise ValueError("could not decode cuts instructions: %s" % str(params))
    if raw:
        # Do not trim samples or select to match tod.det_uid
        return cuts
    # Extract det_uid...
    cuts = cuts.copy(det_uid=tod.det_uid)
    # Trim to match
    cuts = cuts.extract(sample_offset=tod.info.sample_index, nsamps=tod.nsamps)
    return cuts

def get_tod_flags(params, tod=None):
    tod_id = get_tod_id(tod=tod)
    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'], tod_id=tod_id)
    filename = get_product_filename(params, cls=moby2.tod.TODFlags, tod=tod)
    flags_archive = moby2.tod.TODFlagsArchive(filename)
    return flags_archive.get_item(tod_id)

def get_darkDets(params, tod):
    from moby2.analysis.tod_ana.pathologies import darkDets
    filename = get_product_filename(params, cls=darkDets, tod=tod)
    ddObj = darkDets.read_from_path(filename)
    return ddObj

def get_pathologies(params, tod):
    from moby2.analysis.tod_ana.pathologies import Pathologies
    if "paramFile" in params:
        cutParams = moby2.util.MobyDict.from_file(params["paramFile"])
        pathop = cutParams['pathologyParams']
    else: pathop = None
    filename = get_product_filename(params, cls=Pathologies, tod=tod)
    return Pathologies.read_from_path(filename, tod=tod, params=pathop)


def get_flatfield(params, tod_info=None):
    if '_execcfg' in params:
        tod_id = get_tod_id(tod_info=tod_info)
        params = exec_chooser(params['_execcfg'], tod_id=tod_id)

    source = params.get('source', 'column_file')
    columns = params.get('columns', [0,1])
    filename = get_product_filename(params, tod_info=tod_info)

    if is_ascii_columns(source):
        return moby2.detectors.FlatField.from_columns_file(
            filename, fields=['det_uid','cal'],
            columns=columns)
    if source == 'cal_dict':
        return moby2.detectors.FlatField.from_dict(filename)
    if source == 'fits_table':
        return moby2.detectors.FlatField.from_fits_table(filename, boolify=['stable'])

    raise ValueError("Unhandled flatfield source=%s" % source)
    
def get_detector_list(params, tod_info=None):
    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'], tod_id=get_tod_id(tod_info=tod_info))

    source = params.get('format', params.get('source'))

    # Handle easy ones first
    if source == 'list':
        return moby2.detectors.DetectorList(np.array(params['det_uid']))

    # I guess we need the filename
    filename = get_product_filename(params, tod_info=tod_info)

    if is_ascii_columns(source):
        return moby2.detectors.DetectorList.from_columns_file(
            filename, fields=['det_uid'],
            columns=params['columns'])
    if source == 'dict':
        return moby2.detectors.DetectorList.from_dict(filename)
    if source == 'fits':
        return moby2.detectors.DetectorList.from_fits_table(filename)
    if source == 'flatfield':
        db = moby2.util.StructDB.from_fits_table(filename)
        return moby2.detectors.DetectorList(db['det_uid'][db['fiducial']!=0])

    raise ValueError("Unhandled detector list source=%s" % source)


def get_iv_calibration(params, tod_info=None, det_uid=None,
                       instrument=None):
    instrument = _get_instrument(params, tod_info, instrument)
    season, array_name = _get_season_array_name(params, tod_info, instrument)

    iv_source = params.get('source', None)
    if iv_source in ['data', 'runfile']:
        if iv_source == 'data':
            runfile = None
            trace(3, 'Getting IV calibration from data')
        else:
            runfile = params.get('filename')
            trace(3, 'Using user provided IV block (%s)' % runfile)

        if instrument == 'mbac':
            return moby2.detectors.ACTIVCalibration.for_tod(tod_info, runfile=runfile)
        elif instrument == 'actpol':
            return moby2.detectors.IVCalibration.for_tod(tod_info, runfile=runfile)
        raise _UnhandledInstrument(instrument)

    raise ValueError("Unhandled iv_source=%s" % iv_source)
    

def get_calibration(params, tod_info=None, det_uid=None, tod=None):
    """
    Decode params to generate a calibration vector.  For templates you
    need to provide a tod, or tod_info.  To restrict to certain
    detectors, provide tod or pass det_uid.
    """
    if isinstance(params, dict) and '_execcfg' in params:
        params = exec_chooser(params['_execcfg'],
                              tod_id=get_tod_id(tod=tod, tod_info=tod_info))

    # Do our best to get a tod_info together.
    if tod_info is None:
        if hasattr(tod, 'info'):
            tod_info = tod.info
        elif isinstance(tod, basestring):
            tod_info = get_tod_info({'filename': tod})

    # We require det_uid, because this returns a Calibration object,
    # which is enumerated based on det_uid.
    if det_uid is None:
        if hasattr(tod, 'det_uid'):
            det_uid = tod.det_uid
        elif hasattr(tod_info, 'array_data'):
            det_uid = np.arange(tod_info.array_data.ndets)
    if det_uid is None:
        raise ValueError("det_uid could not be initialized from "\
            "tod or tod_info.")
    det_uid = np.asarray(det_uid)
    
    # Calibration object with cal.cal initialized to 1.
    cal = moby2.Calibration(det_uid=det_uid)

    # Translate params... accept list or dict.
    if isinstance(params, dict):
        step = params
    elif isinstance(params, list):
        step = {'type': 'cal_steps',
                'cal_steps': params,
                'name': 'list'}
    else:
        raise ValueError("Cannot decode calibration parameters: %s" % \
            repr(params))

    trace(3, 'Applying cal "%s"' % step.get('name', '<unnamed>'))

    if step.get('type', 'cal_steps') == 'cal_steps':
        # Compile several results.
        for step_params in step['cal_steps']:
            step_cal = get_calibration(step_params, tod_info=tod_info,
                                       det_uid=det_uid, tod=tod)
            cal.cal *= step_cal.cal
        return cal

    if step['type'] == 'constant':
        cal.cal *= step['value']
        return cal

    if step['type'] == 'iv':
        iv_cal = get_iv_calibration(step, tod_info=tod_info,
                                    det_uid=det_uid)
        ok, recal = iv_cal.get_property('dac_to_pW', det_uid)
        cal.cal *= recal
        cal.cal[~ok] = 0.
        return cal

    if step['type'] == 'flatfield':
        ff_cal = get_flatfield(step, tod_info=tod_info)
        ok, recal = ff_cal.get_property('cal', det_uid)
        cal.cal *= recal
        cal.cal[~ok] = 0.
        return cal

    if step['type'] == 'array_data_parameter':
        param_name = step['parameter']
        recal = tod_info.array_data[param_name][det_uid]
        cal.cal *= recal
        return cal

    if step['type'] == 'depot_cal':
        depot = get_depot(step.get('depot', {}))
        filename = get_product_filename(step, cls=moby2.Calibration, tod_info=tod_info)
        #cald = depot.read_object(moby2.Calibration,
        #                         tag=step['tag'], tod_info=tod_info)
        cald = moby2.Calibration.from_dict(filename)
        ok, recal = cald.get_property('cal', det_uid=det_uid) 
        recal[np.isnan(recal)+np.isinf(recal)] = 0.
        cal.cal *= recal
        cal.cal[~ok] = 0.
        return cal

    if step['type'] == 'cal_dict':
        # Load a calgc-style file?
        filename = get_product_filename(step,
                                        cls=moby2.Calibration,
                                        tod=tod,
                                        tod_info=tod_info)
        calgc_cal = moby2.detectors.CalGC.from_dict(filename)
        if len(calgc_cal.cal) > 0:
            ok, recal = calgc_cal.get_property('cal', det_uid)
            cal.cal *= recal
            cal.cal[~ok] = 0.
        else:
            cal.cal[:] = 0
        return cal

    if step['type'] == 'hdf':
        filename = get_product_filename(step,
                                        cls=moby2.Calibration,
                                        tod=tod,
                                        tod_info=tod_info)
        cal_arc = moby2.detectors.CalibrationArchive(filename)
        calg = cal_arc.get_item(tod_info.tod_id)
        ok, recal = calg.get_property('cal', det_uid=det_uid)
        cal.cal *= recal
        cal.cal[~ok] = 0.
        return cal

    if is_ascii_columns(step['type']):
        filename = get_product_filename(step, cls='Calibration',
                                        tod=tod,
                                        tod_info=tod_info)
        columns = step.get('columns', {'det_uid': 0, 'cal': 1})
        caldb = moby2.util.StructDB.from_column_file(
            filename, columns)
        idx = caldb.select_inner({'det_uid': det_uid})
        if 'cal' in columns:
            # If you don't set 'cal': <column>, it acts as a simple veto.
            cal.cal *= caldb['cal'][idx]
        cal.cal[idx<0] = 0. # Kill invalid matches.
        return cal

    if is_ascii_columns(step['type'], 'abs_cal_'): # e.g. 'abscal_column_file'
        if not 'filename' in step:
            if not 'structure' in step:
                step['structure'] = '{cls}/{tag}'
            if not 'cls' in step:
                step['cls'] = 'TODAbsCal'
        filename = get_product_filename(step, tod=tod,
                                        tod_info=tod_info)
        columns = step.get('columns', {'basename': 0, 'cal': 1})
        cal_data = moby2.util.StructDB.from_column_file(filename, columns)
        index = cal_data.select_inner({'basename': [tod_info.basename]})[0]
        if index < 0:
            raise ValueError("Basename %s not found in abscal file %s "\
                "(column %i)" % \
                (tod_info.basename, filename, columns['basename']))
        cal.cal *= cal_data['cal'][index]
        return cal

    if step['type'] == 'per_tod_hdf':
        per_cal = get_pertod_calibration_hdf(params, tod_id=tod_info.tod_id)
        ok, per_cal = per_cal.get_property('cal', det_uid=cal.det_uid)
        cal.cal *= (ok*per_cal)
        return cal

    if step['type'] == 'referenced_cal':
        ref_cal = get_referenced_calibration(params, tod_info=tod_info)
        mask, index = ref_cal.get_index(det_uid)
        cal.cal[~mask] = 0
        cal.cal[mask] *= ref_cal.cal[index[mask]]
        return cal

    if step['type'] == 'remove_readout_filter_gain':
        cal.cal /= tod_info.runfile.ReadoutFilter().gain()
        return cal

    if step['type'] == 'MBAC_abscal':
        calp = get_MBAC_CalParams(tod_info, step['params'])
        cal.cal = calp["cal0"] * calp["cal"] * 1e6
        return cal

    if step['type'] == 'MBAC_transmission':
        transmission = get_MBAC_Transmission(tod_info, step['params'],
                                             pwv=step.get("PWV"))
        cal.cal *= transmission
        return cal

    raise ValueError('Unhandled calibration type "%s"' % step['type'])


def get_pertod_calibration_hdf(params, tod_id=None):
    """Load per-TOD, per-band calibration numbers from an HDF archive.
    params should include the filename (somehow), and the 'dataset'
    name, which defaults to 'abscal'.  The Simple Data Table should
    have columns 'tod_id', 'cal', and possibly 'band_id'.
    """
    import h5py
    if '_execcfg' in params:
        params = exec_chooser(params['_execcfg'], tod_id=tod_id)
    adata = get_array_data_new({'tod_id': tod_id})
    try:
        fcode = adata['code']
    except:
        # remove me when no longer needed.
        fcode = np.array(['f%03i' % n for n in adata['nom_freq']])
    filename = get_product_filename(params)
    dataset_name = params.get('dataset', 'abscal')
    hf = h5py.File(filename, 'r')
    data = np.asarray(hf[dataset_name])
    s = (data['tod_id'] == tod_id)
    idx = s.nonzero()[0]
    cal = moby2.Calibration(det_uid=adata['det_uid'])
    cal.cal[:] = 0.
    for i in idx:
        row = data[i]
        s = np.ones(len(cal.cal), bool)
        if 'band_id' in data.dtype.names:
            s = s * (row['band_id'] == fcode)
            cal.cal[s] = row['cal']
    return cal


def get_referenced_calibration(params, tod_info=None):
    """
    Get calibration with intermediate lookup.  For example, IV and
    biasstep analysis results may be tagged with their acquisition
    times rather than the TOD basename.

    The way we support this is as follows:

    * The caller points to file with "assignments" that maps timestamp
      to some useful "tag".  (The tag could be the ctime of the IV
      file, for example.)
    
    * The caller also points to a "library" of calibration data, and
      the "tag" is used to load one of them.

    E.g., params could be something like:

    params = {
    'name': 'bias_step',
    'type': 'referenced_cal',
    'assignments': {'depot': {'label': 'actpol_shared'},
                     'structure': './BiasStepTimes/intervals_{season}_{array}_{tag}.txt',
                     'tag': '150131',
                     'columns': [0,1,4]},
     'library': {'depot': {'label': 'actpol_egrace1'},
                 'structure': 'calibration/bias_step/{tag}.cal',
                 'type': 'cal_dict',
                 }
    }
    """
    assign_file = get_product_filename(params['assignments'],
                                       tod_info=tod_info)
    columns = params['assignments'].get('columns', [0,1,2])
    t0, t1, tags = moby2.util.ascii.read_columns(assign_file, columns)
    s = (t0 <= tod_info.ctime) * (tod_info.ctime < t1)
    if not s.any():
        raise ValueError("No match for timestamp %i in reference file %s" %\
            (tod_info.ctime, assign_file))
    tag = tags[s.nonzero()[0][0]]
    # Now load the calibration result.  Yes, this is dangerously recursive.
    lib_params = params['library'].copy()
    lib_params['tag'] = tag
    return get_calibration(lib_params, tod_info=tod_info)

def get_MBAC_Transmission(tod_info, paramFile, pwv = None):
    """
    This function was inherited from moby(1) to provide backcompatibility to mbac calibration
    Computes the atmphere transmission, providing also the calibration parameters cal0 and calc
    @return (transmission, cal0, calc) 
    """
    params =  get_MBAC_CalParams(tod_info, paramFile)
    dt = 600
    tauWet = params['tauWet']
    alt0 = params['alt0']*np.pi/180
    pwv0 = params['pwv0']
    if pwv == None:
        ctime = tod_info.ctime
        r = moby2.aux_data.apex.Radiometer()
        pwv = r.get_average(ctime, ctime+dt)[1]
        if (pwv == -1) or np.isnan(pwv):
            pwv = pwv0
            psLib.trace("moby",3,"Using season average PWV: %f" % pwv)
    transmission = np.exp((pwv - pwv0)*tauWet/np.sin(alt0))  # P_det = transmission * P_sky
    return transmission


def get_MBAC_CalParams(tod_info, filename):
    """
    Legacy code from moby(1) to give back compatibility to MBAC
    Retrieve all calibration parameters from a simpleCal parameter file.
    """
    params = {}

    data = []
    f = open(filename)
    for l in f:
        if l[0] != "#":
            data.append(l.split("\n")[0].split())
    f.close()
    data = np.array(data)

    if tod_info.array == "mbac145": ar = "AR1"
    elif tod_info.array == "mbac215": ar = "AR2"
    else: ar = "AR3"
    arSel = np.array(data[:,0]) == ar
    iniSel = np.array(data[:,2], dtype = int) < tod_info.ctime
    endSel = np.array(data[:,3], dtype = int) > tod_info.ctime
    sel = arSel*iniSel*endSel
    if np.any(sel):
        params["cal0"] = float(data[sel,4][0])
        params["pwv0"] = float(data[sel,5][0])
        params["alt0"] = float(data[sel,6][0])
        params["tauWet"] = float(data[sel,11][0])
        params["cal"] = float(data[sel,13][0])
    else:
        raise RuntimeError("No calibration values found for tod %s"%tod_info.basename)
    return params


def get_obs_catalog(params=None):
    if params is None:
        params = moby2.util.get_user_config().get('obs_catalog')
    if params is None:
        raise RuntimeError("To use get_obs_catalog without arguments, set up "\
            "obs_catalog = {...} in your .moby2 file.")
    filename = get_product_filename(params)
    ftype = params.get('type', 'fits')
    if ftype == 'fits':
        return moby2.util.StructDB.from_fits_table(filename)
    if ftype == 'column_file':
        return moby2.util.StructDB.from_column_file(
            filename, params['columns'])
    raise ValueError("Could not decode catalog type '%s'" % ftype)

    
# Someday we should decide to use column_file or columns_file... and not both.
ASCII_COLUMNS_NAMES = ['column_file', 'columns_file', 'ascii', 'ascii_columns']
def is_ascii_columns(format_str, prefix=''):
    return format_str in [prefix+x for x in ASCII_COLUMNS_NAMES]
