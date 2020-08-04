from __future__  import print_function

import moby2
from moby2.analysis import fp_fit
from moby2.analysis import beam_ana

import glob, os, sys, shutil
import numpy as np
DEG = np.pi/180

default_params = {
    'tod_list': 'list.txt',
    'trial_tod': None,
    'map_cfg_tplate': 'map_realign_auto.in.template',
    'map_cfg': '__map_realign_auto.cfg',        #tempfile
    'ptg_file': '__ptg_ptemplate.cfg',  #tempfile
    'map_dir_pattern': 'maps_realign_stage{stage}/',
    'map_name_pattern': '{tod_id}/source_{fcode}_I.fits',
    'shifts_file': 'realign_shifts.txt',
}

def main(args_in):
    """How to use this code -- set up the required inputs, i.e.:

    - map_cfg_tplate: filename for a mapping input file.  The "output" and
      "pointing" blocks will be overwritten!  Best to comment them out
      entirely.

    - base_model: a pointing configuration block to use for the
      "pointing_shifts" part of the pointing configuration.

    - input_template: the detector offsets file to use as the base -- this
      code will generate new offsets based on this one.


    Write the start-point offsets.  (If you know your planet map has the
    signal at position (x,y), then set the offset to (-x, -y)):

      # set shifts for stage 0 to (0,0)
      moby2 fp_realign set-shifts 0 0 0

    Do a trial map (stage 0):

      moby2 fp_realign map 0 --trial

    Measure peak position in stage 0 maps:

      moby2 fp_realign measure 0

    If you are satisfied with that, then great.  If not, update the
    (local) shifts database so as to eliminate the offset:

      moby2 fp_realign measure 0 --set-next-shift --target 0 0

    You can use --target to say where you want the planet to be in the
    map... sometimes 0,0 makes sense but sometimes you might want it to
    appear at the same position as it used to / as it does in other arrays.

    After "set-next-shift", you can run maps for the next stage to verify it worked:

      moby2 fp_realign map 1 --trial
      moby2 fp_realign measure 1

    If satisfied, redo everything without --trial.  To clear the local
    shifts database, run

      moby2 fp_realign reset-shifts 0

    When happy with the maps, write the final shifted template:

      moby2 fp_realign write_template 1 template_out.txt

    """

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config-file', default='fp_realign.cfg')
    sp = parser.add_subparsers(dest='action')
    sap = sp.add_parser('map')
    sap.add_argument('stage', type=int)
    sap.add_argument('--trial', action='store_true')
    sap = sp.add_parser('measure')
    sap.add_argument('stage', type=int)
    sap.add_argument('--trial', action='store_true')
    sap.add_argument('--set-next-shift', action='store_true')
    sap.add_argument('--target', nargs=2, type=float, default=(0,0))
    sap.add_argument('--force', action='store_true')
    sap = sp.add_parser('set-shifts')
    sap.add_argument('stage', type=int)
    sap.add_argument('dx', type=float)
    sap.add_argument('dy', type=float)
    sap.add_argument('--force', action='store_true')
    sap = sp.add_parser('reset-shifts')
    sap.add_argument('stage', type=int)
    sap = sp.add_parser('write-template')
    sap.add_argument('stage', type=int)
    sap.add_argument('output_filename')

    args = parser.parse_args(args_in)
    print(args)

    params = moby2.util.MobyDict.from_file(args.config_file)
    params = params['realign_params']

    if args.action in ['map', 'measure']:
        trial_tod = params.get('trial_tod')
        if not args.trial or trial_tod is None:
            targets = [line.split()[0]  for line in open(params['tod_list'])]
        if args.trial:
            if trial_tod is None:
                targets = targets[len(targets)//2:][:1]
            else:
                targets = [trial_tod]

    if args.action == 'set-shifts':
        shifts = read_shifts(params['shifts_file'])
        if len(shifts) > 0 and args.stage <= max(shifts.keys()):
            if args.force:
                shifts = {k:v for k,v in shits.items()
                          if k < args.stage}
            else:
                parser.error('  Trying to update an old offset.')
        shifts[args.stage] = (args.dx, args.dy)
        write_shifts(shifts, params['shifts_file'])

    elif args.action == 'reset-shifts':
        shifts = read_shifts(params['shifts_file'])
        new_shifts = {k: v for k, v in shifts.items() if k < args.stage}
        write_shifts(new_shifts, params['shifts_file'])

    elif args.action == 'map':
        # What shift to use?
        shifts = read_shifts(params['shifts_file'])
        if args.stage not in shifts:
            parser.error('  No shifts set for stage %i' % args.stage)
        dx, dy = shifts[args.stage]
        print('Updating %s ...' % params['ptg_file'])
        write_shifted_ptg(params['input_template'], (dx, dy), params['ptg_file'])

        ptext = get_pointing_txt(params['map_dir_pattern'].format(stage=args.stage),
                                 params['ptg_file'], params['base_model'])

        # Write new config.
        print('Updating %s ...' % params['map_cfg'])
        with open(params['map_cfg'], 'w') as fout:
            for line in open(params['map_cfg_tplate']):
                fout.write(line)
            fout.write('# and also pointing ...\n')
            fout.write(ptext)

        for tod in targets:
            print(tod)
            exists = False
            for fcode in params['fcodes']:
                info = {'stage': args.stage,
                        'tod_id': tod,
                        'fcode': fcode}
                map_file = (params['map_dir_pattern'] + params['map_name_pattern']).format(**info)
                exists = exists or os.path.exists(map_file)
            if exists:
                print('Skipping %s, exists' % tod)
            else:
                os.system('planetMap %s %s' % (params['map_cfg'], tod))

    elif args.action == 'measure':
        print('Source positions (x, y), in arcminutes:')
        rows = {'all': []}
        for fcode in params['fcodes']:
            print('  fcode=%s' % fcode)
            rows[fcode] = []
            for tod in targets:
                info = {'stage': args.stage,
                        'tod_id': tod,
                        'fcode': fcode}
                map_file = (params['map_dir_pattern'] + params['map_name_pattern']).format(**info)
                if not os.path.exists(map_file):
                    print('    (Not found: %s)' % map_file)
                    continue
                bo = beam_ana.BeamObs(mapfile=map_file)
                bo.fit()
                dx, dy = [bo['gaussian'][k] for k in ['dx', 'dy']]
                print('    %s    %8.3f  %8.3f' % (map_file, dx*60, dy*60))
                rows[fcode].append((dx, dy))
                rows['all'].append((dx, dy))
            if len(rows[fcode]):
                dx0, dy0 = np.mean(rows[fcode], axis=0)
                sx0, sy0 = np.std(rows[fcode], axis=0)
                print('    %s    %8.3f  %8.3f' % ('average', dx0*60, dy0*60))
                print('    %s    %8.3f  %8.3f' % ('stdev  ', sx0*60, sy0*60))
            print()
        print()
        if len(rows['all']):
            print('Super average')
            dx0, dy0 = np.mean(rows['all'], axis=0)
            sx0, sy0 = np.std(rows['all'], axis=0)
            print('    %s    %8.3f  %8.3f' % ('average', dx0*60, dy0*60))
            print('    %s    %8.3f  %8.3f' % ('stdev  ', sx0*60, sy0*60))
        ## Option to set the shifts for next stage.
        if args.set_next_shift:
            print()
            print('Setting shifts for next stage (%i)' % (args.stage+1) +
                  'to target offset (%.4f,%.4f) argmin.' % tuple(args.target))
            shifts = read_shifts(params['shifts_file'])
            if args.stage+1 in shifts:
                if args.force:
                    shifts = {k:v for k,v in shits.items()
                              if k <= args.stage}
                else:
                    parser.error('  We already have a shift for stage %i' % (args.stage+1))
            # We adjust the template to produce planet position target
            # instead of planet position (dx0, dy0).
            shifts[args.stage+1] = [(a-b+c) for a, b, c in zip(
                shifts[args.stage],
                (dx0*60, dy0*60),
                args.target)]
            write_shifts(shifts, params['shifts_file'])

    elif args.action == 'write-template':
        # What shift to use?
        shifts = read_shifts(params['shifts_file'])
        if args.stage not in shifts:
            parser.error('  No shifts set for stage %i' % args.stage)
        dx, dy = shifts[args.stage]
        print('Writing %s based on shifts %.4f %.4f arcmin' % (
            args.output_filename, dx, dy))
        write_shifted_ptg(params['input_template'], (dx, dy), args.output_filename)


# Support functions...

def read_shifts(shifts_file):
    shifts = {} # by stage.
    for line in open(shifts_file):
        casts = (int, float, float)
        stage, dx, dy = [c(w) for c, w in zip(casts, line.split())]
        shifts[stage] = (dx, dy)
    return shifts

def write_shifts(shifts, shifts_file):
    with open('shifts.tmp', 'w') as fout:
        for k, v in sorted(shifts.items()):
            fout.write('%i  %.8e %.8e\n' % (k, v[0], v[1]))
    shutil.move('shifts.tmp', shifts_file)

def get_pointing_txt(map_dir, ptg_template, base_model):
    text = """

### Generated dynamically:
output = {
    # Output for plots and stats
    'prefix': '%s/{tod_id}/',

    'progress_file': 'map.dict',
    #index of channel to plot, or None to suppress.
    'tod_plot_det': 25,
}
pointing = {
    'detector_offsets': {
        'filename': '%s',
        'columns': [('det_uid',0),('x',2),('y',4),('mask',1)],
        'match_bad': 0},
    'pointing_shifts': %s,
}
""" % (map_dir, ptg_template, repr(base_model))
    return text

def write_shifted_ptg(base, shifts_arcmin, output):
    fp = fp_fit.FPFitFile.from_columns_file(base)
    s = fp.ok
    fp.x0[s] += shifts_arcmin[0]/60 * DEG
    fp.y0[s] += shifts_arcmin[1]/60 * DEG
    fp.write(output)
