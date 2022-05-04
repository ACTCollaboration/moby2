.. -*- mode: rst; mode: auto-fill -*-

================
SO Compatibility
================

moby2so
=======

This script is used to index ACT metadata for consumption with
sotodlib context system.

Here's a set of commands that has been used to bundle the complete s19
TOD and metadata set::

    # Write TOD list.
    moby2 obs_catalog | grep s19 > s19_tods.txt

    # Write obsdb.
    moby2so -i s19_tods.txt obsdb

    # Write obsfiledb (this takes a while)
    moby2so -i s19_tods.txt obsfiledb

    # Write detdb
    moby2so -i s19_tods.txt detdb

    # Write abscal
    moby2so -i s19_tods.txt abscal \
            /home/mfh2/depots/actpol_shared//TODAbsCal/abscal_211021.h5

    # Write todofs
    moby2so -i s19_tods.txt pointofs \
            /home/mfh2/depots/actpol_shared//TODOffsets/tod_offsets_200514.txt \
            pointofs_200514.h5
    
    # Write focalplane
    moby2so focalplane focalplane_201016.yaml focalplane.h5

    # Cuts & relcal (this takes a while)
    moby2so cuts_release release_20210616_s19.txt
    
    # Timeconst
    moby2so -i s19_tods.txt timeconst \
            /projects/ACT/mhasse/depots/actpol_mhasse1/TimeConstants/timeconst_201112/s19

    # Write context.yaml
    moby2so -i s19_tods.txt context
