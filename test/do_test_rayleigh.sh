#!/bin/bash
DEST=.
PRE="flotsam0523"
PHASE="Ray"
#PHASE="Iso"
CFG_NADIR="nadir_droplet10micron.cfg"
CFG_POLAR="polar_rayleigh.cfg"
FLOTSAM=../bin/flotsam
ARG="phase_function_file=phase_function_rayleigh.txt"
#ARG="phase_function_file=phase_function_isotropic.txt"
PREFIX="."
PREFIX="/scratch/rd/parr/flotsam"

if [ "$#" = "0" ]
then
    echo "Usage: $0 <run_id>"
    echo "Files to be written:"
    echo "  $DEST/${PRE}_<run_id>_${PHASE}.dat"
    echo "  $DEST/${PRE}_<run_id>_${PHASE}_polar.dat"
    exit
fi

RUN_ID=$1
shift
set -x
$FLOTSAM $ARG $@ $CFG_NADIR  > ${PREFIX}/${PRE}_${RUN_ID}_${PHASE}_out.dat
#$FLOTSAM $ARG $@ $CFG_POLAR  > ${PREFIX}/${PRE}_${RUN_ID}_${PHASE}_polar_out.dat
