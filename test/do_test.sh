#!/bin/bash
DEST=.
PRE="flotsam0522"
PHASE="MieR10"
CFG_NADIR="nadir_droplet10micron.cfg"
CFG_POLAR="polar_droplet10micron.cfg"
FLOTSAM=../bin/flotsam

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
$FLOTSAM $@ $CFG_NADIR  > ${PREFIX}/${PRE}_${RUN_ID}_${PHASE}_out.dat
$FLOTSAM $@ $CFG_POLAR  > ${PREFIX}/${PRE}_${RUN_ID}_${PHASE}_polar_out.dat
