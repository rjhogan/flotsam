#!/bin/bash
#set -x
PHASE="MieR10"
CFG_NADIR="nadir_droplet10micron.cfg"
FLOTSAM=../bin/flotsam

for LAYER_COUNT in 1 2 3 4 5 7 10 14 20 30 40 50 70 100 140 200 300 400 500 700 1000
do
    $FLOTSAM layer_count=$LAYER_COUNT repeat=10 mode=direct  $@ $CFG_NADIR > /dev/null 2> tmp_direct.txt
    $FLOTSAM layer_count=$LAYER_COUNT repeat=1 mode=jacobian $@ $CFG_NADIR > /dev/null 2> tmp_jacobian.txt
    echo $LAYER_COUNT \
	`grep 'time per reflectance calc' tmp_direct.txt | cut -c 43- | awk '{print $1}'` \
	`grep 'time per reflectance calc' tmp_jacobian.txt | cut -c 43- | awk '{print $1}'`
done

