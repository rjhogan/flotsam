#!/bin/bash
set -x
FLOTSAM=../bin/flotsam
CFGFILE=test_phase_function.cfg

for FILE in phase_function*.txt CAMS_benchmark/phase_function*0470nm.txt CAMS_benchmark/phase_function*1024nm.txt
do
    NEWFILE=$(basename $FILE | sed -e 's/phase_function_/pf_/' \
	-e 's/\.txt$/_out\.dat/')
    echo "$FILE -> $NEWFILE"
    $FLOTSAM phase_function_file=$FILE $CFGFILE > $NEWFILE
done
