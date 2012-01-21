#!/bin/bash
#
# how to use
# change to directory with the files you want to run
# and enter:
# $ qsub -V -t 0-9 -d . rtfluka.sh
#
#PBS -N FLUKA_JOB
#
start="$PBS_ARRAYID"
let stop="$start+1"
stop_pad=`printf "%03i\n" $stop`
#
# Init new random number sequence for each calculation.
# This may be a poor fix.
cp $FLUPRO/random.dat ranexample$stop_pad
sed '/RANDOMIZE        1.0/c\RANDOMIZE        1.0 '"${RANDOM}"'.0 \' example.inp &gt; _example.inp
mv _example.inp example.inp
$FLUPRO/flutil/rfluka -N$start -M$stop example -e flukadpm3
