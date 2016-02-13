#!/bin/bash
let My=50
mkdir dispersion_1
for i in `seq 51 70`; do
  echo $i
  # cp file_processing.m datafiles$i
  cd datafiles$i
  cp $My.dat ~/phasefield_lbm/dispersion/code/dispersion_1
  # octave -q file_processing.m $My
  # octave --silent --eval "file_processing(\"$My\")"
  cd ..
  let My=My+5
  # COUNTER=100000
  # while [  $COUNTER -lt 401000 ]; do
  #    echo The counter is $COUNTER
  #    let COUNTER=COUNTER+1000
  # done
done
