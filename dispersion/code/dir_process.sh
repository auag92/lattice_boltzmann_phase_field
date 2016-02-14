#!/bin/bash
let My=50
mkdir dispersion_2
for i in `seq 51 70`; do
  echo $i
  cp file_processing.m datafiles$i
  cd datafiles$i
  octave -q file_processing.m $My
  octave --silent --eval "file_processing(\"$My\")"
  cp $My.dat ~/phasefield_lbm/dispersion/code/dispersion_2
  cd ..
  let My=My+5
  # COUNTER=100000
  # while [  $COUNTER -lt 401000 ]; do
  #    echo The counter is $COUNTER
  #    let COUNTER=COUNTER+1000
  # done
done

# #!/bin/bash
# COUNTER=0
# while [  $COUNTER -lt 10 ]; do
#     echo The counter is $COUNTER
#     let COUNTER=COUNTER+1
# done
# #!/bin/bash
# COUNTER=20
# until [  $COUNTER -lt 10 ]; do
#     echo COUNTER $COUNTER
#     let COUNTER-=1
# done
# #!/bin/bash
# for i in $( ls ); do
#     echo item: $i
# done
