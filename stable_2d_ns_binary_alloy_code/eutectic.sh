#!/bin/bash
#PBS -N apaar_mpi
#PBS -l nodes=4:ppn=32
#cd ~/test_ramanujam/
cd /home/apaar/2d_binary_aniso_mpi/stable_code
NPROCS=`wc -l < $PBS_NODEFILE`
HOSTS=`cat $PBS_NODEFILE | uniq | tr '\n' "," | sed 's|,$||'`
mpirun -np 128 --host $HOSTS ./a.out > output1.log
#mpirun -np 11 grandpot_mp_ternary.out 0 > output1.log
#mpirun -np 15 grandpot_mp.out 0 > output.log
