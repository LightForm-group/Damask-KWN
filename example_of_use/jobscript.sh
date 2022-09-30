#!/bin/bash --login
#$ -cwd
#$ -pe smp.pe 8

source ~/load_DAMASK.sh 
mpirun -n $NSLOTS DAMASK_grid -l load.yaml -g geom.vtr

