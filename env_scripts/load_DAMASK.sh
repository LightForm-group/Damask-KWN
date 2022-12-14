#!/bin/bash

module load tools/env/proxy
module load mpi/intel-18.0/openmpi/4.0.1
module load tools/gcc/cmake/3.13.2
module load apps/binapps/anaconda3/2019.07

export PETSC_DIR=/mnt/eps01-rds/jf01-home01/shared/petsc-3.12.2
export PETSC_ARCH=mkl-opt
export DAMASK_ROOT=~/software/damask-kwn
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1

source $DAMASK_ROOT/env/DAMASK.sh
PATH=$PETSC_DIR/$PETSC_ARCH/bin:$PATH
LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
