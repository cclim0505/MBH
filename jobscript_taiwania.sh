#!/bin/bash
#PBS -N master_test
#PBS -q ctest
#PBS -l select=1:ncpus=20:mpiprocs=20
#PBS -l walltime=00:02:00
#PBS -P MST108285
#PBS -o PBS.log
#PBS -e PBS.err
####PBS -r n
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=2
module load intel/2018_u1
module load mpi/openmpi-3.0.0/intel2018u1
mpirun ./run.out
