#!/bin/bash
#PBS -N Serial_test
#PBS -q serial
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=00:10:00
#PBS -P MST108285
#PBS -o PBS.log
#PBS -e PBS.err
####PBS -r n
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
module load intel/2018_u1
module load mpi/openmpi-3.0.0/intel2018u1
mpirun ./run.out
