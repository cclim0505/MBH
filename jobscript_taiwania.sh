#!/bin/bash
#PBS -N MBH01
#PBS -q ctest
#PBS -l select=1:ncpus=4:mpiprocs=1
#PBS -l walltime=00:30:00
#PBS -P MST108099
#PBS -o PBS.log
#PBS -e PBS.err
####PBS -r n
cd $PBS_O_WORKDIR
module load intel/2018_u1
module load mpi/openmpi-3.0.0/intel2018u1
mpirun ./run.out
