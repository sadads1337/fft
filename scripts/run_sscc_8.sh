#!/bin/bash

module load intel/2019/compilers
module load intel/2019/lib/mkl
module load intel/2019/lib/tbb
module load intel/2019/mpi

#SBATCH --nodes=1
#SBATCH --threads-per-core=1
#SBATCH --time=30:0
#SBATCH --job-name=fft_makarov_8
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --ntasks-per-node=1
#SBATCH --output="fft_makarov_test.txt"

/home/mon.nsu/makarov_i_o/develop/fft/build/projects/Benchmark/Benchmark 8 5.115 65536 1024 7.62939e-06 32
