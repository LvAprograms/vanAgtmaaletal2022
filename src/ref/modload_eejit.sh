#!/bin/sh

module purge
module load opt/all
module load userspace/custom
module load intel-compiler/64/2018.3.222
module load intel-mkl/64/2018.3.222
module load intel-mpi/64/2018.3.222
module load hdf5/1.10.4
#module load Ruby/2.6.2


# - load modules -
#module load gcc	
#module load intel/14.0.1
#module load openmpi/icc-14.0.1
#module load hdf5/1.8.12-intel14.0.1
#module load mkl/11.1.1.106-intel14.0.1
#module load lapack
#module load ruby
#module load matlab 
#module load ddt/5.1

# - export environment variables -
export I2ELVIS_INTEL="-lpthread -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"
export I2ELVIS_HDF5="-L/trinity/opt/apps/gtecton/hdf5fromSource/hdf5-1.10.4/src/.libs -lhdf5"
export I2ELVIS_OPT="-mcmodel=medium -openmp"

export OMP_NUM_THREADS=6

