#
# STM makefile
# 
#    Machine name: Monch (eejit.cscs.ch)
#    Operating system: Linux 
#    CPU Architecture: Intel 10-core Ivy Bridge EP E5-2660v2
#    Intel Compiler: icc
#

# -- Required modules to be loaded --
# module load gcc	
# module load intel/14.0.1
# module load openmpi/icc-14.0.1
# module load hdf5/1.8.12-intel14.0.1
# module load mkl/11.1.1.106-intel14.0.1
# module load lapack
#

# -- Execution sequence --
#
# source modload_eejit.sh
# make 
#
# -- LAB -- if setup<10 --:
# python make_input_lab.py
# sbatch run_lab_eejit.sbatch
#
# -- LS -- if setup>9 --:
# sbatch run_ls_eejit.sbatch
#

# ----------------- USER SHOULD VERIFY SETTINGS BELOW ------------------ #

#
# Select model setup options: 
# LAB (e.g., c.o. mass and momentum; Newtonian viscosity): 
# 0=lab (van Dinther et al., JGR, 2013a), 
# 1=fast setup for sequential data assimilation
# NATURE (e.g., c.o. mass, momentum and energy; P,T,s-dep. viscosity), i.e., large-scale model: 
# 10=nature-like dynw (van Dinther et al., JGR, 2013b), 
# 11=nature-like "offev", incl. off-megathrust events (van Dinther et al., 2014, GRL)
STM_SETUP=11

# Select INItial PRN (for nature-like large scale model):
# 1 = start from t=0 (geometry as defined in init.t3c, use mode_ls_to.t3c),
# 0 = start from prn defined at line in mode_ls_stm.t3c defined by file.t3c
STM_LS_T0 = 1

#
# If the STM code you compile was directly obtained via git, set this value to y
# If HAVE_GIT=y, git version information will be logged in stdout 
HAVE_GIT=n


# ----------------- USER COULD CHANGE SETTINGS BELOW ------------------ #

# Select the intel compiler
PCC = icc -I/trinity/opt/apps/gtecton/hdf5fromSource/hdf5-1.10.4/src

# Define compilation flags
# If need to debug, add -g
HDF5_FLAGS = -DH5_1_8

CFLAGS = -std=c99 -O2 -qopenmp -mcmodel=medium -DSTM_SETUP_DEF=$(STM_SETUP) -DSTM_INIPRN_DEF=$(STM_LS_T0)
CFLAGS += $(HDF5_FLAGS)

# <MKL - REQUIRED DEPENDENCY> Specify name of the MKL libraries - the module system takes care of the paths
MKL_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread 

# <HDF5 - REQUIRED DEPENDENCY>  Specify name of the HDF5 library - the module system takes care of the paths
HDF5_LIBS = -L/trinity/opt/apps/gtecton/hdf5fromSource/hdf5-1.10.4/src/.libs -lhdf5


# ----------------- USER SHOULD NOT EDIT ANYTHING BELOW THIS LINE ------------------ #

# If git info is set to available, append it into the CFLAGS variable
ifeq ($(HAVE_GIT),y)
	CFLAGS += -DHAVE_GIT
endif

LIBS = $(MKL_LIBS) $(HDF5_LIBS)

SRC = 	extrainputlabstm.c \
	gausstm.c \
	hdf5stm.c \
	headstm.c \
	heatstm.c \
	loadstm.c \
	markstm.c \
	movestm.c \
	speclabstm.c \
	speclsstm.c \
	specsdastm.c \


all: releaseinfo i2stm in2stm

releaseinfo:
ifeq ($(HAVE_GIT),y)
	-@echo "Generating git info"
	-@python stm-generate-git-info.py
endif

i2stm: releaseinfo $(SRC) i2stm.o
	$(PCC) $(CFLAGS) -o i2stm i2stm.o $(LIBS)

in2stm: releaseinfo $(SRC) in2stm.o
	$(PCC) $(CFLAGS) -o in2stm in2stm.o $(LIBS)

%.o: %.c $(SRC)
	-@echo "---- Compiling $*.c -----"
	$(PCC) $(CFLAGS) -c $*.c $(INC)

clean:
	rm -rf i2stm i2stm.o
	rm -rf in2stm in2stm.o



