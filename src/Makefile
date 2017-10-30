MACHINE=magus
#MACHINE=clarkp

CC=icc
CFLAGS=-O2 -ipo
#FF=ifort
FFLAGS=

FORT=mpif90
OPTS= -g -O0 -heap-arrays #-openmp -CB -r8

SRC += globals.f90 read_field.f90 tools.f90 slice.f90 statistics.f90 main.f90
MOD += globals.f90 read_field.f90 tools.f90 slice.f90 statistics.f90
MODO = $(addsuffix .o,$(basename $(MOD)))
OBJMOD = $(addsuffix .mod,$(basename $(MOD)))
OBJO  = $(addsuffix .o,$(basename $(SRC)))
EXEC = analyze

EXEC_DIR = ../

LCL_LIB_DIR =/workspace/jmelvin


ifeq ($(MACHINE),athena)

  MKL_LIB = -L$(LAPACK_LIBS) -llapack
  #FFTW=$(LCL_LIB_DIR)/fftw-3.3.3/
  FFTW=$(PFFTW_DIR)/
  FFTW_LIB=-L$(FFTW)/lib -lfftw3
  FFTW_INC=-I$(FFTW)/include

  HDF_LIB = -L$(PHDF5_LIB) -lhdf5_fortran -lhdf5
  HDF_INC = -I$(PHDF5_INC)

#  ESIO_DIR=$(LCL_LIB_DIR)/MOD_LCL/ESIO
  ESIO_DIR=$(LCL_LIB_DIR)/ESIO-release-0.1.9/esio
  ESIO_LIB=-L$(ESIO_DIR)/lib -lesiof -lesio -lz -lm  -Wl,-rpath,$(ESIO_DIR)/lib
  ESIO_INC=-I$(ESIO_DIR)/include/esio -I$(ESIO_DIR)/lib

endif

ifeq ($(MACHINE),lonestar)
  MKL_INC = -I$(TACC_MKL_INC)
  MKL_LIB = -L$(TACC_MKL_LIB) -mkl -lmkl_lapack -lguide -lpthread
  FFTW_LIB = -L$(TACC_FFTW33_LIB)  -lfftw3 -lm #-lfftw3f -lfftw3_omp
  FFTW_INC = -I$(TACC_FFTW33_INC)
  HDF_LIB = -L$(TACC_HDF5_LIB) -lhdf5
  HDF_INC =  -I$(TACC_HDF5_INC)
endif

ifeq ($(MACHINE),CLARK)
#  ulimit -s unlimited  
  MKL_LIB=-L/usr/lib64 -L/usr/lib64/atlas -llapack
  FFTW_LIB=-L$(PFFTW_LIB) -lfftw3_mpi -lfftw3 -lm
  FFTW_INC=-I$(PFFTW_INC)
  HDF_LIB = -L$(PHDF5_LIB) -l:libhdf5_fortran.so.10 -l:libhdf5.so.10
  HDF_INC = -I$(PHDF5_INC)
  ESIO_LIB=-L$(ESIO_DIR)/lib -lesiof -lesio -lz -lm  -Wl,-rpath,$(ESIO_DIR)/lib
  ESIO_INC=-I$(ESIO_DIR)/include/esio -I$(ESIO_DIR)/lib
endif

ifeq ($(MACHINE),magus)
  MKL_LIB = -L$(LAPACK_LIBS) -llapack
  #FFTW=$(LCL_LIB_DIR)/fftw-3.3.3/
  FFTW=$(PFFTW_DIR)/
  FFTW_LIB=-L$(FFTW)/lib -lfftw3
  FFTW_INC=-I$(FFTW)/include

  HDF_LIB = -L$(PHDF5_LIB) -lhdf5_fortran -lhdf5
  HDF_INC = -I$(PHDF5_INC)

  ESIO_DIR=$(LCL_LIB_DIR)/software/esio019
  ESIO_LIB=-L$(ESIO_DIR)/lib -lesiof -lesio -lz -lm  -Wl,-rpath,$(ESIO_DIR)/lib
  ESIO_INC=-I$(ESIO_DIR)/include/esio -I$(ESIO_DIR)/lib
endif


INCLUDE= $(FFTW_INC) $(ESIO_INC) $(HDF_INC)
LIBS= $(MKL_LIB) $(FFTW_LIB) $(ESIO_LIB) $(HDF_LIB)

all: $(EXEC_DIR)$(EXEC)

$(EXEC_DIR)$(EXEC): $(OBJO)
	$(FORT) $(OPTS) $(INCLUDE) $(OBJO) -o $@ $(LDFLAGS) $(LIBS)

%.o: %.f90
	$(FORT) -c  $(OPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(FORT) -c $(OPTS) $(INCLUDE)

%.o: %.c
	$(CC) -c $< $(COPT) $(INCLUDE)

clean:
	rm *.o *~ \#* *.mod $(EXEC)
