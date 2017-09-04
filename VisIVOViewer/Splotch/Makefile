#######################################################################
#  Splotch V4.5                                                       #
#######################################################################

#--------------------------------------- Switch on DataSize
#OPT     +=  -DLONGIDS

#--------------------------------------- Switch on MPI
#OPT	+=  -DUSE_MPI
#OPT	+=  -DUSE_MPIIO

#--------------------------------------- Switch on HDF5

#OPT     +=  -DHDF5
#OPT     +=  -DH5_USE_16_API

#--------------------------------------- Visual Studio Option
#OPT	+=  -DVS

#--------------------------------------- CUDA options
OPT     +=  -DCUDA
OPT     +=  -DNO_WIN_THREAD

#--------------------------------------- Turn on VisIVO stuff
#OPT	+=  -DSPLVISIVO

#--------------------------------------- Select target Computer

#SYSTYPE="SP6"
SYSTYPE="GP"
#SYSTYPE="PLX"
#SYSTYPE="BGP"

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       = mpic++        # sets the C-compiler (default)
else
CC       = g++        # sets the C-compiler (default)
endif
OMP      = -fopenmp

OPTIMIZE = -std=c++98 -pedantic -Wno-long-long -Wfatal-errors -Wextra -Wall -Wstrict-aliasing=2 -Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wold-style-cast -O2 -g    # optimization and warning flags (default)
SUP_INCL = -I. -Icxxsupport -Ic_utils


ifeq (USE_MPIIO,$(findstring USE_MPIIO,$(OPT)))
SUP_INCL += -Impiio-1.0/include/
endif

ifeq ($(SYSTYPE),"SP6")
ifeq (HDF5,$(findstring HDF5,$(OPT)))
HDF5_HOME = /cineca/prod/libraries/hdf5/1.8.4_ser/xl--10.1
LIB_HDF5  = -L$(HDF5_HOME)/lib -lhdf5 -L/cineca/prod/libraries/zlib/1.2.3/xl--10.1/lib/ -lz -L/cineca/prod/libraries/szlib/2.1/xl--10.1/lib/ -lsz
HDF5_INCL = -I$(HDF5_HOME)/include
endif
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  mpCC_r
else
CC       =  xlc++
endif
OPTIMIZE =  -q64 -O3 -qarch=auto -qtune=auto -qinline
LIB_OPT	 =  -bstackpsize:64k -bdatapsize:64k -btextpsize:64k
OMP =
endif

ifeq ($(SYSTYPE),"GP")
CC       =  nvcc -g
CUDA_HOME=/Developer/GPU\ Computing
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  mpicxx -g -I$(CUDA_HOME)/CUDALibraries/common/inc -I$(CUDA_HOME)/C/common/inc -I$(CUDA_HOME)/include
endif
NVCC       =  nvcc -g 
OPTIMIZE = -O2
LIB_OPT  =  -L/usr/local/cuda/lib/ -lcudart
OMP =  
#-Xcompiler -openmp
SUP_INCL += -I$(CUDA_HOME)/CUDALibraries/common/inc -I$(CUDA_HOME)/C/common/inc # -I$(CUDA_HOME)/include  -Icuda
endif

ifeq ($(SYSTYPE),"BGP")
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       = mpixlcxx_r
else
CC       = bgxlC_r
endif
OPTIMIZE = -O3 -qstrict -qarch=450d -qtune=450 # -qipa=inline
LIB_OPT  =
OMP =   -qsmp=omp -qthreaded

endif

ifeq ($(SYSTYPE),"PLX")
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  mpiCC -g 
endif
ifeq (CUDA,$(findstring CUDA,$(OPT)))
NVCC       =  nvcc -g 
# cuda/3.1 :
LIB_OPT  += -Xlinker -L$(NVCC_HOME)/lib64 -lcudart
SUP_INCL += -I$(CUDAUTIL_INC) -I$(NVCC_HOME)/include #-Icuda
# nvcc/2.2 :
#LIB_OPT  += -Xlinker -L$(NVCC_HOME)/lib -lcudart
#SUP_INCL += -I$(CUDASDK_HOME)/common/inc -I$(NVCC_HOME)/include #-Icuda
endif
OPTIMIZE = -O2 -DDEBUG
OMP = -fopenmp
endif

#--------------------------------------- Here we go

OPTIONS = $(OPTIMIZE) $(OPT)

EXEC   = Splotch4.5$(SYSTYPE)

OBJS  =	kernel/transform.o cxxsupport/error_handling.o \
        reader/mesh_reader.o reader/visivo_reader.o \
	cxxsupport/mpi_support.o cxxsupport/paramfile.o cxxsupport/string_utils.o cxxsupport/announce.o cxxsupport/ls_image.o reader/gadget_reader.o \
	reader/millenium_reader.o reader/bin_reader.o reader/bin_reader_mpi.o \
	splotch/splotchutils.o splotch/splotch.o \
	splotch/scenemaker.o splotch/splotch_host.o cxxsupport/walltimer.o c_utils/walltime_c.o \
	booster/mesh_creator.o booster/randomizer.o booster/p_selector.o booster/m_rotation.o

ifeq (HDF5,$(findstring HDF5,$(OPT)))
OBJS += reader/hdf5_reader.o 
endif
ifeq (CUDA,$(findstring CUDA,$(OPT)))
OBJS += cuda/splotch.o cuda/CuPolicy.o cuda/splotch_cuda2.o cuda/deviceQuery.o
endif
ifeq (USE_MPIIO,$(findstring USE_MPIIO,$(OPT)))
LIB_MPIIO = -Lmpiio-1.0/lib -lpartition
endif

INCL   = */*.h Makefile

CPPFLAGS = $(OPTIONS) $(SUP_INCL) $(HDF5_INCL) $(OMP)

CUFLAGS = $(OPTIONS) $(SUP_INCL)

LIBS   = $(LIB_OPT) $(OMP)

.SUFFIXES: .o .cc .cxx .cpp .cu

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cpp.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cu.o:
	$(NVCC) -c $(CUFLAGS) -o "$@" "$<"

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC) $(LIB_MPIIO) $(LIB_HDF5)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)

