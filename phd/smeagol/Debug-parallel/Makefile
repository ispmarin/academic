CC=/home/ispmarin/lib/openmpi-1.3.3-gcc-64/bin/mpicxx
MKLPATH=/home/ispmarin/lib/intel/mkl/lib/em64t

LIB_INTEL = $(MKLPATH)/libmkl_scalapack_lp64.a $(MKLPATH)/libmkl_solver_lp64.a -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -openmp -lpthread

LIB = /home/ispmarin/lib/intel/mkl/lib/em64t/libmkl_scalapack_lp64.a /home/ispmarin/lib/intel/mkl/lib/em64t/libmkl_solver_lp64.a -Wl,--start-group /home/ispmarin/lib/intel/mkl/lib/em64t/libmkl_intel_lp64.a /home/ispmarin/lib/intel/mkl/lib/em64t/libmkl_gnu_thread.a /home/ispmarin/lib/intel/mkl/lib/em64t/libmkl_core.a /home/ispmarin/lib/intel/mkl/lib/em64t/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -fopenmp -lpthread -I/usr/include/mkl -I/usr/include/openmpi 

FILES = ../src-parallel/p_smeagol.cpp ../src-parallel/functions.cpp ../src/aux_functions.cpp ../src/elements.cpp ../src/insert_from_file.cpp ../src/export.cpp 

CFLAGS =  -O3 -msse4.2 -mfpmath=sse -fomit-frame-pointer -march=core2 
CFLAGS_INTEL = -O3 
CFLAGS_DEBUG = -g -Wall -Wstrict-overflow -Wdisabled-optimization -Warray-bounds

INCLUDE = -I/usr/include/mkl -I/usr/include/openmpi


all: source_env debug

intel: $(FILEs)
	$(CC) $(FILES) $(CFLAGS_INTEL) $(LIB_INTEL) $(INCLUDE)

debug: $(FILEs)
	$(CC) $(FILES) $(LIB) $(CFLAGS_DEBUG) $(INCLUDE) -o p_smeagol
	@echo 'Compilation done'
	
build: $(FILEs)
	$(CC) $(FILES) $(LIB) $(CFLAGS) $(INCLUDE) -o p_smeagol
	@echo 'Compilation done'
	
	
source_env:
	/bin/bash /home/ispmarin/env/openmpivars-gcc-64.sh

clean: 
	rm -rf *.o 	
