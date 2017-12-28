# ============================================================================
# Name        : Makefile
# Author      : Ivan Marin
# Version     :
# Copyright   : Proprietary 
# Description : Makefile for carbonate
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
IFORT = ifort
GFORTRAN = gfortran
#FORTRAN_COMPILER = gfortran

COMPILER_OPTIONS_GFORTRAN = -g -fbounds-check -g -Wall -Wuninitialized -O -Wtabs -fbacktrace 
COMPILER_OPTIONS_IFORT = -g  -fpe0 -warn -gen-interfaces -warn interfaces -traceback #-debug extended -check all

LAPACK = -llapack -lblas

PROJECT_DIR = /home/ispmarin/src/lavalwork/carbonate
COMP_DIR = $(PROJECT_DIR)/debug
SRC_DIR = $(PROJECT_DIR)/src
SRC = carbonatetest.f90
SRC_MOD = carbonate.f90
SRC_MOD2 = carbonate_lapack.f90
MODULEOBJ = carbonate.o
MODULEMOD = carbonate.mod
MODULEOBJ2 = carbonate_lapack.o
MODULEMOD2 = carbonate_lapack.mod
BINARYOBJ = carbonatetest.o
OBJECTS = $(MODULEOBJ) $(MODULEOBJ2) $(BINARYOBJ)
BINARY = carb


ifeq ($(FC),gfortran)
	   COMPILER=$(GFORTRAN)
       COMPILER_OPTIONS=$(COMPILER_OPTIONS_GFORTRAN)
     else
       COMPILER=$(IFORT)
	   COMPILER_OPTIONS=$(COMPILER_OPTIONS_IFORT)
     endif


all: $(OBJECTS) 
	$(COMPILER) $(COMPILER_OPTIONS) $(OBJECTS) -o $(COMP_DIR)/$(BINARY) $(LAPACK)
	$(COMPILER) $(COMPILER_OPTIONS) $(SRC_DIR)/carbonate.f -o $(COMP_DIR)/carb-f $(LAPACK)
	#$(COMPILER) $(COMPILER_OPTIONS) $(SRC_DIR)/carbonate_calcium.f -o $(COMP_DIR)/carb-c $(LAPACK)

$(MODULEOBJ) $(MODULEMOD): $(SRC_DIR)/$(SRC_MOD) 
	$(COMPILER) -c $(COMPILER_OPTIONS) $(SRC_DIR)/$(SRC_MOD) $(LAPACK)

$(MODULEOBJ2) $(MODULEMOD2): $(SRC_DIR)/$(SRC_MOD) 
	$(COMPILER) -c $(COMPILER_OPTIONS) $(SRC_DIR)/$(SRC_MOD2) $(LAPACK)

$(BINARYOBJ): $(SRC_DIR)/$(SRC)
	$(COMPILER) -c $(COMPILER_OPTIONS) $(SRC_DIR)/$(SRC)  $(LAPACK)


clean:
	cd $(COMP_DIR)
	rm -f *.mod *.o *__genmod*


