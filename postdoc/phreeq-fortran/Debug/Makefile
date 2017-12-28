# ============================================================================
# Name        : Makefile
# Author      : Ivan Marin
# Version     :
# Copyright   : Proprietary 
# Description : Makefile for phreeqc-fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
IFORT = ifort
GFORTRAN = gfortran
#FORTRAN_COMPILER = gfortran

COMPILER_OPTIONS_GFORTRAN = -g -fbounds-check -g -Wall -Wuninitialized -O -Wtabs -fbacktrace 
COMPILER_OPTIONS_IFORT = -g -check all -fpe0 -warn -gen-interfaces -warn interfaces -traceback #-debug extended 


LIBRARY_PATH_2=/home/ispmarin/src/lib/iphreeqc-2.18.4-6386/lib/
INCLUDE_PATH_2=/home/ispmarin/src/lib/iphreeqc-2.18.4-6386/include/

LIBRARY_PATH_3=/home/ispmarin/src/lib/iphreeqc-3.0.0-7430/lib
INCLUDE_PATH_3=/home/ispmarin/src/lib/iphreeqc-3.0.0-7430/include/

PROJECT_DIR = /home/ispmarin/src/lavalwork/phreeqc-fortran
COMP_DIR = $(PROJECT_DIR)/Debug
SRC_DIR = $(PROJECT_DIR)/src/
SRC = phreeqc-fortran.f90
SRC_MOD = biogeochem.f90
BINARY = phreeqc-fortran 
OBJECTS = biogeochem.o phreeqc-fortran.o

ifeq ($(FC),gfortran)
	   COMPILER=$(GFORTRAN)
       COMPILER_OPTIONS=$(COMPILER_OPTIONS_GFORTRAN)
     else
       COMPILER=$(IFORT)
	   COMPILER_OPTIONS=$(COMPILER_OPTIONS_IFORT)
     endif


all: $(OBJECTS) 
	$(COMPILER) $(COMPILER_OPTIONS) $(OBJECTS) -o $(COMP_DIR)/$(BINARY) -I$(INCLUDE_PATH_3) -L$(LIBRARY_PATH_3) -liphreeqc

biogeochem.o biogeochem.mod: $(SRC_DIR)/$(SRC_MOD) 
	$(COMPILER) -c $(COMPILER_OPTIONS) $(SRC_DIR)/$(SRC_MOD) -I$(INCLUDE_PATH_3) -L$(LIBRARY_PATH_3) -liphreeqc

phreeqc-fortran.o: $(SRC_DIR)/$(SRC)
	$(COMPILER) -c $(COMPILER_OPTIONS) $(SRC_DIR)/$(SRC) -I$(INCLUDE_PATH_3) -L$(LIBRARY_PATH_3) -liphreeqc


clean:
	cd $(COMP_DIR)
	rm -f *.mod *.o 


