
#Nompiler options for generic unix.

F90 = ifort
#LIBS =-L/share/intel/composer_xe_2015.0.090/mkl/lib/mic -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 
F90FLAGS = -O3 -ipo -xHost -fpp -w -r8 -openmp -mkl
#
#
CC = cc
CFLAGS = -O
#-------------------------------------------------------------------------------
#weno_module.o datablock.mod : weno_module.f90
#	$(F90) -c -fpp -r8 weno_module.f90

F90_OBJS = mkl_dfti.o weno_module.o num_recipes.o cosmology.o pm_gravity.o \
           igm_proc.o igm_stat.o readin_out.o weno.o wenopm_main.o

all:	cosmo

cosmo:	$(F90_OBJS)
	$(F90) $(F90FLAGS) $(F90_OBJS) $(LIBS) -o cosmo

.PHONY:clean
clean:
	rm *.o *.mod cosmo

.SUFFIXES: .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<
