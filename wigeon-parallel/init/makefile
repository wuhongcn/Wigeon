
#Nompiler options for generic unix.

F90 = ifort
#LIBS =-L/export/intel/mkl/lib/em64t -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 
F90FLAGS = -O3 -ipo -xHost -fpp -w -r8 -mkl

CC = cc
CFLAGS = -O

F90_OBJS = mkl_dfti.o init_main.o zapp.o mesh.o growth.o power.o numisc.o \
           congauss.o greenfn.o potin.o smooth3d.o 

all:	init

init:	$(F90_OBJS)
	$(F90) $(F90FLAGS) $(F90_OBJS) $(LIBS) -o init 

clean:
	rm init *.o *.mod

.SUFFIXES: .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<


