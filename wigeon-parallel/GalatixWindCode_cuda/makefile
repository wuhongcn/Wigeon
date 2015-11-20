
#Nompiler options for generic unix.

F90=mpiifort
F90FLAGS =-w -r8 -fpp -O3 -ipo -xHost -L /soft/cuda/6.0.37/lib64 -lcudart 
CUDA=nvcc
CUDAFLAGS= -arch sm_35  -use_fast_math -Xptxas -dlcm=cg -Xptxas -v,-abi=no #--ptxas-options=-v#-maxrregcount 32
CUDALIBS=

#
#
CC = icc
CFLAGS = -O
#-------------------------------------------------------------------------------
#weno_module.o datablock.mod : weno_module.f90
#	$(F90) -c -r8 weno_module.f90
F90_OBJS = weno_module_flat.o igm_stat_flat.o  readin_out_flat.o igm_proc_flat.o weno_GPU_final.o weno_flat.o  wenotest_flat.o

all:	sup_flat

sup_flat:	$(F90_OBJS)
	        $(F90) $(F90FLAGS) $(F90_OBJS) $(LIBS) -o sup_flat

.PHONY:clean
clean:
	rm -rf *.o *.mod *.sh.* sup_flat\#* *~

.PHONY:cleanout
cleanout:
	rm *.sh.* prun4000.log

.SUFFIXES: .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<

.SUFFIXES: .cu

.cu.o:
	$(CUDA) $(CUDAFLAGS) -c $<
