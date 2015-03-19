FC=mpif90

#set the compiliation options for the fortran compiler
FCOPTS = -O3 -fopenmp
FCOPTS = -fimplicit-none -Wall -fbounds-check -fbacktrace -fopenmp

#set the location of the libraries
LIB_LIST = -lblas -llapack

all: linear

linear: 
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)linear_empire_vader linear_empire_vader.f90 $(LIB_LIST)

clean:
	rm -rf *.o $(BIN)linear_empire_vader *.mod
