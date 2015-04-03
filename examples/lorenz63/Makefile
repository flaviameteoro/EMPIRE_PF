FC=mpif90
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace


BIN=../../bin/



all: 
	$(FC) $(FCOPTS) -o $(BIN)lorenz63 Lorenz63_empire.f90
	@echo 'Created file executable ' $(BIN)lorenz63
	$(FC) $(FCOPTS) -o $(BIN)lorenz63_v2 Lorenz63_empire_v2.f90
	@echo 'Created file executable ' $(BIN)lorenz63_v2
