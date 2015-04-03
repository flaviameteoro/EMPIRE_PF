FC=mpif90
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace


BIN=../../bin/



all: 
	$(FC) $(FCOPTS) -o $(BIN)lorenz96 Lorenz96_empire.f90
	@echo 'Created file executable ' $(BIN)lorenz96
	$(FC) $(FCOPTS) -o $(BIN)lorenz96_v2 Lorenz96_empire_v2.f90
	@echo 'Created file executable ' $(BIN)lorenz96_v2


