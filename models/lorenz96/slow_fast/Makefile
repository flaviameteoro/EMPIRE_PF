FC=mpif90
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace


BIN=../../../bin/



all: 
	$(FC) $(FCOPTS) -o $(BIN)lorenz96_slow_fast Lorenz96_slow_fast_empire.f90
	@echo 'Created file executable ' $(BIN)lorenz96_slow_fast
	$(FC) $(FCOPTS) -o $(BIN)lorenz96_slow_fast_v2 Lorenz96_slow_fast_empire_v2.f90
	@echo 'Created file executable ' $(BIN)lorenz96_slow_fast_v2

standalone:
	$(FC) $(FCOPTS) -o $(BIN)lorenz96_slow_fast_standalone Lorenz96_slow_fast.f90
	@echo 'Created file executable ' $(BIN)lorenz96_slow_fast_standalone

