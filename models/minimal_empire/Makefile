# makefile for 4DEnVar
FC=mpif90

OBS=../../obs/
MODLOC:=$(OBS)
MODFLAG=-J
BIN=../../bin/

LIB_LIST=-lblas
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
FCOPTS+=$(MODFLAG) $(MODLOC)

OBJECTS_EMPIRE=minimal_empire.o $(OBS)comms.o $(OBS)pf_control.o $(OBS)sizes.o $(OBS)histogram.o $(OBS)timestep_data.o



# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(FC) $(FCOPTS) -c -o $@ $<

all: clean $(OBJECTS_EMPIRE)
	$(FC) $(FCOPTS) -o $(BIN)minimal_empire $(OBJECTS_EMPIRE) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_empire
	@rm -rf *.o

clean:
	rm -rf *.o *.mod
