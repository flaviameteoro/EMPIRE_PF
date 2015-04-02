# makefile for 4DEnVar
FC=mpif90

OBS=../../obs/
MODLOC:=$(OBS)
MODFLAG=-J
BIN=../../bin/

LIB_LIST=-lblas
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
FCOPTS+=$(MODFLAG) $(MODLOC)

OBJECTS_EMPIRE_COMMS=minimal_empire_comms.o $(OBS)comms.o $(OBS)pf_control.o $(OBS)sizes.o $(OBS)histogram.o



# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(FC) $(FCOPTS) -c -o $@ $<

all: clean $(OBJECTS_EMPIRE_COMMS)
	$(FC) $(FCOPTS) -o $(BIN)minimal_empire_comms $(OBJECTS_EMPIRE_COMMS) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_empire_comms
	@rm -rf *.o

clean:
	rm -rf *.o *.mod
