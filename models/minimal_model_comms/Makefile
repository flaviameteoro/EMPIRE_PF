# makefile for 4DEnVar
FC=mpif90

OBS=../../obs/
MODLOC:=$(OBS)
MODFLAG=-J
BIN=../../bin/

LIB_LIST=-lblas
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
FCOPTS+=$(MODFLAG) $(MODLOC)

OBJECTS_COMMS=minimal_model_comms.o
OBJECTS_COMMS_v2=minimal_model_comms_v2.o


# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(FC) $(FCOPTS) -c -o $@ $<

all: clean $(OBJECTS_COMMS) $(OBJECTS_COMMS_v2)
	$(FC) $(FCOPTS) -o $(BIN)minimal_model_comms $(OBJECTS_COMMS) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_model_comms
	$(FC) $(FCOPTS) -o $(BIN)minimal_model_comms_v2 $(OBJECTS_COMMS_v2) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_model_comms_v2
	@rm -rf *.o

clean:
	rm -rf *.o *.mod
