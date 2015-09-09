# makefile for 4DEnVar
FC=mpif90

OBS=../../obs/
MODLOC:=$(OBS)
MODFLAG=-J
BIN=../../bin/

LIB_LIST=-lblas
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
FCOPTS+=$(MODFLAG) $(MODLOC)

OBJECTS_MODEL=minimal_model.o
OBJECTS_MODEL_v2=minimal_model_v2.o
OBJECTS_MODEL_v3=minimal_model_v3.o


# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(FC) $(FCOPTS) -c -o $@ $<

all: clean $(OBJECTS_MODEL) $(OBJECTS_MODEL_v2) $(OBJECTS_MODEL_v3)
	$(FC) $(FCOPTS) -o $(BIN)minimal_model $(OBJECTS_MODEL) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_model
	$(FC) $(FCOPTS) -o $(BIN)minimal_model_v2 $(OBJECTS_MODEL_v2) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_model_v2
	$(FC) $(FCOPTS) -o $(BIN)minimal_model_v3 $(OBJECTS_MODEL_v3) $(LIB_LIST)
	@echo 'Created file executable ' $(BIN)minimal_model_v3
	@rm -rf *.o

clean:
	rm -rf *.o *.mod
