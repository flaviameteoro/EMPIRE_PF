# makefile for 4DEnVar
F90=mpif90
MYFILES=vardata.o 4dEnVar.o objective_function.o objective_gradient.o
LIB_LIST=-lblas
OPTS=-O3 -Wall -fimplicit-none
MODFLAG=-J


OBS=../../obs/
BIN=../../bin/
MODLOC:=$(OBS)
OPTS+=$(MODFLAG) $(MODLOC)

# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(F90) $(OPTS) -c -o $@ $<

all: $(MYFILES)


#sub: $(MYFILES)
#	$(F90) $(MYFILES) $(OPTS) -o subroutine_cg $(LIB_LIST)

clean:
	rm -rf *.o *.mod
