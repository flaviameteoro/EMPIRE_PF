# makefile for 4DEnVar
F90=mpif90
MYFILES=vardata.o 4dEnVar.o 4denvar_fcn.o ../../obs/model_specific.o ../../obs/sizes.o ../../obs/Rdata.o ../../obs/Qdata.o ../../obs/inner_products.o ../../obs/pf_control.o ../../obs/comms.o ../../obs/histogram.o ../../obs/data_io.o ../../obs/gen_rand.o ../../obs/random_d.o
LIB_LIST=-lblas
OPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
MODFLAG=-J


OBS=../../obs/
BIN=../../bin/
MODLOC:=$(OBS)
OPTS+=$(MODFLAG) $(MODLOC)


CGFILES=../optim/CG+/cgsub.o ../optim/CG+/cgfam.o ../optim/CG+/cgsearch.o
LBFGSFILES=../optim/Lbfgsb.3.0/lbfgsb_sub.o ../optim/Lbfgsb.3.0/lbfgs_sub.o ../optim/Lbfgsb.3.0/lbfgsb.o ../optim/Lbfgsb.3.0/linpack.o ../optim/Lbfgsb.3.0/timer.o






# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(F90) $(OPTS) -c -o $@ $<

all: $(MYFILES)
	cd ../../;make
	cd ../optim/CG+;make
	cd ../optim/Lbfgsb.3.0;make sub
	$(F90) $(OPTS) -o EMPIRE_4DENVAR $(MYFILES) $(CGFILES) $(LBFGSFILES) $(LIB_LIST)

#sub: $(MYFILES)
#	$(F90) $(MYFILES) $(OPTS) -o subroutine_cg $(LIB_LIST)

clean:
	rm -rf *.o *.mod
