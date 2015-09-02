# makefile for 4DEnVar
FC=mpif90

OBS=../../obs/
MODLOC:=$(OBS)
MODFLAG=-J
BIN=../../bin/

LIB_LIST=-lblas
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
FCOPTS+=$(MODFLAG) $(MODLOC)

MYFILES=vardata.o 4dEnVar.o 4denvar_fcn.o $(OBS)model_specific.o $(OBS)sizes.o $(OBS)Rdata.o $(OBS)Qdata.o $(OBS)inner_products.o $(OBS)pf_control.o $(OBS)comms.o $(OBS)histogram.o $(OBS)data_io.o $(OBS)gen_rand.o $(OBS)random_d.o $(OBS)timestep_data.o $(OBS)matrix_pf.o $(OBS)generate_pf.o $(OBS)output_mat_tri.o




CGFILES=../optim/CG+/cgsub.o ../optim/CG+/cgfam.o ../optim/CG+/cgsearch.o
LBFGSFILES=../optim/Lbfgsb.3.0/lbfgsb_sub.o ../optim/Lbfgsb.3.0/lbfgs_sub.o ../optim/Lbfgsb.3.0/lbfgsb.o ../optim/Lbfgsb.3.0/linpack.o ../optim/Lbfgsb.3.0/timer.o






# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(FC) $(FCOPTS) -c -o $@ $<

all: clean $(MYFILES)
	cd ../../;make
	cd ../optim/CG+;make
	cd ../optim/Lbfgsb.3.0;make sub
	$(FC) $(FCOPTS) -o $(BIN)empire_4denvar $(MYFILES) $(CGFILES) $(LBFGSFILES) $(LIB_LIST)


clean:
	rm -rf *.o *.mod
