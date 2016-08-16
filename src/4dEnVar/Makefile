# makefile for 4DEnVar
FC=mpif90

OBS=../../obs/
MODLOC:=$(OBS)
MODFLAG=-J
BIN=../../bin/

LIB_LIST=-lblas -llapack
FCOPTS=-O3 -Wall -fimplicit-none -fbounds-check -fbacktrace
FCOPTS+=$(MODFLAG) $(MODLOC)

MYFILES=$(OBS)model_specific.o $(OBS)sizes.o $(OBS)Rdata.o $(OBS)Qdata.o $(OBS)inner_products.o $(OBS)pf_control.o $(OBS)comms.o $(OBS)output_empire.o $(OBS)histogram.o $(OBS)data_io.o $(OBS)gen_rand.o $(OBS)random_d.o $(OBS)timestep_data.o $(OBS)matrix_pf.o $(OBS)generate_pf.o $(OBS)output_mat_tri.o $(OBS)output_rmse.o $(OBS)output_variance.o $(OBS)ziggurat.o var_data.o fourdenvardata.o 4dEnVar.o 4denvar_fcn.o 




CGFILES=$(OBS)cgsub.o $(OBS)cgfam.o $(OBS)cgsearch.o
LBFGSFILES=$(OBS)lbfgsb_sub.o $(OBS)lbfgs_sub.o $(OBS)lbfgsb.o $(OBS)linpack.o $(OBS)timer.o






# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f90
	$(FC) $(FCOPTS) -c -o $@ $<

all: clean EMPIRE CG Lbfgsb $(MYFILES)
	$(FC) $(FCOPTS) -o $(BIN)empire_4denvar $(MYFILES) $(CGFILES) $(LBFGSFILES) $(LIB_LIST)

EMPIRE:
	cd ../../;make EMPIRE

CG:
	cd ../optim/CG+;make

Lbfgsb:
	cd ../optim/Lbfgsb.3.0;make sub

clean:
	rm -rf *.o *.mod
