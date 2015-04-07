#modify the following options:
#set FC to be the fortran compiler (mpi wrapper generally)
FC=mpif90

#set the compiliation options for the fortran compiler
FCOPTS = -O3 -fopenmp
FCOPTS = -fimplicit-none -Wall -fbounds-check -fbacktrace -fopenmp

#set the location of the libraries
LIB_LIST = -L$(METISDIR) -l$(METISLIB) -lblas -llapack

#set the option for placing module files
MODFLAG=-J


all: EMPIRE ALLTESTS TEST_R TEST_Q TEST_HQHTR #TEST_H

current_dir = $(shell pwd)/
OBS=$(current_dir)obs/
BIN=$(current_dir)bin/
MODLOC:=$(OBS)
SR_FILTS=$(current_dir)src/filters/
SR_UTILS=$(current_dir)src/utils/
SR_CONTS=$(current_dir)src/controllers/
SR_DATAS=$(current_dir)src/data/
SR_TESTS=$(current_dir)src/tests/
SR_OPERS=$(current_dir)src/operations/
OBJSQ= sizes.o pf_couple.o Qdata.o Rdata.o equivalent_weights_filter.o comms.o gen_rand.o random_d.o proposal_filter.o histogram.o pf_control.o data_io.o model_specific.o operator_wrappers.o quicksort.o resample.o diagnostics.o perturb_particle.o update_state.o genQ.o sir_filter.o stochastic_model.o tests.o letkf_analysis.o deterministic_model.o inner_products.o trajectories.o
OBJS=$(addprefix $(OBS),$(OBJSQ))
FCOPTS+=$(MODFLAG) $(MODLOC)

$(OBS)trajectories.o: $(SR_UTILS)trajectories.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)trajectories.f90 -o $@

$(OBS)histogram.o: $(SR_UTILS)histogram.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)histogram.f90 -o $@

$(OBS)stochastic_model.o: $(SR_FILTS)stochastic_model.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)stochastic_model.f90 -o $@

$(OBS)deterministic_model.o: $(SR_FILTS)deterministic_model.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)deterministic_model.f90 -o $@

$(OBS)sir_filter.o: $(SR_FILTS)sir_filter.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)sir_filter.f90 -o $@

$(OBS)genQ.o: $(SR_UTILS)genQ.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)genQ.f90 -o $@

$(OBS)Qdata.o: $(SR_DATAS)Qdata.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_DATAS)Qdata.f90 -o $@

$(OBS)Rdata.o: $(SR_DATAS)Rdata.f90
	$(FC) $(FCOPTS) -c $(SR_DATAS)Rdata.f90 -o $@

$(OBS)inner_products.o: $(SR_OPERS)inner_products.f90
	$(FC) $(FCOPTS) -c $(SR_OPERS)inner_products.f90 -o $@

$(OBS)model_specific.o: model_specific.f90 $(OBS)sizes.o 
	$(FC) $(FCOPTS) -c model_specific.f90 -o $@

$(OBS)operator_wrappers.o: $(SR_OPERS)operator_wrappers.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)operator_wrappers.f90 -o $@

$(OBS)pf_control.o: $(SR_CONTS)pf_control.f90 $(OBS)sizes.o $(OBS)histogram.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)pf_control.f90 -o $@

$(OBS)perturb_particle.o: $(SR_OPERS)perturb_particle.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)perturb_particle.f90 -o $@

$(OBS)update_state.o: $(SR_OPERS)update_state.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)update_state.f90 -o $@

$(OBS)data_io.o: $(SR_UTILS)data_io.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)data_io.f90 -o $@

$(OBS)proposal_filter.o: $(SR_FILTS)proposal_filter.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)proposal_filter.f90 -o $@

$(OBS)sizes.o: $(SR_CONTS)sizes.f90
	$(FC) $(FCOPTS) -c $(SR_CONTS)sizes.f90 -o $@

$(OBS)comms.o: $(SR_UTILS)comms.f90 $(OBS)sizes.o $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)comms.f90 -o $@

$(OBS)equivalent_weights_filter.o: $(SR_FILTS)equivalent_weights_filter.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)equivalent_weights_filter.f90 -o $@

$(OBS)random_d.o: $(SR_UTILS)random_d.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)random_d.f90 -o $@

$(OBS)resample.o: $(SR_OPERS)resample.f90 $(OBS)pf_control.o $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)resample.f90 -o $@

$(OBS)diagnostics.o: $(SR_UTILS)diagnostics.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)diagnostics.f90 -o $@

$(OBS)quicksort.o: $(SR_UTILS)quicksort.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)quicksort.f90 -o $@

$(OBS)pf_couple.o: $(SR_CONTS)pf_couple.f90 $(OBS)comms.o $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)pf_couple.f90 -o $@

$(OBS)alltests.o: $(SR_TESTS)alltests.f90 $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)alltests.f90 -o $@

$(OBS)test_h.o: $(SR_TESTS)test_h.f90 $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_h.f90 -o $@

$(OBS)test_r.o: $(SR_TESTS)test_r.f90 $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_r.f90 -o $@

$(OBS)tests.o: $(SR_TESTS)tests.f90 $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)tests.f90 -o $@

$(OBS)test_q.o: $(SR_TESTS)test_q.f90 $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_q.f90 -o $@

$(OBS)test_hqhtr.o: $(SR_TESTS)test_hqhtr.f90 $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_hqhtr.f90 -o $@

$(OBS)gen_rand.o: $(SR_OPERS)gen_rand.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)gen_rand.f90 -o $@


#ENKF SECTION:
$(OBS)letkf_analysis.o: $(SR_FILTS)letkf_analysis.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)letkf_analysis.f90 -o $@







EMPIRE: $(OBJS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)empire $(OBJS) $(LIB_LIST)


OBJS_ALLTEST = $(shell echo $(OBJS) | sed 's/pf_couple/alltests/g') 
ALLTESTS: $(OBJS_ALLTEST)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)alltests $(OBJS_ALLTEST) $(LIB_LIST)

OBJS_TEST_H = $(shell echo $(OBJS) | sed 's/pf_couple/test_h/g') 
TEST_H: $(OBJS_TEST_H)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_h $(OBJS_TEST_H) $(LIB_LIST)

OBJS_TEST_R = $(shell echo $(OBJS) | sed 's/pf_couple/test_r/g') 
TEST_R: $(OBJS_TEST_R)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_r $(OBJS_TEST_R) $(LIB_LIST)

OBJS_TEST_Q = $(shell echo $(OBJS) | sed 's/pf_couple/test_q/g') 
TEST_Q: $(OBJS_TEST_Q)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_q $(OBJS_TEST_Q) $(LIB_LIST)

OBJS_TEST_HQHTR = $(shell echo $(OBJS) | sed 's/pf_couple/test_hqhtr/g') 
TEST_HQHTR: $(OBJS_TEST_HQHTR)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_hqhtr $(OBJS_TEST_HQHTR) $(LIB_LIST)


docs: doc_latex

doc_html: FORCE
	doxygen .Doxyfile

models: FORCE
	cd models;make -e

export
lorenz63:
	cd models/lorenz63;make -e

lorenz96:
	cd models/lorenz96;make -e

lorenz: lorenz63 lorenz96

linear:
	cd models/linear;make -e

4dEnVar:
	cd src/4dEnVar;make -e


doc_latex: doc_html
	cd doc/latex && make

FORCE:


clean:
	rm -f obs/* bin/*

minimal:
	cd examples/;make -e minimal
