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

-include Makefile.in

default: EMPIRE ALLTESTS TEST_R TEST_Q TEST_HQHTR 

all: EMPIRE ALLTESTS TEST_R TEST_Q TEST_HQHTR models 4dEnVar

current_dir = $(shell pwd)/
OBS=$(current_dir)obs/
BIN=$(current_dir)bin/
MODLOC:=$(OBS)
SR_FILTS=$(current_dir)src/filters/
SR_SMOOT=$(current_dir)src/smoothers/
SR_UTILS=$(current_dir)src/utils/
SR_CONTS=$(current_dir)src/controllers/
SR_USERS=$(current_dir)src/user/
SR_USER_MDL=$(current_dir)src/user/model/
SR_TESTS=$(current_dir)src/tests/
SR_OPERS=$(current_dir)src/operations/
SR_4DENVAR=$(current_dir)src/4dEnVar/
SR_VAR=$(current_dir)src/var/
SR_CG=$(current_dir)src/optim/CG+/
CGOPTS = $(FCOPTS) #$(shell echo $(FCOPTS) | sed 's/-\<[-a-zA-Z0-9]*implicit[-a-zA-Z0-9]*\>//g')
SR_LBFGSB=$(current_dir)src/optim/Lbfgsb.3.0/

CGFILES_local=cgsub.o cgfam.o cgsearch.o
LBFGSFILES_local=lbfgsb_sub.o lbfgs_sub.o lbfgsb.o linpack.o timer.o

OBJSQ= compile_options.o $(CGFILES_local) $(LBFGSFILES_local) comm_version.o output_empire.o timestep_data.o sizes.o empire_main.o Bdata.o Qdata.o Rdata.o equivalent_weights_filter.o comms.o var_data.o ziggurat.o gen_rand.o random_d.o proposal_filter.o histogram.o allocate_pf.o pf_control.o letks.o matrix_pf.o data_io.o model_specific.o operator_wrappers.o quicksort.o resample.o diagnostics.o perturb_particle.o update_state.o genQ.o sir_filter.o stochastic_model.o tests.o letkf_analysis.o deterministic_model.o inner_products.o trajectories.o user_perturb_particle.o generate_pf.o output_mat_tri.o equivalent_weights_filter_zhu.o lambertw.o randperm.o user_initialise_mpi.o loc_function.o phalf_etkf.o phalf.o threedvar_data.o three_d_var_all_particles.o threedvar_fcn.o three_d_var.o fcn.o output_spatial_rmse.o output_variance.o output_forecast.o output_ens_rmse.o model_as_subroutine_data.o model_as_subroutine_return.o model_as_subroutine_start.o model_as_subroutine_initialise.o
OBJS=$(addprefix $(OBS),$(OBJSQ))
FCOPTS+=$(MODFLAG) $(MODLOC)

$(OBS)cgsub.o: $(SR_CG)cgsub.f90
	$(FC) $(FCOPTS) -c $(SR_CG)cgsub.f90 -o $@

$(OBS)cgfam.o: $(SR_CG)cgfam.f
	$(FC) $(CGOPTS) -c $(SR_CG)cgfam.f -o $@

$(OBS)cgsearch.o: $(SR_CG)cgsearch.f
	$(FC) $(CGOPTS) -c $(SR_CG)cgsearch.f -o $@

$(OBS)lbfgsb_sub.o: $(SR_LBFGSB)lbfgsb_sub.f90
	$(FC) $(FCOPTS) -c $(SR_LBFGSB)lbfgsb_sub.f90 -o $@

$(OBS)lbfgs_sub.o: $(SR_LBFGSB)lbfgs_sub.f90
	$(FC) $(FCOPTS) -c $(SR_LBFGSB)lbfgs_sub.f90 -o $@

$(OBS)lbfgsb.o: $(SR_LBFGSB)lbfgsb.f
	$(FC) $(FCOPTS) -c $(SR_LBFGSB)lbfgsb.f -o $@

$(OBS)linpack.o: $(SR_LBFGSB)linpack.f
	$(FC) $(FCOPTS) -c $(SR_LBFGSB)linpack.f -o $@

$(OBS)timer.o: $(SR_LBFGSB)timer.f
	$(FC) $(FCOPTS) -c $(SR_LBFGSB)timer.f -o $@

$(OBS)output_mat_tri.o: $(SR_UTILS)output_mat_tri.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)output_mat_tri.f90 -o $@

$(OBS)output_spatial_rmse.o: $(SR_UTILS)output_spatial_rmse.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)output_spatial_rmse.f90 -o $@

$(OBS)output_variance.o: $(SR_UTILS)output_variance.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)output_variance.f90 -o $@

$(OBS)output_forecast.o: $(SR_UTILS)output_forecast.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)output_forecast.f90 -o $@

$(OBS)output_ens_rmse.o: $(SR_UTILS)output_ens_rmse.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)output_ens_rmse.f90 -o $@

$(OBS)matrix_pf.o: $(SR_UTILS)matrix_pf.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)matrix_pf.f90 -o $@

$(OBS)lambertw.o: $(SR_UTILS)lambertw.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)lambertw.f90 -o $@

$(OBS)randperm.o: $(SR_UTILS)randperm.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)randperm.f90 -o $@

$(OBS)loc_function.o: $(SR_UTILS)loc_function.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)loc_function.f90 -o $@

$(OBS)generate_pf.o: $(SR_UTILS)generate_pf.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)generate_pf.f90 -o $@

$(OBS)trajectories.o: $(SR_UTILS)trajectories.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)trajectories.f90 -o $@

$(OBS)histogram.o: $(SR_UTILS)histogram.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)histogram.f90 -o $@

$(OBS)threedvar_data.o: $(SR_VAR)threedvar_data.f90
	$(FC) $(FCOPTS) -c $(SR_VAR)threedvar_data.f90 -o $@

$(OBS)threedvar_fcn.o: $(SR_VAR)threedvar_fcn.f90
	$(FC) $(FCOPTS) -c $(SR_VAR)threedvar_fcn.f90 -o $@

$(OBS)three_d_var.o: $(SR_VAR)three_d_var.f90
	$(FC) $(FCOPTS) -c $(SR_VAR)three_d_var.f90 -o $@

$(OBS)three_d_var_all_particles.o: $(SR_VAR)three_d_var_all_particles.f90
	$(FC) $(FCOPTS) -c $(SR_VAR)three_d_var_all_particles.f90 -o $@

$(OBS)var_data.o: $(SR_4DENVAR)var_data.f90
	$(FC) $(FCOPTS) -c $(SR_4DENVAR)var_data.f90 -o $@

$(OBS)fcn.o: $(SR_VAR)fcn.f90
	$(FC) $(FCOPTS) -c $(SR_VAR)fcn.f90 -o $@

$(OBS)stochastic_model.o: $(SR_FILTS)stochastic_model.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)stochastic_model.f90 -o $@

$(OBS)deterministic_model.o: $(SR_FILTS)deterministic_model.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)deterministic_model.f90 -o $@

$(OBS)sir_filter.o: $(SR_FILTS)sir_filter.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)sir_filter.f90 -o $@

$(OBS)genQ.o: $(SR_UTILS)genQ.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)genQ.f90 -o $@

$(OBS)Bdata.o: $(SR_USERS)Bdata.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_USERS)Bdata.f90 -o $@

$(OBS)Qdata.o: $(SR_USERS)Qdata.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_USERS)Qdata.f90 -o $@

$(OBS)Rdata.o: $(SR_USERS)Rdata.f90
	$(FC) $(FCOPTS) -c $(SR_USERS)Rdata.f90 -o $@

$(OBS)user_perturb_particle.o: $(SR_USERS)user_perturb_particle.f90
	$(FC) $(FCOPTS) -c $(SR_USERS)user_perturb_particle.f90 -o $@

$(OBS)user_initialise_mpi.o: $(SR_USERS)user_initialise_mpi.f90
	$(FC) $(FCOPTS) -c $(SR_USERS)user_initialise_mpi.f90 -o $@

$(OBS)inner_products.o: $(SR_OPERS)inner_products.f90
	$(FC) $(FCOPTS) -c $(SR_OPERS)inner_products.f90 -o $@

$(OBS)phalf_etkf.o: $(SR_OPERS)phalf_etkf.f90
	$(FC) $(FCOPTS) -c $(SR_OPERS)phalf_etkf.f90 -o $@

$(OBS)phalf.o: $(SR_OPERS)phalf.f90
	$(FC) $(FCOPTS) -c $(SR_OPERS)phalf.f90 -o $@

$(OBS)comm_version.o: comm_version.f90
	$(FC) $(FCOPTS) -c comm_version.f90 -o $@

$(OBS)model_specific.o: model_specific.f90 $(OBS)sizes.o 
	$(FC) $(FCOPTS) -c model_specific.f90 -o $@

$(OBS)operator_wrappers.o: $(SR_OPERS)operator_wrappers.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)operator_wrappers.f90 -o $@

$(OBS)compile_options.o: $(SR_CONTS)compile_options.f90
	$(FC) $(FCOPTS) -c $(SR_CONTS)compile_options.f90 -o $@

$(OBS)pf_control.o: $(SR_CONTS)pf_control.f90 $(OBS)var_data.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)pf_control.f90 -o $@

$(OBS)perturb_particle.o: $(SR_OPERS)perturb_particle.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)perturb_particle.f90 -o $@

$(OBS)update_state.o: $(SR_OPERS)update_state.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)update_state.f90 -o $@

$(OBS)data_io.o: $(SR_UTILS)data_io.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)data_io.f90 -o $@

$(OBS)allocate_pf.o: $(SR_UTILS)allocate_pf.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)allocate_pf.f90 -o $@

$(OBS)proposal_filter.o: $(SR_FILTS)proposal_filter.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)proposal_filter.f90 -o $@

$(OBS)sizes.o: $(SR_CONTS)sizes.f90
	$(FC) $(FCOPTS) -c $(SR_CONTS)sizes.f90 -o $@

$(OBS)timestep_data.o: $(SR_CONTS)timestep_data.f90
	$(FC) $(FCOPTS) -c $(SR_CONTS)timestep_data.f90 -o $@

$(OBS)output_empire.o: $(SR_CONTS)output_empire.f90
	$(FC) $(FCOPTS) -c $(SR_CONTS)output_empire.f90 -o $@

$(OBS)comms.o: $(SR_UTILS)comms.f90 $(OBS)sizes.o $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)comms.f90 -o $@

$(OBS)equivalent_weights_filter.o: $(SR_FILTS)equivalent_weights_filter.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)equivalent_weights_filter.f90 -o $@

$(OBS)equivalent_weights_filter_zhu.o: $(SR_FILTS)equivalent_weights_filter_zhu.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)equivalent_weights_filter_zhu.f90 -o $@

$(OBS)random_d.o: $(SR_UTILS)random_d.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)random_d.f90 -o $@

$(OBS)ziggurat.o: $(SR_UTILS)ziggurat.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)ziggurat.f90 -o $@

$(OBS)resample.o: $(SR_OPERS)resample.f90 $(OBS)pf_control.o $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)resample.f90 -o $@

$(OBS)diagnostics.o: $(SR_UTILS)diagnostics.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)diagnostics.f90 -o $@

$(OBS)quicksort.o: $(SR_UTILS)quicksort.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)quicksort.f90 -o $@

$(OBS)empire_main.o: $(SR_CONTS)empire_main.f90 $(OBS)comms.o $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)empire_main.f90 -o $@

$(OBS)letks_test.o: $(SR_CONTS)letks_test.f90 $(OBS)comms.o $(OBS)pf_control.o $(OBS)letks.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)letks_test.f90 -o $@

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

$(OBS)letks.o: $(SR_SMOOT)letks.f90
	$(FC) $(FCOPTS) -c $(SR_SMOOT)letks.f90 -o $@


#model as subroutine bits
$(OBS)model_as_subroutine_data.o: $(SR_USER_MDL)model_as_subroutine_data.f90
	$(FC) $(FCOPTS) -c $(SR_USER_MDL)model_as_subroutine_data.f90 -o $@

$(OBS)model_as_subroutine_initialise.o: $(SR_USER_MDL)model_as_subroutine_initialise.f90
	$(FC) $(FCOPTS) -c $(SR_USER_MDL)model_as_subroutine_initialise.f90 -o $@

$(OBS)model_as_subroutine_start.o: $(SR_USER_MDL)model_as_subroutine_start.f90
	$(FC) $(FCOPTS) -c $(SR_USER_MDL)model_as_subroutine_start.f90 -o $@

$(OBS)model_as_subroutine_return.o: $(SR_USER_MDL)model_as_subroutine_return.f90
	$(FC) $(FCOPTS) -c $(SR_USER_MDL)model_as_subroutine_return.f90 -o $@





EMPIRE: $(OBJS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)empire $(OBJS) $(LIB_LIST)

OBJS_LETKS_TEST = $(shell echo $(OBJS) | sed 's/empire_main/letks_test/g') 
LETKS_TEST: $(OBJS_LETKS_TEST)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)letks_test $(OBJS_LETKS_TEST) $(LIB_LIST)


OBJS_ALLTEST = $(shell echo $(OBJS) | sed 's/empire_main/alltests/g') 
ALLTESTS: $(OBJS_ALLTEST)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)alltests $(OBJS_ALLTEST) $(LIB_LIST)

OBJS_TEST_H = $(shell echo $(OBJS) | sed 's/empire_main/test_h/g') 
TEST_H: $(OBJS_TEST_H)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_h $(OBJS_TEST_H) $(LIB_LIST)

OBJS_TEST_R = $(shell echo $(OBJS) | sed 's/empire_main/test_r/g') 
TEST_R: $(OBJS_TEST_R)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_r $(OBJS_TEST_R) $(LIB_LIST)

OBJS_TEST_Q = $(shell echo $(OBJS) | sed 's/empire_main/test_q/g') 
TEST_Q: $(OBJS_TEST_Q)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_q $(OBJS_TEST_Q) $(LIB_LIST)

OBJS_TEST_HQHTR = $(shell echo $(OBJS) | sed 's/empire_main/test_hqhtr/g') 
TEST_HQHTR: $(OBJS_TEST_HQHTR)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN)test_hqhtr $(OBJS_TEST_HQHTR) $(LIB_LIST)


docs: doc_latex

doc_html: FORCE
	sed -i "s/PROJECT_NUMBER         =.*/PROJECT_NUMBER         = `git describe --abbrev=0 --tags`/g" .Doxyfile
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
	rm -f obs/*

veryclean: clean
	rm -f bin/*

minimal: EMPIRE
	cd models/minimal_empire;make -e
	cd models/minimal_empire_comms;make -e
	cd models/minimal_model;make -e
	cd models/minimal_model_comms;make -e

v1:
	sed -i 's/comm_version=.*/comm_version=1/1' comm_version.f90

v2:
	sed -i 's/comm_version=.*/comm_version=2/1' comm_version.f90

v3:
	sed -i 's/comm_version=.*/comm_version=3/1' comm_version.f90

v4:
	sed -i 's/comm_version=.*/comm_version=4/1' comm_version.f90

v5:
	sed -i 's/comm_version=.*/comm_version=5/1' comm_version.f90

linear_identity: clean v1 linear
	cp model_specific.f90 model_specific.f90_backup
	cp examples/linear_identity/model_specific_linear_identity.f90 model_specific.f90
	-mv $(BIN)empire $(BIN)empire_temp
	make EMPIRE
	mv $(BIN)empire $(BIN)empire_linear_identity
	-mv $(BIN)empire_temp $(BIN)empire
	mv model_specific.f90_backup model_specific.f90
