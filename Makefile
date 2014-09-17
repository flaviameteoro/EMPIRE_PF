FC=mpif90


all: pf_couple

FCOPTS = -O3 -fopenmp

LIB_LIST = -L$(METISDIR) -l$(METISLIB) -lblas


OBS=obs/
MODFLAG=-J
MODLOC:=$(OBS)
SR_FILTS=src/filters/
SR_UTILS=src/utils/
SR_CONTS=src/controlers/
SR_DATAS=src/data/
SR_TESTS=src/test/
SR_OPERS=src/operations/
OBJSQ= sizes.o pf_couple.o Qdata.o Rdata.o equivalent_weights_step.o comms.o gen_rand.o random_d.o proposal_filter.o pf_control.o data_io.o model_specific.o operator_wrappers.o quicksort.o resample.o diagnostics.o perturb_particle.o genQ.o sir_filter.o
OBJS=$(addprefix $(OBS),$(OBJSQ))
FCOPTS+=$(MODFLAG) $(MODLOC)

$(OBS)sir_filter.o: $(SR_FILTS)sir_filter.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)sir_filter.f90 -o $@

$(OBS)genQ.o: $(SR_UTILS)genQ.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)genQ.f90 -o $@

$(OBS)Qdata.o: $(SR_DATAS)Qdata.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_DATAS)Qdata.f90 -o $@

$(OBS)Rdata.o: $(SR_DATAS)Rdata.f90
	$(FC) $(FCOPTS) -c $(SR_DATAS)Rdata.f90 -o $@

$(OBS)model_specific.o: model_specific.f90 $(OBS)sizes.o 
	$(FC) $(FCOPTS) -c model_specific.f90 -o $@

$(OBS)operator_wrappers.o: $(SR_OPERS)operator_wrappers.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)operator_wrappers.f90 -o $@

$(OBS)pf_control.o: $(SR_CONTS)pf_control.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)pf_control.f90 -o $@

$(OBS)perturb_particle.o: $(SR_OPERS)perturb_particle.f90 $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)perturb_particle.f90 -o $@

$(OBS)data_io.o: $(SR_UTILS)data_io.f90 $(OBS)pf_control.o $(OBS)sizes.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)data_io.f90 -o $@

$(OBS)proposal_filter.o: $(SR_FILTS)proposal_filter.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)proposal_filter.f90 -o $@

$(OBS)sizes.o: $(SR_CONTS)sizes.f90
	$(FC) $(FCOPTS) -c $(SR_CONTS)sizes.f90 -o $@

$(OBS)comms.o: $(SR_UTILS)comms.f90 $(OBS)sizes.o $(OBS)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)comms.f90 -o $@

$(OBS)equivalent_weights_step.o: $(SR_FILTS)equivalent_weights_step.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)equivalent_weights_step.f90 -o $@

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

$(OBS)gen_rand.o: $(SR_OPERS)gen_rand.f90 $(OBS)random_d.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)gen_rand.f90 -o $@

pf_couple: $(OBJS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o pf_couple $(OBJS) $(LIB_LIST)


#OBJS2 = $(shell echo $(OBJS) | sed 's/pf_couple/getdata/g') 

#getdata: $(OBJS2)
#	$(FC) $(FCOPTS) $(LOADOPTS) -o getdata $(OBJS2) $(LIB_LIST)

clean:
	rm -f obs/* pf_couple

