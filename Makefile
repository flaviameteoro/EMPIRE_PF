#module swap PrgEnv-cray PrgEnv-pgi; module load xt-libnagfl
all: pf_couple

#detect the programming environment and set the compiler options accordingly
ProEnv = $(shell `eval `/opt/modules/3.2.6.6/bin/modulecmd bash list`` 2>&1 | grep 'PrgEnv' | cut -d - -f2 | cut -d / -f1)

ifeq ($(ProEnv),cray)
FCOPTS = -s real64
else
ifeq ($(ProEnv),pgi)
FCOPTS=  -r8 -O2 -Kieee -fastsse -Mbounds -traceback -Ktrap=fp 
FCOPTS = -fast -Mipa=fast,inline -mp -Mprefetch -O3
else
ifeq ($(ProEnv),gnu)
#FCOPTS = -fimplicit-none -Wall -fbounds-check -fopenmp 
FCOPTS = -O3 -fopenmp
else
$(error Programming environment not correct. Detected as: $(ProEnv))
endif
endif
endif



FC=ftn
#FCOPTS=  -r8 -O2 -Kieee -fastsse
#FCOPTS= -r8 -I ../PrimitiveEquationPF 
LOAD=mpif90
LOADOPTS=
EA20DIR=/home/n02/n02/pbrowne/hsl/ea20/hsl_ea20-1.0.0/src/

#OBJS= pf_couple.o hadcm3_config.o nudge_data.o equal_weight_filter.o quicksort.o comms.o gen_rand.o random_d.o extra.o proposal_filter.o pf_control.o data_io.o model_specific.o operator_wrappers.o
OBJS= common90.o ddeps90.o hsl_ma87d.o extra.o hadcm3_config.o pf_couple.o Qdata.o Rdata.o hsl_ea20d_new.o equal_weight_filter.o comms.o gen_rand.o random_d.o proposal_filter.o pf_control.o data_io.o model_specific.o operator_wrappers.o quicksort.o resample.o diagnostics.o perturb_particle.o genQ.o


MA87_LOC = /home/n02/n02/pbrowne/hsl/ma87/hsl_ma87-2.1.0/src/
METISDIR = /home/n02/n02/pbrowne/metis/metis-4.0.3/
METISLIB = metis
LIB_LIST = -L$(METISDIR) -l$(METISLIB)



common90.o: $(MA87_LOC)common90.f90
	$(FC) $(FCOPTS) -c $(MA87_LOC)common90.f90

ddeps90.o: $(MA87_LOC)ddeps90.f90
	$(FC) $(FCOPTS) -c $(MA87_LOC)ddeps90.f90

hsl_ma87d.o: $(MA87_LOC)hsl_ma87d.f90 ddeps90.o common90.o
	$(FC) $(FCOPTS) -c $(MA87_LOC)hsl_ma87d.f90

genQ.o: genQ.f90
	$(FC) $(FCOPTS) -c genQ.f90

Qdata.o: Qdata.f90 extra.o
	$(FC) $(FCOPTS) -c Qdata.f90

Qevd.o: Qevd.f90 extra.o
	$(FC) $(FCOPTS) -c Qevd.f90

Rdata.o: Rdata.f90
	$(FC) $(FCOPTS) -c Rdata.f90

hsl_ea20d_new.o: $(EA20DIR)hsl_ea20d_new.f90
	$(FC) $(FCOPTS) -c $(EA20DIR)hsl_ea20d_new.f90

hadcm3_config.o: hadcm3_config.f90 
	$(FC) $(FCOPTS) -c hadcm3_config.f90 

model_specific.o: model_specific.f90 extra.o hadcm3_config.o
	$(FC) $(FCOPTS) -c model_specific.f90 

operator_wrappers.o: operator_wrappers.f90 pf_control.o extra.o
	$(FC) $(FCOPTS) -c operator_wrappers.f90 

pf_control.o: pf_control.f90 extra.o
	$(FC) $(FCOPTS) -c pf_control.f90

perturb_particle.o: perturb_particle.f90 extra.o
	$(FC) $(FCOPTS) -c perturb_particle.f90

cake:
	cake

data_io.o: data_io.f90 pf_control.o extra.o
	$(FC) $(FCOPTS) -c data_io.f90

proposal_filter.o: proposal_filter.f90 random_d.o
	$(FC) $(FCOPTS) -c proposal_filter.f90

extra.o: extra.f90
	$(FC) $(FCOPTS) -c extra.f90 

comms.o: comms.f90 extra.o pf_control.o
	$(FC) $(FCOPTS) -c comms.f90 

equal_weight_filter.o: equal_weight_filter.f90 random_d.o
	$(FC) $(FCOPTS) -c equal_weight_filter.f90 

random_d.o: random_d.f90
	$(FC) $(FCOPTS) -c random_d.f90 

resample.o: resample.f90 pf_control.o random_d.o
	$(FC) $(FCOPTS) -c resample.f90

diagnostics.o: diagnostics.f90 pf_control.o extra.o
	$(FC) $(FCOPTS) -c diagnostics.f90

quicksort.o: quicksort.f90
	$(FC) $(FCOPTS) -c quicksort.f90

pf_couple.o: pf_couple.f90 comms.o pf_control.o
	$(FC) $(FCOPTS) -c pf_couple.f90 

test_Q.o: test_Q.f90
	$(FC) $(FCOPTS) -c test_Q.f90

gen_rand.o: gen_rand.f90 random_d.o
	$(FC) $(FCOPTS) -c gen_rand.f90


pf_couple: $(OBJS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o pf_couple $(OBJS) $(LIB_LIST)


#OBJS2 = $(shell echo $(OBJS) | sed 's/pf_couple/getdata/g') 

#getdata: $(OBJS2)
#	$(FC) $(FCOPTS) $(LOADOPTS) -o getdata $(OBJS2) $(LIB_LIST)

clean:
	rm -f *.o *.oo *.mod pf_couple

