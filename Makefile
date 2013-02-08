#module swap PrgEnv-cray PrgEnv-pgi; module load xt-libnagfl
all: pf_couple

#detect the programming environment and set the compiler options accordingly
ProEnv = $(shell `eval `/opt/modules/3.2.6.6/bin/modulecmd bash list`` 2>&1 | grep 'PrgEnv' | cut -d - -f2 | cut -d / -f1)

ifeq ($(ProEnv),cray)
FCOPTS = -d I -m 0 -R b -h omp
else
ifeq ($(ProEnv),pgi)
FCOPTS=  -r8 -O2 -Kieee -fastsse -Mbounds -traceback -Ktrap=fp 
else
ifeq ($(ProEnv),gnu)
FCOPTS = -fimplicit-none -Wall -fbounds-check #-fopenmp 
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

OBJS= pf_couple.o hadcm3_config.o nudge_data.o equal_weight_filter.o kb05d.o comms.o gen_rand.o random_d.o extra.o proposal_filter.o pf_control.o data_io.o model_specific.o operator_wrappers.o

#nag90_nagVersion.o

hadcm3_config.o: hadcm3_config.f90 
	$(FC) $(FCOPTS) -c hadcm3_config.f90 

model_specific.o: model_specific.f90 
	$(FC) $(FCOPTS) -c model_specific.f90 

operator_wrappers.o: operator_wrappers.f90 
	$(FC) $(FCOPTS) -c operator_wrappers.f90 

pf_control.o: pf_control.f90
	$(FC) $(FCOPTS) -c pf_control.f90

data_io.o: data_io.f90
	$(FC) $(FCOPTS) -c data_io.f90

proposal_filter.o: proposal_filter.f90 random_d.o
	$(FC) $(FCOPTS) -c proposal_filter.f90

extra.o: extra.f90 hadcm3_config.o
	$(FC) $(FCOPTS) -c extra.f90 

comms.o: comms.f90 extra.o
	$(FC) $(FCOPTS) -c comms.f90 

enssmm.o: enssmm.f90 extra.o
	$(FC) $(FCOPTS) -c enssmm.f90 

equal_weight_filter.o: equal_weight_filter.f90 random_d.o
	$(FC) $(FCOPTS) -c equal_weight_filter.f90 

nudge_data.o: nudge_data.f90
	$(FC) $(FCOPTS) -c nudge_data.f90 

random_d.o: random_d.f90
	$(FC) $(FCOPTS) -c random_d.f90 

#nag90_nagVersion.o: nag90_nagVersion.f90
#	$(FC) $(FCOPTS) -c nag90_nagVersion.f90 

#nag90.o: nag90.f90
#	$(FC) $(FCOPTS) -c nag90.f90 

analyse.o: analyse.f90 extra.o
	$(FC) $(FCOPTS) -c analyse.f90 

#fullQ.o: fullQ.f90 extra.o minresModule.o
#	$(FC) $(FCOPTS) -c fullQ.f90 

minresModule.o: minresModule.f90 minresDataModule.o
	$(FC) $(FCOPTS) -c minresModule.f90

minresDataModule.o: minresDataModule.f90
	$(FC) $(FCOPTS) -c minresDataModule.f90

statistics.o: statistics.f90 extra.o
	$(FC) $(FCOPTS) -c statistics.f90 

pf_couple.o: pf_couple.f90 comms.o pf_control.o
	$(FC) $(FCOPTS) -c pf_couple.f90 

gen_rand.o: gen_rand.f90 random_d.o
	$(FC) $(FCOPTS) -c gen_rand.f90

kb05d.o: /home/n02/n02/pbrowne/hsl/kb05/kb05d-1.0.0/kb05d.f
	$(FC) $(FCOPTS) -c /home/n02/n02/pbrowne/hsl/kb05/kb05d-1.0.0/kb05d.f



pf_couple: $(OBJS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o pf_couple $(OBJS)



clean:
	rm -f *.o *.mod pf_couple
