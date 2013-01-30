#module swap PrgEnv-cray PrgEnv-pgi; module load xt-libnagfl
all: pf_couple

FC=ftn
FCOPTS=  -r8 -O2 -Kieee -fastsse
#FCOPTS= -r8 -I ../PrimitiveEquationPF 
LOAD=mpif90
LOADOPTS=

OBJS= pf_couple.o hadcm3_config.o nudge_data.o comms.o nag90.o extra.o statistics.o fullQ.o analyse.o minresModule.o minresDataModule.o enssmm.o 

#nag90_nagVersion.o

hadcm3_config.o: hadcm3_config.f90 
	$(FC) $(FCOPTS) -c hadcm3_config.f90 

extra.o: extra.f90 hadcm3_config.o
	$(FC) $(FCOPTS) -c extra.f90 

comms.o: comms.f90 extra.o
	$(FC) $(FCOPTS) -c comms.f90 

enssmm.o: enssmm.f90 extra.o
	$(FC) $(FCOPTS) -c enssmm.f90 

nudge_data.o: nudge_data.f90
	$(FC) $(FCOPTS) -c nudge_data.f90 

#nag90_nagVersion.o: nag90_nagVersion.f90
#	$(FC) $(FCOPTS) -c nag90_nagVersion.f90 

nag90.o: nag90.f90
	$(FC) $(FCOPTS) -c nag90.f90 

analyse.o: analyse.f90 fullQ.o extra.o
	$(FC) $(FCOPTS) -c analyse.f90 

fullQ.o: fullQ.f90 extra.o minresModule.o
	$(FC) $(FCOPTS) -c fullQ.f90 

minresModule.o: minresModule.f90 minresDataModule.o
	$(FC) $(FCOPTS) -c minresModule.f90

minresDataModule.o: minresDataModule.f90
	$(FC) $(FCOPTS) -c minresDataModule.f90

statistics.o: statistics.f90 nag90.o extra.o
	$(FC) $(FCOPTS) -c statistics.f90 

pf_couple.o: pf_couple.f90 comms.o
	$(FC) $(FCOPTS) -c pf_couple.f90 

pf_couple: $(OBJS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o pf_couple $(OBJS)



clean:
	rm -f *.o *.mod pf_couple
