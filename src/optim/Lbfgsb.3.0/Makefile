FC = mpif90

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

DRIVER1_90 = driver1.f90
DRIVER2_90 = driver2.f90
DRIVER3_90 = driver3.f90

LBFGSB  = lbfgsb.f
LINPACK = linpack.f
BLAS    = -lblas
TIMER   = timer.f

LIB_LIST = $(BLAS)

all :  lbfgsb_90_1 lbfgsb_90_2 lbfgsb_90_3 sub

OBS = lbfgsb_sub.o lbfgs_sub.o lbfgsb.o linpack.o timer.o call.o objective_function.o objective_gradient.o fcn.o

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

sub: $(OBS)
	$(FC) $(OBS) $(FFLAGS) -o subroutine_lbfgs $(LIB_LIST)

lbfgsb_90_1 : $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(TIMER) -o x.lbfgsb_90_1 $(LIB_LIST)

lbfgsb_90_2 : $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(TIMER) -o x.lbfgsb_90_2 $(LIB_LIST)

lbfgsb_90_3 : $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(TIMER) -o x.lbfgsb_90_3 $(LIB_LIST)

clean:
	rm -rf *.o
