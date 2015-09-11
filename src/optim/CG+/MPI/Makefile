# makefile for CG code
F77=mpif90
F90=mpif90
FILES = main.o cgfam.o cgsearch.o timer.o fcn.o
MYFILES=call.o cgsub.o cgfam.o cgsearch.o timer.o fcn.o objective_function.o objective_gradient.o
LIB_LIST=-lblas
OPTS=-Wall -fimplicit-none -fbounds-check -fbacktrace

all: sub

# let's use the special macros $@ and $^, which are the
# left and right sides of the :, respectively, 
# to make the overall compilation rule more general
%.o: %.f
	$(F77) $(OPTS) -c -o $@ $<

%.o: %.f90
	$(F90) $(OPTS) -c -o $@ $<


main: $(FILES)
	$(F77) $(FILES) $(OPTS) -o main $(LIB_LIST)

sub: $(MYFILES)
	$(F90) $(MYFILES) $(OPTS) -o subroutine_cg $(LIB_LIST)

clean:
	rm -rf *.o
