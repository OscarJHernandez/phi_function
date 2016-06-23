# This is the makefile program for compiling the Eigen Value solver

# The compiler
FC = gfortran -g -fbounds-check -O2
#FC = ifort -g -parallel -traceback -O1 -heap-arrays
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
FCFLAGS += -Wall -fbounds-check

# libraries needed for linking
LDFLAGS = -llapack

# List of executables to be built within the package
PROGRAMS =  quadrature.f special_functions.f90 coulcc36-f90.f bessel.f90 phiFuncParams.f90 phiFunction.f90 matrixElements.f90 quadpack.f90 qDependent_PionInFlight_MEC.f90 qDependent_Seagull_MEC.f90

# The List of objects that will be built
_OBJS=  quadrature.o special_functions.o coulcc36-f90.o bessel.o phiFuncParams.o phiFunction.o matrixElements.o quadpack.o qDependent_PionInFlight_MEC.o qDependent_Seagull_MEC.o
OBJS = $(patsubst %,Objects/%,$(_OBJS))
#PROG= gen_H_FB
PROG= test

# Instructions for building all Objects needed:
Objects/%.o: Programs/%.f90
	$(FC)  -o  $@  -c  $<

Objects/%.o: Programs/%.f
	$(FC)  -o  $@  -c  $<

# Instructions for building the final executable program
$(PROG): $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o $(PROG) $(PROG).f90 $(LDFLAGS)

clean:
	rm -f Objects/*.o *.mod *.MOD
	rm -f Programs/*.o *.mod *.MOD
	
