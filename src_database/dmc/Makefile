#Generic makefile for fortran90
FC := ifort -r8 
LIBS := $(PES_LIB) -mkl=sequential
FFLAGS := -O $(PES_MOD)

vpath %.f90 src

OBJ := \
pes_shell.o

OBJDMC := \
dmc_proc.o  dmc.o

# Targets:
.PHONY : clean 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 
 
dmc.x  :  $(OBJ) $(OBJDMC)
	$(FC) -o $@  $^ $(FFLAGS) $(LIBS)
clean :
	rm *.o *.mod 
