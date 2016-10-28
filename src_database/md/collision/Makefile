FC       = ifort
FFLAGS   = -r8 -O
LIBS     = -mkl=sequential

all : dynam.x

dynam.x : pes_shell.o fmd_util.o dynam.o
	$(FC) $^ -o dynam.x $(FFLAGS) $(LIBS)

%.o : %.f90
	$(FC) -c $(FFLAGS) $<
