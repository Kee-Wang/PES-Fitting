FC       = ifort
FFLAGS   = -r8 -O
LIBS     = -mkl=sequential

all : dyn.x

dyn.x : pes_shell.o fmd_util.o dynam-oct09.o
	$(FC) $^ -o dyn.x $(FFLAGS) $(LIBS)

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

