FC       = ifort
FFLAGS   = -r8 -O -static
LIBS     = -mkl=sequential

all : dmd.x

dmd.x : pes_util.o  pes_shell.o  fmd_util.o  xdmd.o
	$(FC) $^ -o dmd.x $(FFLAGS) $(LIBS)

%.o : %.f90
	$(FC) -c $(FFLAGS) $<
