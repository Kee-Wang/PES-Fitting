#compiles pes_fit program
FC = ifort -r8 -O
FFLAGS  = $(PES_MOD)
LIBS    = $(PES_LIB) -mkl=sequential

SRC=pes_shell.f90  opt.f90  xopt.f90
OBJ=$(SRC:.f90=.o)

xopt.x : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
	rm -rf *.o *.mod

