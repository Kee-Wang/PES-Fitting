FC =ifort -r8  -O # gfortran -fdefault-real-8 -fdefault-double-8 #ifort -r8 -O
FFLAGS  = -I ../pes_shell/lib/mod
LIBS    = -L../pes_shell/lib/pes-xyz -lpes -lpx -mkl=sequential


SRC1=fit.f90
OBJ1=$(SRC1:.f90=.o)


1_fit.x : $(OBJ1)
	$(FC) $(FFLAGS) -o fit.x $(OBJ1) $(LIBS)
	rm -rf *.o *.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
	rm -rf *.o *.mod
