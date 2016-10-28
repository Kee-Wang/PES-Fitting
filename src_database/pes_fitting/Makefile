#compiles pes_fit and getpot program
FC = ifort -r8 -O
FFLAGS  = $(PES_MOD)
LIBS    = $(PES_LIB) -mkl=sequential

#SRC1=fit.f90
#OBJ1=$(SRC1:.f90=.o)

SRC2=pes_shell.f90 getpot.f90
OBJ2=$(SRC2:.f90=.o)

SRC3=pes_shell.f90 readE.f90
OBJ3=$(SRC3:.f90=.o)

#fit.x : $(OBJ1)
#	$(FC) $(FFLAGS) -o fit.x $(OBJ1) $(LIBS)

getpot.x : $(OBJ2)
	$(FC) $(FFLAGS) -o getpot.x $(OBJ2) $(LIBS)

readE.x : $(OBJ3)
	$(FC) $(FFLAGS) -o  readE.x $(OBJ3) $(LIBS)

#%.o : %.f90
#	$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
	rm -rf *.o *.mod
