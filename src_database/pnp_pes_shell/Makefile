#compiles pes_fit and getpot program
FC = ifort -r8  -O
FFLAGS  = $(PES_MOD)
LIBS    = $(PES_LIB)  -mkl=sequential


SRC=co2_coordinate.f90 \
co2potlongrange.f90 \
constants.f90 \
pot_monomer_mod.f90 \
pot_monomer.f90 \
shell_water.f90 \
pes_shell.f90 \
xxralpha.f90 \
update.f90 \

OBJ=$(SRC:.f90=.o)

update.x : $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o update.x
	echo 'update sucessfull'

%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
	rm -rf *.o *.mod
