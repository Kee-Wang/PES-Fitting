
FC := ifort 
FFLAGS :=  -static $(PES_MOD)
OPTS := -r8 -O 
LIBS   := $(PES_LIB)  -mkl=sequential


OBJ := co2_coordinate.o co2potlongrange.o constants.o pot_monomer_mod.o pot_monomer.o shell_water.o pes_shell.o user.CO2H2O.o  memo.vscf.4.9.EMORY.o  user.vscf.4.9.0.o \
react.vscf.eckart.o\
eckart.vscf.4.9.0.o contr.vscf.4.9.0.o lanczos.vscf.4.9.0.o molpro.vscf.4.9.0.o vib.vscf.4.9.0.o \
curve.vscf.4.9.0.o  normal.vscf.4.9.0.o diag.vscf.4.9.0.o  projec.vscf.4.9.0.o  rot.vscf.4.9.0.o \
drive.vscf.4.9.0.o  prop.vscf.4.9.0.o \
blas.vscf.3.4.o  lib.vscf.3.4.o lap.vscf.3.4.o \
matinv.vscf.2.3.o  timit.vscf.2.3.o

%.o : %.f
	$(FC) -c $< $(OPTS) $(FFLAGS)

%.o : %.f90
	$(FC) -c $< $(OPTS) $(FFLAGS)

#%.o : %.o
#	@$(FC) -c $< $(OPTS)

.PHONY : clean

all : mm.x

mm.x : $(OBJ) 
	@$(FC) $^ -o $@ $(FFLAGS) $(OPTS) $(LIBDIR) $(LIB) $(LIBS) 

clean :
	-rm -f $(OBJ) 
