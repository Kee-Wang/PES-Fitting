# Compile MULTIMODE.5.1.4 fortran77 codes

# missing react.vscf.eckart.5.1.3.f

FC := ifort
FFLAGS := -r8 -O -static $(PES_MOD)
BLD = /usr/bin/ar -rcv
LIBS   := $(PES_LIB)  -mkl=sequential

OBJ := co2_coordinate.o co2potlongrange.o constants.o pot_monomer_mod.o pot_monomer.o shell_water.o pes_shell.o user.CO2H2O.o\
 blas.vscf.3.4.o\
 contr.vscf.5.1.1.o\
 curve.vscf.5.1.4.MVF.o\
 diag.vscf.5.1.3.o\
 eckart.vscf.5.1.1.o\
 lanczos.vscf.5.1.4.MVF.o\
 lap.vscf.3.4.o\
 lib.vscf.3.4.o\
 matinv.vscf.2.3.o\
 drive.vscf.5.1.4.MVF.f\
memo.vscf.5.1.3.o\
 molpro.vscf.5.0.o\
 normal.vscf.4.9.0.o\
 projec.vscf.4.9.1.o\
prop.vscf.5.1.3.F\
psplot.o\
 react.vscf.eckart.5.1.3.o\
 rot.vscf.5.1.4.MVF.o\
 spectrum.o\
 timit.vscf.4.9.1.o\
 user.vscf.4.9.2.o\
 vib.vscf.5.1.4.MVF.o

%.o : %.f
	$(FC) -c $< $(FFLAGS) $(OPTS)

%.o : %.f90
	$(FC) -c $< $(FFLAGS) $(OPTS)

.PHONY : clean

all : mm.x 

mm.x : $(OBJ)
	@$(FC) $^ -o $@ $(FFLAGS) $(OPTS) $(LIBDIR) $(LIB) $(LIBS)

clean :
	-rm -f $(OBJ) libmm.a 
