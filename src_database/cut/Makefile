#compiles pes_fit and getpot program

iinput=rlx_tau.f90
output=a.x

#----------------------------------
FC = ifort -r8  -O
LIBS    = $(PES_LIB)  -mkl=sequential
PAATH= ~/co2h2o/4_fit_pes/14_v2b_pes/pnp_pes_shell/
FFLAGS  = $(PES_MOD) -I $(PAATH)
SRC=  \
$(PAATH)co2_coordinate.o \
$(PAATH)co2potlongrange.o \
$(PAATH)constants.o \
$(PAATH)pot_monomer_mod.o \
$(PAATH)pot_monomer.o \
$(PAATH)shell_water.o \
$(PAATH)pes_shell.o \
$(PAATH)xxralpha.o


%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<

OBJ=$(iinput:.f90=.o) $(SRC) 


$(output) : $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o $(output)
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod
