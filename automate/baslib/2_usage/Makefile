#compiles pes_fit and getpot program

iinput=spt.f90
output=a.x

#----------------------------------
FC = ifort -r8  -O
LIBS    = $(PES_LIB)  -mkl=sequential
PAATH= ../pes_shell/
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
