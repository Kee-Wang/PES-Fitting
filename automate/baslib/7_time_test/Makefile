
#=====================User define area===================
input1=opt.f90 1_grad_num.f90 
input2=2_pes.f90
output1=1_g_n.x
output2=2_p.x
#========================================================

FC = ifort -r8 -O
LIBS = $(PES_LIB) -mkl=sequential
PAATH= ../pes_shell/
MMOD= $(PES_MOD)  -I $(PAATH)

SRC1 =\
$(PAATH)co2_coordinate.o \
$(PAATH)co2potlongrange.o \
$(PAATH)constants.o \
$(PAATH)pot_monomer_mod.o \
$(PAATH)pot_monomer.o \
$(PAATH)shell_water.o \
$(PAATH)xxralpha.o\
$(PAATH)pes_shell.o


OBJECT1 := $(SRC1) $(input1:.f90=.o) 
OBJECT2 := $(SRC1) $(input2:.f90=.o)


$(output1) : $(OBJECT1)
	$(FC) -o $@ $^ $(LIBS)
	rm -f *.o *.mod 

$(output2) : $(OBJECT2)
	$(FC) -o $@ $^ $(LIBS) 
	rm -f *.o *.mod 

%.o: %.f90
	$(FC) -c $< $(MMOD)
clean:
	rm -f *.o *.mod *.x	
