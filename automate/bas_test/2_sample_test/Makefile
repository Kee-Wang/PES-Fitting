#Please enter your input files here:
INPUT =spt.f90#../pes_shell/test/update.o# spt.f90

#Name of the excutable:
EXE = spt.x

#Please enter the relative path from excutable file and 'pes_shell':
PATH1 = ../pes_shell/#Do not leave whitespace after `/`

#==========SRC code for test==========
SRC=\
$(PATH1)co2_coordinate.o \
$(PATH1)co2potlongrange.o \
$(PATH1)h2o_constants.o \
$(PATH1)h2o_pot_monomer_mod.o \
$(PATH1)h2o_pot_monomer.o \
$(PATH1)h2o_shell_water.o \
$(PATH1)pes_shell.o \
$(PATH1)xxralpha.o
#=====================================

#Compiler and flags
FC = ifort -r8  -O
#Mod and flags
MOD  = -I $(PATH1)lib/mod -I $(PATH1) 
#Lib and flags
LIBS  = -L$(PATH1)lib/pes-xyz -lpes -lpx  -mkl=sequential

#Compile .f90 into .o file according to the rule.
OBJ1=$(INPUT:.f90=.o) 
OBJ2=$(SRC) #$(SRC:.f90=.o)
OBJ = $(OBJ2) $(OBJ1) #Combine two sets of objects into one OBJ

#The rule for conversion
%.o : %.f90
	$(FC) -c  $(MOD) $(LIBS) $<

#The rule to link all objs and generate executable. Default command `make`
$(EXE) :  $(OBJ)
	$(FC) $(MOD) $(OBJ) $(LIBS) -o $(EXE)


#Clear all process files with command `make clean`
#Notice, DO NOT clean mod files in `pes_shell` bececause it will be called from other Makefiles.
clean:
	rm -rf *.o *.mod
