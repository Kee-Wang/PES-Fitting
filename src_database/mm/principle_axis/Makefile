FC := ifort -r8 
#parallel 
LIBS := -static -mkl=sequential
FFLAGS := -O

vpath %.f90 src

OBJ := \
util.o\
pa_main.o

# Targets:
.PHONY : clean 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 
%.o : %.f
	$(FC) -c $(FFLAGS) $< 

all : pa.x 

pa.x :  $(OBJ) 
	$(FC) -o $@  $^ $(FFLAGS) $(LIBS)


clean :
	rm *.o *.mod 
