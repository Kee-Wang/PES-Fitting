FC := pgfortran
FFLAGS := -r8 -O


vpath %.f90 src

OBJ  := \
blas.vscf.3.4.o \
contr.vscf.5.1.1.o\
curve.vscf.5.1.4.o\
diag.vscf.5.2.o\
drive.vscf.5.2.new.o\
eckart.vscf.5.1.1.o\
lanczos.vscf.5.2.o\
lap.vscf.3.4.o\
lib.vscf.3.4.o\
matinv.vscf.2.3.o\
memo.vscf.5.1.3.o\
molpro.vscf.5.0.o\
normal.vscf.4.9.0.o\
projec.vscf.4.9.1.o\
prop.vscf.5.2.o\
psplot.o\
react.vscf.eckart.5.1.3.o\
rot.vscf.5.2.o\
spectrum.o\
timit.vscf.4.9.1.o\
user.CO2.o\
user.vscf.4.9.2.o\
vib.vscf.5.2.o

# Targets:
.PHONY : clean 

%.o : %.f
	$(FC) -c $(FFLAGS) $< 
 
all : mm.x

mm.x   :  $(OBJ)
	$(FC) -o $@  $^ $(FFLAGS)

clean :
	rm *.o
