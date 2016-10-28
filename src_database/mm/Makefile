FC := ifort
FFLAGS := -r8 -O

vpath %.f90 src

OBJ  := \
user.H2CS.5.1.3.o \
memo.vscf.5.1.3.o  user.vscf.4.9.2.o  react.vscf.eckart.5.1.3.o \
eckart.vscf.5.1.1.o  contr.vscf.5.1.1.o  lanczos.vscf.5.1.4.o  molpro.vscf.5.0.o \
vib.vscf.5.1.4.new.o  curve.vscf.5.1.4.o  normal.vscf.4.9.0.o  diag.vscf.5.1.3.o \
projec.vscf.4.9.1.o  rot.vscf.5.1.4.o  drive.vscf.5.1.4.new.o  prop.vscf.5.1.3.o \
blas.vscf.3.4.o  lib.vscf.3.4.o  lap.vscf.3.4.o  matinv.vscf.2.3.o  psplot.o \
timit.vscf.4.9.1.o  spectrum.o

# Targets:
.PHONY : clean 

%.o : %.f
	$(FC) -c $(FFLAGS) $< 
 
all : mm.x

mm.x   :  $(OBJ)
	$(FC) -o $@  $^ $(FFLAGS)

clean :
	rm *.o
