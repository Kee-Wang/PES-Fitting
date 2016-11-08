SYSTEM := $(shell uname)
ifeq ($(SYSTEM),SunOS)
 FC = /development/WS12/SUNWspro/bin/f90
 FFLAGS = -u -fast -fsimple=1 -xarch=native64 -xlic_lib=sunperf \
   -moddir=$(HOME)/mod
 BLD = /usr/xpg4/bin/ar rcv
endif
ifeq ($(SYSTEM),Linux)
 FC = ifort
# FFLAGS = -O -ip-no-inlining -u -module ../mod -vec-report0
 FFLAGS = -O -module ../mod #  -openmp -parallel
 BLD = /usr/bin/ar -rcv
# with earlier versions of ifort, may have to do this:
# FFLAGS = -O -Wl,-m elf_i386 -static-libcxa -module $(HOME)/mod
# BLD = /usr/bin/ar -rcv
# For the Pathscale compiler, do this:
# FC = pathf90
# FFLAGS = $(ACML) -module $(HOME)/mod
# BLD = /usr/bin/ar -rcv
endif
ifeq ($(SYSTEM),OSF1)
 FC = f90
 FFLAGS = -u -module $(HOME)/mod
 BLD = ar rcv
endif

# CAUTION: We have specified very few of the dependencies.  It is
# advised to do `make realclean ; make lib' in order to re-make the
# libraries.  A proper Makefile is yet to be written.

# Set objlib<foo> to the list of files to be included in lib<foo>.a
objlibpes = pes0.o \
  pes_x1.o pes_y1.o pes_z1.o pes_u1.o pes_x2.o pes_x1y1.o pes_y2.o \
  pes_x1z1.o pes_y1z1.o pes_z2.o pes_x1u1.o pes_y1u1.o pes_z1u1.o \
  pes_x3.o pes_x2y1.o pes_x1y2.o pes_y3.o pes_x2z1.o pes_x1y1z1.o \
  pes_y2z1.o pes_x1z2.o pes_y1z2.o pes_z3.o \
  pes_x2u1.o pes_x1y1u1.o pes_y2u1.o pes_x1z1u1.o pes_y1z1u1.o \
  pes_x4.o pes_x3y1.o pes_x2y2.o pes_x1y3.o \
  pes_x3z1.o pes_x2y1z1.o pes_x1y2z1.o pes_y3z1.o \
  pes_x2z2.o pes_x1y1z2.o pes_y2z2.o pes_x1z3.o pes_y1z3.o \
  pes_x3u1.o pes_x2y1u1.o pes_x1y2u1.o \
  pes_x2z1u1.o pes_x1y1z1u1.o pes_y2z1u1.o \
  pes_x5.o pes_x4y1.o pes_x3y2.o pes_x2y3.o pes_x4z1.o \
  pes_x3y1z1.o pes_x2y2z1.o pes_x1y3z1.o \
  pes_x3z2.o pes_x2y1z2.o pes_x1y2z2.o pes_y3z2.o \
  pes_x2z3.o pes_x1y1z3.o pes_y2z3.o \
  pes_x4u1.o pes_x3y1u1.o pes_x2y2u1.o \
  pes_x3z1u1.o pes_x2y1z1u1.o pes_x1y2z1u1.o \
  pes_x6.o pes_x5y1.o pes_x4y2.o pes_x3y3.o \
  pes_x5z1.o pes_x4y1z1.o pes_x3y2z1.o pes_x2y3z1.o \
  pes_x4z2.o pes_x3y1z2.o pes_x2y2z2.o pes_x1y3z2.o \
  pes_x3z3.o pes_x2y1z3.o pes_x1y2z3.o pes_y3z3.o \
  pes_x5u1.o pes_x4y1u1.o pes_x3y2u1.o pes_x4z1u1.o pes_x3y1z1u1.o \
  pes_x2y2z1u1.o \
  pes_x7.o pes_x6y1.o pes_x5y2.o pes_x4y3.o \
  pes_x5y1z1.o pes_x4y2z1.o pes_x3y3z1.o \
  pes_x4y1z2.o pes_x3y2z2.o pes_x2y3z2.o \
  pes_x3y1z3.o pes_x2y2z3.o pes_x1y3z3.o \
  pes_x6u1.o pes_x5y1u1.o pes_x4y2u1.o \
  pes_x5z1u1.o pes_x4y1z1u1.o pes_x3y2z1u1.o \
  pes_x7y1.o pes_x6y2.o pes_x5y3.o pes_x4y4.o \
  pes_x5y2z1.o pes_x4y3z1.o \
  pes_x4y2z2.o pes_x3y3z2.o pes_x3y2z3.o pes_x2y3z3.o \
  pes_x6y1u1.o pes_x6z1u1.o pes_x5y1z1u1.o pes_x4y2z1u1.o \
  pes_x7y2.o pes_x6y3.o pes_x5y4.o pes_x5y3z1.o pes_x4y4z1.o pes_x3y3z3.o \
  pes_x7y3.o pes_x6y3z1.o pes_x4y4z2.o pes_x9y4.o \
  pes1_xyz.o pes1c_xyz.o pes1c.o \
  pes2.o pes.o
objlibpx = px.o

objects = $(objlibpes) $(objlibpx)
libraries = libpes.a libpx.a

#-----------------------------------------------------------------------
# Rules.

all : $(objects) $(libraries)
obj : $(objects)
lib : $(libraries)

.PHONY : all obj lib clean realclean
.SUFFIXES : $(nil)
FORCE :

# Pattern rule to make a %.o file from a %.f90 file.
%.o : %.f90
	$(FC) -c $(FFLAGS) $<

# Here are the rules to make the libraries.
libpes.a : $(objlibpes) ../inv/libinv.a
	\rm -f $@
	cp ../inv/libinv.a $(@F)
	$(BLD) $(@F) $(objlibpes)
libpx.a : $(objlibpx)
	\rm -f $@
	$(BLD) $(@F) $(objlibpx)

# targets 'clean' and 'realclean' clean up the directory.
clean : FORCE
	rm -f *.o *~ core junk
realclean : clean
	rm -f *.a *.mod
# DO NOT DELETE
