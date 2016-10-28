FC = ifort -r8 -O
LIBS = -mkl=sequential

%.o: %.f90
	$(FC) -c $<

OBJFIT := bemsa.o fit.o

OBJGET := bemsa.o dbemsa.o pes_shell.o getpot.o

fit.x : $(OBJFIT) 
	$(FC) -o $@ $^ $(LIBS)

getpot.x : $(OBJGET)
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *.mod *.x	
