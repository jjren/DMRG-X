FFTINLCUDE = /usr/include
FFTLIBS = /usr/lib/x86_64-linux-gnu/
FFTLIBFLAGS = -lfftw3 -lfftw3_threads -lm
MKLROOT=/opt/intel/intel2016/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64
mklinc1=$(MKLROOT)/include
MKLFLAGS= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm



F90 = ifort

.SUFFIXES: .f90

%.o : %.f90
	$(F90) -c -g -I/$(FFTINCLUDE) -I/$(mklinc) -I/$(mklinc1) $<

object= corrfuncdft.o

corrfuncdft : $(object)
	$(F90) -o $@ $^ -L$(FFTLIBS) $(FFTLIBFLAGS) -L$(MKLLIB) $(MKLFLAGS)

clean:
	rm -f *.o 
