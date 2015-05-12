MKLROOT=/export/home/mahe/intel/composer_xe_2013.2.146/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64 
mklinc1=$(MKLROOT)/include


FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm

FC=mpiifort

FCCOMPILEOPTS= -O2

.SUFFIXES: .f90 .f

%.o : %.f90
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
%.o : %.f
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<

object = kinds_mod.o communicate.o exit_mod.o variables.o \
	   ppp_term.o contructquanta.o mathlib.o davidson.o \
	   symmetry.o InitialGuess.o  GetHDiag.o \
	   hamiltonian.o infinit_MPS.o infinit_initmat.o \
       max_overlap.o\
	   davidson_wrapper.o\
	   meanfield.o\
	   loadbalance.o\
	   count.o\
	   main.o\
	   excitedbasis.o\
	   onesitematrix.o op.o readinput.o Renormalization.o \
	   store_operator.o system_big.o \
	   enviro_big.o finit_MPS.o sweep.o \
	   selectstates.o\
	   splitsvdL.o splitsvdR.o\
	   C2_copy.o\
	   transmoment.o fullmat.o\
	   

DMRG-X : $(object)
	$(FC) -o $@ $^ -L$(MKLLIB) $(FCCFLAG)

clean:
	rm -f *.o *.mod DMRG-X
