MKLROOT=/opt/intel/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64 
mklinc1=$(MKLROOT)/include


FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm

FC=mpif90
FCCOMPILEOPTS= -g -debug

.SUFFIXES: .f90 .f

%.o : %.f90
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
%.o : %.f
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<

object = kinds_mod.o communicate.o exit_mod.o variables.o checkinfo.o \
	   ppp_term.o contructquanta.o mathlib.o module_sparse.o davidson.o \
	   symmetry.o InitialGuess.o excitedbasis.o Renormalization.o \
	   stateOverlap.o coefftosparse.o hamiltonian.o GetHDiag.o \
	   op.o onesitematrix.o  splitsvd_direct.o\
	   hamiltonian.o infinit_MPS.o  \
	   loadbalance.o  \
	     system_big.o readinput.o  \
	    store_operator.o   \
	   enviro_big.o finit_MPS.o sweep.o \
	   selectstates.o\
	    meanfield.o C2_copy.o\
	  transmoment.o bondord.o analysis.o infinit_initmat.o count.o main.o   \

DMRG-X : $(object)
	$(FC) -o $@ $^ -L$(MKLLIB) $(FCCFLAG)

clean:
	rm -f *.o *.mod DMRG-X
