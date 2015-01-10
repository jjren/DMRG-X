MKLROOT=/opt/intel/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=/opt/intel/mkl/include/intel64/lp64 
mklinc1=/opt/intel/mkl/include

FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm

FC=mpif90

DMRG-X:$(object) 
	$(FC) -o DMRG-X $(object) -I$(mklinc) -I$(mklinc1)-L$(MKLLIB) $(FCCFLAG)
#op.o: 
	$(FC) -c op.f90 -I$(mklinc) -L$(MKLLIB) $(FCCFLAG)
#Renormalization.o: 
#	$(FC) -c Renormalization.f90 -I$(mklinc) -L$(MKLLIB) $(FCCFLAG)
#InitialGuess.o: 
#	$(FC) -c InitialGuess.f90 -I$(mklinc) -L$(MKLLIB) $(FCCFLAG)
fullmat.o: 
	$(FC) -c fullmat.f90 -I$(mklinc) -L$(MKLLIB) $(FCCFLAG)

object = contructquanta.o davidson.o davidson_wrapper.o GetHDiag.o \
	   hamiltonian.o infinit_MPS.o infinit_smallL.o infinit_smallR.o \
	   InitialGuess.o loadbalance.o main.o \
	   mathlib.o onesitematrix.o op.o PPP_term.o readinput.o Renormalization.o \
	   spin_reversal.o store_operator.o system_bigL.o system_bigR.o variables.o \
	   fullmat.o halfmat.o
