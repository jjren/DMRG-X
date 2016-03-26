MKLROOT=/opt/intel/intel2016/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64 
mklinc1=$(MKLROOT)/include
JDLIB=./lib
#MPILIB=/opt/intel/parallel_studio_xe_2016/impi/5.1.1.109/lib64
MPILIB=
MPIINC=

FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm

FC=mpiifort
FCCOMPILEOPTS= -g -debug
#FCCOMPILEOPTS= -O3

.SUFFIXES: .f90 .f

%.o : %.f90
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
%.o : %.f
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
# 
object = kinds_mod.o communicate.o exit_mod.o variables.o SpMatTrans.o basisindex.o checkinfo.o \
	   ppp_term.o CreatFCIDUMP.o onesitematrix.o contructquanta.o mathlib.o module_sparse.o checkmem.o \
	   symmetry.o system_big.o InitialGuess.o coefftosparse.o GetHDiag.o pre_perturbation.o\
	   op.o davidson.o masterdiag.o  perturbation.o \
	   noise.o Renormalization.o splitsvd_direct.o\
	   hamiltonian.o infinit_MPS.o  \
	   loadbalance.o readinput.o store_operator.o   \
	   enviro_big.o finit_MPS.o checkmat.o sweep.o \
	   selectstates.o excitedbasis.o \
	   meanfield.o C2_copy.o\
	   transmoment.o bondord.o localspinnew.o analysis.o infinit_initmat.o count.o free_DMRG.o main.o  \

DMRG-X : $(object)
	$(FC) -o $@ $^ -L$(JDLIB) -ljadamilu -L$(MPILIB) -L$(MKLLIB) $(FCCFLAG)

clean:
	rm -f *.o *.mod DMRG-X

	

