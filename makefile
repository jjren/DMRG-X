MKLROOT=/opt/intel/intel2016/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64 
mklinc1=$(MKLROOT)/include
JDLIB=./lib
#MPILIB=/opt/intel/parallel_studio_xe_2016/impi/5.1.1.109/lib64
MPILIB=
MPIINC=
LANCLIB=/home/jjren/Code/planso.jl/external/PLAN/liblanso.a

FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm 

FC=mpiifort
FCCOMPILEOPTS= -g -debug -traceback 
#FCCOMPILEOPTS= -O3

.SUFFIXES: .f90 .f

%.o : %.f90
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
%.o : %.f
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
# 
object = BLASmodified.o kinds_mod.o communicate.o exit_mod.o variables.o SpMatTrans.o basisindex.o checkinfo.o \
	   ppp_term.o CreatFCIDUMP.o onesitematrix.o contructquanta.o mathlib.o module_sparse.o checkmem.o \
	   symmetry.o system_big.o InitialGuess.o coefftosparse.o op.o GetHDiag.o pre_perturbation.o\
	   davidson.o masterdiag.o  perturbation.o \
	   noise.o Renormalization.o splitsvd_direct.o\
	   hamiltonian.o infinit_MPS.o  \
	   loadbalance.o peierls.o readintegral.o  store_operator.o   \
	   enviro_big.o finit_MPS.o sweep.o \
	   selectstates.o excitedbasis.o \
	   meanfield.o C2_copy.o\
	   opexpec.o corrfunc.o bondord.o localspin.o transmoment.o  infinit_initmat.o free_DMRG.o free_program.o \
	   lanczos.o opc.o dyn_prop.o analysis.o readinput.o main.o  \

DMRG-X : $(object)
	$(FC) -g -traceback -o $@ $^ $(LANCLIB) -L$(JDLIB) -ljadamilu -L$(MPILIB) -L$(MKLLIB) $(FCCFLAG) 

clean:
	rm -f *.o *.mod DMRG-X

	

