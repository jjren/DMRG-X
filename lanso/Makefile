
FC = ifort
FFLAGS = -O4
NOOPT = -O1

AR=ar
RANLIB=ranlib

F77FLAGS = ${FFLAGS}

MKLROOT=/opt/intel/intel2016/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64 
mklinc1=$(MKLROOT)/include
LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm 
#
D4OBJ=	landr.o lanso.o purge.o ritvec.o startv.o
D4SRC=	landr.f lanso.f purge.f ritvec.f startv.f
#
# eXtra supporting objects:
#
#	MACHAR is written by W. C. Cody,
#	TQL2 and PYTHAG are part of EISPACK.
#
# Note that TQLB needs PYTHAG also.
#
X_OBJ=	machar.o tql2.o pythag.o newtql2.o 
X_SRC=	machar.f tql2.f pythag.f newtql2.f90 

#
# common objects
#
C_OBJ=	${X_OBJ} dsort2.o tqlb.o ortbnd.o random.o store.o stpone.o
C_SRC=	${X_SRC} dsort2.f tqlb.f ortbnd.f random.f store.f stpone.f

#
# library liblanso.a
#
lib: ${C_OBJ} ${D4OBJ} datx.o
	$(AR) ruv liblanso.a ${C_OBJ} ${D4OBJ} datx.o
	${RANLIB} liblanso.a
#

clean:
	rm -f *.o *.mod liblanso.a

#
machar.o: machar.f
	${FC} ${NOOPT} -c machar.f
.f.o:
	${FC} ${FFLAGS} -c $<
%.o : %.f90
	$(FC) -c $(FFLAGS) -I/$(mklinc) -I/$(mklinc1) $<
#
