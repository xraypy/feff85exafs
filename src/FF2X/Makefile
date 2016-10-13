

FORTRAN = gfortran
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib
RM	= rm

objects = exconv.o feffdt.o fdtarr.o fdthea.o ff2afs.o ff2chi.o ff2xmu.o	\
ffmod6.o fprime.o ff2gen.o rdfbin.o reff2x.o xscorr.o

archives = ../json-fortran/libjsonfortran.a ../DEBYE/libfeffdw.a	\
../COMMON/libfeffcom.a ../PAR/libfeffpar.a ../MATH/libfeffmath.a	\
../JSON/libfeffjson.a


all:	ff2x

ff2x:	$(objects) $(archives)
	$(FORTRAN) $(objects) $(archives) -o $@

%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.o ff2x

.PHONEY:	clean all
