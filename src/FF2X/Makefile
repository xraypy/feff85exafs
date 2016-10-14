

FORTRAN ?= gfortran
CFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM	?= rm

COPY    ?= cp -v

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

fdtarr.o:	fdtarr.f ../HEADERS/const.h ../HEADERS/dim.h
fdthea.o:	fdthea.f ../HEADERS/const.h ../HEADERS/dim.h
fdtxdi.o:	fdtxdi.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
feffdt.o:	feffdt.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
ff2afs.o:	ff2afs.f ../HEADERS/const.h ../HEADERS/dim.h
ff2chi.o:	ff2chi.f ../HEADERS/const.h ../HEADERS/dim.h
ff2gen.o:	ff2gen.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
ff2xmu.o:	ff2xmu.f ../HEADERS/const.h ../HEADERS/dim.h
ffmod6.o:	ffmod6.f                    ../HEADERS/dim.h ../HEADERS/parallel.h
fprime.o:	fprime.f ../HEADERS/const.h ../HEADERS/dim.h
reff2x.o:	reff2x.f ../HEADERS/const.h ../HEADERS/dim.h
xscorr.o:	xscorr.f ../HEADERS/const.h ../HEADERS/dim.h


clean:
	$(RM) *.o ff2x

install:
	$(COPY) ff2x $(PREFIX)/bin

.PHONEY:	clean install all
