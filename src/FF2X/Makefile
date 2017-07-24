

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM	?= rm

COPY    ?= cp -v
ARCHV   ?= .a

objects = exconv.o feffdt.o fdtarr.o fdthea.o ff2afs.o ff2chi.o ff2xmu.o	\
ffmod6.o fprime.o ff2gen.o rdfbin.o reff2x.o xscorr.o

archives = ../json-fortran/libjsonfortran$(ARCHV)	\
../DEBYE/libfeffdw$(ARCHV) ../COMMON/libfeffcom$(ARCHV)	\
../PAR/libfeffpar$(ARCHV) ../MATH/libfeffmath$(ARCHV)	\
../JSON/libfeffjson$(ARCHV)


all:	otherdirs feff8l_ff2x$(EXEXT)

feff8l_ff2x$(EXEXT):	$(objects) $(archives)
	$(FORTRAN) $(objects) $(archives) -o $@

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<

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
	$(RM) *.o feff8l_ff2x$(EXEXT)

install:
	$(COPY) feff8l_ff2x$(EXEXT) $(BINDIR)

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran$(ARCHV)
	$(MAKE) -C ../DEBYE        libfeffdw$(ARCHV)
	$(MAKE) -C ../COMMON       libfeffcom$(ARCHV)
	$(MAKE) -C ../PAR          libfeffpar$(ARCHV)
	$(MAKE) -C ../MATH         libfeffmath$(ARCHV)
	$(MAKE) -C ../JSON         libfeffjson$(ARCHV)

.PHONEY:	clean install all otherdirs
