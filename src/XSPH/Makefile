

PREFIX  ?= /usr/local


FORTRAN ?= gfortran
FCFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM	?= rm

COPY    ?= cp -v

objects = acoef.o ffmod2.o getedg.o phmesh.o rdgrid.o rholat.o szlz.o	\
wrxsph.o xsect.o axafs.o fmssz.o phase.o radint.o rexsph.o rholsz.o	\
wphase.o xmult.o xsph.o

archives = ../json-fortran/libjsonfortran.a ../ATOM/libfeffatom.a	\
../FOVRG/libfeffpha.a ../FMS/libfefffms.a ../POT/libfeffint.a		\
../EXCH/libfeffexch.a ../COMMON/libfeffcom.a ../PAR/libfeffpar.a	\
../MATH/libfeffmath.a ../JSON/libfeffjson.a

otherobj = ../FF2X/xscorr.o ../POT/grids.o ../MATH/terpc.o

all:	otherdirs xsph

xsph:	$(objects) $(archives) $(otherobj)
	$(FORTRAN) $(objects) $(archives) $(otherobj) -o $@

../FF2X/xscorr.o:	../FF2X/xscorr.f ../HEADERS/dim.h ../HEADERS/const.h
../POT/grids.o:		../POT/grids.f                    ../HEADERS/const.h
../MATH/terpc.o:	../MATH/terpc.f

acoef.o:	acoef.f   ../HEADERS/dim.h
axafs.o:	axafs.f   ../HEADERS/dim.h ../HEADERS/const.h
ffmod2.o:	ffmod2.f  ../HEADERS/dim.h                     ../HEADERS/parallel.h
fmssz.o:	fmssz.f   ../HEADERS/dim.h
getedg.o:	getedg.f  ../HEADERS/dim.h
phase.o:	phase.f   ../HEADERS/dim.h ../HEADERS/const.h
phmesh.o:	phmesh.f  ../HEADERS/dim.h ../HEADERS/const.h
phmesh2.o:	phmesh2.f ../HEADERS/dim.h ../HEADERS/const.h
radint.o:	radint.f  ../HEADERS/dim.h
rdgrid.o:	rdgrid.f                   ../HEADERS/const.h
rexsph.o:	rexsph.f  ../HEADERS/dim.h ../HEADERS/const.h
rholat.o:	rholat.f  ../HEADERS/dim.h ../HEADERS/const.h
rholsz.o:	rholsz.f  ../HEADERS/dim.h ../HEADERS/const.h
szlz.o:		szlz.f    ../HEADERS/dim.h ../HEADERS/const.h
wphase.o:	wphase.f  ../HEADERS/dim.h
wrxsph.o:	wrxsph.f  ../HEADERS/dim.h
xmult.o:	xmult.f                    ../HEADERS/const.h
xsect.o:	xsect.f   ../HEADERS/dim.h ../HEADERS/const.h
xsph.o:		xsph.f    ../HEADERS/dim.h ../HEADERS/const.h



%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.o xsph

install:
	$(COPY) xsph $(PREFIX)/bin

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran.a
	$(MAKE) -C ../ATOM         libfeffatom.a
	$(MAKE) -C ../FOVRG        libfeffpha.a
	$(MAKE) -C ../FMS          libfefffms.a
	$(MAKE) -C ../POT          libfeffint.a
	$(MAKE) -C ../EXCH         libfeffexch.a
	$(MAKE) -C ../COMMON       libfeffcom.a
	$(MAKE) -C ../PAR          libfeffpar.a
	$(MAKE) -C ../MATH         libfeffmath.a
	$(MAKE) -C ../JSON         libfeffjson.a

.PHONEY:	clean install all otherdirs
