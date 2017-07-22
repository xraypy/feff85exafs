
PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM	?= rm
ARCHV   ?= .a
EXEXT   ?=

COPY    ?= cp -v

objects = acoef.o ffmod2.o getedg.o phmesh.o rdgrid.o rholat.o szlz.o	\
wrxsph.o xsect.o axafs.o fmssz.o phase.o radint.o rexsph.o rholsz.o	\
wphase.o xmult.o xsph.o

archives = ../json-fortran/libjsonfortran$(ARCHV)		\
../ATOM/libfeffatom$(ARCHV) ../FOVRG/libfeffpha$(ARCHV)		\
../FMS/libfefffms$(ARCHV) ../POT/libfeffint$(ARCHV)		\
../EXCH/libfeffexch$(ARCHV) ../COMMON/libfeffcom$(ARCHV)	\
../PAR/libfeffpar$(ARCHV) ../MATH/libfeffmath$(ARCHV)		\
../JSON/libfeffjson$(ARCHV)

otherobj = ../FF2X/xscorr.o ../POT/grids.o ../MATH/terpc.o

all:	otherdirs feff8l_xsph$(EXEXT)

feff8l_xsph$(EXEXT):	$(objects) $(archives) $(otherobj)
	$(FORTRAN) $(FLINKARGS) $(objects) $(archives) $(otherobj) -o $@

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
	$(RM) *.o feff8l_xsph$(EXEXT)

install:
	$(COPY) feff8l_xsph$(EXEXT) $(BINDIR)

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran$(ARCHV)
	$(MAKE) -C ../ATOM         libfeffatom$(ARCHV)
	$(MAKE) -C ../FOVRG        libfeffpha$(ARCHV)
	$(MAKE) -C ../FMS          libfefffms$(ARCHV)
	$(MAKE) -C ../POT          libfeffint$(ARCHV)
	$(MAKE) -C ../EXCH         libfeffexch$(ARCHV)
	$(MAKE) -C ../COMMON       libfeffcom$(ARCHV)
	$(MAKE) -C ../PAR          libfeffpar$(ARCHV)
	$(MAKE) -C ../MATH         libfeffmath$(ARCHV)
	$(MAKE) -C ../JSON         libfeffjson$(ARCHV)

.PHONEY:	clean install all otherdirs
