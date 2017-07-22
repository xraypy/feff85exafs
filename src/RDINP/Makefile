

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
ARCHV   ?= .a
EXEXT   ?=

COPY    ?= cp -v

objects = ffsort.o iniall.o mkptz.o rdinp_l.o rdline.o setedg.o	wrtall.o wrtjsn.o
archives = ../json-fortran/libjsonfortran$(ARCHV)	\
../COMMON/libfeffcom$(ARCHV) ../PAR/libfeffpar$(ARCHV)	\
../MATH/libfeffmath$(ARCHV) ../JSON/libfeffjson$(ARCHV)

all:	otherdirs feff8l_rdinp$(EXEXT)

feff8l_rdinp$(EXEXT):	$(objects) $(archives)
	$(FORTRAN) $(FLINKARGS) $(objects) $(archives) -o $@

ffsort.o:	ffsort.f  allinp.h ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
iniall.o:	iniall.f  allinp.h                    ../HEADERS/dim.h
mkptz.o:	mkptz.f            ../HEADERS/const.h
rdinp_l.o:	rdinp_l.f allinp.h ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
wrtall.o:	wrtall.f  allinp.h                    ../HEADERS/dim.h                   ../HEADERS/parallel.h
wrtjsn.o:	wrtjsn.f  allinp.h                    ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.o feff8l_rdinp$(EXEXT)

install:
	$(COPY) feff8l_rdinp$(EXEXT) $(BINDIR)

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran$(ARCHV)
	$(MAKE) -C ../COMMON       libfeffcom$(ARCHV)
	$(MAKE) -C ../PAR          libfeffpar$(ARCHV)
	$(MAKE) -C ../MATH         libfeffmath$(ARCHV)
	$(MAKE) -C ../JSON         libfeffjson$(ARCHV)

.PHONEY:	clean install all otherdirs
