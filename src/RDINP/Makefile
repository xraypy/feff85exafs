

PREFIX  ?= /usr/local

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

COPY    ?= cp -v

objects = ffsort.o iniall.o mkptz.o rdinp_l.o rdline.o setedg.o	wrtall.o wrtjsn.o
archives = ../json-fortran/libjsonfortran.a ../COMMON/libfeffcom.a ../PAR/libfeffpar.a ../MATH/libfeffmath.a ../JSON/libfeffjson.a

all:	otherdirs rdinp

rdinp:	$(objects) $(archives)
	$(FORTRAN) $(objects) $(archives) -o $@

ffsort.o:	ffsort.f  allinp.h ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
iniall.o:	iniall.f  allinp.h                    ../HEADERS/dim.h
mkptz.o:	mkptz.f            ../HEADERS/const.h
rdinp_l.o:	rdinp_l.f allinp.h ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
wrtall.o:	wrtall.f  allinp.h                    ../HEADERS/dim.h                   ../HEADERS/parallel.h
wrtjsn.o:	wrtjsn.f  allinp.h                    ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.o rdinp

install:
	$(COPY) rdinp $(PREFIX)/bin

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran.a
	$(MAKE) -C ../COMMON       libfeffcom.a
	$(MAKE) -C ../PAR          libfeffpar.a
	$(MAKE) -C ../MATH         libfeffmath.a
	$(MAKE) -C ../JSON         libfeffjson.a

.PHONEY:	clean install all otherdirs
