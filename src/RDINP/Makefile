

PREFIX  ?= /usr/local

FORTRAN ?= gfortran
CFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

COPY    ?= cp -v

objects = ffsort.o iniall.o mkptz.o rdinp_l.o rdline.o setedg.o	wrtall.o wrtjsn.o
archives = ../json-fortran/libjsonfortran.a ../COMMON/libfeffcom.a ../PAR/libfeffpar.a ../MATH/libfeffmath.a ../JSON/libfeffjson.a

rdinp:	$(objects) $(archives)
	$(FORTRAN) $(objects) $(archives) -o $@

ffsort.o:	ffsort.f  allinp.h ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
iniall.o:	iniall.f  allinp.h                    ../HEADERS/dim.h
mkptz.o:	mkptz.f            ../HEADERS/const.h
rdinp_l.o:	rdinp_l.f allinp.h ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
wrtall.o:	wrtall.f  allinp.h                    ../HEADERS/dim.h                   ../HEADERS/parallel.h
wrtjsn.o:	wrtjsn.f  allinp.h                    ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h

%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


all:	rdinp

clean:
	$(RM) *.o rdinp

install:
	$(COPY) rdinp $(PREFIX)/bin

.PHONEY:	clean install all
