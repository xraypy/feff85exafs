
PREFIX  ?= /usr/local

FORTRAN ?= gfortran
SHARED  ?= -shared
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

CC      ?= gcc
CCFLAGS ?= -c -g -fPIC

COPY    ?= cp -v

all:	otherdirs libfeffgenfmt.a libonepath.so libfeffpath.so genfmt

libobj = genfmt_prep.o fmtrxi.o mmtr.o mmtrxi.o import.o rot3i.o sclmz.o setlam.o snlm.o xstar.o trig.o

libfeffgenfmt.a:	$(libobj)
	$(AR) $(ARFLAGS) libfeffgenfmt.a $(libobj)
	$(RANLIB) libfeffgenfmt.a

shobjects = ../FF2X/fdtarr.o ../FF2X/fdthea.o ../FF2X/fdtxdi.o		\
../RDINP/mkptz.o ../COMMON/wlog.o ../COMMON/chopen.o			\
../COMMON/getxk.o ../COMMON/head.o ../COMMON/padlib.o			\
../COMMON/pertab.o ../COMMON/pijump.o ../COMMON/rdxsph.o		\
../COMMON/setgam.o ../COMMON/setkap.o ../COMMON/str.o ../MATH/bcoef.o	\
../MATH/cwig3j.o ../MATH/dist.o ../MATH/polint.o ../MATH/rotwig.o	\
../MATH/terp.o ../MATH/trap.o ../PAR/sequential.o


../FF2X/fdtarr.o:	../FF2X/fdtarr.f    ../HEADERS/const.h ../HEADERS/dim.h
../FF2X/fdthea.o:	../FF2X/fdthea.f    ../HEADERS/const.h ../HEADERS/dim.h
../FF2X/fdtxdi.o:	../FF2X/fdtxdi.f    ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
../RDINP/mkptz.o:	../RDINP/mkptz.f    ../HEADERS/const.h
../COMMON/wlog.o:	../COMMON/wlog.f                                        ../HEADERS/parallel.h
../COMMON/head.o:	../COMMON/head.f    ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
../COMMON/padlib.o:	../COMMON/padlib.f                                      ../COMMON/padlib.h
../COMMON/pijump.o:	../COMMON/pijump.f  ../HEADERS/const.h
../COMMON/rdxsph.o:	../COMMON/rdxsph.f                     ../HEADERS/dim.h
../MATH/bcoef.o:	../MATH/bcoef.f                        ../HEADERS/dim.h
../PAR/sequential.o:	../PAR/sequential.f                                     ../HEADERS/parallel.h

libonepath.so:	onepath.o pathgeom.o $(libobj) $(shobjects)
	$(FORTRAN) $(SHARED) onepath.o pathgeom.o $(libobj) $(shobjects) -o $@

libfeffpath.so:	feffpath.o libonepath.so
	$(FORTRAN) $(SHARED) feffpath.o libonepath.so -o $@


binobj = ffmod5.o genfmt.o rdpath.o regenf.o
archives = ../json-fortran/libjsonfortran.a ../DEBYE/libfeffdw.a	\
../COMMON/libfeffcom.a ../PAR/libfeffpar.a ../MATH/libfeffmath.a	\
../JSON/libfeffjson.a



genfmt:	$(binobj) libfeffgenfmt.a $(archives)
	$(FORTRAN) $(binobj) libfeffgenfmt.a $(archives) -o $@


feffpath.o:	feffpath.c feffpath.h
ffmod5.o:	ffmod5.f       ../HEADERS/dim.h                    ../HEADERS/parallel.h
fmtrxi.o:	fmtrxi.f       ../HEADERS/const.h ../HEADERS/dim.h
genfmt.o:	genfmt.f       ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
genfmt_prep.o:	genfmt_prep.f  ../HEADERS/const.h ../HEADERS/dim.h
import.o:	import.f       ../HEADERS/const.h ../HEADERS/dim.h
mmtr.o:		mmtr.f         ../HEADERS/const.h ../HEADERS/dim.h
mmtrxi.o:	mmtrxi.f       ../HEADERS/const.h ../HEADERS/dim.h
onepath.o:	onepath.f      ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h
pathgeom.o:	pathgeom.f     ../HEADERS/const.h ../HEADERS/dim.h
rdpath.o:	rdpath.f       ../HEADERS/const.h ../HEADERS/dim.h
read_onepath.o:	read_onepath.f ../HEADERS/const.h ../HEADERS/dim.h
regenf.o:	regenf.f       ../HEADERS/const.h ../HEADERS/dim.h
rot3i.o:	rot3i.f                           ../HEADERS/dim.h
sclmz.o:	sclmz.f        ../HEADERS/const.h ../HEADERS/dim.h
setlam.o:	setlam.f       ../HEADERS/const.h ../HEADERS/dim.h
snlm.o:		snlm.f                            ../HEADERS/dim.h
xstar.o:	xstar.f        ../HEADERS/const.h ../HEADERS/dim.h



%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<

%.o:	%.c
	$(CC) $(CCFLAGS) -o $@ $<


clean:
	$(RM) *.so *.o *.a *.so genfmt

INSTARCH = @echo "Not installing .a files"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffgenfmt.a $(PREFIX)/lib
endif

install:
	$(COPY) genfmt         $(PREFIX)/bin
	$(COPY) libonepath.so  $(PREFIX)/lib
	$(COPY) libfeffpath.so $(PREFIX)/lib
	$(COPY) feffpath.h     $(PREFIX)/include
	$(INSTARCH)

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran.a
	$(MAKE) -C ../DEBYE        libfeffdw.a
	$(MAKE) -C ../COMMON       libfeffcom.a
	$(MAKE) -C ../PAR          libfeffpar.a
	$(MAKE) -C ../MATH         libfeffmath.a
	$(MAKE) -C ../JSON         libfeffjson.a

.PHONEY:	clean install all otherdirs
