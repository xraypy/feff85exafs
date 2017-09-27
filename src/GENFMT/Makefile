
PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin
LIBDIR  ?= $(PREFIX)/lib
INCDIR  ?= $(PREFIX)/include

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

CLINK = gcc -O2 -fPIC -shared
CLINKARGS = -lgfortran

ARCHV = .a
SHOBJ = .so
EXEXT =
FCFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FLINKARGS =


ifeq ($(OS),Windows_NT)
    UNAME = Windows_NT
    SHOBJ = .dll
    EXEXT = .exe
    FCFLAGS  = -c -O3 -ffree-line-length-none -g -Wall
    CLINK = gcc -O2 -shared
    CLINKARGS = -Wl,--add-stdcall-alias -lgcc -lgfortran -lm -mwindows
else
    UNAME := $(shell uname -s)
    ifeq ($(UNAME),Darwin)
        SHOBJ = .dylib
        FLINKARGS =  -headerpad_max_install_names
        FORTRAN = gfortran -arch x86_64
	CC      = gcc -O2 -arch x86_64
	CLINK   = gcc -O2 -arch x86_64 -fPIC -shared
        CLINKARGS = -L/usr/local/gfortran/lib/ -lgfortran
    endif
endif


all:	otherdirs libfeffgenfmt$(ARCHV) libonepath$(SHOBJ) libfeff8lpath$(SHOBJ) feff8l_genfmt$(EXEXT)

libobj = genfmt_prep.o fmtrxi.o mmtr.o mmtrxi.o import.o rot3i.o sclmz.o setlam.o snlm.o xstar.o trig.o

libfeffgenfmt$(ARCHV):	$(libobj)
	$(AR) $(ARFLAGS) libfeffgenfmt$(ARCHV) $(libobj)
	$(RANLIB) libfeffgenfmt$(ARCHV)

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

libonepath$(SHOBJ):	onepath.o pathgeom.o $(libobj) $(shobjects)
	$(FORTRAN) $(SHARED) onepath.o pathgeom.o $(libobj) $(shobjects) -o $@

libfeff8lpath$(SHOBJ):	feffpath.o onepath.o pathgeom.o $(libobj) $(shobjects)
	$(CLINK) -o $@ feffpath.o onepath.o pathgeom.o $(libobj) $(shobjects)  $(CLINKARGS)


binobj = ffmod5.o genfmt.o rdpath.o regenf.o
archives = ../json-fortran/libjsonfortran$(ARCHV)	\
../DEBYE/libfeffdw$(ARCHV) ../COMMON/libfeffcom$(ARCHV)	\
../PAR/libfeffpar$(ARCHV) ../MATH/libfeffmath$(ARCHV)	\
../JSON/libfeffjson$(ARCHV)


feff8l_genfmt$(EXEXT):	$(binobj) libfeffgenfmt$(ARCHV) $(archives)
	$(FORTRAN) $(FLINKARGS) $(binobj) libfeffgenfmt$(ARCHV) $(archives) -o $@

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
	$(RM) *$(SHOBJ) *.o *$(ARCHV) feff8l_genfmt$(EXEXT)

INSTARCH = @echo "Not installing $(ARCHV) files"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffgenfmt$(ARCHV) $(PREFIX)/lib
endif

install:
	$(COPY) feff8l_genfmt$(EXEXT)      $(BINDIR)
	$(COPY) libonepath$(SHOBJ)  $(LIBDIR)
	$(COPY) libfeff8lpath$(SHOBJ) $(LIBDIR)
	$(COPY) feffpath.h          $(INCDIR)
	$(INSTARCH)

otherdirs:
	$(MAKE) -C ../json-fortran libjsonfortran$(ARCHV)
	$(MAKE) -C ../DEBYE        libfeffdw$(ARCHV)
	$(MAKE) -C ../COMMON       libfeffcom$(ARCHV)
	$(MAKE) -C ../PAR          libfeffpar$(ARCHV)
	$(MAKE) -C ../MATH         libfeffmath$(ARCHV)
	$(MAKE) -C ../JSON         libfeffjson$(ARCHV)

.PHONEY:	clean install all otherdirs
