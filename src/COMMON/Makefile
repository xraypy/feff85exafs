

PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
COPY    ?= cp -v

ARCHV   ?= .a

objects = chopen.o fixvar.o head.o itoken.o pertab.o rdcmt.o rdxsph.o	\
str.o xx.o fixdsp.o getorb.o iniptz.o nxtunt.o pijump.o rdhead.o	\
setgam.o str2dp.o fixdsx.o getxk.o isnum.o padlib.o qsorti.o rdpot.o	\
setkap.o wlog.o


all:	libfeffcom$(ARCHV)

libfeffcom$(ARCHV):	$(objects) padlib.h
	$(AR) $(ARFLAGS) libfeffcom$(ARCHV) $(objects)
	$(RANLIB) libfeffcom$(ARCHV)

fixdsp.o:	fixdsp.f ../HEADERS/dim.h ../HEADERS/const.h
fixdsx.o:	fixdsx.f ../HEADERS/dim.h ../HEADERS/const.h
fixvar.o:	fixvar.f ../HEADERS/dim.h ../HEADERS/const.h
head.o:		head.f   ../HEADERS/dim.h ../HEADERS/const.h ../HEADERS/vers.h
padlib.o:	padlib.f padlib.h
pijump.o:	pijump.f ../HEADERS/const.h
rdpot.o:	rdpot.f  ../HEADERS/dim.h
rdxsph.o:	rdxsph.f ../HEADERS/dim.h
wlog.o:		wlog.f   ../HEADERS/parallel.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.a *.o

INSTARCH = @echo "Nothing to install in COMMON"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffcom$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
