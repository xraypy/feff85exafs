

PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
COPY    ?= cp -v
ARCHV   ?= .a

objects = bcoef.o besjn.o conv.o csomm.o cwig3j.o determ.o exjlnl.o	\
phamp.o quartc.o sdist.o somm2.o terp.o trap.o besjh.o bjnser.o		\
cpl0.o csomm2.o czeros.o dist.o lu.o polint.o rotwig.o somm.o strap.o	\
terpc.o

all:	libfeffmath$(ARCHV)

libfeffmath$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) libfeffmath$(ARCHV) $(objects)
	$(RANLIB) libfeffmath$(ARCHV)

bcoef.o:	bcoef.f ../HEADERS/dim.h
besjh.o:	besjh.f ../HEADERS/dim.h
confv.o:	conv.f  ../HEADERS/dim.h ../HEADERS/const.h
phamp.o:	phamp.f ../HEADERS/dim.h ../HEADERS/const.h


%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<


clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in MATH"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffmath$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
