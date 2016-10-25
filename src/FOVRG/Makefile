

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

objects = aprdec.o dfovrg.o dsordc.o intout.o nucdec.o potdvp.o		\
solin.o wfirdc.o yzktec.o aprdep.o diff.o inmuac.o muatcc.o ortdac.o	\
potex.o solout.o yzkrdc.o


all:	libfeffpha$(ARCHV)

libfeffpha$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) libfeffpha$(ARCHV) $(objects)
	$(RANLIB) libfeffpha$(ARCHV)

dfovrg.o:	dfovrg.f ../HEADERS/dim.h ../HEADERS/const.h
diff.o:		diff.f   ../HEADERS/dim.h
dsordc.o:	dsordc.f ../HEADERS/dim.h
inmuac.o:	inmuac.f ../HEADERS/dim.h
intout.o:	intout.f ../HEADERS/dim.h ../HEADERS/const.h
muatcc.o:	muatcc.f ../HEADERS/dim.h
nucdec.o:	nucdec.f ../HEADERS/dim.h
ortdac.o:	ortdac.f ../HEADERS/dim.h
potdvp.o:	potdvp.f ../HEADERS/dim.h
potex.o:	potex.f  ../HEADERS/dim.h
solin.o:	solin.f  ../HEADERS/dim.h ../HEADERS/const.h
solout.o:	solout.f ../HEADERS/dim.h ../HEADERS/const.h
wfirdc.o:	wfirdc.f ../HEADERS/dim.h
yzkrdc.o:	yzkrdc.f ../HEADERS/dim.h
yzktec.o:	yzktec.f ../HEADERS/dim.h



%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<


clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in FOVRG"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffpha$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
