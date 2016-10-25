
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

objects = akeato.o cofcon.o etotal.o fpf0.o lagdat.o nucdev.o	\
potslw.o soldir.o wfirdf.o aprdev.o dentfa.o fdmocc.o inmuat.o	\
messer.o ortdat.o s02at.o tabrat.o yzkrdf.o bkmrdf.o dsordf.o	\
fdrirk.o intdir.o muatco.o potrdf.o scfdat.o vlda.o yzkteg.o


deps = ../COMMON/wlog.o ../COMMON/getorb.o ../COMMON/str.o		\
../PAR/sequential.o ../MATH/determ.o ../MATH/somm.o ../MATH/somm2.o	\
../MATH/cwig3j.o ../EXCH/vbh.o ../EXCH/edp.o


all:	libfeffatom$(ARCHV)

libfeffatom$(ARCHV):	$(objects) $(deps)
	$(AR) $(ARFLAGS) libfeffatom$(ARCHV) $(objects) $(deps)
	$(RANLIB) libfeffatom$(ARCHV)


../COMMON/wlog.o:	../COMMON/wlog.f
../COMMON/getorb.o:	../COMMON/getorb.f
../COMMON/str.o:	../COMMON/str.f
../PAR/sequential.o:	../PAR/sequential.f
../MATH/determ.o:	../MATH/determ.f
../MATH/somm.o:		../MATH/somm.f
../MATH/somm2.o:	../MATH/somm2.f
../MATH/cwig3j.o:	../MATH/cwig3j.f
../EXCH/vbh.o:		../EXCH/vbh.f
../EXCH/edp.o:		../EXCH/edp.f

fpf0.o:		fpf0.f   ../HEADERS/dim.h ../HEADERS/const.h ../HEADERS/parallel.h
scfdat.o:	scfdat.f ../HEADERS/dim.h                    ../HEADERS/parallel.h
vlda.o:		vlda.f                    ../HEADERS/const.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in ATOM"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffatom$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
