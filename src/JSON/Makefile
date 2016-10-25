
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

objects = bailout.o read_geom.o read_atoms.o read_global.o read_pot.o	\
read_fpf0.o json_xsect.o read_xsect.o json_nnnn.o

all:	libfeffjson$(ARCHV)

libfeffjson$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) libfeffjson$(ARCHV) $(objects)
	$(RANLIB) libfeffjson$(ARCHV)

json_nnnn.o:	json_nnnn.f ../HEADERS/vers.h ../HEADERS/dim.h ../HEADERS/const.h
json_read_libpotph.o:	json_read_libpotph.f ../HEADERS/dim.h ../HEADERS/const.h
json_xsect.o:	json_xsect.f  ../HEADERS/vers.h ../HEADERS/dim.h ../HEADERS/parallel.h ../RDINP/allinp.h
read_atoms.o:	read_atoms.f  ../HEADERS/dim.h ../HEADERS/const.h
read_fpf0.o:	read_fpf0.f   ../HEADERS/dim.h ../HEADERS/const.h
read_geom.o:	read_geom.f   ../HEADERS/dim.h ../HEADERS/const.h
read_global.o:	read_global.f ../HEADERS/dim.h ../HEADERS/const.h
read_pot.o:	read_pot.f    ../HEADERS/dim.h ../HEADERS/const.h
read_xsect.o:	read_xsect.f  ../HEADERS/dim.h ../HEADERS/const.h



%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in JSON"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffjson$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
