
PREFIX  ?= /usr/local

F90      ?= gfortran
F90FLAGS ?= -std=f2008 -c -O2 -fbacktrace -g -Wall -Wextra -Wno-maybe-uninitialized -pedantic -fPIC
    # several "Unused dummy argument" warnings remain: padio, iomod, and atomicpotio

## apparently these are the flags for ifort
#F90FLAGS  = -O2 -warn -stand f08 -diag-disable 7601 -traceback -module ${CURDIR}

AR       ?= ar
ARFLAGS  ?= rvc
RANLIB   ?= ranlib
RM       ?= rm
COPY     ?= cp -v

FORTRAN ?= gfortran
SHARED  ?= -shared
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC


archives = ../JSON/libfeffjson.a ../json-fortran/libjsonfortran.a	\
../COMMON/libfeffcom.a ../PAR/libfeffpar.a ../MATH/libfeffmath.a

objects = oca_dimsmod.o oca_atoms_inp.o oca_geometry_inp.o		\
oca_global_inp.o oca_opcons_inp.o oca_errormod.o oca_padio.o		\
oca_iofiles.o oca_iomod.o oca_kinds.o oca_mtdp.o oca_atomicpotio.o	\
oca_potential_inp.o oca_errorfile.o oca_constants.o oca_par.o eps.o	\
getelement.o bwords_nc.o

all:	otherdirs libfeffopconsat.a opconsat eps2exc libfeffloss.so

opconsat:	opconsat.o e2emisc.o $(objects) $(archives)
	$(FORTRAN) -o $@ $< e2emisc.o $(objects) $(archives)

eps2exc:	eps2exc.o e2emisc.o rdloss.o ../MATH/libfeffmath.a ../COMMON/libfeffcom.a
	$(FORTRAN) -o $@ $< e2emisc.o rdloss.o ../MATH/libfeffmath.a ../COMMON/libfeffcom.a


%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<

%.o:	%.f90
	$(F90) $(FCFLAGS) -J${CURDIR} -o $@ $<

libfeffloss.o:	libfeffloss.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h
loss.o:		loss.f        ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h

libfeffopconsat.a:	$(objects)
	$(AR) $(ARFLAGS) $@ $(objects)
	$(RANLIB) $@

libfeffloss.so:	libfeffloss.o e2emisc.o eps.o libfeffopconsat.a
	$(FORTRAN) $(SHARED) libfeffloss.o e2emisc.o eps.o libfeffopconsat.a -o $@

otherdirs:
	$(MAKE) -C ../COMMON       libfeffcom.a
	$(MAKE) -C ../PAR          libfeffpar.a
	$(MAKE) -C ../MATH         libfeffmath.a
	$(MAKE) -C ../JSON         libfeffjson.a
	$(MAKE) -C ../json-fortran libjsonfortran.a

clean:
	$(RM) *.a *.o *.mod opconsat eps2exc

INSTARCH = @echo "Not installing .a files"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffopconsat.a $(PREFIX)/lib
endif

install:
	$(COPY) opconsat         $(PREFIX)/bin
	$(COPY) eps2exc          $(PREFIX)/bin
	$(COPY) libfeffloss.so   $(PREFIX)/lib
	$(INSTARCH)

.PHONEY:	clean install all otherdirs
