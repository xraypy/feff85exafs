
PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin
LIBDIR  ?= $(PREFIX)/lib
INCDIR  ?= $(PREFIX)/include

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
ARCHV    ?= .a
SHOBJ    ?= .so
EXEXT    ?=

FORTRAN ?= gfortran
SHARED  ?= -shared
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC


archives = ../JSON/libfeffjson$(ARCHV)					\
../json-fortran/libjsonfortran$(ARCHV) ../COMMON/libfeffcom$(ARCHV)	\
../PAR/libfeffpar$(ARCHV) ../MATH/libfeffmath$(ARCHV)

objects = oca_dimsmod.o oca_atoms_inp.o oca_geometry_inp.o		\
oca_global_inp.o oca_opcons_inp.o oca_errormod.o oca_padio.o		\
oca_iofiles.o oca_iomod.o oca_kinds.o oca_mtdp.o oca_atomicpotio.o	\
oca_potential_inp.o oca_errorfile.o oca_constants.o oca_par.o eps.o	\
getelement.o bwords_nc.o

all:	otherdirs libfeffopconsat$(ARCHV) opconsat$(EXEXT) eps2exc$(EXEXT) libfeffloss$(SHOBJ)

opconsat$(EXEXT):	opconsat.o e2emisc.o $(objects) $(archives)
	$(FORTRAN) -o $@ $< e2emisc.o $(objects) $(archives)

eps2exc$(EXEXT):	eps2exc.o e2emisc.o rdloss.o ../MATH/libfeffmath$(ARCHV) ../COMMON/libfeffcom$(ARCHV)
	$(FORTRAN) -o $@ $< e2emisc.o rdloss.o ../MATH/libfeffmath$(ARCHV) ../COMMON/libfeffcom$(ARCHV)


%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<

%.o:	%.f90
	$(F90) $(FCFLAGS) -J${CURDIR} -o $@ $<

libfeffloss.o:	libfeffloss.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h

loss.o:		loss.f        ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h
loss$(EXEXT):	loss.o libfeffloss.o e2emisc.o libfeffopconsat$(ARCHV)
	$(FORTRAN) -o $@  libfeffloss.o e2emisc.o libfeffopconsat$(ARCHV) $< ../POT/inipotph.o ../JSON/json_read_libpotph.o $(archives)

libfeffopconsat$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) $@ $(objects)
	$(RANLIB) $@

libfeffloss$(SHOBJ):	libfeffloss.o e2emisc.o eps.o libfeffopconsat$(ARCHV)
	$(FORTRAN) $(SHARED) libfeffloss.o e2emisc.o eps.o libfeffopconsat$(ARCHV) -o $@

otherdirs:
	$(MAKE) -C ../COMMON       libfeffcom$(ARCHV)
	$(MAKE) -C ../PAR          libfeffpar$(ARCHV)
	$(MAKE) -C ../MATH         libfeffmath$(ARCHV)
	$(MAKE) -C ../JSON         libfeffjson$(ARCHV)
	$(MAKE) -C ../json-fortran libjsonfortran$(ARCHV)

clean:
	$(RM) *$(ARCHV) *.o *.mod opconsat$(EXEXT) eps2exc$(EXEXT)

INSTARCH = @echo "Not installing $(ARCHV) files"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffopconsat$(ARCHV) $(PREFIX)/lib
endif

install:
	$(COPY) opconsat$(EXEXT)    $(BINDIR)
	$(COPY) eps2exc$(EXEXT)     $(BINDIR)
	$(COPY) libfeffloss$(SHOBJ) $(LIBDIR)
	$(INSTARCH)

.PHONEY:	clean install all otherdirs
