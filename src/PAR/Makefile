
PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
ARCHV   ?= .a

objects = sequential.o
#objects = parallel.o

all:	libfeffpar$(ARCHV)

libfeffpar$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) libfeffpar$(ARCHV) $(objects)
	$(RANLIB) libfeffpar$(ARCHV)

sequential.o:	sequential.f ../HEADERS/parallel.h
#parallel.o:	parallel.f   ../HEADERS/parallel.h mpif.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<


clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in PAR"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffpar$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
