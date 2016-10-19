
PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

objects = sequential.o
#objects = parallel.o

all:	libfeffpar.a

libfeffpar.a:	$(objects)
	$(AR) $(ARFLAGS) libfeffpar.a $(objects)
	$(RANLIB) libfeffpar.a

sequential.o:	sequential.f ../HEADERS/parallel.h
#parallel.o:	parallel.f   ../HEADERS/parallel.h mpif.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<


clean:
	$(RM) *.a *.o

INSTARCH = @echo "Nothing to install in PAR"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffpar.a $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
