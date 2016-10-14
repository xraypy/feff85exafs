
FORTRAN ?= gfortran
CFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
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
	$(FORTRAN) $(CFLAGS) -o $@ $<

#%.os:	%.f
#	$(FORTRAN) $(CFLAGS) -o $@ $<



clean:
	$(RM) *.a *.o

.PHONEY:	clean all
