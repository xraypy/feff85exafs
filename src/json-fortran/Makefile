
PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

F90      ?= gfortran
F90FLAGS ?= -std=f2008 -c -O2 -fbacktrace -g -Wall -Wextra -Wno-maybe-uninitialized -pedantic -fPIC
AR       ?= ar
ARFLAGS  ?= rvc
RANLIB   ?= ranlib
RM       ?= rm
COPY     ?= cp -v
ARCHV    ?= .a

all:	json_module.mod json_module.o libjsonfortran$(ARCHV)

json_module.mod json_module.o:	json_module.f90
	$(F90) $(F90FLAGS) -o $@ $<

libjsonfortran$(ARCHV):
	$(AR) $(ARFLAGS) libjsonfortran$(ARCHV) json_module.o
	$(RANLIB) libjsonfortran$(ARCHV)

clean:
	$(RM) *$(ARCHV) *.o *.mod

install:
	$(COPY) libjsonfortran$(ARCHV) $(LIBDIR)
	$(COPY) json_module.mod        $(LIBDIR)

.PHONEY:	clean install all
