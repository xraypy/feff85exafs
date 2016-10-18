

F90      ?= gfortran
F90FLAGS ?= -std=f2008 -c -O2 -fbacktrace -g -Wall -Wextra -Wno-maybe-uninitialized -pedantic -fPIC
AR       ?= ar
ARFLAGS  ?= rvc
RANLIB   ?= ranlib
RM       ?= rm
COPY     ?= cp -v

all:	json_module.mod json_module.o libjsonfortran.a

json_module.mod json_module.o:	json_module.f90
	$(F90) $(F90FLAGS) -o $@ $<

libjsonfortran.a:
	$(AR) $(ARFLAGS) libjsonfortran.a json_module.o
	$(RANLIB) libjsonfortran.a

clean:
	$(RM) *.a *.o *.mod

install:
	$(COPY) libjsonfortran.a $(PREFIX)/lib
	$(COPY) json_module.mod  $(PREFIX)/lib

.PHONEY:	clean install all
