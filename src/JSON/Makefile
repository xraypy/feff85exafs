
FORTRAN ?= gfortran
CFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

objects = bailout.o read_geom.o read_atoms.o read_global.o read_pot.o	\
read_fpf0.o json_xsect.o read_xsect.o json_nnnn.o

all:	libfeffcom.a

libfeffjson.a:	$(objects)
	$(AR) $(ARFLAGS) libfeffjson.a $(objects)
	$(RANLIB) libfeffjson.a

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
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.a *.o

.PHONEY:	clean all
