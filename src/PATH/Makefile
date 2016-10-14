

FORTRAN ?= gfortran
CFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

COPY    ?= cp -v

objects = ccrit.o heap.o mcrith.o mpprmd.o mrb.o paths.o phash.o	\
repath.o timrep.o ffmod4.o ipack.o mcritk.o mpprmp.o outcrt.o		\
pathsd.o prcrit.o sortix.o
archives = ../json-fortran/libjsonfortran.a ../COMMON/libfeffcom.a ../PAR/libfeffpar.a ../MATH/libfeffmath.a ../JSON/libfeffjson.a

pathfinder:	$(objects)
	$(FORTRAN) $(objects) $(archives) -o $@


ccrit.o:	ccrit.f  ../HEADERS/const.h ../HEADERS/dim.h
ffmod4.o:	ffmod4.f ../HEADERS/dim.h   ../HEADERS/const.h ../HEADERS/parallel.h
mcrith.o:	mcrith.f ../HEADERS/const.h ../HEADERS/dim.h
mcritk.o:	mcritk.f ../HEADERS/const.h ../HEADERS/dim.h
mpprmd.o:	mpprmd.f ../HEADERS/dim.h
mpprmp.o:	mpprmp.f ../HEADERS/dim.h
mrb.o:		mrb.f    ../HEADERS/dim.h
outcrt.o:	outcrt.f ../HEADERS/const.h ../HEADERS/dim.h
paths.o:	paths.f  ../HEADERS/dim.h
pathsd.o:	pathsd.f ../HEADERS/const.h ../HEADERS/dim.h
phash.o:	phash.f  ../HEADERS/dim.h
prcrit.o:	prcrit.f ../HEADERS/const.h ../HEADERS/dim.h
repath.o:	repath.f ../HEADERS/const.h ../HEADERS/dim.h
timrep.o:	timrep.f ../HEADERS/dim.h


%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


all:	pathfinder

clean:
	$(RM) *.o pathfinder

install:
	$(COPY) pathfinder $(PREFIX)/bin

.PHONEY:	clean install all
