

FORTRAN = gfortran
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib
RM      = rm

objects = bcoef.o besjn.o conv.o csomm.o cwig3j.o determ.o exjlnl.o	\
phamp.o quartc.o sdist.o somm2.o terp.o trap.o besjh.o bjnser.o		\
cpl0.o csomm2.o czeros.o dist.o lu.o polint.o rotwig.o somm.o strap.o	\
terpc.o

libfeffmath.a:	$(objects)
	$(AR) $(ARFLAGS) libfeffmath.a $(objects)
	$(RANLIB) libfeffmath.a

bcoef.o:	bcoef.f ../HEADERS/dim.h
besjh.o:	besjh.f ../HEADERS/dim.h
confv.o:	conv.f  ../HEADERS/dim.h ../HEADERS/const.h
phamp.o:	phamp.f ../HEADERS/dim.h ../HEADERS/const.h


%.o:	%.f
	$(FORTRAN) $(CFLAGS) -o $@ $<


all:	libfeffmath.a

clean:
	$(RM) *.a *.o

.PHONEY:	clean all
