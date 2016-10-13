

FORTRAN = gfortran
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib
RM      = rm

objects = aprdec.o dfovrg.o dsordc.o intout.o nucdec.o potdvp.o		\
solin.o wfirdc.o yzktec.o aprdep.o diff.o inmuac.o muatcc.o ortdac.o	\
potex.o solout.o yzkrdc.o


libfeffpha.a:	$(objects)
	$(AR) $(ARFLAGS) libfeffpha.a $(objects)
	$(RANLIB) libfeffpha.a

dfovrg.o:	dfovrg.f ../HEADERS/dim.h ../HEADERS/const.h
diff.o:		diff.f   ../HEADERS/dim.h
dsordc.o:	dsordc.f ../HEADERS/dim.h
inmuac.o:	inmuac.f ../HEADERS/dim.h
intout.o:	intout.f ../HEADERS/dim.h ../HEADERS/const.h
muatcc.o:	muatcc.f ../HEADERS/dim.h
nucdec.o:	nucdec.f ../HEADERS/dim.h
ortdac.o:	ortdac.f ../HEADERS/dim.h
potdvp.o:	potdvp.f ../HEADERS/dim.h
potex.o:	potex.f  ../HEADERS/dim.h
solin.o:	solin.f  ../HEADERS/dim.h ../HEADERS/const.h
solout.o:	solout.f ../HEADERS/dim.h ../HEADERS/const.h
wfirdc.o:	wfirdc.f ../HEADERS/dim.h
yzkrdc.o:	yzkrdc.f ../HEADERS/dim.h
yzktec.o:	yzktec.f ../HEADERS/dim.h



%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


all:	libfeffpha.a

clean:
	$(RM) *.a *.o

.PHONEY:	clean all
