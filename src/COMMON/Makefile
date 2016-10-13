

FORTRAN = gfortran
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib

objects = chopen.o fixvar.o head.o itoken.o pertab.o rdcmt.o rdxsph.o	\
str.o xx.o fixdsp.o getorb.o iniptz.o nxtunt.o pijump.o rdhead.o	\
setgam.o str2dp.o fixdsx.o getxk.o isnum.o padlib.o qsorti.o rdpot.o	\
setkap.o wlog.o


libfeffcom.a:	$(objects) padlib.h
	$(AR) $(ARFLAGS) libfeffcom.a $(objects)
	$(RANLIB) libfeffcom.a

fixdsp.o:	fixdsp.f ../HEADERS/dim.h ../HEADERS/const.h
fixdsx.o:	fixdsx.f ../HEADERS/dim.h ../HEADERS/const.h
fixvar.o:	fixvar.f ../HEADERS/dim.h ../HEADERS/const.h
head.o:		head.f   ../HEADERS/dim.h ../HEADERS/const.h ../HEADERS/vers.h
padlib.o:	padlib.f padlib.h
pijump.o:	pijump.f ../HEADERS/const.h
rdpot.o:	rdpot.f  ../HEADERS/dim.h
rdxsph.o:	rdxsph.f ../HEADERS/dim.h
wlog.o:		wlog.f   ../HEADERS/parallel.h

%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


all:	libfeffcom.a

clean:
	rm *.a *.o

.PHONEY:	clean all
