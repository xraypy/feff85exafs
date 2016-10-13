

FORTRAN = gfortran
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib
RM      = rm

objects = csigma.o csigz.o cubic.o edp.o ffq.o fndsng.o imhl.o mpse.o	\
qsorti.o quinn.o rhl.o rhlbp.o vbh.o xcpot.o

libfeffexch.a:	$(objects)
	$(AR) $(ARFLAGS) libfeffexch.a $(objects)
	$(RANLIB) libfeffexch.a

csigma.o:	csigma.f ../HEADERS/const.h
csigz.o:	csigz.f  ../HEADERS/const.h
edp.o:		edp.f    ../HEADERS/const.h
imhl.o:		imhl.f   ../HEADERS/const.h
mpse.o:		mpse.f   ../HEADERS/const.h   ../HEADERS/dim.h
quinn.o:	quinn.f  ../HEADERS/const.h
rhl.o:		rhl.f    ../HEADERS/const.h
rhlbp.o:	rhlbp.f  ../HEADERS/const.h
xcpot.o:	xcpot.o  ../HEADERS/const.h   ../HEADERS/dim.h


%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


all:	libfeffexch.a

clean:
	$(RM) *.a *.o

.PHONEY:	clean all
