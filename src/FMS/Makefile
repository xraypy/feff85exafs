

FORTRAN = gfortran
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib
RM      = rm

objects = fmsie.o fmspack.o ggbi.o gggm.o gglu.o ggrm.o ggtf.o xstaff.o yprep.o


libfefffms.a:	$(objects)
	$(AR) $(ARFLAGS) libfefffms.a $(objects)
	$(RANLIB) libfefffms.a

ffmod3.o:	ffmod3.f  ../HEADERS/dim.h ../HEADERS/parallel.h
fmsie.o:	fmsie.f   ../HEADERS/dim.h
fmspack.o:	fmspack.f ../HEADERS/dim.h ../HEADERS/parallel.h xparam.h
fmstot.o:	fmstot.f  ../HEADERS/dim.h ../HEADERS/const.h ../HEADERS/parallel.h
ggbi.o:		ggbi.f    ../HEADERS/dim.h xparam.h
gggm.o:		gggm.f    ../HEADERS/dim.h xparam.h
gglu.o:		gglu.f    ../HEADERS/dim.h xparam.h
ggrm.o:		ggrm.f    ../HEADERS/dim.h xparam.h
ggtf.o:		ggtf.f    ../HEADERS/dim.h xparam.h
reafms.o:	reafms.f  ../HEADERS/dim.h ../HEADERS/const.h
xprep.o:	xprep.f   ../HEADERS/dim.h ../HEADERS/parallel.h xparam.h
xstaff.o:	xstaff.f  ../HEADERS/dim.h xparam.h
yprep.o:	yprep.f   ../HEADERS/dim.h xparam.h


%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


all:	libfefffms.a

clean:
	$(RM) *.a *.o

.PHONEY:	clean all
