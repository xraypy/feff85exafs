

PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
COPY    ?= cp -v
ARCHV   ?= .a

objects = fmsie.o fmspack.o ggbi.o gggm.o gglu.o ggrm.o ggtf.o xstaff.o yprep.o


libfefffms$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) libfefffms$(ARCHV) $(objects)
	$(RANLIB) libfefffms$(ARCHV)

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
	$(FORTRAN) $(FCFLAGS) -o $@ $<


all:	libfefffms$(ARCHV)

clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in FMS"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfefffms$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
