

PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
COPY    ?= cp -v
ARCHV   ?= .a

objects = csigma.o csigz.o cubic.o edp.o ffq.o fndsng.o imhl.o mpse.o	\
qsorti.o quinn.o rhl.o rhlbp.o vbh.o xcpot.o

all:	libfeffexch$(ARCHV)

libfeffexch$(ARCHV):	$(objects)
	$(AR) $(ARFLAGS) libfeffexch$(ARCHV) $(objects)
	$(RANLIB) libfeffexch$(ARCHV)

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
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *$(ARCHV) *.o

INSTARCH = @echo "Nothing to install in EXCH"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffexch$(ARCHV) $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
