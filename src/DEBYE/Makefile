

PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
COPY    ?= cp -v

objects = sigcl.o  sigm3.o  sigms.o  sigrem.o  sigte3.o


deps = ../COMMON/wlog.o ../COMMON/str.o ../COMMON/itoken.o	\
../COMMON/pertab.o ../COMMON/chopen.o ../MATH/dist.o		\
../PAR/sequential.o


all:	libfeffdw.a

libfeffdw.a:	$(objects) $(deps)
	$(AR) $(ARFLAGS) libfeffdw.a $(objects) $(deps)
	$(RANLIB) libfeffdw.a


../COMMON/wlog.o:	../COMMON/wlog.f
../COMMON/str.o:	../COMMON/str.f
../COMMON/itoken.o:	../COMMON/itoken.f
../COMMON/pertab.o:	../COMMON/pertab.f
../COMMON/chopen.o:	../COMMON/chopen.f
../MATH/dist.o:		../MATH/dist.f
../PAR/sequential.o:	../PAR/sequential.f

sigcl.o:	sigcl.f  ../HEADERS/const.h
sigms.o:	sigms.f  ../HEADERS/const.h
sigrem.o:	sigrem.f ../HEADERS/parallel.h

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<


clean:
	$(RM) *.a *.o

INSTARCH = @echo "Nothing to install in DEBYE"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffdw.a $(LIBDIR)
endif

install:
	$(INSTARCH)

.PHONEY:	clean all install
