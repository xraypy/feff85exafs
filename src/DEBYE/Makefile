

FORTRAN ?= gfortran
CFLAGS  ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm

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
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.a *.o

.PHONEY:	clean all
