

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

FORTRAN ?= gfortran
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM	?= rm

COPY    ?= cp -v
ARCHV   ?= .a


objects= feff.o rdinp.o chopen.o mkptz.o setgam.o terp.o potph.o rpotph.o \
	 rdhead.o atom.o indata.o getorb.o dirac.o inouh.o inth.o potsl.o \
	 somm.o totale.o cdsld.o ykdir.o renorm.o potslw.o ovrlp.o sumax.o \
	 frnrm.o istprm.o sidx.o istval.o fermi.o head.o wpot.o phmesh.o \
	 fixvar.o phase.o besjn.o bjnser.o xcpot.o rhl.o imhl.o ffq.o \
	 cubic.o quinn.o edp.o fovrg.o diff.o intpol.o wphase.o prcrit.o \
	 rphbin.o cpl0.o paths.o sortix.o ccrit.o mrb.o mcrith.o mcritk.o \
	 ipack.o heap.o pathsd.o timrep.o mpprmp.o phash.o outcrt.o \
	 strap.o mpprmd.o genfmt.o snlm.o rdpath.o setlam.o rot3i.o mmtr.o \
	 sclmz.o mmtrxi.o fmtrxi.o trap.o pijump.o ff2chi.o sigms.o dist.o \
	 sdist.o xx.o getxk.o dalp.o fpot.o exchan.o exchee.o pertab.o	\
	 genfmt_misc.o padlib.o strutils.o echo.o get_inpfile.o iread.o compat.o misc_io.o


archives =

all:	feff6l$(EXEXT)

feff6l$(EXEXT):	$(objects)
	$(FORTRAN) $(FLINKARGS) $(objects) -o $@

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<

clean:
	$(RM) *.o feff6l$(EXEXT)

install:
	$(COPY) feff6l$(EXEXT) $(BINDIR)

otherdirs:

.PHONEY:	clean install all otherdirs
