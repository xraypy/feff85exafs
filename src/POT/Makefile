

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin
LIBDIR  ?= $(PREFIX)/lib
INCDIR  ?= $(PREFIX)/include

FORTRAN ?= gfortran
SHARED  ?= -shared
FCFLAGS ?= -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   ?= -I../json-fortran -J../json-fortran
AR      ?= ar
ARFLAGS ?= rvc
RANLIB  ?= ranlib
RM      ?= rm
ARCHV   ?= .a
SHOBJ   ?= .so
EXEXT   ?=

CC      ?= gcc
CCFLAGS ?= -c -g -fPIC

COPY    ?= cp -v

RDINPobj = ../RDINP/mkptz.o

COMMONobj = ../COMMON/chopen.o ../COMMON/fixdsp.o ../COMMON/fixdsx.o	\
../COMMON/fixvar.o ../COMMON/getorb.o ../COMMON/getxk.o			\
../COMMON/head.o ../COMMON/padlib.o ../COMMON/setkap.o			\
../COMMON/rdcmt.o ../COMMON/str.o ../COMMON/xx.o

ATOMobj = ../ATOM/akeato.o ../ATOM/aprdev.o ../ATOM/cofcon.o		\
../ATOM/bkmrdf.o ../ATOM/dentfa.o ../ATOM/dsordf.o			\
../ATOM/etotal.o ../ATOM/fdmocc.o ../ATOM/fdrirk.o			\
../ATOM/inmuat.o ../ATOM/intdir.o ../ATOM/lagdat.o			\
../ATOM/messer.o ../ATOM/muatco.o ../ATOM/nucdev.o			\
../ATOM/ortdat.o ../ATOM/potrdf.o ../ATOM/potslw.o			\
../ATOM/scfdat.o ../ATOM/soldir.o ../ATOM/s02at.o			\
../ATOM/tabrat.o ../ATOM/vlda.o ../ATOM/wfirdf.o ../ATOM/yzkteg.o	\
../ATOM/yzkrdf.o

MATHobj = ../MATH/bcoef.o ../MATH/besjh.o ../MATH/besjn.o		\
../MATH/bjnser.o ../MATH/csomm.o ../MATH/csomm2.o ../MATH/cwig3j.o	\
../MATH/czeros.o ../MATH/determ.o ../MATH/dist.o ../MATH/exjlnl.o	\
../MATH/lu.o ../MATH/phamp.o ../MATH/polint.o ../MATH/rotwig.o		\
../MATH/somm.o ../MATH/somm2.o ../MATH/terp.o ../MATH/terpc.o


FMSobj = ../FMS/fmspack.o ../FMS/fmsie.o ../FMS/gglu.o ../FMS/ggbi.o	\
../FMS/ggrm.o ../FMS/gggm.o ../FMS/ggtf.o ../FMS/yprep.o		\
../FMS/xstaff.o

FF2Xobj = ../FF2X/xscorr.o

EXCHobj = ../EXCH/cubic.o ../EXCH/csigma.o ../EXCH/csigz.o	\
../EXCH/edp.o ../EXCH/ffq.o ../EXCH/fndsng.o ../EXCH/imhl.o	\
../EXCH/qsorti.o ../EXCH/quinn.o ../EXCH/rhl.o ../EXCH/rhlbp.o	\
../EXCH/vbh.o ../EXCH/xcpot.o

FOVRGobj = ../FOVRG/dfovrg.o ../FOVRG/inmuac.o ../FOVRG/diff.o		\
../FOVRG/wfirdc.o ../FOVRG/potex.o ../FOVRG/solin.o ../FOVRG/solout.o	\
../FOVRG/muatcc.o ../FOVRG/nucdec.o ../FOVRG/potdvp.o			\
../FOVRG/yzkrdc.o ../FOVRG/yzktec.o ../FOVRG/aprdep.o			\
../FOVRG/aprdec.o ../FOVRG/intout.o

PARobj = ../PAR/sequential.o

XSPHobj = ../XSPH/acoef.o ../XSPH/axafs.o ../XSPH/fmssz.o		\
../XSPH/phase.o ../XSPH/phmesh.o ../XSPH/radint.o			\
../XSPH/rholsz.o ../XSPH/szlz.o ../XSPH/wphase.o ../XSPH/wrxsph.o	\
../XSPH/xmult.o ../XSPH/xsect.o ../XSPH/xsph.o


objects = libpotph.o inipotph.o ffsort_potph.o		\
wlog_potph.o pot.o afolp.o inipot.o moveh.o wpot.o ovrlp.o broydn.o	\
ff2g.o grids.o corval.o wrpot.o scmt.o scmtmp.o sumax.o rholie.o	\
frnrm.o coulom.o istprm.o movrlp.o ovp2mt.o fermi.o sidx.o
# libpotph_xsect.o

archives = ../POT/libfeffint$(ARCHV) ../ATOM/libfeffatom$(ARCHV)			\
../FOVRG/libfeffpha$(ARCHV) ../COMMON/libfeffcom$(ARCHV) ../PAR/libfeffpar$(ARCHV)	\
../MATH/libfeffmath$(ARCHV) ../JSON/libfeffjson$(ARCHV) ../EXCH/libfeffexch$(ARCHV)	\
../FMS/libfefffms$(ARCHV) ../json-fortran/libjsonfortran$(ARCHV)
# ../OPCONSAT/libfeffloss$(ARCHV)


binobj = ffmod1.o pot.o scmtmp.o afolp.o broydn.o corval.o coulom.o	\
ff2g.o frnrm.o grids.o inipot.o istval.o moveh.o ovrlp.o reapot.o	\
rholie.o scmt.o sumax.o wpot.o wrpot.o

shobjects = $(objects) $(RDINPobj) $(COMMONobj) $(ATOMobj) $(MATHobj) $(FMSobj) $(FF2Xobj) $(EXCHobj) $(FOVRGobj) $(PARobj) $(XSPHobj)

JSONobj = ../JSON/json_read_libpotph.o ../JSON/json_xsect.o ../JSON/bailout.o ../json-fortran/libjsonfortran$(ARCHV)



all:	otherdirs libfeffint$(ARCHV) libpotph$(SHOBJ) libfeff8lpotph$(SHOBJ) feff8l_pot$(EXEXT)

libfeffint$(ARCHV):	istprm.o movrlp.o ovp2mt.o fermi.o sidx.o
	$(AR) $(ARFLAGS) libfeffint$(ARCHV) istprm.o movrlp.o ovp2mt.o fermi.o sidx.o
	$(RANLIB) libfeffint$(ARCHV)

libpotph$(SHOBJ):	$(shobjects) $(JSONobj)
	$(FORTRAN) $(SHARED) $(shobjects) $(JSONobj) -o $@

libfeff8lpotph$(SHOBJ):	feffphases.o $(shobjects) $(JSONobj)
	$(FORTRAN) $(SHARED) feffphases.o $(shobjects) $(JSONobj) -o $@


feff8l_pot$(EXEXT):	$(binobj) $(archives)
	$(FORTRAN) $(FLINKARGS) $(binobj) $(archives) -o $@


../COMMON/fixdsp.o:	../COMMON/fixdsp.f  ../HEADERS/const.h ../HEADERS/dim.h
../COMMON/fixdsx.o:	../COMMON/fixdsx.f  ../HEADERS/const.h ../HEADERS/dim.h
../COMMON/fixvar.o:	../COMMON/fixvar.f  ../HEADERS/const.h ../HEADERS/dim.h
../COMMON/padlib.o:	../COMMON/padlib.f                                      ../COMMON/padlib.h
../COMMON/pijump.o:	../COMMON/pijump.f  ../HEADERS/const.h
../COMMON/rdxsph.o:	../COMMON/rdxsph.f                     ../HEADERS/dim.h

../ATOM/scfdat.o:	../ATOM/scfdat.f ../HEADERS/dim.h                    ../HEADERS/parallel.h
../ATOM/vlda.o:		../ATOM/vlda.f                    ../HEADERS/const.h

../MATH/bcoef.o:	../MATH/bcoef.f ../HEADERS/dim.h
../MATH/besjh.o:	../MATH/besjh.f ../HEADERS/dim.h
../MATH/phamp.o:	../MATH/phamp.f ../HEADERS/dim.h ../HEADERS/const.h

../FMS/fmsie.o:		../FMS/fmsie.f   ../HEADERS/dim.h
../FMS/fmspack.o:	../FMS/fmspack.f ../HEADERS/dim.h ../HEADERS/parallel.h ../FMS/xparam.h
../FMS/ggbi.o:		../FMS/ggbi.f    ../HEADERS/dim.h                       ../FMS/xparam.h
../FMS/gggm.o:		../FMS/gggm.f    ../HEADERS/dim.h                       ../FMS/xparam.h
../FMS/gglu.o:		../FMS/gglu.f    ../HEADERS/dim.h                       ../FMS/xparam.h
../FMS/ggrm.o:		../FMS/ggrm.f    ../HEADERS/dim.h                       ../FMS/xparam.h
../FMS/ggtf.o:		../FMS/ggtf.f    ../HEADERS/dim.h                       ../FMS/xparam.h
../FMS/xstaff.o:	../FMS/xstaff.f  ../HEADERS/dim.h                       ../FMS/xparam.h
../FMS/yprep.o:		../FMS/yprep.f   ../HEADERS/dim.h                       ../FMS/xparam.h

../FF2X/xscorr.o:	../FF2X/xscorr.f ../HEADERS/const.h ../HEADERS/dim.h

../EXCH/csigma.o:	../EXCH/csigma.f ../HEADERS/const.h
../EXCH/csigz.o:	../EXCH/csigz.f  ../HEADERS/const.h
../EXCH/edp.o:		../EXCH/edp.f    ../HEADERS/const.h
../EXCH/imhl.o:		../EXCH/imhl.f   ../HEADERS/const.h
../EXCH/quinn.o:	../EXCH/quinn.f  ../HEADERS/const.h
../EXCH/rhl.o:		../EXCH/rhl.f    ../HEADERS/const.h
../EXCH/rhlbp.o:	../EXCH/rhlbp.f  ../HEADERS/const.h
../EXCH/xcpot.o:	../EXCH/xcpot.f  ../HEADERS/const.h ../HEADERS/dim.h

../FOVRG/dfovrg.o:	../FOVRG/dfovrg.f ../HEADERS/dim.h ../HEADERS/const.h
../FOVRG/diff.o:	../FOVRG/diff.f   ../HEADERS/dim.h
../FOVRG/inmuac.o:	../FOVRG/inmuac.f ../HEADERS/dim.h
../FOVRG/intout.o:	../FOVRG/intout.f ../HEADERS/dim.h ../HEADERS/const.h
../FOVRG/muatcc.o:	../FOVRG/muatcc.f ../HEADERS/dim.h
../FOVRG/nucdec.o:	../FOVRG/nucdec.f ../HEADERS/dim.h
../FOVRG/potdvp.o:	../FOVRG/potdvp.f ../HEADERS/dim.h
../FOVRG/potex.o:	../FOVRG/potex.f  ../HEADERS/dim.h
../FOVRG/solin.o:	../FOVRG/solin.f  ../HEADERS/dim.h ../HEADERS/const.h
../FOVRG/solout.o:	../FOVRG/solout.f ../HEADERS/dim.h ../HEADERS/const.h
../FOVRG/wfirdc.o:	../FOVRG/wfirdc.f ../HEADERS/dim.h
../FOVRG/yzkrdc.o:	../FOVRG/yzkrdc.f ../HEADERS/dim.h
../FOVRG/yzktec.o:	../FOVRG/yzktec.f ../HEADERS/dim.h

../PAR/sequential.o:	../PAR/sequential.f ../HEADERS/parallel.h
#../PAR/parallel.o:	../PAR/parallel.f   ../HEADERS/parallel.h mpif.h

../XSPH/acoef.o:	../XSPH/acoef.f   ../HEADERS/dim.h
../XSPH/axafs.o:	../XSPH/axafs.f   ../HEADERS/dim.h ../HEADERS/const.h
../XSPH/fmssz.o:	../XSPH/fmssz.f   ../HEADERS/dim.h
../XSPH/phase.o:	../XSPH/phase.f   ../HEADERS/dim.h ../HEADERS/const.h
../XSPH/phmesh.o:	../XSPH/phmesh.f  ../HEADERS/dim.h ../HEADERS/const.h
../XSPH/radint.o:	../XSPH/radint.f  ../HEADERS/dim.h
../XSPH/rholsz.o:	../XSPH/rholsz.f  ../HEADERS/dim.h ../HEADERS/const.h
../XSPH/szlz.o:		../XSPH/szlz.f    ../HEADERS/dim.h ../HEADERS/const.h
../XSPH/wphase.o:	../XSPH/wphase.f  ../HEADERS/dim.h
../XSPH/wrxsph.o:	../XSPH/wrxsph.f  ../HEADERS/dim.h
../XSPH/xmult.o:	../XSPH/xmult.f                    ../HEADERS/const.h
../XSPH/xsect.o:	../XSPH/xsect.f   ../HEADERS/dim.h ../HEADERS/const.h
../XSPH/xsph.o:		../XSPH/xsph.f    ../HEADERS/dim.h ../HEADERS/const.h


afolp.o:	afolp.f    ../HEADERS/const.h ../HEADERS/dim.h
broydn.o:	broydn.f   ../HEADERS/const.h ../HEADERS/dim.h
corval.o:	corval.f   ../HEADERS/const.h ../HEADERS/dim.h
coulom.o:	coulom.f   ../HEADERS/const.h ../HEADERS/dim.h
fermi.o:	fermi.f    ../HEADERS/const.h
ff2g.o:		ff2g.f     ../HEADERS/const.h ../HEADERS/dim.h
ffmod1.o:	ffmod1.f                      ../HEADERS/dim.h ../HEADERS/parallel.h
ffsort_potph.o:	ffsort_potph.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/vers.h ../HEADERS/parallel.h
frnrm.o:	frnrm.f                       ../HEADERS/dim.h
grids.o:	grids.f    ../HEADERS/const.h
inipot.o:	inipot.f                      ../HEADERS/dim.h
inipotph.o:	inipotph.f ../HEADERS/const.h ../HEADERS/dim.h
istprm.o:	istprm.f   ../HEADERS/const.h ../HEADERS/dim.h
istval.o:	istval.f                      ../HEADERS/dim.h
libpotph.o:	libpotph.f ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h
libpotph_xsect.o:	libpotph_xsect.f      ../HEADERS/dim.h
moveh.o:	moveh.f    ../HEADERS/dim.h
movrlp.o:	movrlp.f   ../HEADERS/const.h ../HEADERS/dim.h
ovp2mt.o:	ovp2mt.f   ../HEADERS/const.h ../HEADERS/dim.h
ovrlp.o:	ovrlp.f    ../HEADERS/const.h ../HEADERS/dim.h
pot.o:		pot.f      ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h
reapot.o:	reapot.f   ../HEADERS/const.h ../HEADERS/dim.h
rholie.o:	rholie.f   ../HEADERS/const.h ../HEADERS/dim.h
scmt.o:		scmt.f     ../HEADERS/const.h ../HEADERS/dim.h
scmtmp.o:	scmtmp.f   ../HEADERS/const.h ../HEADERS/dim.h ../HEADERS/parallel.h
wpot.o:		wpot.f     ../HEADERS/const.h ../HEADERS/dim.h
wrpot.o:	wrpot.f                       ../HEADERS/dim.h


%.o:	%.f
	$(FORTRAN) $(FCFLAGS) $(FJSON) -o $@ $<

%.o:	%.c
	$(CC) $(CCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *$(SHOBJ) *.o *$(ARCHV) *$(SHOBJ) feff8l_pot$(EXEXT)

INSTARCH = @echo "Not installing $(ARCHV) files"
ifdef INSTALL_ARCHIVES
INSTARCH = $(COPY) libfeffint$(ARCHV) $(PREFIX)/lib
endif

install:
	$(COPY) feff8l_pot$(EXEXT)           $(BINDIR)
	$(COPY) libpotph$(SHOBJ)      $(LIBDIR)
	$(COPY) libfeff8lpotph$(SHOBJ) $(LIBDIR)
	$(COPY) feffphases.h          $(INCDIR)
	$(INSTARCH)

archives = ../POT/libfeffint$(ARCHV) ../ATOM/libfeffatom$(ARCHV)	\
../FOVRG/libfeffpha$(ARCHV) ../COMMON/libfeffcom$(ARCHV)		\
../PAR/libfeffpar$(ARCHV) ../MATH/libfeffmath$(ARCHV)			\
../JSON/libfeffjson$(ARCHV) ../EXCH/libfeffexch$(ARCHV)			\
../FMS/libfefffms$(ARCHV) ../json-fortran/libjsonfortran$(ARCHV)
# ../OPCONSAT/libfeffloss$(ARCHV)

otherdirs:
	$(MAKE) -C ../POT          libfeffint$(ARCHV)
	$(MAKE) -C ../ATOM         libfeffatom$(ARCHV)
	$(MAKE) -C ../FOVRG        libfeffpha$(ARCHV)
	$(MAKE) -C ../COMMON       libfeffcom$(ARCHV)
	$(MAKE) -C ../PAR          libfeffpar$(ARCHV)
	$(MAKE) -C ../MATH         libfeffmath$(ARCHV)
	$(MAKE) -C ../JSON         libfeffjson$(ARCHV)
	$(MAKE) -C ../EXCH         libfeffexch$(ARCHV)
	$(MAKE) -C ../FMS          libfefffms$(ARCHV)
	$(MAKE) -C ../json-fortran libjsonfortran$(ARCHV)
#	$(MAKE) -C ../OPCONSAT     libfeffloss$(ARCHV)

.PHONEY:	clean install all otherdirs
