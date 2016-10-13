

FORTRAN = gfortran
SHARED  = -shared
CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
AR      = ar
ARFLAGS = rc
RANLIB  = ranlib
RM      = rm

CC      = gcc
CCFLAGS = -c -g -fPIC


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

archives = ../POT/libfeffint.a ../ATOM/libfeffatom.a			\
../FOVRG/libfeffpha.a ../COMMON/libfeffcom.a ../PAR/libfeffpar.a	\
../MATH/libfeffmath.a ../JSON/libfeffjson.a ../EXCH/libfeffexch.a	\
../FMS/libfefffms.a ../json-fortran/libjsonfortran.a

binobj = ffmod1.o pot.o scmtmp.o afolp.o broydn.o corval.o coulom.o	\
ff2g.o frnrm.o grids.o inipot.o istval.o moveh.o ovrlp.o reapot.o	\
rholie.o scmt.o sumax.o wpot.o wrpot.o

shobjects = $(objects) $(RDINPobj) $(COMMONobj) $(ATOMobj) $(MATHobj) $(FMSobj) $(FF2Xobj) $(EXCHobj) $(FOVRGobj) $(PARobj) $(XSPHobj)

JSONobj = ../JSON/json_read_libpotph.o ../JSON/json_xsect.o ../JSON/bailout.o



all:	libfeffint.a libpotph.so libfeffphases.so pot

libfeffint.a:	istprm.o movrlp.o ovp2mt.o fermi.o sidx.o
	$(AR) $(ARFLAGS) libfeffint.a istprm.o movrlp.o ovp2mt.o fermi.o sidx.o
	$(RANLIB) libfeffint.a

libpotph.so:	$(shobjects) $(JSONobj)
	$(FORTRAN) $(SHARED) $(shobjects) $(JSONobj) -o $@

libfeffphases.so:	feffphases.o libpotph.so
	$(FORTRAN) $(SHARED) feffphases.o libpotph.so -o $@


pot:	$(binobj) $(archives)
	$(FORTRAN) $(binobj) $(archives) -o $@


%.o:	%.f
	$(FORTRAN) $(CFLAGS) $(FJSON) -o $@ $<

%.o:	%.c
	$(CC) $(CCFLAGS) $(FJSON) -o $@ $<


clean:
	$(RM) *.o *.a *.so pot

.PHONEY:	clean all
