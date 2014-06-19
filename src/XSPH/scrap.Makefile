# directory contents: cross section and phase shifts calculation
# main routine: ffmod2.f

# customize compile command for different machines:
#   SGI
#F77 = f77 -trapuv -fullwarn -C -c
 F77 = f77 -trapuv  -C
# fast
#F77 = f77 -Ofast -LNO:opt=0 -c
#   Linux/g77
# F77 = g77 -Wall -O2

LIC = ../HEADERS/license.h
SRC = ffmod2.f axafs.f  phase.f phmesh2.f phmesh.f radint.f wphase.f wrxsph.f \
 xmult.f rexsph.f xsect.f xsph.f szlz.f acoef.f fmssz.f rholsz.f rholat.f \
 getedg.f rdgrid.f 

SRCX = $(SRC) ../FF2X/xscorr.f  ../POT/grids.f

OBJ = ${SRCX:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h
PDATA = ../HEADERS/pdata.h
#LFLAG = -ltd -lpha -lint -lexch -lfms -lcom -lmath -lpar 
LFLAG = -lpha -lint -lexch -lfms -lcom -lmath -lpar 

TARGET = ../../bin/Seq/xsph
TARMPI = ../../bin/MPI/xsph

#suffix rules
.f.o :
	$(F77) -c $*.f

# default target (first in list)
default:  $(TARGET)
$(TARGET): $(OBJ) ../LIB/libcom.a ../LIB/libmath.a ../MATH/lu.o \
 ../LIB/libfms.a ../LIB/libtd.a \
  ../LIB/libatom.a ../LIB/libexch.a ../LIB/libpha.a ../LIB/libint.a
	$(F77) $(OBJ) ../MATH/lu.o -o $(TARGET)  -L../LIB $(LFLAG)

mpi:$(TARMPI)
$(TARMPI): $(OBJ) ../LIB/libint.a ../LIB/libcom.a ../LIB/libmath.a \
 ../MATH/lu.o ../LIB/libfms.a ../LIB/libtd.a \
  ../LIB/libatom.a ../LIB/libexch.a ../LIB/libpha.a 
	$(F77) $(OBJ) ../MATH/lu.o -o $(TARMPI)  -L../LIB $(LFLAG) -lmpi

clean:
	rm -f *.o

src:
	cat $(LIC) $(SRCX) ../FOVRG/fovrg_tot.f \
 ../POT/ist_tot.f  ../EXCH/exch_tot.f ../FMS/fmslib_tot.f \
 ../MATH/math_tot.f ../COMMON/common_tot.f ../TDLDA/td_tot.f \
 ../MATH/lu.f ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../xsph_tot.f
	rm inc.f

mono:
	cat $(LIC) $(SRC) > inc.f
	perl ../Utility/Uninclude -b inc.f > ../xsph_tot.f
	rm inc.f

# dependencies:
 ffmod2.o: $(DIMH)
 axafs.o: $(CONSTH) $(DIMH)
 phase.o: $(CONSTH) $(DIMH)
 phmesh.o: $(CONSTH) $(DIMH)
 wphase.o: $(DIMH)
 wrxsph.o: $(DIMH)
 radint.o: $(DIMH)
 xmult.o: $(CONSTH)
 xsect.o: $(CONSTH) $(DIMH) 
 xsph.o: $(CONSTH) $(DIMH)
# szlz staff
 acoef.o:  $(DIMH)
 fmssz.o: $(DIMH)
 rholsz.o: $(CONSTH) $(DIMH)
 rholat.o: $(CONSTH) $(DIMH)
 szlz.o: $(CONSTH) $(DIMH)

