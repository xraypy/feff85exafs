# directory contents: self-consistent density and potential calculation
# main routine: ffmod1.f

# customize compile command for different machines:
#   SGI with a lot of warnings and initialization check (trapuv)
#F77 = f77 -trapuv -fullwarn -C -c
 F77 = f77 -trapuv  -C 
#  SGI fast
#F77 = f77 -Ofast -LNO:opt=0 -c
#   Linux/g77
# F77 = g77 -Wall -O2

# define max number of processors for parallel execution (require MPI library)
# use 'make clean' if you edited line below to change dimensions properly
SUBMPI = 's/Maxprocs = 1/Maxprocs = 128/'

# output executable for sequential and parallel code
TARGET = ../../bin/Seq/pot
TARMPI = ../../bin/MPI/pot

# other macro definitions
SRCIST = istprm.f movrlp.f ovp2mt.f fermi.f sidx.f

SRC = ffmod1.f afolp.f broydn.f corval.f coulom.f ff2g.f frnrm.f \
grids.f inipot.f istval.f moveh.f ovrlp.f reapot.f rholie.f \
scmt.f sumax.f wpot.f wrpot.f $(SRCIST)

LIC = ../HEADERS/license.h 

OBJ = ${SRC:.f=.o} 
OBJIST = ${SRCIST:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h
PDATA = ../HEADERS/pdata.h
OBJEXT = ../MATH/lu.o 
LIB = ../LIB/libint.a

# how to make .o file from .f file (suffix rule)
.f.o :
	$(F77) -c $*.f

# Main targets: default, src, mono, mpi, srcmpi, clean

# compile sequential code
default:  $(TARGET) $(LIB)
$(TARGET): $(OBJ) pot.o scmtmp.o $(OBJEXT) ../LIB/libcom.a ../LIB/libmath.a \
  ../LIB/libatom.a ../LIB/libexch.a ../LIB/libpha.a ../LIB/libfms.a
	$(F77) $(OBJ) pot.o scmtmp.o $(OBJEXT) -o $(TARGET)  -L../LIB \
 -latom -lfms  -lpha -lcom -lmath  -lexch -lpar

# compile parallel code (output is ffmod1, not mpmod1)
mpi:  $(TARMPI) $(LIB)
 $(TARMPI): $(OBJ) pot_p.o scmtmp_p.o \
    $(OBJEXT) ../LIB/libcom.a ../LIB/libmath.a \
  ../LIB/libatom.a ../LIB/libexch.a ../LIB/libpha.a ../LIB/libfms.a
	$(F77) $(OBJ) pot_p.o scmtmp_p.o $(OBJEXT) -o $(TARMPI)  -L../LIB \
 -latom -lfms  -lpha -lcom -lmath  -lexch -lpar -lmpi

# make a library  of 'istprm' routines (parallel or sequential)
$(LIB): $(OBJIST)
	$(lbtool) $(LIB) $(OBJIST)

# instructions how to make pot_p.f and scmtmp_p.f for parallel code
pot_p.f: pot.f
	sed $(SUBMPI) pot.f > pot_p.f
scmtmp_p.f: scmtmp.f
	sed $(SUBMPI) scmtmp.f > scmtmp_p.f
	
# use before switching between parallel and sequential codes
clean:
	rm -f *.o *_tot.f *_p.f

# concatenate files for sequential modular code
src:
	cat $(SRCIST) > inc.f
	perl ../Utility/Uninclude -b inc.f > ist_tot.f
	rm inc.f
	cat $(LIC) $(SRC) pot.f scmtmp.f \
 ../ATOM/atom_tot.f ../FMS/fmslib_tot.f \
 ../FOVRG/fovrg_tot.f ../EXCH/exch_tot.f  ../MATH/lu.f \
 ../MATH/math_tot.f ../COMMON/common_tot.f \
 ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../pot_tot.f
	rm inc.f

# concatenate files for parallel modular code
srcmpi: pot_p.f scmtmp_p.f
	cat $(SRCIST) > inc.f
	perl ../Utility/Uninclude -b inc.f > ist_tot.f
	rm inc.f
	cat $(LIC) $(SRC) pot_p.f scmtmp_p.f \
 ../ATOM/atom_tot.f ../FMS/fmslib_tot.f \
 ../FOVRG/fovrg_tot.f ../EXCH/exch_tot.f  ../MATH/lu.f \
 ../MATH/math_tot.f ../COMMON/common_tot.f \
 ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../pot_tot.f
	rm inc.f

# concatenate files in this directory for sequential monolithic code
mono:
	cat $(SRC) pot.f scmtmp.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../pot_tot.f
	rm inc.f

# dependencies (for default and mpi):
 ffmod1.o: $(VERS)
 afolp.o: $(CONSTH) $(DIMH)
 broydn.o: $(CONSTH) $(DIMH)
 corval.o: $(CONSTH) $(DIMH)
 coulom.o: $(CONSTH) $(DIMH)
 fermi.o: $(CONSTH)
 ff2g.o: $(CONSTH) $(DIMH)
 frnrm.o: $(DIMH)
 grids.o: $(CONSTH)
 inipot.o: $(DIMH)
 istprm.o: $(CONSTH) $(DIMH)
 istval.o: $(DIMH)
 moveh.o: $(DIMH)
 movrlp.o: $(CONSTH) $(DIMH)
 ovp2mt.o: $(CONSTH) $(DIMH)
 ovrlp.o: $(CONSTH) $(DIMH)
 pot.o: $(CONSTH) $(DIMH)
 rholie.o: $(CONSTH) $(DIMH)
 scmt.o: $(CONSTH) $(DIMH)
 scmtmp.o: $(CONSTH) $(DIMH)
 wpot.o: $(CONSTH) $(DIMH)
 wrpot.o: $(DIMH)
