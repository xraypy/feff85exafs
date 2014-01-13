# directory content: final calculations for various spectroscopies
# main program: ffmod6.f

# things to customize for different machines:
# compile command 
# Generic
# F77  = f77 -O1
#   SGI
#F77 = f77 -trapuv -fullwarn -C -c
 F77 = f77 -trapuv  -C

# fast
#F77 = f77 -Ofast -LNO:opt=0 -c
#   Linux/g77
# F77 = g77 -Wall -O2

SRC = ffmod6.f exconv.f feffdt.f ff2afs.f ff2chi.f ff2gen.f ff2xmu.f \
fprime.f rdfbin.f reff2x.f xscorr.f

LIC = ../HEADERS/license.h

OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h

TARGET = ../../bin/Seq/ff2x
TARMPI = ../../bin/MPI/ff2x

#suffix rules
.f.o :
	$(F77) -c $*.f

# default target (first in list)
default:  $(TARGET)
 $(TARGET): $(OBJ)  ../LIB/libdw.a ../LIB/libcom.a ../LIB/libmath.a
	$(F77) $(OBJ) -o $(TARGET)  -L../LIB -ldw -lcom -lmath -lpar

mpi:  $(TARMPI)
# output is always ffmod6, and not mpmod6
 $(TARMPI): $(OBJ)  ../LIB/libdw.a ../LIB/libcom.a ../LIB/libmath.a
	$(F77) $(OBJ) -o $(TARMPI)  -L../LIB -ldw -lcom -lmath -lpar -lmpi

clean:
	rm -f *.o

src:
	cat $(LIC) $(SRC) ../DEBYE/debye_tot.f \
 ../COMMON/common_tot.f  ../MATH/math_tot.f \
 ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../ff2x_tot.f
	rm inc.f

mono:
	cat $(LIC) $(SRC) > inc.f
	perl ../Utility/Uninclude -b inc.f > ../ff2x_tot.f
	rm inc.f


# dependencies: need updating
 ffmod6.o:  $(DIMH)
 reff2x.o:  $(DIMH) $(CONSTH)
