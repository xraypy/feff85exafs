# directory contents: single input file (feff.inp) reader
# main routines: rdinp.f and ffsort.f

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

LIC = ../HEADERS/license.h
SRC = rdinp.f  ffsort.f iniall.f  mkptz.f  rdline.f  setedg.f \
       wrtall.f
SRCLIGHT = rdinp_l.f  ffsort.f iniall.f  mkptz.f  rdline.f  setedg.f \
       wrtall_l.f

OBJ = ${SRCLIGHT:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h
IORH = ../RDINP/allinp.h
PARH = ../HEADERS/parallel.h

TARGET = ../../bin/Seq/rdinp
TARMPI = ../../bin/MPI/rdinp
LIB = ../LIB/libcom.a ../LIB/libmath.a 
LIBF77 =  -L../LIB  -lcom -lmath -lpar

# suffix rules
.f.o :
	$(F77) -c $*.f

# default target (first in list)
default:  $(TARGET)

 $(TARGET): $(OBJ) $(LIB)
	$(F77) $(OBJ) -o $(TARGET)  $(LIBF77)

mpi:  $(TARMPI)
 $(TARMPI): $(OBJ) $(LIB)
	$(F77) $(OBJ) -o $(TARMPI)  $(LIBF77) -lmpi

clean:
	rm -f *.o 

src: 
	cat $(LIC) $(SRC) \
 ../MATH/math_tot.f ../COMMON/common_tot.f \
 ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../rdinp_tot.f
	rm inc.f

mono:
	cat $(LIC) $(SRC) > inc.f
	perl ../Utility/Uninclude -b inc.f > ../rdinp_tot.f
	rm inc.f

light:
	cat $(LIC) $(SRCLIGHT) > inc.f
	perl ../Utility/Uninclude -b inc.f > ../rdinp_tot.f
	rm inc.f

# dependencies:
 rdinp.o:  $(CONSTH) $(DIMH) $(VERS) $(IORH) $(PARH)
 ffsort.o: $(CONSTH) $(DIMH) $(PARH)
 iniall.o: $(DIMH) $(IORH)
 mkptz.o:  $(CONSTH)
 wrtall.o: $(DIMH) $(IORH) $(PARH)

