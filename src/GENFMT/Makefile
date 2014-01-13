# directory contents: scattering F-matrix multiplication for each path
# main routine: ffmod5.f

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
SRC = ffmod5.f fmtrxi.f genfmt.f mmtr.f mmtrxi.f rdpath.f regenf.f \
rot3i.f sclmz.f setlam.f snlm.f xstar.f 

OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h
PDATA = pdata.h

TARGET = ../../bin/Seq/genfmt
TARMPI = ../../bin/MPI/genfmt

# suffix rules
.f.o :
	$(F77) -c $*.f

# default target (first in list)
default:  $(TARGET)
 $(TARGET): $(OBJ)  ../LIB/libcom.a ../LIB/libmath.a
	$(F77) $(OBJ) -o $(TARGET)  -L../LIB -lcom -lmath -lpar

mpi:  $(TARMPI)
 $(TARMPI): $(OBJ)  ../LIB/libcom.a ../LIB/libmath.a
	$(F77) $(OBJ) -o $(TARMPI)  -L../LIB -lcom -lmath -lpar -lmpi

clean:
	rm -f *.o

src:
	cat $(LIC) $(SRC) \
 ../MATH/math_tot.f ../COMMON/common_tot.f \
 ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../genfmt_tot.f
	rm inc.f

mono:
	cat $(LIC) $(SRC) > inc.f
	perl ../Utility/Uninclude -b inc.f > ../genfmt_tot.f
	rm inc.f

# dependencies:
 ffmod5.o: $(DIMH)
 fmtrxi.o: $(DIMH) $(PDATA) nlm.h lambda.h clmz.h fmatrx.h rotmat.h
 genfmt.o: $(CONSTH) $(DIMH) $(VERS) $(PDATA) \
            nlm.h lambda.h clmz.h fmatrx.h rotmat.h 
 mmtr.o: $(CONSTH) $(DIMH) $(PDATA) rotmat.h
 mmtrxi.o: $(CONSTH) $(DIMH) $(PDATA) \
            nlm.h lambda.h clmz.h fmatrx.h rotmat.h
 rdpath.o: $(CONSTH) $(DIMH) $(PDATA)
 regenf.o: $(DIMH)
 rot3i.o: $(DIMH) $(PDATA) rotmat.h
 sclmz.o:  $(CONSTH) $(DIMH) clmz.h
 setlam.o: $(CONSTH) $(DIMH) $(PDATA) lambda.h
 snlm.o: $(DIMH) nlm.h
 xstar.o:  $(CONSTH) $(DIMH) $(PDATA)

