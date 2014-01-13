# directory contents: pathfinder for multiple scattering expansion
# main routine: ffmod4.f

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
SRC = ffmod4.f ccrit.f heap.f ipack.f mcrith.f mcritk.f mpprmd.f \
 mpprmp.f mrb.f outcrt.f repath.f paths.f pathsd.f phash.f prcrit.f \
 sortix.f timrep.f

OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h

TARGET = ../../bin/Seq/path
TARMPI = ../../bin/MPI/path

#suffix rules
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
	perl ../Utility/Uninclude -b inc.f > ../path_tot.f
	rm inc.f

mono:
	cat $(LIC) $(SRC) > inc.f
	perl ../Utility/Uninclude -b inc.f > ../path_tot.f
	rm inc.f

# dependencies:
 ffmod4.o: $(CONSTH)
 ccrit.o: $(CONSTH) $(DIMH)
 mcrith.o: $(CONSTH) $(DIMH)
 mcritk.o: $(CONSTH) $(DIMH)
 mpprmd.o: $(DIMH)
 mpprmp.o: $(DIMH)
 mrb.o: $(DIMH)
 outcrt.o: $(CONSTH) $(DIMH)
 repath.o: $(DIMH)
 paths.o: $(DIMH)
 pathsd.o: $(CONSTH) $(DIMH)
 phash.o: $(DIMH)
 prcrit.o: $(CONSTH) $(DIMH)
 timrep.o: $(DIMH)
 
