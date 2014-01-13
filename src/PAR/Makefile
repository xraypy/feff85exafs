# directory contents: routines for parallel code

#  customize compile command for different machines:
#   SGI
#F77 = f77 -trapuv -fullwarn -C
 F77 = f77 -trapuv  -C
# fast
#F77 = f77 -Ofast -LNO:opt=0
#   Linux/g77
# F77 = g77 -Wall -O2

#Macros
LIC = license.h
SRC = sequential.f
SRCMP = parallel.f

OBJ = ${SRC:.f=.o} 
OBJMP = ${SRCMP:.f=.o} 

LIB = ../LIB/libpar.a
LIBMP = ../LIB/lmppar.a
#lmppar.a currently is never created (see below)

# suffix rules
.f.o :
	$(F77) -c $*.f

# Main targets
# default target (first in list)
default: $(LIB) sequential.o 
 $(LIB): $(OBJ)
	$(lbtool) $(LIB)  $(OBJ)

# for parallel code (note: library is written in libpar.a, not lmppar.a)
mpi: $(LIBMP) parallel.o 
 $(LIBMP): $(OBJMP)
	$(lbtool) $(LIB)  $(OBJMP)
#	$(lbtool) $(LIBMP)  $(OBJMP)

clean:
	 rm -f *.o par_tot.f  $(LIB)
src:
	cat $(LIC) $(SRC) > par_inc.f
	perl ../Utility/Uninclude -b par_inc.f > par_tot.f
	rm par_inc.f
srcmpi:
	cat $(LIC) $(SRCMP) > par_inc.f
	perl ../Utility/Uninclude -b par_inc.f > par_tot.f
	rm par_inc.f

