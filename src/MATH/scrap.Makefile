# directory contents: mathematical routines for bessel fuctions,
#   integration, interpolation, 3j symbols, rotation matrices, etc..

# things to customize for different machines:
# compile command 
# Generic
# F77  = f77 -O1
#   SGI
#F77 = f77 -trapuv -fullwarn -C
 F77 = f77 -trapuv  -C

# fast
#F77 = f77 -Ofast -LNO:opt=0
#   Linux/g77
# F77 = g77 -Wall -O2

LIC =  ../HEADERS/license.h
SRC =  bcoef.f besjn.f besjh.f bjnser.f conv.f cpl0.f csomm.f csomm2.f \
 cwig3j.f determ.f  dist.f  rotwig.f \
 phamp.f exjlnl.f polint.f sdist.f somm.f somm2.f strap.f terp.f terpc.f trap.f czeros.f

OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
LIB = ../LIB/libmath.a

# default target (first in list)
default: $(LIB) lu.o 

.f.o :
	$(F77) -c $*.f

# dependencies:
 conv.o:  $(CONSTH) $(DIMH)
 phamp.o:  $(CONSTH) $(DIMH)

 $(LIB): $(OBJ)
	$(lbtool) $(LIB)  $(OBJ)

clean:
	rm -f *.o *_tot.f
src:
	cat $(LIC) $(SRC) > math_inc.f
	perl ../Utility/Uninclude -b math_inc.f > math_tot.f
	rm math_inc.f

mono:
	cat $(LIC) $(SRC) lu.f > math_inc.f
	perl ../Utility/Uninclude -b math_inc.f > math_tot.f
	rm math_inc.f
