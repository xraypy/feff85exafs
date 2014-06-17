# directory contents: self-energy (energy dep exchage-correlation) routines
# main routine: xcpot.f

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

SRC = xcpot.f cubic.f edp.f ffq.f imhl.f quinn.f rhl.f rhlbp.f vbh.f csigma.f csigz.f \
fndsng.f 

LIC = ../HEADERS/license.h

OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
LIB = ../LIB/libexch.a

# default target (first in list)
default: $(LIB) 

.f.o :
	$(F77) -c $*.f

# dependencies:
 xcpot.o:  $(CONSTH) $(DIMH)
 edp.o:  $(CONSTH)
 imhl.o:  $(CONSTH)
 quinn.o:  $(CONSTH)
 rhl.o:  $(CONSTH)
 rhlbp.o:  $(CONSTH)
 $(LIB): $(OBJ)
	$(lbtool) $(LIB)  $(OBJ)

clean:
	rm -f *.o  *_tot.f
src:
	cat $(LIC) $(SRC) > common_inc.f
	perl ../Utility/Uninclude -b common_inc.f > exch_tot.f
	rm common_inc.f

