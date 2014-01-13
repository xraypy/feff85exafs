# directory contents: calculation of Debye-Waller factors
#   sigms.f: correlated Debye model
#   sigrem.f: recursion and equation of motion methods

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

SRC =  sigms.f sigrem.f sigte3.f sigm3.f sigcl.f

LIC =  ../HEADERS/license.h

OBJ = ${SRC:.f=.o} 

CONSTH = ../HEADERS/const.h
DW = dwpar.h
LIB = ../LIB/libdw.a

# default target (first in list)
default: $(LIB) 

.f.o :
	$(F77) -c $*.f

# dependencies:
 sigms.o:  $(CONSTH)
 sigrem.o: $(DW)
 $(LIB): $(OBJ)
	$(lbtool) $(LIB)  $(OBJ)

clean:
	rm -f *.o *_tot.f

src:
	cat $(LIC) $(SRC) > common_inc.f
	perl ../Utility/Uninclude -b common_inc.f > debye_tot.f
	rm common_inc.f
