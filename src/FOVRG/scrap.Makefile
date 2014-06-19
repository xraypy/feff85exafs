# directory contents: Dirac equation solver for complex potential
# main routine: dfovrg.f

# things to customize for different machines:
# compile command 
# Generic
# F77  = f77 -O1
#   SGI
#F77 = f77 -trapuv -fullwarn -C
#F77 = pgf77 -C
 F77 = f77 -C

# fast
#F77 = f77 -Ofast -LNO:opt=0
#   Linux/g77
# F77 = g77 -Wall -O2

LIC = ../HEADERS/license.h
SRC = aprdec.f aprdep.f dfovrg.f diff.f dsordc.f inmuac.f intout.f \
 muatcc.f nucdec.f ortdac.f potdvp.f potex.f solin.f solout.f wfirdc.f \
 yzkrdc.f yzktec.f  


OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
LIB = ../LIB/libpha.a

# default target (first in list)
default: $(LIB) 

.f.o :
	$(F77) -c $*.f

# dependencies:
 dfovrg.o:  $(CONSTH) $(DIMH)
 diff.o:  $(DIMH)
 dsordc.o:  $(DIMH)
 inmuac.o:  $(DIMH)
 intout.o:  $(CONSTH) $(DIMH)
 muatcc.o:  $(DIMH)
 nucdec.o:  $(DIMH)
 ortdac.o:  $(DIMH)
 potdvp.o:  $(DIMH)
 potex.o:  $(DIMH)
 solin.o:  $(CONSTH) $(DIMH)
 solout.o:  $(CONSTH) $(DIMH)
 wfirdc.o:  $(DIMH)
 yzkrdc.o:  $(DIMH)
 yzktec.o:  $(DIMH)
 $(LIB): $(OBJ)
	$(lbtool) $(LIB)  $(OBJ)

clean:
	rm -f *.o  *_tot.f
src:
	cat $(LIC) $(SRC) > common_inc.f
	perl ../Utility/Uninclude -b common_inc.f > fovrg_tot.f
	rm common_inc.f

