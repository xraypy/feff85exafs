# directory contents: routines that used repeatedly in several modules
#   work with files: chopen.f nxtunt.f 
#   work with strings: str.f wlog.f head.f
#   packed ascii data (PAD) routines: padlib.f
#   read PAD formatted files: rdpot.f rdxsph.f
#   simple functions: xx.f getorb.f getxk.f pijump.f pertab.f
#   interpolation on fine grid: fixdsp.f fixdsx.f fixvar.f

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

# notice that padlib.f has it's own license, and therefore is listed last
SRC = chopen.f fixdsp.f fixdsx.f fixvar.f getorb.f getxk.f head.f \
itoken.f nxtunt.f pertab.f pijump.f rdhead.f rdpot.f rdxsph.f \
setkap.f str.f str2dp.f isnum.f wlog.f xx.f padlib.f rdcmt.f setgam.f \
iniptz.f qsorti.f #KJ added iniptz  1-06

LIC = ../HEADERS/license.h

OBJ = ${SRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h
LIB = ../LIB/libcom.a

# default target (first in list)
default: $(LIB) 

.f.o :
	$(F77) -c $*.f -lmath

# dependencies:
 fixdsp.o:  $(CONSTH) $(DIMH)
 fixdsx.o:  $(CONSTH) $(DIMH)
 fixvar.o:  $(CONSTH) $(DIMH)
 head.o:  $(CONSTH) $(DIMH) $(VERS)
 rdpot.o: $(DIMH)
 rdxsph.o: $(DIMH)
 pijump.o: $(CONSTH)
 padlib.o: padlib.h
 $(LIB): $(OBJ)
	$(lbtool) $(LIB)  $(OBJ)

clean:
	rm -f *.o *_tot.f
src:
	cat $(LIC) $(SRC) > common_inc.f
	perl ../Utility/Uninclude -b common_inc.f > common_tot.f
	rm common_inc.f

