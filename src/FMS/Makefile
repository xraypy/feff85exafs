# directory contents: full multiple scattering routines
# main routines: ffmod3.f and fmsie.f

# things to customize for different machines:
# compile command 
# Generic
# F77  = f77 -O1
#   SGI
#F77 = f77 -trapuv -fullwarn -C -c
 F77 = f77 

# fast
#F77 = f77 -Ofast -LNO:opt=0 -c
#   Linux/g77
# F77 = g77 -Wall -O2

LIC =  ../HEADERS/license.h 
LSRC = fmsie.f fmspack.f gglu.f ggbi.f ggrm.f gggm.f ggtf.f yprep.f xstaff.f
SRC =  ffmod3.f fmspack.f gglu.f ggbi.f ggrm.f gggm.f ggtf.f \
       fmstot.f reafms.f xprep.f xstaff.f 
MONO = license.h fmsie.f fmspack.f gglu.f ggbi.f ggrm.f gggm.f ggtf.f \
       fmstot.f reafms.f xprep.f xstaff.f yprep.f 
LIGHT = license.h fmsie.f fmspack.f gglu.f ggbi.f ggrm.f gggm.f ggtf.f \
         xstaff.f yprep.f

OBJ = ${SRC:.f=.o} 
LOBJ = ${LSRC:.f=.o} 

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h
VERS = ../HEADERS/vers.h
XPAR = ../FMS/xparam.h

LIBFMS = ../LIB/libfms.a
TARGET = ../../bin/Seq/fms
TARMPI = ../../bin/MPI/fms
LIB = ../LIB/libdw.a ../LIB/libcom.a ../LIB/libmath.a
LIBF77 =  -L../LIB  -ldw -lcom -lmath

#suffix rules
.f.o :
	$(F77) -c $*.f

# default target (first in list)
default: $(LIBFMS)

 $(TARGET): $(OBJ) $(LIB)
	$(F77) $(OBJ) ../MATH/lu.o  ../LIB/libpar.a $(LIBF77) -o $(TARGET) 
 $(LIBFMS): $(LOBJ)
	$(lbtool) ../LIB/libfms.a $(LOBJ)

mpi:  $(LIBFMS) $(TARMPI)
# note: ouput is written in ffmod3 (not TARMPI - mpmod3) 
$(TARMPI): $(OBJ) $(LIB)
	$(F77) $(OBJ) ../MATH/lu.o ../LIB/libpar.a $(LIBF77) -lmpi -o $(TARMPI) 

clean:
	rm -f *.o  *_tot.f

src:
	cat $(LSRC) >inc.f
	perl ../Utility/Uninclude -b inc.f > fmslib_tot.f
	rm inc.f
	cat $(LIC) $(SRC) ../MATH/lu.f ../DEBYE/debye_tot.f \
 ../MATH/math_tot.f ../COMMON/common_tot.f \
 ../PAR/par_tot.f > inc.f
	perl ../Utility/Uninclude -b inc.f > ../fms_tot.f
	rm inc.f

mono:
	cat $(MONO) >inc.f
	perl ../Utility/Uninclude -b inc.f > ../fms_tot.f
	rm inc.f
light:
	cat $(LIGHT) >inc.f
	perl ../Utility/Uninclude -b inc.f > ../fms_tot.f
	rm inc.f

# dependencies:
 ffmod3.o: $(DIMH)
 fmsie.o:  $(DIMH)
 fmspack.o:  $(XPAR) $(DIMH)
 gglu.o:     $(XPAR) $(DIMH)
 ggbi.o:     $(XPAR) $(DIMH)
 ggrm.o:     $(XPAR) $(DIMH)
 gggm.o:     $(XPAR) $(DIMH)
 ggtf.o:     $(XPAR) $(DIMH)
 fmstot.o:  $(CONSTH) $(DIMH)
 reafms.o:  $(CONSTH) $(DIMH)
 xprep.o:  $(XPAR) $(DIMH)
 yprep.o:  $(XPAR) $(DIMH)
