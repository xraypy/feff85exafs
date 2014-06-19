# directory contents: single configuration Dirac-Fock atomic code
# main subroutine: scfdat.f

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

# SRC + ../MATH/determ.f is a self-contained code
SRC = scfdat.f akeato.f  fdmocc.f  intdir.f  nucdev.f  potrdf.f  \
yzkrdf.f aprdev.f  dentfa.f  fdrirk.f  lagdat.f  ortdat.f  potslw.f \
tabrat.f  yzkteg.f bkmrdf.f  dsordf.f  fpf0.f    messer.f \
s02at.f   vlda.f cofcon.f  etotal.f  inmuat.f  muatco.f  soldir.f  wfirdf.f

CRIGHT = ../HEADERS/license.h 

AUX = ../COMMON/getorb.f ../COMMON/wlog.f ../COMMON/str.f ../MATH/cwig3j.f \
 ../MATH/somm.f ../MATH/determ.f ../EXCH/vbh.f ../EXCH/edp.f ../PAR/sequential.f

OBJ = ${SRC:.f=.o}  
LIB = ../LIB/libatom.a

DIMH = ../HEADERS/dim.h
CONSTH = ../HEADERS/const.h

# default target 
default: $(LIB)

.f.o :
	$(F77) -c $*.f

# dependencies:
 scfdat.o:  $(DIMH)
 vlda.o:  $(CONSTH)
 fpf0.o:  $(CONSTH) $(DIMH)
 atom.f: $(SRC)  $(CONSTH) $(DIMH)
	cat $(CRIGHT) $(SRC) > atom.f
 $(LIB): $(OBJ)
	$(lbtool) ../LIB/libatom.a $(OBJ)

clean:
	rm -f *.o  *_tot.f

src:
	cat $(CRIGHT) $(SRC) > common_inc.f
	perl ../Utility/Uninclude -b common_inc.f > atom_tot.f
	rm common_inc.f

atom:
	cat atom.f $(CRIGHT) $(SRC) $(AUX) > atom_temp.f
	perl ../Utility/Uninclude -b atom_temp.f > attot.f
	rm atom_temp.f

