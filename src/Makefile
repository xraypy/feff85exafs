#Use commands:
#  make        - compiles sequential modular code (see script 'feff')
#  make mpi    - compiles parallel modular code (see script 'feffmpi')
#                edit POT/Makefile to change number of processors
#  make src    - generate sequential modular code  (files below)
#               ff2x_tot.f    genfmt_tot.f  path_tot.f    rdinp_tot.f
#               fms_tot.f     ldos_tot.f    pot_tot.f     xsph_tot.f
#  make srcmpi - generate parallel modular code  (file above)
#  make mono   - generate sequential monolithic code (feff_tot.f)
#  make clean  - to switch between parallel and sequential version,
#                also to clean directories from *.o files
# Makefile system written by A. Ankudinov 04/2001

# define global compiler below (use 'make clean' if compiler options changed)
# notice that this global value overwrites F77 in each subdirectory
# g77
#F77 = g77 -Wall -O2


# ifort
#F77 = g77 -g
F77 = ifort

# gfortran
F77 = gfortran -O3 -ffree-line-length-none -finit-local-zero

# Mac xlf
#F77 = xlf -qextern=trap

#Linux ifort
#F77 = ifort -O3

#fraangelico
#MPIFlags = -qextern=trap
#lbtool = libtool -static -o

#bernini
#MPIFlags = -O3 -xK -tpp7 -unroll:1
#lbtool = ar ru

#rodin
#MPIFlags = -static-libcxa -O3 -xK -tpp7 -unroll:1
MPIFlags = -fastsse
lbtool = ar ru

MPIF77 = mpif77 $(MPIFlags)

parmessg = "\n\n\n\tmodular src is in ../mod/MPI/\n\tPlease cd there and compile\n\twith the appropriate script."
messg = "\n\n\n\tmodular src is in ../mod/Seq/\n\tPlease cd there and compile\n\twith the appropriate script."
monmessg = "\n\n\n\tsrc is in ../mod/MONO/\n\tPlease cd there and compile\n\twith the appropriate script."
smessg1 = "\n\tBinaries are located in ../bin/"
smessg2 = "\n\tYou can use the script named\n\t\t"
smessg3 = "\n\tlocated in ../bin/ to run feff"

MODS = ff2x_tot.f  fms_tot.f  genfmt_tot.f  path_tot.f  pot_tot.f \
	rdinp_tot.f xsph_tot.f 

#MONO = PAR/par_tot.f COMMON/common_tot.f MATH/math_tot.f ATOM/atom_tot.f \
#    FOVRG/fovrg_tot.f EXCH/exch_tot.f DEBYE/debye_tot.f rdinp_tot.f pot_tot.f \
#    ldos_tot.f xsph_tot.f  fms_tot.f path_tot.f genfmt_tot.f ff2x_tot.f
MONO = rdinp_tot.f pot_tot.f xsph_tot.f \
    fms_tot.f path_tot.f genfmt_tot.f ff2x_tot.f \
    DEBYE/debye_tot.f ATOM/atom_tot.f  FOVRG/fovrg_tot.f  EXCH/exch_tot.f \
    MATH/math_tot.f  COMMON/common_tot.f PAR/par_tot.f
    
LIGHT = rdinp_tot.f pot_tot.f xsph_tot.f \
    fms_tot.f path_tot.f genfmt_tot.f ff2x_tot.f \
    DEBYE/debye_tot.f ATOM/atom_tot.f  FOVRG/fovrg_tot.f  EXCH/exch_tot.f \
    MATH/math_tot.f  COMMON/common_tot.f PAR/par_tot.f 


# default target (first in list)
default:  
#  library
	if grep -q "$(F77)" comp.last ; then true; else make -s clean; fi;
	cd PAR; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd COMMON ; make  "F77 = $(F77)" "lbtool = $(lbtool)"
	cd MATH ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd ATOM ;  make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd FOVRG ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd EXCH ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd DEBYE ; make "F77 = $(F77)" "lbtool = $(lbtool)"
# library + module
	cd FMS ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd POT ; make "F77 = $(F77)" "lbtool = $(lbtool)"
# module
	cd RDINP ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd FF2X ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd GENFMT ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd PATH ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	cd XSPH ; make "F77 = $(F77)" "lbtool = $(lbtool)"
	. Utility/MkFeffScript seq
	@echo  $(F77) > comp.last
	@echo -e $(smessg1)Seq$(smessg2)feff$(smessg3)

clean:
#  library
	cd PAR ; make clean
	cd ATOM ;  make clean
	cd COMMON ; make clean
	cd DEBYE ; make clean
	cd EXCH ; make clean
	cd MATH ; make clean
	cd FOVRG ; make clean
# library + module 
	cd FMS ; make clean
	cd POT ; make clean
# module
	cd FF2X ; make clean
	cd GENFMT ; make clean
	cd PATH ; make clean
	cd RDINP ; make clean
	cd XSPH ; make clean

mpi: 
#  library
	if grep -q "$(MPIF77)" comp.last ; then true; else make -s clean; fi;
	cd PAR ; make  mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd ATOM ;  make "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd COMMON ; make "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd DEBYE ; make "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd EXCH ; make "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd MATH ; make "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd FOVRG ; make "F77 = $(MPIF77)" "lbtool = $(lbtool)"
# library + module
	cd FMS ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd POT ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
# module
	cd FF2X ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd GENFMT ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd PATH ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd RDINP ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	cd XSPH ; make mpi "F77 = $(MPIF77)" "lbtool = $(lbtool)"
	echo "`hostname` cpu=2" > ../bin/MPI/mynodes 
	. Utility/MkFeffScript mpi
	@echo $(MPIF77) > comp.last
	@echo -e $(smessg1)MPI$(smessg2)feffmpi$(smessg3)

srcmpi: 
#  library
	cd PAR ; make srcmpi
	cd ATOM ;  make src
	cd COMMON ; make src
	cd DEBYE ; make src
	cd EXCH ; make src
	cd MATH ; make src
	cd FOVRG ; make src
# library + module
	cd FMS ; make src
	cd POT ; make srcmpi
# module
	cd FF2X ; make src
	cd GENFMT ; make src
	cd PATH ; make src
	cd RDINP ; make src
	cd XSPH ; make src
	cp $(MODS) ../mod/MPI/
	rm $(MODS)
	echo "`hostname` cpu=2" > ../bin/MPI/mynodes 
	. Utility/MkFeffScript mpi
	@echo -e $(parmessg)
mono: 
#  library
	cd PAR ; make src
	cd ATOM ;  make src
	cd COMMON ; make src
	cd DEBYE ; make src
	cd EXCH ; make src
	cd MATH ; make mono
	cd FOVRG ; make src
# library + module
	cd FMS ; make mono
	cd POT ; make mono
# module
	cd FF2X ; make  mono
	cd GENFMT ; make mono
	cd PATH ; make mono
	cd RDINP ; make mono
	cd XSPH ; make mono
	cat HEADERS/cright.h HEADERS/license.h HEADERS/feff.f \
 $(MONO) > feff_temp.f
	sed -f Utility/subpro feff_temp.f > feff_tot.f
	rm feff_temp.f
	mv feff_tot.f ../mod/MONO/
	rm *tot*
	. Utility/MkFeffScript
	@echo -e $(monmessg)

src:
#  library
	cd PAR ; make src
	cd ATOM ;  make src
	cd COMMON ; make src
	cd DEBYE ; make src
	cd EXCH ; make src
	cd MATH ; make mono
	cd FOVRG ; make src
# library + module
	cd FMS ; make light
	cd POT ; make mono
# module
	cd FF2X ; make  mono
	cd GENFMT ; make mono
	cd PATH ; make mono
	cd RDINP ; make light
	cd XSPH ; make mono
	cat HEADERS/cright.h HEADERS/license.h HEADERS/feff.f \
 $(LIGHT) > feff_temp.f
	sed -f Utility/subpro feff_temp.f > feff_tot.f
	rm feff_temp.f
	mv feff_tot.f ../mod/MONO/feff85L.f
	rm *tot*
	. Utility/MkFeffScript
	@echo -e $(lightmessg)
