
export PREFIX  = /usr/local
export BASEDIR = ${CURDIR}
export BINDIR  = $(PREFIX)/bin		# installation location for programs
export LIBDIR  = $(PREFIX)/lib		# installation location for libraries
export INCDIR  = $(PREFIX)/include	# installation location for include files

###########################################################################################
# gcc (tested on Ubuntu linux with gcc 5.4.0)                                             #
###########################################################################################
export FORTRAN  = gfortran	## compile Feff's Fortran, including the F90 bits
export SHARED   = -shared
export FCFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
export FJSON    = -I$(BASEDIR)/src/json-fortran -J$(BASEDIR)/src/json-fortran

export CC       = gcc	 	## compile the C wrappers for phases and paths
export CCFLAGS  = -c -g -fPIC

export F90      = gfortran	## compile json-fortran, which is F2008
export F90FLAGS = -std=f2008 -c -O2 -fbacktrace -g -Wall -Wextra -Wno-maybe-uninitialized -pedantic -fPIC
###########################################################################################

###########################################################################################
## gnu tools on linux
export AR      = ar
export ARFLAGS = rvc
export RANLIB  = ranlib
export RM      = rm -f
export COPY    = cp -v
###########################################################################################

###########################################################################################
## file extensions, shared object libraries, object archives, and executables
export SHOBJ = .so
export ARCHV = .a
export EXEXT =
###########################################################################################


all:
	$(MAKE) -C src/PAR
	$(MAKE) -C src/COMMON
	$(MAKE) -C src/json-fortran
	$(MAKE) -C src/JSON
	$(MAKE) -C src/MATH
	$(MAKE) -C src/ATOM
	$(MAKE) -C src/DEBYE
	$(MAKE) -C src/EXCH
	$(MAKE) -C src/FOVRG
	$(MAKE) -C src/FMS
	$(MAKE) -C src/RDINP
#	$(MAKE) -C src/OPCONSAT
	$(MAKE) -C src/POT
	$(MAKE) -C src/XSPH
	$(MAKE) -C src/PATH
	$(MAKE) -C src/FF2X
	$(MAKE) -C src/GENFMT

install:
	$(MAKE) -C src/ATOM   install
	$(MAKE) -C src/COMMON install
	$(MAKE) -C src/DEBYE  install
	$(MAKE) -C src/EXCH   install
	$(MAKE) -C src/FF2X   install
	$(MAKE) -C src/FMS    install
	$(MAKE) -C src/FOVRG  install
	$(MAKE) -C src/GENFMT install
	$(MAKE) -C src/JSON   install
	$(MAKE) -C src/MATH   install
#	$(MAKE) -C src/OPCONSAT install
	$(MAKE) -C src/PAR    install
	$(MAKE) -C src/PATH   install
	$(MAKE) -C src/POT    install
	$(MAKE) -C src/RDINP  install
	$(MAKE) -C src/XSPH   install
	$(MAKE) -C src/json-fortran install

clean:
	$(MAKE) -C src/ATOM   clean
	$(MAKE) -C src/COMMON clean
	$(MAKE) -C src/DEBYE  clean
	$(MAKE) -C src/EXCH   clean
	$(MAKE) -C src/FF2X   clean
	$(MAKE) -C src/FMS    clean
	$(MAKE) -C src/FOVRG  clean
	$(MAKE) -C src/GENFMT clean
	$(MAKE) -C src/JSON   clean
	$(MAKE) -C src/MATH   clean
#	$(MAKE) -C src/OPCONSAT clean
	$(MAKE) -C src/PAR    clean
	$(MAKE) -C src/PATH   clean
	$(MAKE) -C src/POT    clean
	$(MAKE) -C src/RDINP  clean
	$(MAKE) -C src/XSPH   clean
	$(MAKE) -C src/json-fortran clean


.PHONEY: 	all install clean
