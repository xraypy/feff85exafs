export PREFIX  = /usr/local
export PREFIX  = ${HOME}/local

export BASEDIR = ${CURDIR}
export BINDIR  = $(PREFIX)/bin		# installation location for programs
export LIBDIR  = $(PREFIX)/lib		# installation location for libraries
export INCDIR  = $(PREFIX)/include	# installation location for include files

###########################################################################################
# gcc (tested on Ubuntu linux with gcc 5.4.0, darwin with gcc 6.3)                        #
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
export MAKEDIR = mkdir -p
###########################################################################################

###########################################################################################
## file extensions for static object archives, shared object libraries, and executables
export ARCHV = .a
export SHOBJ = .so
export EXEXT =

ifeq ($(OS),Windows_NT)
    UNAME = Windows_NT
    SHOBJ = .dll
    EXEXT = .exe
    FCFLAGS  = -c -O3 -ffree-line-length-none -g -Wall
else
    UNAME := $(shell uname -s)
    ifeq ($(UNAME),Darwin)
        SHOBJ = .dylib
    endif
endif
###########################################################################################


all:
	$(MAKE) -C src/PAR -j4
	$(MAKE) -C src/COMMON -j4
	$(MAKE) -C src/json-fortran -j1
	$(MAKE) -C src/JSON -j4
	$(MAKE) -C src/MATH -j4
	$(MAKE) -C src/ATOM -j4
	$(MAKE) -C src/DEBYE -j4
	$(MAKE) -C src/EXCH -j4
	$(MAKE) -C src/FOVRG -j4
	$(MAKE) -C src/FMS -j4
	$(MAKE) -C src/RDINP -j4
#	$(MAKE) -C src/OPCONSAT -j4
	$(MAKE) -C src/POT -j4
	$(MAKE) -C src/XSPH -j4
	$(MAKE) -C src/PATH -j4
	$(MAKE) -C src/FF2X -j4
	$(MAKE) -C src/GENFMT -j4
	$(MAKE) -C src/feff6l -j4

install:
	$(MAKEDIR) $(BINDIR) $(LIBDIR) $(INCDIR)
	$(MAKE) -C src/ATOM   -j4 install
	$(MAKE) -C src/COMMON -j4 install
	$(MAKE) -C src/DEBYE  -j4 install
	$(MAKE) -C src/EXCH   -j4 install
	$(MAKE) -C src/FF2X   -j4 install
	$(MAKE) -C src/FMS    -j4 install
	$(MAKE) -C src/FOVRG  -j4 install
	$(MAKE) -C src/GENFMT -j4 install
	$(MAKE) -C src/JSON   -j4 install
	$(MAKE) -C src/MATH   -j4 install
#	$(MAKE) -C src/OPCONSAT -j4 install
	$(MAKE) -C src/PAR    -j4 install
	$(MAKE) -C src/PATH   -j4 install
	$(MAKE) -C src/POT    -j4 install
	$(MAKE) -C src/RDINP  -j4 install
	$(MAKE) -C src/XSPH   -j4 install
	$(MAKE) -C src/json-fortran -j1 install
	$(MAKE) -C src/feff6l -j4 install

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
	$(MAKE) -C src/feff6l clean

.PHONEY: 	all install clean
