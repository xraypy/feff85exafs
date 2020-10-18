export PREFIX  = /usr/local
export PREFIX =  ${CURDIR}/../local_install

export BASEDIR = ${CURDIR}
# installation location for programs, libraries, and include files:
export BINDIR  = $(PREFIX)/bin
export LIBDIR  = $(PREFIX)/lib
export INCDIR  = $(PREFIX)/include

###########################################################################################
# gcc (tested on Ubuntu linux with gcc 5.4.0, darwin with gcc 6.3)                        #
###########################################################################################

## compile Feff's Fortran, including the F90 bits
export FORTRAN  = gfortran
export SHARED   = -shared
export FCFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC -fno-automatic
export FJSON    = -I${CURDIR}/json-fortran -J${CURDIR}/json-fortran
export FJSON    = -I${CURDIR}/json-fortran
export FLINKARGS =

## compile the C wrappers for phases and paths
export CC       = gcc
export CCFLAGS  = -c -g -fPIC

## compile json-fortran, which is F2008
export F90      = gfortran
export F90FLAGS = -std=f2008 -c -O2 -fbacktrace -g -Wall -Wextra -Wno-maybe-uninitialized -pedantic -fPIC
###########################################################################################

###########################################################################################
## gnu tools on linux
export AR      = ar
export ARFLAGS = rvc
export RANLIB  = ranlib
export RM      = rm -f
export COPY    = cp
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
        FLINKARGS =  -headerpad_max_install_names
    endif
endif
###########################################################################################

all:
	$(MAKE) -C PAR -j4
	$(MAKE) -C COMMON -j4
	$(MAKE) -C json-fortran -j1
	$(MAKE) -C JSON -j4
	$(MAKE) -C MATH -j4
	$(MAKE) -C ATOM -j4
	$(MAKE) -C DEBYE -j4
	$(MAKE) -C EXCH -j4
	$(MAKE) -C FOVRG -j4
	$(MAKE) -C FMS -j4
	$(MAKE) -C RDINP -j4
#	$(MAKE) -C OPCONSAT -j4
	$(MAKE) -C POT -j4
	$(MAKE) -C XSPH -j4
	$(MAKE) -C PATH -j4
	$(MAKE) -C FF2X -j4
	$(MAKE) -C GENFMT -j4
	$(MAKE) -C feff6l -j4

install: all
	$(MAKE) -C ATOM   -j4 install
	$(MAKE) -C COMMON -j4 install
	$(MAKE) -C DEBYE  -j4 install
	$(MAKE) -C EXCH   -j4 install
	$(MAKE) -C FF2X   -j4 install
	$(MAKE) -C FMS    -j4 install
	$(MAKE) -C FOVRG  -j4 install
	$(MAKE) -C GENFMT -j4 install
	$(MAKE) -C JSON   -j4 install
	$(MAKE) -C MATH   -j4 install
#	$(MAKE) -C OPCONSAT -j4 install
	$(MAKE) -C PAR    -j4 install
	$(MAKE) -C PATH   -j4 install
	$(MAKE) -C POT    -j4 install
	$(MAKE) -C RDINP  -j4 install
	$(MAKE) -C XSPH   -j4 install
	$(MAKE) -C json-fortran -j1 install
	$(MAKE) -C feff6l -j4 install


clean:
	$(MAKE) -C ATOM   clean
	$(MAKE) -C COMMON clean
	$(MAKE) -C DEBYE  clean
	$(MAKE) -C EXCH   clean
	$(MAKE) -C FF2X   clean
	$(MAKE) -C FMS    clean
	$(MAKE) -C FOVRG  clean
	$(MAKE) -C GENFMT clean
	$(MAKE) -C JSON   clean
	$(MAKE) -C MATH   clean
#	$(MAKE) -C OPCONSAT clean
	$(MAKE) -C PAR    clean
	$(MAKE) -C PATH   clean
	$(MAKE) -C POT    clean
	$(MAKE) -C RDINP  clean
	$(MAKE) -C XSPH   clean
	$(MAKE) -C json-fortran clean
##	$(MAKE) -C feff6l clean

.PHONEY: 	all install clean
