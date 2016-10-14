

export FORTRAN = gfortran
export SHARED  = -shared
export CFLAGS  = -c -O3 -ffree-line-length-none -g -Wall -fPIC
export FJSON   = -I/home/bruce/git/feff85exafs/src/json-fortran -J/home/bruce/git/feff85exafs/src/json-fortran
export AR      = ar
export ARFLAGS = rvc
export RANLIB  = ranlib
export RM      = rm

export CC      = gcc
export CCFLAGS = -c -g -fPIC

export COPY    = cp -v

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
	$(MAKE) -C src/XSPH
	$(MAKE) -C src/POT
	$(MAKE) -C src/PATH
	$(MAKE) -C src/FF2X
	$(MAKE) -C src/GENFMT


.PHONEY: 	all
