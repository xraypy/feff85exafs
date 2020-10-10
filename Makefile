# Makefile for Feff85L

export PREFIX = ${CURDIR}/local_install

# installation location for programs, libraries, and include files:
export BINDIR  = $(PREFIX)/bin
export LIBDIR  = $(PREFIX)/lib
export INCDIR  = $(PREFIX)/include
export MAKEDIR = mkdir -p
export COPY    = cp
export REMOVE  = rm -rf

all:
	cd src && $(MAKE) all

install:
	$(MAKEDIR) $(BINDIR) $(LIBDIR) $(INCDIR)
	$(COPY) bin/feff8l $(BINDIR)
	cd src && $(MAKE) install

clean:
	cd src && $(MAKE) clean

realclean: clean
	$(REMOVE) $(BINDIR) $(LIBDIR) $(INCDIR)

test: install
	cd tests && python run_tests.py

.PHONY:	all install clean test
#
