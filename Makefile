# Makefile for Feff85L

export PREFIX = ${CURDIR}/local_install

# installation location for programs, libraries, and include files:
export BINDIR  = $(PREFIX)/bin
export LIBDIR  = $(PREFIX)/lib
export INCDIR  = $(PREFIX)/include
export MAKEDIR = mkdir -p
export COPY    = cp
export REMOVE  = rm -rf

ifeq ($(OS),Windows_NT)
    MAKEDIST = sh ./makedist.sh
    DISTDIR = win64
    CMDEXT = bat
else
    MAKEDIST = sh ./makedist.sh
    DISTDIR = linux64
    CMDEXT = sh
    UNAME := $(shell uname -s)
    ifeq ($(UNAME),Darwin)
        DISTDIR = darwin64
    endif
endif

all:
	cd src && $(MAKE) all

install:
	$(MAKEDIR) $(BINDIR) $(LIBDIR) $(INCDIR)
	$(COPY) bin/feff8l.py $(BINDIR)/.
	$(COPY) bin/feff8l.$(CMDEXT) $(BINDIR)/.
	cd src && $(MAKE) install

clean:
	cd src && $(MAKE) clean

realclean: clean
	$(REMOVE) $(BINDIR) $(LIBDIR) $(INCDIR)

test: install
	cd tests && python run_tests.py

dist: install
	cd dist/$(DISTDIR) && $(MAKEDIST) $(PREFIX)

.PHONY:	all install clean test
#
