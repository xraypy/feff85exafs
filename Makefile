
export PREFIX =  ${CURDIR}/local_install

# installation location for programs, libraries, and include files:
export BINDIR  = $(PREFIX)/bin
export LIBDIR  = $(PREFIX)/lib
export INCDIR  = $(PREFIX)/include

export MAKEDIR = mkdir -p
export COPY    = cp

all:
	cd src && $(MAKE) all

install:
	$(MAKEDIR) $(BINDIR) $(LIBDIR) $(INCDIR)
	$(COPY) bin/feff8l $(BINDIR)
	cd src && $(MAKE) install

clean:
	cd src && $(MAKE) clean

.PHONEY: 	all install clean
