PREFIX  = /usr/local

FORTRAN = gfortran
FCFLAGS = -c -g

all:	makepotph makepath

makepotph:	makepotph.o
	$(FORTRAN) -o $@ $< -lpotph -ljsonfortran

makepath:	makepath.o
	$(FORTRAN) -o $@ $< -lonepath

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<

clean:
	rm -f makepotph.o makepotph makepath.o makepath

.PHONEY:	all clean
