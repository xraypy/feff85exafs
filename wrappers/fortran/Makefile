
## on linux, /usr/local/lib should be in the load path already.  elsewhere, you may need to set this
#LDPATH = -I/usr/local/lib
LDPATH = 

#########################################################################################
# if using this option, run the programs as 					        #
#   LD_LIBRARY_PATH="../../src/GENFMT:../../src/POT:../../src/json-fortran" ./makepotph #
#   LD_LIBRARY_PATH="../../src/GENFMT:../../src/POT:../../src/json-fortran" ./makepath  #
#########################################################################################
#LDPATH = -L../../src/GENFMT -L../../src/POT -L../../src/json-fortran


FORTRAN = gfortran
FCFLAGS = -c -g

all:	makepotph makepath

makepotph:	makepotph.o
	$(FORTRAN) -o $@ $< $(LDPATH) -lpotph -ljsonfortran

makepath:	makepath.o
	$(FORTRAN) -o $@ $< $(LDPATH) -lonepath

%.o:	%.f
	$(FORTRAN) $(FCFLAGS) -o $@ $<

clean:
	rm -f makepotph.o makepotph makepath.o makepath

.PHONEY:	all clean
