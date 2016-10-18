
## on linux, /usr/local/lib should be in the load path already.  elsewhere, you may need to set this
#LDPATH = -I/usr/local/lib
LDPATH = 

#########################################################################################
# if using this option, run the programs as 					        #
#   LD_LIBRARY_PATH="../../src/GENFMT:../../src/POT:../../src/json-fortran" ./makepotph #
#   LD_LIBRARY_PATH="../../src/GENFMT:../../src/POT:../../src/json-fortran" ./makepath  #
#########################################################################################
#LDPATH = -L../../src/GENFMT -L../../src/POT -L../../src/json-fortran


CC = gcc
CCFLAGS = -c -g

all:	pherr makephases errors makepath

makephases:	makephases.o
	$(CC) -o $@ $< $(LDPATH) -lfeffphases -lpotph

pherr:	pherr.o
	$(CC) -o $@ $< $(LDPATH) -lfeffphases -lpotph

makepath:	makepath.o
	$(CC) -o $@ $< $(LDPATH) -lfeffpath -lonepath

errors:		errors.o
	$(CC) -o $@ $< $(LDPATH) -lfeffpath -lonepath

%.o:	%.c
	$(CC) $(CCFLAGS) -o $@ $<

clean:
	rm -f makephases.o makephases makepath.o makepath errors.o errors pherr.o pherr

.PHONEY:	all clean
