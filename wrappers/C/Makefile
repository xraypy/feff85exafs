
PREFIX  = /usr/local

CC = gcc
CCFLAGS = -c -g

all:	pherr makephases errors makepath

makephases:	makephases.o
	$(CC) -o $@ $< -lfeffphases -lpotph

pherr:	pherr.o
	$(CC) -o $@ $< -lfeffphases -lpotph

makepath:	makepath.o
	$(CC) -o $@ $< -lfeffpath -lonepath

errors:		errors.o
	$(CC) -o $@ $< -lfeffpath -lonepath

%.o:	%.c
	$(CC) $(CCFLAGS) -o $@ $<

clean:
	rm -f makephases.o makephases makepath.o makepath errors.o errors pherr.o pherr

.PHONEY:	all clean
