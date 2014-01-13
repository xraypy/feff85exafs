# feff8 8.2 with MPI makefile
FFLAGC  = 
opt     = -fast
## This line is needed for use prefered LAM lib's, 11/21/00
LIBS    =  -L/usr/local/lampgi/lib -lmpi++ -lmpi -limpi -lmpi -limpi -ltstdio -ltrillium -largs -lt 
#LIBS    =  -L/usr/local/lampgi/lib -llammpio -llamf77mpi -lmpi -llam -lm
input   = 
CFLAGC  = -fast

## In order to work with pgi and LAM compiler -Msecond_underscore is needed 11/21/00
FCOMPL  = pgf77 -Msecond_underscore -I/usr/local/lampgi/include
#FCOMPL  = /usr/local/lam/bin/mpif77 -Msecond_underscore -I/usr/local/lampgi/include

CCOMPL  = cc 
# list of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .f .c .s .o

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) $<

.s.o:
	as $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $<

rdinp_src   = rdinp_tot.f
ffmod1_src  = pot_tot.f
ldos_src    = ldos_tot.f
ffmod2_src  = xsph_tot.f
ffmod3_src  = fms_tot.f
ffmod4_src  = path_tot.f
ffmod5_src  = genfmt_tot.f
ffmod6_src  = ff2x_tot.f

ALL:    $(rdinp_src) $(ffmod1_src) $(ldos_src) $(ffmod2_src) $(ffmod3_src) \
        $(ffmod4_src) $(ffmod5_src) $(ffmod6_src)
	$(FCOMPL) $(FFLAGC) -o ../bin/rdinp $(rdinp_src) $(LIBS)
	$(FCOMPL) $(FFLAGC) -o ../bin/ffmod1 $(ffmod1_src) $(LIBS)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ldos $(ldos_src) $(LIBS)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod2 $(ffmod2_src) $(LIBS)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod3 $(ffmod3_src) $(LIBS)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod4 $(ffmod4_src) $(LIBS)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod5 $(ffmod5_src) $(LIBS)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod6 $(ffmod6_src) $(LIBS)

rdinp: $(rdinp_src)
	$(FCOMPL) $(FFLAGC) -o ../bin/rdinp $(rdinp_src) $(LIBS)

ffmod1: $(ffmod1_src)
	$(FCOMPL) $(FFLAGC) -o ../bin/ffmod1 $(ffmod1_src) $(LIBS)

ldos:   $(ldos_src)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ldos $(ldos_src) $(LIBS)

ffmod2: $(ffmod2_src)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod2 $(ffmod2_src) $(LIBS)

ffmod3: $(ffmod3_src)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod3 $(ffmod3_src) $(LIBS)

ffmod4: $(ffmod4_src)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod4 $(ffmod4_src) $(LIBS)

ffmod5: $(ffmod5_src)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod5 $(ffmod5_src) $(LIBS)

ffmod6: $(ffmod6_src)
	$(FCOMPL) $(FFLAGC)  -o ../bin/ffmod6 $(ffmod6_src) $(LIBS)

clean:
	@rm -f *.o ../bin/*

