
# Compilation warnings (gfortran on linux)

 1. In `Math/`: 

           gfortran -o polint.o -c polint.f
           polint.f:30.72:

                       if (den.eq.0) pause 'failure in polint'                     
                                                                        1
           Warning: Deleted feature: PAUSE statement at (1)

 2. In `PAR/`: missing `mpif.h`.  Symlink to same from fdmnes allowed
    compilation to finish.  That seems unlikely....


 3. It seems as though the NOHOLE option won't work.  rdinp allows it,
    but XSPH/xsph.f wants, when nohole.eq.2, some data files (see
    lines 209 to 224) that are nowhere written in f85e

# General confusion

 1. What is `libtd.a`?

 2. EXCH/xcpot.f lines 379, 386 write to fort.39

 3. EXCH/xcpot.f line 390 write to fort.38

# For Josh

 1. In file `EXCH/mpse.f`, `mpse.f` does not compile with rank
    mismatches at lines 37 and 57.  Variable `edens` is defined as
    rank 2 and called as rank 1.  How is this to be fixed?
	
 
