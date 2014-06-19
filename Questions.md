
# Compilation warnings (gfortran on linux)

 1. In `Math/`: 

           gfortran -o polint.o -c polint.f
           polint.f:30.72:

                       if (den.eq.0) pause 'failure in polint'                     
                                                                        1
           Warning: Deleted feature: PAUSE statement at (1)

 2. In `PAR/`: missing `mpif.h`.  Symlink to same from fdmnes allowed
    compilation to finish.  That seems unlikely....


# General confusion

 1. What is `libtd.a`?

# For Josh

 1. In file `EXCH/mpse.f`, `mpse.f` does not compile with rank
    mismatches at lines 37 and 57.  Variable `edens` is defined as
    rank 2 and called as rank 1.  How is this to be fixed?
	
 
