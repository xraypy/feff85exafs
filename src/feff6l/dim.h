c{vers.h -*-fortran-*-
       integer nphx, npotx, nfrx, novrx, natx
       integer ltot, nrptx, nex, lamtot, mtot, ntot
       integer legtot, npatx
       parameter (nphx  = 7)     !max number of unique potentials (potph)
       parameter (npotx = nphx)  !max number of unique potentials (genfmt, paths)
       parameter (nfrx  = nphx)  !max number of free atom types
       parameter (novrx =   8)   !max number of overlap shells
       parameter (natx  = 500)   !max number of atoms in problem
       parameter (ltot  = 24)    !max number of ang mom (arrays 1:ltot+1)
       parameter (nrptx = 250)   !Loucks r grid used through overlap
       parameter (nex   = 100)   !Number of energy points genfmt, etc.
       
       parameter (lamtot=  15)   !Max number of distinct lambda's for genfmt
                                 !15 handles iord 2 and exact ss
       parameter (mtot=4, ntot=2) !vary mmax and nmax independently
       parameter (legtot=9)     !matches path finder, used in GENFMT
       parameter (npatx = 8)    !max number of path atoms, used in path
                                !finder, NOT in genfmt
c}
