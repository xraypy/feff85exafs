c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
