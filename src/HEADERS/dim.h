C -*-fortran-*-
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      integer nclusx
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      integer nclxtd
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      integer nspx
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      integer natx
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      integer nattx
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      integer lx
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      integer nphx
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      integer ltot
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      integer nrptx
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      integer nex
      parameter (nex = 150)
c      Max number of distinct lambda values for genfmt
c      15 handles iord 2 and exact ss
      integer lamtot
      parameter (lamtot=15)
c      vary mmax and nmax independently
      integer mtot, ntot
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      integer npatx
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      integer legtot
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      integer novrx
      parameter (novrx=8)
c      max number of header lines
      integer nheadx
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      integer MxPole
      parameter (MxPole=1000)
