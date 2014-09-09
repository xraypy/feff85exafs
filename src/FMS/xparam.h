c -*- fortran -*-
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
