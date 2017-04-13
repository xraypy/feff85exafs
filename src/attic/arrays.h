c  Notes:                                        -*-fortran-*-
c        nat	number of atoms in problem
c        nph	number of unique potentials
c        ihole	hole code of absorbing atom
c        iph=0 for central atom

c  Specific atom input data
c     given specific atom, which unique pot?
      dimension iphat(natx)
c     cartesian coords of specific atom
      dimension rat(3,natx)

c  Unique potential input data
c     iatph(0:nphx)  - given unique pot, which atom is model?
c     (0 if none specified for this unique pot)
      dimension iatph(0:nphx)
c     xnatph(0:nphx) - given unique pot, how many atoms are there 
c     of this type? (used for interstitial calc)
      dimension xnatph(0:nphx)
c     labels for user convienence
      character*6 potlbl(0:nphx)
c     overlap factor for rmt calculation
      dimension folp(0:nphx)
c     number of overlap shells for unique pot
      dimension novr(0:nphx)
c     unique pot for this overlap shell
      dimension iphovr(novrx,0:nphx)
c     number of atoms in overlap shell
      dimension nnovr(novrx,0:nphx)
c     r for overlap shell
      dimension rovr(novrx,0:nphx)

c  Free atom data
c     ionicity, input
      dimension xion(0:nphx)
c     atomic number, input
      dimension iz(0:nphx)

c  ATOM output
c  Note that ATOM output is dimensioned 251, all other r grid
c  data is set to nrptx, currently 250
c     density*4*pi
      dimension rho(251,0:nphx)
c     coulomb potential
      dimension vcoul(251,0:nphx)

c  Overlap calculation results
c     overlapped density*4*pi
      dimension edens(251,0:nphx)
c     overlapped coul pot
      dimension vclap(251,0:nphx)
c     overlapped total potential
      dimension vtot (251,0:nphx)

c  Muffin tin calculation results
c     r mesh index just inside rmt
      dimension imt(0:nphx)
c     r mesh index just inside rnorman
      dimension inrm(0:nphx)
c     muffin tin radius
      dimension rmt(0:nphx)
c     norman radius
      dimension rnrm(0:nphx)

c  PHASE output
c     interstitial energy ref
      complex*16 eref(nex)
c     phase shifts
      complex*16 ph(nex,ltot+1,0:nphx)
c     number of ang mom levels
      dimension lmax(0:nphx)
