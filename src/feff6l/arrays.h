c     Notes:
c        nat	number of atoms in problem
c        nph	number of unique potentials
c        nfr	number of unique free atoms
c        ihole	hole code of absorbing atom
c        iph=0 for central atom
c        ifr=0 for central atom

c  Specific atom input data
c     given specific atom, which unique pot?
      dimension iphat(natx)

c     cartesian coords of specific atom
      dimension rat(3,natx)

c Unique potential input data
c     given unique pot, which atom is model?
c     (0 if none specified for this unique pot)
      dimension iatph(0:nphx)

c     given unique pot, which free atom?
      dimension ifrph(0:nphx)

c     given unique pot, how many atoms are there
c     of this type? (used for interstitial calc)
      dimension xnatph(0:nphx)

c     label for user convienence
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

c Free atom data
c     ionicity, input  ,  atomic number, input
      dimension ion(0:nfrx),  iz(0:nfrx)

c ATOM output
c     Note that ATOM output is dimensioned 251, all other r grid
c     data is set to nrptx, currently 250

c     density*4*pi, coulomb potential
      dimension rho(251,0:nfrx), vcoul(251,0:nfrx)

c     Overlap calculation results
c     overlapped density*4*pi, 	overlapped coul pot
      dimension edens(nrptx,0:nphx), vclap(nrptx,0:nphx)
c     overlapped total potential
      dimension vtot(nrptx,0:nphx)

c     Muffin tin calculation results
c     r mesh index just inside rmt, r mesh index just inside rnorman
      dimension imt(0:nphx),  inrm(0:nphx)
c     muffin tin radius, norman radius
      dimension rmt(0:nphx), rnrm(0:nphx)

c     PHASE output
c     interstitial energy ref, phase shifts
      complex*16 eref(nex), ph(nex,ltot+1,0:nphx)
c     number of ang mom levels
      dimension lmax(0:nphx)
