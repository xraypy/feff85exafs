c     Notes:
c        nat	number of atoms in problem
c        nph	number of unique potentials
c        nfr	number of unique free atoms
c        ihole	hole code of absorbing atom
c        iph=0 for central atom
c        ifr=0 for central atom

c     Specific atom input data
      dimension iphat(natx)	!given specific atom, which unique pot?
      dimension rat(3,natx)	!cartesian coords of specific atom

c     Unique potential input data
      dimension iatph(0:nphx)	!given unique pot, which atom is model?
				!(0 if none specified for this unique pot)
      dimension ifrph(0:nphx)	!given unique pot, which free atom?
      dimension xnatph(0:nphx)	!given unique pot, how many atoms are there 
				!of this type? (used for interstitial calc)
      character*6 potlbl(0:nphx)	!label for user convienence

      dimension folp(0:nphx)	!overlap factor for rmt calculation
      dimension novr(0:nphx)	!number of overlap shells for unique pot
      dimension iphovr(novrx,0:nphx)	!unique pot for this overlap shell
      dimension nnovr(novrx,0:nphx)	!number of atoms in overlap shell
      dimension rovr(novrx,0:nphx)	!r for overlap shell

c     Free atom data
      dimension ion(0:nfrx)	!ionicity, input
      dimension iz(0:nfrx)	!atomic number, input

c     ATOM output
c     Note that ATOM output is dimensioned 251, all other r grid
c     data is set to nrptx, currently 250
      dimension rho(251,0:nfrx)		!density*4*pi
      dimension vcoul(251,0:nfrx)	!coulomb potential

c     Overlap calculation results
      dimension edens(nrptx,0:nphx)	!overlapped density*4*pi
      dimension vclap(nrptx,0:nphx) 	!overlapped coul pot
      dimension vtot (nrptx,0:nphx)	!overlapped total potential

c     Muffin tin calculation results
      dimension imt(0:nphx)	!r mesh index just inside rmt
      dimension inrm(0:nphx)	!r mesh index just inside rnorman
      dimension rmt(0:nphx)	!muffin tin radius
      dimension rnrm(0:nphx)	!norman radius

c     PHASE output
      complex*16 eref(nex)		!interstitial energy ref
      complex*16 ph(nex,ltot+1,0:nphx)	!phase shifts
      dimension lmax(0:nphx)		!number of ang mom levels
