      subroutine genfmt_prep(phpad, ispin,
c     arguments for rdxsph
     &       ne, ne1, ne3, npot, ihole, rnrmav,
     &       xmu, edge, ik0, ixc, rs, vint,
     &       em, eref2, iz, potlbl, ph4, rkk2, lmax, lmaxp1,
c     arguments for setkap
     &       kinit, linit, ilinit,
c     argument for snlm (also a return)
     &       xnlm,
c     things set here
     &       eref, ph, xk, ck, ckmag, xkr,
     &       nsp, ll, npath, ntotal, nused, xportx)

c+---------------------------------------------------------------------
c     "Based on or developed using Distribution: FEFF8.5L
c      Copyright (c) [2013] University of Washington"
c
C  See ../HEADERS/license.h for full llicense information
c+---------------------------------------------------------------------
c  abstract out the initialization parts of genfmt.  this can then be
c  dropped into genfmt for normal use or be used as part of a single-path
c  library
c+----------------------------------------------------------------------
c     phpad  - specify path to phase.pad     (character*256)
c  Energy grid information
c     em     - complex energy grid
c     eref   - V_int + i*gamach/2 + self-energy correction
c     ne     - total number of points in complex energy grid
c     ne1    - number of points on main horizontal axis
c     ne2    - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3    - number of points on auxilary horizontal axis (need for f')
c     xmu    - Fermi energy
c     edge   - x-ray frequency for final state at Fermi level
c     ik0    - grid point index at Fermi level
c  Potential type information
c     npot   - number of potential types
c     iz     - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax   - max orb momentum for each potential type
c     ihole  - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph     - complex scattering phase shifts
c     rkk    - complex multipole matrix elements
c+----------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      complex*16 ph(nex,-ltot:ltot,0:nphx), eref(nex), em(nex)
      complex*16 ph4(nex,-ltot:ltot, nspx, 0:nphx), eref2(nex,nspx)
      complex*16 ck(nex)
      complex*16 rkk2(nex,8,nspx)
      double precision xk(nex), ckmag(nex), xkr(nex)
      double precision rnrmav, xmu, edge
      dimension xnlm(ltot+1,mtot+1)
      integer lmax(nex,0:nphx), iz(0:nphx)
      integer ne, ne1, ne3, npot, ihole, ik0, kinit, linit, ll
      integer nsp, ispin
      character*6  potlbl(0:nphx)

      character*2 atsym
      external atsym

      character*256 phpad

c     Read phase calculation input
c      print *, 'calling rdxsph'
      call rdxsph (phpad,
     1     ne, ne1, ne3, npot, ihole, rnrmav, xmu, edge, ik0,
     2     ixc, rs, vint,
     3     em, eref2, iz, potlbl, ph4, rkk2, lmax, lmaxp1)
      call setkap (ihole, kinit, linit)
      ilinit = linit + 1

c
c     need to sum over spin-up and -down for |ispin|=1 (fix later)
      nsp = 1
      if (ispin.eq.1) nsp = nspx

      if (nsp.eq.1) then
c        for ispin=2 the variables already written into is=1 positions
         is = 1
         do 10 ie = 1, ne
            eref(ie) = eref2(ie,is)
 10      continue
         do 20 iph = 0, npot
            do 22 ie = 1, ne
               do 24 il = -lmax(ie, iph), lmax(ie, iph)
                  ph(ie,il, iph) = ph4(ie, il, is, iph)
 24            continue
 22         continue
 20      continue
      else
c        average over two spin direction
         do 12 ie = 1, ne
            eref(ie) = (eref2(ie,1) + eref2(ie,nsp)) /2
 12      continue
c        !KJ  12     eref(ie) = (eref2(ie,1) + eref2(ie,2)) /2
         do 30 iph = 0, npot
            do 32 ie = 1, ne
               do 34 il = -lmax(ie, iph), lmax(ie, iph)
                  ph(ie,il, iph) = (ph4(ie, il, 1,iph) + 
     &                   ph4(ie, il, nsp,iph)) /2
 34            continue
 32         continue
 30      continue
c        !KJ  22     ph(ie,il, iph) =(ph4(ie, il, 1,iph) + ph4(ie, il, 2,iph)) /2
      endif

c     Set nlm factors in common /nlm/ for use later
      call snlm (ltot+1, mtot+1, xnlm)

      do 380 i = 0, npot
         if (potlbl(i).eq.' ') potlbl(i)  = atsym(iz(i))
         if (potlbl(i).eq.' ') potlbl(i)  = 'null'
 380  continue 

c     Make xk and ck array for later use
      do 850  ie = 1, ne
c        real momentum (k)
         xk(ie) = getxk (dble(em(ie)) - edge)
c        complex momentum (p)
         ck(ie) = sqrt (2*(em(ie) - eref(ie)))
         ckmag(ie) = abs(ck(ie))
         xkr(ie) = real(xk(ie))
 850  continue

c     Central atom phase shifts
      ll = linit+1
      if (kinit.lt.0) ll = -ll

      npath  = 0
      ntotal = 0
      nused  = 0
      xportx = -1


      return
      end
