      subroutine fdtarr(ne, reff, lzero, achi, phchi, caps, xk, ck,
     &       col1, col2, col3, col4, col5, col6, col7)

c+---------------------------------------------------------------------
c     "Based on or developed using Distribution: FEFF8.5L
c      Copyright (c) [2013] University of Washington"
c
C  See ../HEADERS/license.h for full llicense information
c+---------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      parameter (eps = 1.0d-16)
      real xk(nex), reff
      complex caps(nex), ck(nex)
      real achi(nex), phchi(nex)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
      complex*16 cchi, cfms
      integer ne, lzero, ie

c     Make the feff.dat stuff and write it to feff.dat
c     Also write out for inspection to fort.66
c     note that dimag takes complex*16 argument, aimag takes
c     single precision complex argument.  Stuff from feff.pad
c     is single precision, cchi is complex*16

      do 450  ie = 1, ne
c        Consider chi in the standard XAFS form.  Use R = rtot/2.
         cchi = achi(ie) * exp (coni*phchi(ie))
         xlam = 1.0e10
         if (abs(aimag(ck(ie))) .gt. eps) xlam = 1/aimag(ck(ie))
         redfac = exp (-2 * aimag (caps(ie)))
         cdelt = 2*dble(caps(ie))
         cfms = cchi * xk(ie) * reff**2 *
     1          exp(2*reff/xlam) / redfac
         if (abs(cchi) .lt. eps)  then
            phff = 0
         else
            phff = atan2 (dimag(cchi), dble(cchi))
         endif
c        remove 2 pi jumps in phases
         if (ie .gt. 1)  then
            call pijump (phff, phffo)
            call pijump (cdelt, cdelto)
         endif
         phffo = phff
         cdelto = cdelt


c        columns
c        . 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c        . 2 central atom phase shift (real part),
c        . 3 magnitude of feff,
c        . 4 phase of feff,
c        . 5 absorbing atom reduction factor,
c        . 6 mean free path = 1/(Im (p))
c        . 7 real part of local momentum p

         col1(ie) = xk(ie)/bohr
         col2(ie) = cdelt + lzero*pi
         col3(ie) = abs(cfms) * bohr
         col4(ie) = phff - cdelt - lzero*pi
         col5(ie) = redfac
         col6(ie) = xlam * bohr
         col7(ie) = dble(ck(ie))/bohr
 450  continue
      return
      end
