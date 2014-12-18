      subroutine import(ne1, nsp, ik0, reff, deg, ckmag, em, eref2,
     &        cchi, xportx, crit)

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      complex*16  ckp
      dimension   ffmag(nex), ckmag(nex)
      complex*16  eref(nex), em(nex), eref2(nex,nspx), cchi(nex)
      complex*16  ck(nex)

      external trap

c     Make importance factor, deg*(integral (|chi|*d|p|))
c     make ffmag (|chi|)
c     xport   importance factor
      do 10  ie = 1, ne1
c        this gets ck(ie) set correctly; see lines 257 and 272 in genfmt.f
         ck(ie) = sqrt (2* (em(ie) - eref2(ie,1)))
         if (nsp.eq.2) then
            eref(ie) = (eref2(ie,1) + eref2(ie,nsp)) /2
c           !KJ eref(ie) = (eref2(ie,1) + eref2(ie,2)) /2
            ck(ie) = sqrt (2* (em(ie) - eref(ie)))
         endif
         ckp = ck(ie)
         xlam0 = dimag(ck(ie)) - dimag(ckp)
         ffmag(ie) = abs( cchi(ie) * exp(2*reff*xlam0) )
 10   continue

c     integrate from edge (ik0) to ne
      nemax = ne1 - ik0 + 1
      call trap (ckmag(ik0), ffmag(ik0), nemax, xport)
      xport = abs(deg*xport)
      if (xportx.le.0)  xportx = xport
      crit = 100 * xport / xportx
c     use line below to disable importance factor (e.g. for dichroism)
c     crit = crit0+1

      return
      end
