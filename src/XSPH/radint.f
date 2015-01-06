      subroutine radint (ifl, mult, bf, kinit, dgc0, dpc0, ikap, p, q,
     1   pn, qn, ri, dx, ilast, iold, xrc, xnc, xrcold, xncold, xirf)
c     performs radial integration for multipole matrix element
c     or central atom absorption depending on flag 'ifl'.
      implicit double precision (a-h, o-z)

c     INPUT
c     ifl - number corresponds to the calling order in xsect.f
c         - 1 - calculate matrix element (rkk)
c         - -1 - calculate matrix element (rkk) in nonrelativistic limit
c         - 2 - calculate cross section (xsec)
c             cross term needed for spin-dependent potential only
c         - 3 - cross term (xsec) with irregular part for current kappa
c         - 4 - cross term (xsec) with regular part for current kappa
c     mult - specifies multipole transition
c     bf - bessel functions for x-ray k-vector for l=0,1,2
c     kinit - initial kappa
c     dgc0,dpc0 - large (small) dirac components for initial orbital 
c     ikap  = final state kappa
c     p,q   Dirac components for regular (R) final state solution
c     pn,qn  Dirac components for irregular(N) final state solution
c     ri,dx - radial grid
c     ilast - last integration point
c     iold  - 0 - do nothing to xrcold, xncold (ic3=0 case)
c             1 - store intermediate results in xrcold, xncold (ic3=1)
c             2 - use intermediate results in xrcold, xncold (ic3=1)
c
c     OUTPUT
c     xrcold,xncold - coupling to regular (R) and irregular(N) solutions
c                     both output and input
c     xirf  - value of the radial integral

      include '../HEADERS/dim.h'
      dimension ri(nrptx), dgc0(nrptx), dpc0(nrptx)
      dimension bf(0:2, nrptx)
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx)
c     storage for calculation of cross term (SPIN 1 only)
      complex*16 xrcold(nrptx) , xncold(nrptx)
      complex*16  xirf, temp

c     local staff
      complex*16  xm(4)
      complex*16 xrc(nrptx), xnc(nrptx)
      complex*16 coni
      parameter (coni = (0.d0, 1.d0))

      temp = (1.,0.)

      linit = kinit
      if (kinit.lt.0) linit = - kinit - 1
      lfin = ikap
      if (ikap.lt.0) lfin = - ikap - 1
c     set multipliers  from Grant,Advan.Phys.,v.19,747(1970) eq. 6.30,
c     using Messiah's "Q.M." appendices to calculate 9j,3j symbols
      if (ifl.lt.0) then
        ji2 = 2*abs(kinit)-1
        jf2 = 2*abs(ikap)-1
        if (mult.eq.0 .or. mult.eq.2) then
           ll = 1
           if (mult.eq.2) ll = 2
           ll2 = 2*ll
           temp = sqrt(dble((ji2+1)*(jf2+1))) *cwig3j(jf2,ll2,ji2,1,0,2)
c          sign of temp is (-)**(j+1/2): compare eq. 6.2 and 6.30 
c          of Grant, Adv. Phys. 19, 747 (1970).
           temp = temp * (-1)**(abs(ikap))
           ls = ll-1
           xm(1) = temp * (ll2+1) *coni**ls *(2*ls+1) *
     1     cwig3j(ls,1,ll,0,0,1) * cwig3j(ls,1,ll,0,1,1)
           ls = ll+1
           xm(3) = 0
c          xm(3) = temp * (ll2+1) *coni**ls *(2*ls+1) *
c    1     cwig3j(ls,1,ll,0,0,1) * cwig3j(ls,1,ll,0,1,1)
        else
c          if (mult.eq.1) then
           stop 'not set up for M1 transition in nonrelativistic limit'
        endif
      elseif (mult.eq.0) then
        call xmult( ikap, kinit, 0, 1, xm(1), xm(2))
        call xmult( ikap, kinit, 2, 1, xm(3), xm(4))
      else
        xm(3) = 0
        xm(4) = 0
        if (mult.eq.2) then
          call xmult( ikap, kinit, 1, 2, xm(1), xm(2))
        else
c         mult=1 - M1 transition
          call xmult( ikap, kinit, 1, 1, xm(1), xm(2))
        endif
      endif

c     radial integrals depending on case
      ia = abs(ifl)
      is = ifl /ia
      if (ia.eq.1) then
c       single radial integral for rkk - reduced matrix elements
c       xirf = <f |p| i> relativistic version of dipole m.e.
        do 10  i = 1, ilast
          xnc(i) = 0.0d0
          if (is.gt.0) then
           call xrci(mult,xm,dgc0(i),dpc0(i),p(i),q(i),bf(0,i),xrc(i))
          else
c          nonrelativistic case 
           if (mult.eq.0) then
             temp = xm(1)*bf(0,i)+ xm(3)*bf(2,i)
           elseif (mult.eq.2) then
             temp = xm(1)*bf(1,i)
           endif
           temp = temp *coni
           xrc(i) = ri(i) * (dgc0(i)*p(i) + dpc0(i)*q(i)) *temp
c          xrc(i) = ri(i) * (dgc0(i)*p(i) ) *temp
          endif

c         store xrc if needed
          if (iold.eq.1) xrcold(i) = xrc(i)
  10    continue
        xirf=lfin+linit+2
        if (mult.gt.0) xirf = xirf + 1
        call csomm (ri, xrc, xnc, dx, xirf, 0, ilast)
      else
c       need to perform double radial integral in all cases below
        if (ia.eq.2) then
c         combine regular(kdif) and irregular(kdif) solution into
c         the central atom absorption coefficient xsec (mu = dimag(xsec))
c         thus for real energy dimag(xsec)=xsnorm
          do 20  i = 1, ilast
           if (is.gt.0) then
           call xrci(mult,xm,dgc0(i),dpc0(i),pn(i),qn(i),bf(0,i),xnc(i))
           call xrci(mult,xm,dgc0(i),dpc0(i),p(i),q(i),bf(0,i),xrc(i))
           else
c            nonrelativistic case 
             if (mult.eq.0) then
               temp = xm(1)*bf(0,i)+ xm(3)*bf(2,i)
             elseif (mult.eq.2) then
               temp = xm(1)*bf(1,i)
             endif
             temp = temp*coni
             xrc(i) = ri(i) * (dgc0(i)*p(i) + dpc0(i)*q(i)) *temp
             xnc(i) = ri(i) * (dgc0(i)*pn(i) + dpc0(i)*qn(i)) *temp
c            xrc(i) = ri(i) * (dgc0(i)*p(i) ) *temp
c            xnc(i) = ri(i) * (dgc0(i)*pn(i) ) *temp
           endif
c           store irregular contribution for later use
            if (iold.eq.1) xncold(i) = xnc(i)
  20      continue
        elseif (ifl.eq.3 .and. iold.eq.2) then
c         combine regular(k1) and irregular (kdif) solutions into the
c         central atom absorption coefficient xsec (mu = dimag(xsec))
c         nonzero only for |ispin=1| and same angular momenta in k1,kdif
          do 30  i = 1, ilast
            xrc(i)= xrcold(i)
           call xrci(mult,xm,dgc0(i),dpc0(i),pn(i),qn(i),bf(0,i),xnc(i))
  30      continue
        elseif(ifl.eq.4 .and. iold.eq.2) then
c         combine regular(kdif) and irregular (k1) solutions into the
c         central atom absorption coefficient xsec (mu = dimag(xsec))
c         nonzero only for |ispin=1| and same angular momenta in k1,kdif
          do 40  i = 1, ilast
            call xrci( mult,xm,dgc0(i),dpc0(i),p(i),q(i),bf(0,i),xrc(i))
            xnc(i) = xncold(i)
  40      continue
        endif

c       same staff for all double integrals
        if ((iold.eq.0.and.ia.eq.2) .or. (ifl.gt.2.and.iold.eq.2)) then
c          do radial integration for r'>r first
c          power of xrc near zero
           lpwr = lfin + linit +2
c          factor 2 since integral(r<r')=integral(r>r')
           xirf = 2 * xrc(1) * ri(1) /(lpwr+1)
           xnc(1) = xnc(1) * xirf
           do 70 i = 2, ilast
             xirf = xirf + (xrc(i-1)+xrc(i)) * (ri(i)-ri(i-1))
             xnc(i) = xnc(i) * xirf
  70       continue
           do 80 i = 1,ilast
  80       xrc(i) = 0
           xirf = lpwr+1+linit+1-lfin
c          ready for second integral over r from 0 to \infty
           call csomm (ri, xrc, xnc, dx, xirf, 0, ilast)
        endif
      endif

      return
      end

      subroutine xrci( mult, xm, dgc0, dpc0, p, q, bf, value)
c     r-dependent multipole matrix element (before r-integration)
      implicit double precision (a-h, o-z)
      complex*16 xm(4), p, q, value
      dimension bf(0:2)

      if (mult.eq.0) then
c       el. dipole transition with both j0 and j2 contributions
        value = dgc0*q* (xm(2)*bf(0) + xm(4)*bf(2)) +
     1         dpc0*p* (xm(1)*bf(0) + xm(3)*bf(2))
       else
         value = (xm(2)*dgc0*q+xm(1)*dpc0*p) * bf(1)
       endif

      return
      end
