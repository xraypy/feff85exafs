      subroutine bjnser (x, l, jl, nl, ifl)

c-----------------------------------------------------------------------
c
c     subroutine: bjnser (x,l,jl,nl,ifl)
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c
c     arguments:
c       x = argument of jl and nl
c       l = l value calculated (no offset)
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c       ifl = 0 return both jl and nl
c             1 return jl only
c             2 return nl only
c
c     notes:  jl and nl are calculated by a series
c             expansion according to 10.1.2 and 10.1.3
c             in abramowitz and stegun (ninth printing),
c             page 437
c
c             error msgs written with PRINT statements.
c
c     first coded by r. c. albers on 26 jan 83
c
c     version 2
c
c     last modified: 27 jan 83 by r. c. albers
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      complex*16 x,u,ux,del,pj,pn
      complex*16 jl,nl

      parameter (niter = 20, tol = 1.e-15)

       if (l .lt. 0) then
          call echo('Error in bjnser: l<0')
       elseif (dble(x).lt. 0.) then
          call echo('Error in bjnser: x<0')
       else 
          lp1 = l+1
          u = x**2 / 2

c     make djl = 1 * 3 * 5 * ... * (2*l+1),
c          dnl = 1 * 3 * 5 * ... * (2*l-1)
          djl = 1
          fac = -1
          do 50 il = 1, lp1
             fac = fac + 2
             djl = fac * djl
 50       continue
          dnl = djl / (2*l+1)


          if (ifl .eq. 2)   goto 90
c      make jl
c     pj is term in { } in 10.1.2, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
          pj = 1
          nf = 1
          nfac = 2*l + 3
          den = nfac
          sgn = -1
          ux = u
          do 60 il = 1, niter
             del = sgn*ux / den
             pj = pj + del
             trel = abs (del / pj)
             if (trel .le. tol)  goto 80
             sgn = -sgn
             ux = u*ux
             nf = nf+1
             nfac = nfac+2
             den = nf * nfac * den
 60       continue
          call echo('Error in bjsner: jl does not converge')
 80       jl = pj * (x**l) / djl

 90       if (ifl.eq.1) return
c     make nl
c      pn is term in { } in 10.1.3, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
          pn = 1
          nf = 1
          nfac = 1 - 2*l
          den = nfac
          sgn = -1
          ux = u
          do 100  il = 1, niter
             del = sgn * ux / den
             pn = pn + del
             trel = abs (del / pn)
             if (trel .le. tol) goto 120
             sgn = -sgn
             ux = u*ux
             nf = nf+1
             nfac = nfac+2
             den = nf * nfac * den
 100      continue
          call echo('Error in bjnser: nl does not converge')
 120      continue 
          nl = -pn * dnl / (x**lp1)
       endif

      return
      end
