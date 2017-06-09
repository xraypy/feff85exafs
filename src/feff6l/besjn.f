      subroutine besjn (x, jl, nl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 30 (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this array nl = abramowitz yl.
c       jl and nl must be dimensioned 
c            complex*16 jl(ltot+2), nl(ltot+2), with ltot defined in 
c            dim.h.
c
c     notes:  jl and nl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c
c-----------------------------------------------------------------------

       implicit double precision (a-h, o-z)
       include 'dim.h'
       
       complex*16 x
       complex*16 jl(ltot+2), nl(ltot+2)
       complex*16 cjl(ltot+2), sjl(ltot+2), cnl(ltot+2), snl(ltot+2)
       
       complex*16 xjl,xnl,asx,acx
       complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11
       double precision xcut, xcut1, xcut2
       
       parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)

      if (dble(x) .le. 0)  call fstop(' at besjn: '//
     $      'Re(x) is <=0')

      lmaxp1 = ltot+2

      if (dble(x) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 il = 1,lmaxp1
            l = il-1
            ifl = 0
            call bjnser (x,l,xjl,xnl,ifl)
            jl(il) = xjl
            nl(il) = xnl
   10    continue

      elseif (dble(x) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmaxp1-1,xjl,xnl,1)
         jl(lmaxp1) = xjl

         call bjnser (x,lmaxp1-2,xjl,xnl,1)
         jl(lmaxp1-1) = xjl

         if (dble(x) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi*xi
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            tlxp1 = 2*lp1 -3
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lx = 3,lmaxp1
            lp1 = lmaxp1+1-lx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(lp1) = tlxp3 * jl(lp1+1) / x  -  jl(lp1+2)
   60    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3.d0*xi3 - xi
         sjl(4) = 15.d0*xi4 - 6.d0*xi2
         sjl(5) = 105.d0*xi5 - 45.d0*xi3 + xi
         sjl(6) = 945.d0*xi6 - 420.d0*xi4 + 15.d0*xi2
         sjl(7) = 10395.d0*xi7 - 4725.d0*xi5 + 210.d0*xi3 - xi
         sjl(8) = 135135.d0*xi8 - 62370.d0*xi6 + 3150.d0*xi4 -
     $        28.d0*xi2
         sjl(9) = 2027025.d0*xi9 - 945945.d0*xi7 + 51975.d0*xi5 -
     1        630.d0*xi3 + xi
         sjl(10) = 34459425.d0*xi10 - 16216200.d0*xi8 + 
     1        945945.d0*xi6 - 13860.d0*xi4 + 45.d0*xi2
         sjl(11) = 654729075.d0*xi11 - 310134825.d0*xi9 + 
     1        18918900.d0*xi7 - 315315.d0*xi5 + 1485.d0*xi3 - xi

         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3.d0*xi2
         cjl(4) = -15.d0*xi3 + xi
         cjl(5) = -105.d0*xi4 + 10.d0*xi2
         cjl(6) = -945.d0*xi5 + 105.d0*xi3 - xi
         cjl(7) = -10395.d0*xi6 + 1260.d0*xi4 - 21.d0*xi2
         cjl(8) = -135135.d0*xi7 + 17325.d0*xi5 - 378.d0*xi3 + xi
         cjl(9) = -2027025.d0*xi8 + 270270.d0*xi6 - 6930.d0*xi4 +
     $        36.d0*xi2
         cjl(10) = -34459425.d0*xi9 + 4729725.d0*xi7 - 135135.d0*xi5 +
     1        990.d0*xi3 - xi
         cjl(11) = -654729075.d0*xi10 + 91891800.d0*xi8 - 
     1        2837835.d0*xi6 + 25740.d0*xi4 - 55.d0*xi2
         do 80 ie = 1,11
            snl(ie) = cjl(ie)
            cnl(ie) = -sjl(ie)
   80    continue
         do 90 lp1 = 12,lmaxp1
            tlxp1    = 2*lp1-3
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
            snl(lp1) = tlxp1*xi*snl(lp1-1)-snl(lp1-2)
            cnl(lp1) = tlxp1*xi*cnl(lp1-1)-cnl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         do 110 lp1 = 1,lmaxp1
            jl(lp1) = asx*sjl(lp1)+acx*cjl(lp1)
            nl(lp1) = asx*snl(lp1)+acx*cnl(lp1)
  110    continue
      endif

      return
      end
