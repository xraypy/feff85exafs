      subroutine polint( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
c
c Use of a pause is a poor idea, not just because it is deprecated, but
c because it's normal use is probably not what is desired here.  The most
c common intepretation of "pause" is to wait for a carriage return as
c acknowledgment.  In this case, it will result in a divide-by-zero segfault
c with an unhelpful error message.
c I have used a functionally equivalent, write/read replacement, but it's
c still the wrong thing to do. -BR
c
c            if (den.eq.0) pause 'failure in polint'
            if (den.eq.0) then
               write( *, * ) 'failure in polint'
               read( *, * ) 
            end if
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
