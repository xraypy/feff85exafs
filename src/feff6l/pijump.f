      subroutine pijump (ph, old)
      implicit double precision (a-h, o-z)

c     removes jumps of 2*pi in phases

c     ph = current value of phase (may be modified on output, but
c          only by multiples of 2*pi)
c     old = previous value of phase

      include 'const.h'
      parameter (twopi = 2 * pi)
      dimension xph(3)

      xph(1) = ph - old
      jump   =  (abs(xph(1))+ pi) / twopi
      xph(2) = xph(1) - jump*twopi
      xph(3) = xph(1) + jump*twopi


      xphmin = min (abs(xph(1)), abs(xph(2)), abs(xph(3)))
      isave = 1
      do 10  i = 1, 3
         if (abs (xphmin - abs(xph(i))) .le. 0.01d0)  isave = i
   10 continue

      ph = old + xph(isave)

      return
      end
