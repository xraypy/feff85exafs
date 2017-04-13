      subroutine phash (npat, ipat, rx, ry, rz, dhash)
c     hashes a path into double precision real dhash

      include 'dim.h'
      double precision dhash
      dimension rx(npatx), ry(npatx), rz(npatx), ipat(npatx+1)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      double precision xx

      parameter (iscale = 1000)
      parameter (factor = 16.12345678)

c     Hashing scheme: Assume about 15 significant digits in a double 
c     precision number.  This is 53 bit mantissa and 11 bits for sign 
c     and exponent, vax g_floating and probably most other machines.
c     With max of 9 legs, 47**9 = 1.12e15, so with a number less than 
c     47, we can use all these digits, scaling each leg's data by 
c     47**(j-1).  Actually, since our numbers can go up to about 10,000,
c     we should keep total number < 1.0e11, 17**9 = 1.18e11, which means
c     a factor a bit less than 17.  Choose 16.12345678, a non-integer,
c     to help avoid hash collisions.

c     iscale and 'int' below are to strip off trailing digits, which
c     may contain roundoff errors

      dhash = 0
      do 210  j = 1, npat
         xx = factor**(j-1)
         dhash = dhash + xx * (nint(rx(j)*iscale) +
     1               nint(ry(j)*iscale)*0.894375 +
     2               nint(rz(j)*iscale)*0.573498)
  210 continue
      do 220  j = 1, npat
         xx = factor**(j-1)
         dhash = dhash + xx * ipot(ipat(j))
  220 continue
      dhash = dhash + npat * 40 000 000

      return
      end
