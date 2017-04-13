c-*-fortran-*-
      double precision pi, one, zero, third, raddeg, fa
      double precision bohr, ryd, hart, alpinv, alphfs
      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1.d0, zero = 00.d0)
      parameter (third = one/3.d0)
      parameter (raddeg = 180.d0 / pi)
      complex*16 coni
      parameter (coni = (0.d0, 1.d0))
c     kf = fa/rs with fa = (9*pi/4)**third
c     see Ashcroft and Mermin, p 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.52917721067d0)
      parameter (ryd  = 13.60569301d0)
      parameter (hart = 2.0d0 * ryd)
c     fine structure alpha
      parameter (alpinv = 137.035999139d0)
      parameter (alphfs = 1.d0 / alpinv)
