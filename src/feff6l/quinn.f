      subroutine quinn (x, rs, wp, ef, ei)
      implicit double precision (a-h, o-z)

c     input  x, rs, wp, ef
c     output ei

c***********************************************************************
c
c     quinn: calculates low energy gamma (approx. proportional to e**2)
c             formula taken from john j. quinn, phys. rev. 126,
c             1453 (1962); equation (7).
c             a cut-off is set up at quinn's cutoff + ef = ekc; it is a
c             rounded inverted step function (a fermi function)
c             theta = 1/( 1 + exp((e-ekc)/gam)) )
c             where the rounding factor gam is set to be about 0.3 ekc.
c     modified by j. rehr (oct 1991) based on coding of r. albers
c     subroutines quinn.f and quinnc.f
c
c     variables:
c        x  = p/pf
c        rs = ws density parameter
c        ei = imaginary self energy
c        pfqryd = quinn's prefactor in atomic-rydberg units
c        wkc = quinn's plasmon threshold
c
c***********************************************************************

      include 'const.h'

      parameter (alphaq = 1/ fa)

c     calculate quinn prefactor in atomin Hartree units
      pisqrt = sqrt(pi)
      pfq = pisqrt / (32 * (alphaq*rs)**1.5)
      temp1 = atan (sqrt (pi / (alphaq*rs)))
      temp2 = sqrt(alphaq*rs/pi) / (1 + alphaq*rs/pi)
      pfq = pfq * (temp1 + temp2)

c     calculate quinn cutoff
c     wkc = quinn's plasmon threshold
c     wkc is cut-off of quinn, pr126, 1453, 1962, eq. (11)
c     in formulae below wp=omegap/ef
      wkc = (sqrt(1+wp) - 1)**2
      wkc = (1 + (6./5.) * wkc / wp**2) * wp * ef

c     we add fermi energy to get correct energy for
c     plasma excitations to turn on
      ekc = wkc + ef

c     calculate gamma
c     gamryd = 2 * (pfqryd/x) * (x**2-1)**2
      gam = (pfq/x) * (x**2-1)**2

c     put in fermi function cutoff
      eabs = ef * x**2
      arg = (eabs-ekc) / (0.3*ekc)
      f = 0
      if (arg .lt. 80)  f = 1 / (1 + exp(arg))

      ei = -gam * f / 2

      return
      end
