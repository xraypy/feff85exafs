      subroutine fixvar (rmt, edens, vtot,
     1                   vint, rhoint, nr, dx, x0, ri,
     2                   vtotph, rhoph)

      implicit double precision (a-h, o-z)

      include 'dim.h'
      include 'const.h'

      dimension edens(nrptx), vtot (nrptx)
      dimension vtotph(nr), rhoph(nr)
      dimension ri(nr)

c     PHASE needs
c     vtot = total potential including gs xcorr, no r**2
c     edens = rho, charge density, no factor of 4*pi, no r**2
c     From overlapping, vtot = potential only, ok as is
c                       edens = density*4*pi, so fix this here.

c     If new grid is different from old one, be sure to interpolate
c     somehow...

c     Only values inside the muffin tin are used, except that XCPOT
c     (in PHASE) uses values at imt+1 and requires these to be the
c     interstitial values.  So set the last part of the arrays to
c     interstitial values...

      imt = ii(rmt)

      do 190  i = 1, imt
         vtotph(i) = vtot(i)
         rhoph(i) = edens(i)/(4*pi)
  190 continue
      do 200  i = imt+1, nrptx
         vtotph(i) = vint
         rhoph(i) = rhoint/(4*pi)
  200 continue

      return
      end
