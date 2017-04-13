      subroutine fermi (rhoint, vint, xmu, rs, xf)

      implicit double precision (a-h, o-z)

      include 'const.h'

c     calculate fermi level of the system (mu) according to formula
c     mu=vcoulomb(interstitial)+vxc(interstitial)+kf(interstitial)^2
c     formula  2.13 in lee and beni, phys. rev. b15,2862(1977)

c     note that vint includes both coulomb and ground state
c     exchange-correlation potentials

c     den is the interstitial density
c     rs is the density parameter
c     xf is the interstital fermi momentum
c     xmu is the fermi level in rydbergs

      den = rhoint / (4*pi)
      rs = (3 / (4*pi*den)) ** third
      xf = fa / rs
      xmu = vint + xf**2

      return
      end
