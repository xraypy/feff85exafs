      subroutine sigte3 (iz1,iz2, sig2, alphat, thetad, reff, sig1,sig3)
c     single scattering only.

c     input: sig2
c     iz1, iz2 are iz at central atom and neighbor
c     alphat coeef of thermal expansion at high T
c     reff

c     output: sig1 sig3
      implicit double precision (a-h, o-z)
      real reff

c     con=hbar**2/(kB*amu)*10**20   in ang**2 units
c     hbar, amu, kb updated 2017
      parameter (hbar = 1.054571800d-34)
      parameter (amu = 1.660539040d-27)
      parameter (xkb = 1.38064852d-23)
      parameter (con = 48.50875019927435d0)

      ami=atwtd(iz1)*amu
      amj=atwtd(iz2)*amu

c     reduced mass
      xmu = 1 / (1/ami + 1/amj)
c     Einstein frequency
      omega = (2 * xkb * thetad) / (3 * hbar)
      xks = xmu * omega**2
      xk3 = xks**2 * reff * alphat / (3 * xkb)
      sig02 = hbar * omega / xks
      sig1 = -3 * (xk3 / xks) * sig2
      sig3 = 2 - (4.0/3.0) * (sig02 / sig2)**2
      sig3 = sig3 * (sig1 * sig2)

      return
      end
