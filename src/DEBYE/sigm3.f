      subroutine sigm3(sig1, sig2, sig3, tk, alphat, thetae)

c     using correlated Einstein-model with a morse potential
c     Nguyen Van Hung & J.J.Rehr Phys. Rev. B 56 , 43

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'

      real sig02,  sig01, z
c     dimension alphat=[1/anstroems]
      parameter(three=3.d0)
      parameter(four=4.d0 )
      parameter(fourthird= four/three)
      parameter(threequarter= three/four)

      alphat= alphat * bohr
      z     = real(exp(- thetae/tk))
      sig02 = real((1-z)/ (1+z) * sig2)
      sig01 = real(alphat) * sig02 * real(threequarter)
      sig1  = sig01 * sig2 / sig02
      sig3  = (2- fourthird * (sig02/sig2) **2)* sig1 * sig2

      return
      end
