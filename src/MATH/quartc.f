ccccccccccccccccccccccccccccccccccccccccccccccccccc
c      PROGRAM Quartic
      SUBROUTINE Quartic(q)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Finds the roots of the quartic polynomial   c
c             a*q**4 + b*q**2 + c*q + d             c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE

c     Variables
      COMPLEX*16 a, b, c, d, A1, B1, D1, D2, P, QPlus, QMin, q(4), F, G

c     PARAMETERS
      DOUBLE PRECISION One, Two, Four, Eight, Twelve, Twnt7, Sevnt2,
     $     Tto1O3, Tto2O3, Root6, amp
c     Integer-double Parameters
      PARAMETER(One = 1.d0, Two = 2.d0, Four = 4.d0, Eight = 8.d0)
      PARAMETER(Twelve = 12.d0, Twnt7 = 27.d0, Sevnt2 = 72.d0)
c     Irrational Parameters - Tto1O3 = 2**(1/3)
c                           - Tto1O3 = 2**(2/3)
c                           - Root6  = Sqrt(6)
      PARAMETER(Tto1O3 = 1.259921049894873, Tto2O3 = 1.587401051968199)
      PARAMETER(Root6 = 2.449489742783178)

      a = q(1)
      b = q(2)
      c = q(3)
      d = q(4)

c      a = 1.d0
c      b = 2.d0
c      c = 3.d0
c      d = 4.d0

      F  = b**2+ Twelve*a*d
      G  = Two*b**3 + Twnt7*a*c**2 - Sevnt2*a*b*d

      A1 = (G + SQRT(-Four*F**3 + G**2))**(1.d0/3.d0)
      B1 = Two*Tto1O3*F

      P  = SQRT( (-Four*b + B1/A1 + Tto2O3*A1)/a )

      D1 = Eight*b + B1/A1 + Tto2O3*A1
      D2 = Twelve*Root6*c/P

      QPlus = SQRT( -(D1 + D2)/a )
      QMin  = SQRT( -(D1 - D2)/a )

      amp = 1.d0/(Two*Root6)

      q(1) = amp*(P - QPlus)
      q(2) = amp*(P + QPlus)
      q(3) = -amp*(P + QMin)
      q(4) = amp*(-P + QMin)

c      PRINT*, q(1)
c      PRINT*, q(2)
c      PRINT*, q(3)
c      PRINT*, q(4)

      RETURN
      END
