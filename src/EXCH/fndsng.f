      SUBROUTINE FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
c     Josh Kas
c     This subroutine finds the singularities in the integrands of eq. 13
c     in
c     Single-particle Spectrum of the Degenerate Electron Gas
c     II. Numerical Results for Electrons Coupled to Plasmons
c     Phys. kondens. Materie, Bd. 6 (1967)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     In practice this routine only solves for the singularities of one
c     of the three integrands, then checks to see that the singularity
c     is within the limits of integration, and throws out singularities
c     that are not.
c     In units of the Fermi energy the equations to solve are:
c     1) +/- k*q**3 + 2*(3*k**2 - E - 2/3)*q**2 +/- 4*k*(k**2 - E)*q +
c                     [(k**2 - E)**2 - Wp**2] = 0 
c     2) q**4 + 4/3*q**2 + Wp**2 - (1 - E)**2 = 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Input:
c     Limit1 - Lower limit of integration
c     Limit2 - Upper limit of integration
c     CPar   - Array of complex parameters passed to function
c              CPar(1) = ck/kFermi
c              CPar(2) = Energy/EFermi + i*Gamma/EFermi
c     DPPar  - Array of double precision parameters passed to function
c              DPPar(1) = Wp/EFermi
c              DPPar(2) = Gamma/EFermi
c              DPPar(3) = Energy/EFermi
c              DPPar(4) = xeg (gap energy)
c     iFcn   - Integer denoting which function is the integrand
c              iFcn = 1: solve eqs 1 and 2 for q
c              iFcn = 2: solve eq 1 for q      
      COMPLEX*16 Limit1, Limit2, CPar(10)
      DOUBLE PRECISION DPPar(10)
      INTEGER iFcn
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     XSing  - Array of singularities
c     NSing  - Number of singularities
      COMPLEX*16 XSing(20)
      INTEGER NSing
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     Coef   - Coefficients of q**n of eq. to solve.
c              eq = Coef(1)*q**n + Coef(2)*q**(n-1)...
c     Sol(4) - Array of solutions to the equation.
c     XSing2 - Temp XSing      
c     Test   - Used to test solution of equation.
c     Zero   - Tolerance for testing solution to eqs.
c     NSol   - Number of solutions to eq.
c     Order  - Used to order singulaties from smallest to largest
      COMPLEX*16 Coef(4), Sol(4), XSing2(4)
      DOUBLE PRECISION Test, Zero
      INTEGER NSol, Order(100)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Loop variables
      INTEGER i1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialization
      NSing = 0
      Zero=1.d-4

c     Solve eq 1 for q with + sign
      Coef(1) = 4.d0*CPar(1)
      Coef(2) = 2.d0*(3.d0*CPar(1)**2 - DPPar(3) - 2.d0/3.d0)
      Coef(3) = 4.d0*CPar(1)*(CPar(1)**2 - DPPar(3))
      Coef(4) = (CPar(1)**2 - DPPar(3))**2 - DPPar(1)**2
      
      CALL CCubic(Coef, Sol, NSol)
      
c     Test solutions. If Sol is a solution and it is real
c     and it lies between Limit1 and Limit2, add it to list
c     of singularities.         
      DO i1 = 1, NSol
         Test = ABS((CPar(1)+Sol(i1))**2 - DPPar(3) +
     &        SQRT(Sol(i1)**4 + 4.d0/3.d0*Sol(i1)**2 + DPPar(1)**2))
         IF(Test.lt.Zero) THEN               
            IF((DBLE(Sol(i1)).ge.DBLE(Limit1)).and.
     &           (DBLE(Sol(i1)).le.DBLE(Limit2)).and.
     &           (ABS(DIMAG(Sol(i1))).le.Zero)) THEN               
               NSing = NSing + 1
               XSing(NSing) = DBLE(Sol(i1))
            END IF
         END IF
      END DO
      
c     Now solve eq. 1 for q with - sign
      Coef(1) = -Coef(1)
      Coef(3) = -Coef(3)
      
      CALL CCubic(Coef, Sol, NSol)
      
c     Test solutions as before.
      DO i1 = 1, NSol
         Test = ABS((CPar(1)-Sol(i1))**2 - DPPar(3) -
     &        SQRT(Sol(i1)**4 + 4.d0/3.d0*Sol(i1)**2 + DPPar(1)**2))
         IF(Test.lt.Zero) THEN
            IF((DBLE(Sol(i1)).ge.DBLE(Limit1)).and.
     &           (DBLE(Sol(i1)).le.DBLE(Limit2)).and.
     &           (ABS(DIMAG(Sol(i1))).le.Zero)) THEN
               NSing = NSing + 1
               XSing(NSing) = DBLE(Sol(i1))
            END IF
         END IF
      END DO

c     If iFcn = 1 (Solving for singularities of r1(q))
      IF(iFcn.eq.1) THEN         
c        Solve eq. 2 for q
         Coef(1) = 1.d0
         Coef(2) = 4.d0/3.d0
         Coef(3) = DPPar(1)**2

         CALL CQdrtc(Coef,Sol,NSol)
         DO i1 = 1, NSol
            XSing2(2*i1-1) =  DBLE(SQRT(Sol(i1)))
            XSing2(2*i1)   = -DBLE(SQRT(Sol(i1)))
         END DO

c        Test Solutions
         DO i1 = 1, 2*NSol
            IF((DBLE(XSing2(i1)).ge.DBLE(Limit1)).and.
     &           (DBLE(XSing2(i1)).le.DBLE(Limit2)).and.
     &              (ABS(DIMAG(Sol(i1))).le.Zero)) THEN
               NSing = NSing + 1
               XSing(NSing) = XSing2(i1)
            END IF
         END DO
      END IF
      
c     Sort singularities
      CALL QSORTI(Order,NSing,Xsing)
      DO i1 = 1, NSing
         XSing2(i1) = XSing(i1)
      END DO
      DO i1 = 1, NSing
         XSing(i1) = XSing2(Order(i1))
      END DO

      RETURN
      END
