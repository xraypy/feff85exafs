! Subroutine getomi finds poles and weights for representing epsInv.
! Written by J. Kas, 2006
      SUBROUTINE getomi(wgrid, epsInv, NPts, NPoles, omi, gi, Deltai,
     &     eps0)
      IMPLICIT NONE
!     Input variables
!     wGrid(NPts)  - grid of omega
!     epsInv(NPts) - loss function
!     NPoles       - Number of poles to represent epsInv
      
!     Output
!     omi          - omega position of poles
!     gi           - pole strengths
!     Deltai       - energry region about omi where
      INTEGER NPoles, NPts
      DOUBLE PRECISION wGrid(NPts), epsInv(NPts), omi(NPoles),
     &     gi(NPoles), Deltai(NPoles), eps0

!     Local variables
      INTEGER maxPts, i1, i2
      PARAMETER (maxPts = 50000)
      DOUBLE PRECISION M0, M0i, M1i, MM1i, M1Old, M1New, omNew, omOld,
     $       pi, Tmp1(maxPts), Tmp2(maxPts), IntGrl, a, b, c, frac

c     INTEGER MaxPoles, iMax, 
c     DOUBLE PRECISION egap
      PARAMETER (pi = 3.14159 26535 89793 23846 26433d0)
      EXTERNAL IntGrl
!     Find zeroth moment
!      frac = (1.d0-0.01d0)
      frac = 0.999d0
      M0 = 0.5d0*epsInv(1)*wGrid(1)
      DO i1 = 2, NPts
         M0 = M0 + 0.5d0*(epsInv(i1)+epsInv(i1-1))*
     &        (wGrid(i1)-wGrid(i1-1))
         Tmp1(i1) = epsInv(i1)/wGrid(i1)
         Tmp2(i1) = epsInv(i1)*wGrid(i1)
      END DO
!     Find gap energy
!      DO i1 = 1, NPts
      
!     Now find energy regions such that the integral of epsinv in each region is
!     equal.
      omOld = 0.d0
      omNew = 0.d0
      M1Old = 0.d0
      M1New = 0.d0
      b     = 0.d0
      ! Save the first 10 poles for low energies.
      DO i1 = 1, NPoles
         M0i = 0.5d0*epsInv(1)*wGrid(1)
         M1New = 0.5d0*epsInv(1)*wGrid(1)**2
         DO i2 = 2, NPts
            IF(M0i.ge.(DBLE(i1)*frac*M0/DBLE(NPoles))) THEN
               omOld = omNew
               omNew = wGrid(i2-1)
               IF(i2.eq.1) THEN
                  M0i = 0.d0
               ELSEIF(i2.eq.2) THEN
                  M0i = 0.5d0*epsInv(1)*wGrid(1)
               ELSE
                  M0i = M0i - 0.5d0*(epsInv(i2-1)+ epsInv(i2-2))*
     $              (wGrid(i2-1)-wGrid(i2-2))
               END IF
               GOTO 5
            ELSE
               M0i = M0i + 0.5d0*(epsInv(i2)+epsInv(i2-1))*
     $              (wGrid(i2)-wGrid(i2-1))
               M1New = M1New + 0.5d0*(epsInv(i2)*wGrid(i2)+ epsInv(i2-1)
     $              *wGrid(i2-1))*(wGrid(i2)-wGrid(i2-1))
            END IF
         END DO
 5       CONTINUE
 
!        The region has been found up to the wGrid accuracy. Now interpolate
!        to find omi
         IF(i2.eq.2) THEN
            a = epsInv(i2-1)/omNew
         ELSE
            a = (epsInv(i2-1) - epsInv(i2-2))/(omNew - wGrid(i2-2))
            b = epsInv(i2-2)
         END IF
         a = MAX(a, 1.d-6)
         c = 2.d0/a*(M0i - frac*M0*(DBLE(i1)/DBLE(NPoles)))

         omNew = wGrid(i2-2) - b/a + SQRT((b/a)**2 - c)
         IF((omNew.lt.wGrid(i2-2)).or.(omNew.gt.wGrid(i2-1))) THEN
            omNew = wGrid(i2-2) - b/a - SQRT((b/a)**2 - c)
         END IF

         IF(.FALSE.) THEN
            omOld = DBLE(i1-1)*500.d0/DBLE(NPoles)
            omNew = DBLE(i1)*500.d0/DBLE(NPoles) 
         END IF
!        Now set omi and g s.t. the inverse and first moments are preserved for
!        this region.
         MM1i = 0.d0
         M1i  = 0.d0
         IF(omOld.eq.0.d0) THEN
            MM1i = 0.5d0*epsInv(1)
            M1i  = 0.5d0*epsInv(1)*wGrid(1)**2
            omOld = wGrid(1)
         END IF
         MM1i = MM1i + IntGrl(omOld, omNew, wGrid, Tmp1, NPts)
         M1i  = M1i + IntGrl(omOld, omNew, wGrid, Tmp2, NPts)
         gi(i1)  = 2.d0/pi*MM1i
         omi(i1) = SQRT(M1i/MM1i)
         Deltai(i1) = omNew-omOld
         M1Old  = M1New
         !PRINT '(a27,i6,20f20.10)', 'i1, omi, omOld, omNew, a, b', i1,
         !&    omi(i1), omOld, omNew, a, b
      END DO

!     Shift and scale the poles to set inverse moment to 1-1/eps0, without changing
!     first moment. (eps0 is input by user)
      CALL getdom(gi,omi,NPoles,eps0)
      
      RETURN
      END

      DOUBLE PRECISION FUNCTION IntGrl(a, b, x, y, n)
      INTEGER n, i, iStart, iLast
      DOUBLE PRECISION x(n), y(n), dx, a, b, 
     &     ya, yb
c      DOUBLE PRECISION aEnd, bEnd,
      
      iStart = -1
      iLast = -1
      IntGrl = 0.d0
      dx = x(2) - x(1)
      
      DO i = 1, n
         IF(x(i).ge.a.and.x(i).lt.b) THEN
            iStart = i
            goto 5
         ELSE
            iStart = -1
         END IF
      END DO
 5    CONTINUE
      
      DO i = n, 1, -1
         IF(x(i).le.b.and.x(i).gt.a) THEN
            iLast = i
            goto 10
         ELSE
            iLast = -1
         END IF
      END DO
 10   CONTINUE

      IF (iStart.lt.0.and.iLast.lt.0) THEN
         CALL terp(x, y, n, 1, a, ya)
         CALL terp(x, y, n, 1, b, yb)
         IntGrl = IntGrl + 0.5d0*(yb+ya)*(b-a)
      ELSE
         DO i = iStart, iLast-1
            dx=x(i+1)-x(i)
            IntGrl = IntGrl + 0.5d0*(y(i)+y(i+1))*dx
         END DO
         
         CALL terp(x, y, n, 1, a, ya)
         CALL terp(x, y, n, 1, b, yb)
         
         IntGrl = IntGrl + 0.5d0*(ya+y(iStart))*(x(iStart) - a)
         IntGrl = IntGrl + 0.5d0*(yb+y(iLast))*(b-x(iLast))
      END IF
      RETURN
      END

! This will change the pole representation of the loss function without
! changing the f-sum rules. Written by J. Kas 2006
      SUBROUTINE getdom(gi,omi,NPoles,eps0)
      IMPLICIT NONE

!     Input/output
      INTEGER NPoles
      DOUBLE PRECISION gi(NPoles), omi(NPoles), eps0

!     Local
      INTEGER i1, i2, MxIter
      DOUBLE PRECISION MM1, tol, dom, dom1, dom2, sumg1, sumg2
      PARAMETER(MxIter = 10000, tol = 1.d-6)
      IF(eps0.lt.-1.5d0) RETURN
      IF(eps0.lt.0.d0) THEN
         MM1 = 1.d0
      ELSE
         MM1   = 1.d0 - 1.d0/eps0
      END IF
      dom1  = -omi(1)*0.999d0
      dom2  = omi(1)*0.999d0
!     Bracket dom between dom1 and dom2
      DO i1 = 1, MxIter
         sumg1 = 0.d0
         sumg2 = 0.d0
         DO i2 = 1, NPoles
            sumg1 = sumg1 + gi(i2)*(omi(i2)/(omi(i2)-dom1))**2
            sumg2 = sumg2 + gi(i2)*(omi(i2)/(omi(i2)-dom2))**2
         END DO
         IF((sumg1.lt.MM1).and.(sumg2.gt.MM1)) GOTO 5
         !IF(sumg1.gt.MM1) dom1 = dom1*2.d0
         !IF(sumg2.lt.MM1) dom2 = dom2/2.d0
         PRINT*, 'Error in getdom.f.'
         PRINT*, 'sumg1, sumg2', i1, sumg1, sumg2, MM1, dom1, dom2
         STOP
      END DO
 5    CONTINUE                  ! dom bracketed

      dom = 0.5d0*(dom1+dom2)
      DO i1 = 1, MxIter
         sumg1 = 0.d0
         DO i2 = 1, NPoles
            sumg1 = sumg1 + gi(i2)*(omi(i2)/(omi(i2)-dom))**2
         END DO
         IF(abs((sumg1-MM1)/MM1).lt.tol) GOTO 10
         IF((sumg1-MM1).gt.0.d0) THEN
            dom2 = dom
         ELSE
            dom1 = dom
         END IF
         dom = 0.5d0*(dom1+dom2)
      END DO
 10   CONTINUE

      DO i1 = 1, NPoles
         gi(i1)  = gi(i1)*(omi(i1)/(omi(i1)-dom))**2
         omi(i1) = omi(i1)-dom
      END DO
      RETURN
      END
      
      
