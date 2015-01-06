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
      INTEGER maxPts, MaxPoles, iMax, i1, i2
      PARAMETER (maxPts = 50000)
      DOUBLE PRECISION M0, M0i, M1i, MM1i, M1Old, M1New, omNew, omOld,
     $     pi, Tmp1(maxPts), Tmp2(maxPts), IntGrl, a, b, c, frac, egap
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
      subroutine terp (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polint( x(k), y(k), m+1, x0, y0, dy)

      return
      end


      subroutine polint( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.
                                                                                                                                                             
      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)
                                                                                                                                                             
      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) pause 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue
                                                                                                                                                             
      return
      end


      function locat (x, n, xx)
      integer  u, m, n
      double precision x, xx(n)
                                                                                                                                                             
c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)
                                                                                                                                                             
      locat = 0
      u = n+1
                                                                                                                                                             
   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif
                                                                                                                                                             
      return
      end

      DOUBLE PRECISION FUNCTION IntGrl(a, b, x, y, n)
      INTEGER n, i, iStart, iLast
      DOUBLE PRECISION x(n), y(n), dx, a, b, aEnd, bEnd,
     &     ya, yb

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
      SUBROUTINE rdloss(Uloss,Dat,NPts)
      IMPLICIT NONE
c     Variables passed
      INTEGER Uloss, NPts
      DOUBLE PRECISION Dat(50000,2)

      
c     local variables
      INTEGER nterp,n,i
      DOUBLE PRECISION dum
      CHARACTER(1) comment, ch
      CHARACTER(200) Line
      comment='#*!c'
      nterp=1
      n=0

c     Read data into Dat
 200  DO i=1,50000
c     Read past comments
         CALL rdcmt(Uloss,comment)
         READ(Uloss,*,end=250) Dat(i,1), Dat(i,2)
         n=i
      END DO      
 250  CONTINUE
      NPts=n
      
300   CLOSE(Uloss)
      RETURN 
      END
      SUBROUTINE rdcmt(iUnt,Cmt)
      INTEGER iUnt, i1
      CHARACTER(300) line
      CHARACTER(4) Cmt
      CHARACTER TmpCmt(4), ch
      LOGICAL CmtLin

      CmtLin = .true.
      DO i1 = 1, 4
         TmpCmt(i1) = Cmt(i1:i1)
      END DO
 5    CONTINUE
      READ(iUnt,*,END=10) ch
      DO i1 = 1, 4
         IF(ch.eq.TmpCmt(i1)) goto 5
      END DO
      
      BACKSPACE(iUnt)
      
 10   CONTINUE
      
      RETURN
      END
! Written by J. Kas, 2006
! eps2exc takes loss function from loss.dat and gives poles and weights
! corresponding to the many pole model detailed in
! Phys. Rev. B 76, 195116 (2007). 
! First, regions are chosen 
      PROGRAM eps2exc
      DOUBLE PRECISION Dat(50000,2), g(10000), gamma, omi(10000),
     &     Delta(10000),eps0, sumrl, xNElec, csumrl
      INTEGER i,i1,ios1,ios2,imax, NPoles,NPts,NFile
      CHARACTER*80 infl, outfl, comment
      CHARACTER ch
      EXTERNAL IntGrl
      NPts = 1000
      infl='loss.dat'
      outfl='exc.dat'
      PRINT*, '# Enter number of poles:'
      READ*, NPoles
      ! This input can be used to correct the sumrules
      !PRINT*, '# Enter eps^-1 sumrule and N_el:'
      !READ*, sumrl, xNElec
      sumrl = 1.d0
      xNElec = 1.d0
      PRINT*, 'Is this a metal? (y/n)'
      READ*, ch
      IF(ch.EQ.'y'.OR.ch.EQ.'Y') THEN
         eps0 = -1.d0
      ELSE
         PRINT*, 'Would you like to set the dielectric constant? (y/n)'
         READ*, ch
         IF(ch.EQ.'n'.OR.ch.EQ.'N') THEN
            eps0 = -2.d0
         ELSE
       
            ! This input can be used to correct the dielectric constant,
            ! which is related to the inverse moment.
            ! Use eps0 = -2 to ignore this correction.
            PRINT*, '# Enter dielectric constant: '
            READ*, eps0
         END IF
      END IF 
      PRINT*, eps0
      gamma = 0.01
      csumrl= xNElec/sumrl
      
      OPEN(unit=12,file=infl,status='old')
      OPEN(unit=13,file=outfl,status='replace')
      CALL rdloss(12,Dat,NPts)

      DO i1 = 1, NPts
         Dat(i1,2) = Dat(i1,2)*csumrl
      END DO
      ! getomi finds poles and weights
      CALL getomi(Dat(1,1), Dat(1,2), NPts, NPoles, omi, g, Delta, eps0)
      WRITE(13,'(A33,I4,A6)') '# Loss function represented with ',
     $     NPoles, ' poles'
      WRITE(13,'(A23,f8.4)') '# Dielectric constant: ', eps0
      DO i = 1, NPoles
         WRITE(13,'(4f20.10)'), omi(i), omi(i)*gamma, g(i), Delta(i)
      END DO
      
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
         PRINT*, 'sumg1, sumg2', sumg1, sumg2, MM1, dom1, dom2
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
      
