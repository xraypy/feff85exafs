! Written by J. Kas, 2006
! eps2exc takes loss function from loss.dat and gives poles and weights
! corresponding to the many pole model detailed in
! Phys. Rev. B 76, 195116 (2007). 
! First, regions are chosen 
      PROGRAM eps2exc
      DOUBLE PRECISION Dat(50000,2), g(10000), gamma, omi(10000),
     &     Delta(10000),eps0, sumrl, xNElec, csumrl
      INTEGER i,i1,NPoles,NPts
c     INTEGER ios1,ios2,imax,NFile
      CHARACTER*80 infl, outfl
c     CHARACTER*80 comment
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
