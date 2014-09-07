      SUBROUTINE mpse(edens,jIntrs,E,iPl,Mu)
      IMPLICIT NONE
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
c     Josh Kas - This subroutine calculates the many pole self
c     energy at each energy pt and at a few r pts
!     and saves it in mpse.bin to be used later by xcpot
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input Variables:
c     densty - density for all unique potentials
c     Mu     - Fermi Energy
      DOUBLE PRECISION edens(251,0:nphx)
      DOUBLE PRECISION Mu, E, MaxDens, MinDens
c     jIntrs - index of the last point before intersitial level
      INTEGER jIntrs, iPl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parameters:
      INTEGER NRPts
      PARAMETER (NRPts = 12)
      
c     Local Variables:
      INTEGER iR, iph, iPole
      INTEGER ios, RsMax, RsMin, Rs, RsInt
      DOUBLE PRECISION WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole),
     &     Densty, dR, delrHL(NRPts), deliHL(NRPts),
     &     Rs1(NRPts)
c      DOUBLE PRECISION MinDns
      COMPLEX*16 ZRnrm, ZTemp

c     Read exc.dat
      OPEN(file='exc.dat', unit=47, status='old',iostat=ios)
      CALL chopen (ios, 'exc.dat', 'ffmod2(xsect)')
      DO iPole = 1, MxPole
         CALL rdcmt(47,'#*cC')
         READ(47,*,END=5) WpCorr(iPole), Gamma(iPole), AmpFac(iPole)
         Gamma(iPole)  = Gamma(iPole)/hart
         WpCorr(iPole) = (WpCorr(iPole)/hart) /
     &        SQRT(3.d0/((3 / (4*pi*edens(jIntrs+1,1))) ** third)**3)
      END DO
 5    CONTINUE
      WpCorr(iPole) = -1.d30
      CLOSE(47)
      
c     find the minimum and maximum densities and calculate RsMin, RsMax
      MaxDens=-1.d30
      MinDens=1.d30
      DO iR = 1, 251
         DO iph = 0, nphx
            Densty=edens(iR,iph)
            IF(MaxDens.lt.Densty) MaxDens=Densty
            IF(MinDens.gt.Densty) MinDens=Densty
         END DO
      END DO

c     Calculate RsMax(MinDens), RsMin(MaxDens)
      RsMax = int(MIN( 10.d0, (3 / (4*pi*MinDens)) ** third ))
      RsMin = int(MAX( 0.001d0, (3 / (4*pi*MaxDens)) ** third ))
      RsInt = int((3 / (4*pi*edens(jIntrs+1,1))) ** third )
      
c     Calculate Sigma on a grid of NRPts points from RsMin to RsMax
c     and write to mpse.bin
      OPEN(file='mpse.bin', unit=47, status='replace',iostat=ios)
      CALL chopen (ios, 'mpse.bin', 'ffmod2(xsect)')
      dR = (RsMax - RsMin)/(NRPts - 1)
      DO iR = 1, NRPts
         delrHL(iR) = 0.d0
         deliHL(iR) = 0.d0
         Rs = int(RsMin + (iR-1)*dR)
         IF(iPl.gt.1) THEN
               IF((iPl.eq.2).or.(iR.eq.NRPts)) THEN                  
c                  CALL CSigZ(E,Mu,Rs,delrHL(iR),deliHL(iR),ZTemp,
c     &                 WpCorr,Gamma,AmpFac)
                  CALL CSigZ(E,Mu,Rs,delrHL(iR),deliHL(iR),ZTemp,
     &                 WpCorr,AmpFac)
               ELSE IF(iPl.eq.3) THEN
                  delrHL(iR) = 0.d0
                  deliHL(iR) = 0.d0
               ELSE
                  delrHL(iR) = delrHL(NRPts)
               END IF
               IF(iR.eq.NRPts) ZRnrm = ZTemp
            ELSE
c               CALL CSigma(E,Mu,Rs1(iR),delrHL(iR),deliHL(iR),WpCorr,
c     &              Gamma,AmpFac)
               CALL CSigma(E,Mu,Rs1(iR),delrHL(iR),deliHL(iR),WpCorr,
     &              AmpFac)
            END IF  
      
c           Write Sigma(Rs,E) to mpse.bin
            WRITE(47,'(6F20.10)') E, Rs, delrHL(iR), deliHL(iR),
     &           DBLE(ZTemp), DIMAG(ZTemp)
         END DO
         CLOSE(47)
         
         RETURN
         END
