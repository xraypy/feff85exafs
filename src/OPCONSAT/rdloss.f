      SUBROUTINE rdloss(Uloss,Dat,NPts)
      IMPLICIT NONE
c     Variables passed
      INTEGER Uloss, NPts
      DOUBLE PRECISION Dat(50000,2)

      
c     local variables
      INTEGER nterp,n,i
c      DOUBLE PRECISION dum
      CHARACTER(LEN=4) comment
c      CHARACTER ch
c      CHARACTER(200) Line
      comment='#*!c'
      nterp=1
      n=0

c     Read data into Dat
      DO i=1,50000
c     Read past comments
         CALL rdcmt(Uloss,comment)
         READ(Uloss,*,end=250) Dat(i,1), Dat(i,2)
         n=i
      END DO      
 250  CONTINUE
      NPts=n
      
      CLOSE(Uloss)
      RETURN 
      END
