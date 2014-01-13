      PROGRAM test
      INCLUDE '../../HEADERS/dim.h'
      INTEGER iprint, ispec, ne, ne1, ik0, ne3, i1, igrid
      DOUBLE PRECISION edge, emu, vi0, gamach, xkmax, xkstep, vixan, ecv
      COMPLEX*16 em1(nex), em2(nex)
      iprint = 0
      ispec = 1
      ne = 0
      ne1 = 0
      ik0 = 0
      ne3 = 0
      edge = 0.d0
      emu = 10.d0
      vi0 = 0.d0
      gamach = 1.d0
      xkmax = 10.1d0
      xkstep = 1.05d0
      vixan = 1.d0
      ecv = 0.d0
      igrid = 1
      
c      CALL phmesh(iprint, ispec, edge, emu, vi0, gamach, ecv,
c     &                  xkmax, xkstep, vixan, ne, ne1, em1, ik0, ne3)

      CALL phmesh2(iprint, ispec, edge, emu, vi0, gamach,
     &     xkmax, xkstep, vixan, ne, ne1, em2, ik0, ne3, 1)
c      DO i1 = 1, ne
c         print*, i1, SQRT(2*ABS(DBLE(em1(i1)-edge)))/0.529,
c     &        SQRT(2*ABS(DBLE(em2(i1)-edge)))/0.529
c      END DO
      PRINT*, 'next'
      DO i1 = 1, ne
c     &        DIMAG(em1(i1)), DIMAG(em2(i1))
         print*, i1, DBLE(em2(i1))*27.2, DIMAG(em2(i1))*27.2
c     &        DIMAG(em1(i1)), DIMAG(em2(i1))
      END DO

      END
