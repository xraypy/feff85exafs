      subroutine inipot (dgc, dpc, edenvl, vvalgs, xnmues)
c     initialize values of arrays to zero
      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      parameter (zero=0.0d0)

      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension edenvl(251,0:nphx), vvalgs (251,0:nphx)
      dimension xnmues(0:lx,0:nphx)

      do iph  = 0,nphx+1
         do iorb = 1,30
            do i = 1,251
               dgc(i,iorb,iph) = zero
            enddo
         enddo
      enddo

      do iph  = 0,nphx+1
         do iorb = 1,30
            do i = 1,251
               dpc(i,iorb,iph) = zero
            enddo
         enddo
      enddo
      do iph = 0, nphx
         do i = 1, 251
            edenvl(i, iph) = zero
         enddo
      enddo
      do iph = 0, nphx
         do i = 1, 251
            vvalgs(i, iph) = zero
         enddo
      enddo
      do iph = 0, nphx
         do ll = 0, lx
            xnmues (ll, iph) = zero
         enddo
      enddo
      return
      end
