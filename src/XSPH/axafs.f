      subroutine axafs(em, emu, xsec,ne1,ik0)
c     extract axafs from xsec
c     written by a.l.ankudinov Dec. 1998

c     the file axafs.dat (format as in xmu.dat) will be written if
c     you use PRINT 0 1 0 0 0 0 (ipr2 > 0), and ran the second module.

c     the code draws a parabola using least mean square method
c     through xsec(i) * ee (i)**xn 
c     the weight for each point i, is defined as (ee(i)-E_F)**mm*
c     (ee(i+1)- ee(i-1)), where the last multiplier is used since the 
c     grid is not regular in energy.
c     E_F - energy that corresponds to Fermi level.

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

      complex*16 em(nex), xsec(nex)
      dimension ee(nex), xmu(nex), wt(nex)
      dimension xx(0:4), yy(0:2), xm(3,3)

c     empirically I found that the best curve is drawn if xn=0 and mm=1
c     alex ankudinov, january 1999.
      xn = 0
      mm = 1
      np = ne1 - ik0
      ef = emu

      do ie = 1, np
         ee(ie) = dble(em(ik0+ie)-em(ik0)) +emu
         xmu(ie) = dimag(xsec(ik0+ie)) * ee(ie)**xn
      enddo
      do ie = 1, np
         if (ie.eq.1) then
            wt(ie) = (ee(ie+1)-ef) * (abs(ee(ie)-ef))**mm
         elseif (ie.eq.np) then
            wt(ie) = (ee(ie)-ee(ie-1)) * (ee(ie)-ef)**mm
         else
          wt(ie) = (ee(ie+1)-ee(ie-1)) * (ee(ie)-ef)**mm
       endif
      enddo
      do i = 0, 4
         xx(i) = 0
      enddo
      do i = 0, 2
         yy(i) = 0
      enddo

      do ie = 1, np
         do i = 0,4
            xx(i) = xx(i) + wt(ie)*ee(ie)**i
         enddo
         do i = 0,2
            yy(i) = yy(i) + wt(ie)*xmu(ie)*ee(ie)**i
         enddo
      enddo

      do i=1,3
         do j=1,3
            xm(i,j) = xx(i+j-2)
         enddo
      enddo
      denom = determ (xm, 3, 3)

      do i=1,3
         do j=1,3
            xm(i,j) = xx(i+j-2)
         enddo
      enddo
      do i=1,3
         xm(i,1) = yy (i-1)
      enddo
      aa = determ (xm,3,3)
      aa = aa / denom

      do i=1,3
         do j=1,3
            xm(i,j) = xx(i+j-2)
         enddo
      enddo
      do i=1,3
         xm(i,2) = yy (i-1)
      enddo
      bb = determ (xm,3,3)
      bb = bb / denom

      do i=1,3
         do j=1,3
            xm(i,j) = xx(i+j-2)
         enddo
      enddo
      do i=1,3
         xm(i,3) = yy (i-1)
      enddo
      cc = determ (xm,3,3)
      cc = cc / denom

c     find normalization at edge+100 eV
      eee = ee(1) + 100/hart
      xnorm = (aa+bb*eee+cc*eee**2) / eee**xn

      open (unit=1,file='axafs.dat', status='unknown')
      write (1,*) '# File contains AXAFS. See manual for details.'
      write (1,*)
     1 '#--------------------------------------------------------------'
      write(1,*) '#  e, e(wrt edge), k,',
     1           ' mu_at=(1+chi_at)*mu0_at, mu0_at, chi_at @#'
      do ie = 1, np
        xmu(ie) = dimag(xsec(ie+ik0))
        xmu0 = (aa+bb*ee(ie)+cc*ee(ie)**2) / ee(ie)**xn
        chiat = (xmu(ie) - xmu0) / xmu0
        eee = ee(ie) -ef
        if (eee.ge.0.d0) then
           xk = sqrt(2*eee) /bohr
        else
           xk = -sqrt(-2*eee) /bohr
        endif
        write (1, 410) ee(ie)*hart, (ee(ie)-emu)*hart, xk,
     1              xmu(ie)/xnorm, xmu0/xnorm, chiat
 410    format (1x, 2f11.3, f8.3, 1p, 3e13.5)
      enddo
      close (unit=1)

      return
      end
         

