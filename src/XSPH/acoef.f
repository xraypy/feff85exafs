      subroutine acoef(ispin, amat)
c     performs the sum of the product of 4 3j symbols
c       ispin - type of spin calculation
c       amat  - matrix to calculate density via
c             \mu=\mu_at*(1- Im \sum_kp,kpp rkk(kp)*rkk(kpp)*
c                         \sum_m1,m2 bmat(kp,kpp,m1,m2)*G_lp,m2;lpp,m1 )

      implicit double precision (a-h,o-z)
      include '../HEADERS/dim.h'
      integer ispin

      real  amat(-lx:lx,2,2, 3,0:lx)
      real t3j(0:lx,0:lx,0:1), operls (0:1, 0:1, 3)
      real xms, xml, xmj

      external cwig3j

      do i5 =  0,lx
         do i4 =  1,3 
            do i3 =  1,2
               do i2 =  1,2
                  do i1 = -lx,lx
                     amat(i1,i2,i3,i4,i5)=0
                  enddo
               enddo
            enddo
         enddo
      enddo
      ms = 1
      if (ispin.lt.0) ms=0
      print*, ' Spin = ', 2*ms-1

      do ml = -lx, lx
        mj = 2*ml + (2*ms-1)
        xmj = 0.5e0*mj
        mj = -mj
c       mj is conserved for all operators of interst (s_z, l_z, t_z)
c       tabulate necessary Clebsh-Gordon coefficients
        do lp = 0,lx
           do jp = 0,lx
              do mp = 0,1
                 lp2 = 2*lp
                 jp2 = 2*jp+1
                 mp2 = 2*mp-1
                 t3j(lp,jp,mp) = (-1)**lp *sqrt(jp2+1.e0) * 
     1                real( cwig3j ( 1, jp2, lp2, mp2, mj, 2) )
              enddo
           enddo
        enddo

        do lpp = 0,lx
           do m1=0,1
              do m2=0,1
                 do iop=1,3
                    operls(m1,m2,iop) = 0
                    if (m1.eq.m2) then
                       xms =  m1 - 0.5e0
                       xml = xmj-xms
                       if (abs(ml+ms-m1).le.lpp) then
                          if (ispin.eq.0) then
c                 occupation numbers N_l, N_j- , N_j+
                             operls(m1,m2,iop) = 2
                          elseif (iop.eq.1) then
c                 s_z operator in ls basis
                             operls(m1,m2,iop) = xms
                          elseif (iop.eq.2 .and. abs(ispin).eq.1) then
c                 l_z operator
                             operls(m1,m2,iop) = xml
                          elseif (iop.eq.2 .and. abs(ispin).eq.2) then
c                 unit operator for occupation numbers
                             operls(m1,m2,iop) = 1
                          elseif (iop.eq.3 .and. abs(ispin).eq.1) then
c                 t_z operator
                             operls(m1,m2,iop) = xms*2*(3*xml**2-
     $                            lpp*(lpp+1))/(2*lpp+3) /(2*lpp-1)
                          elseif (iop.eq.3 .and. abs(ispin).eq.2) then
c                 occupation number for j=l+1/2
                             operls(m1,m2,iop) = t3j(lpp,lpp,m1)**2
                          endif
                       endif
                    else
c             t_z only has nonzero off diagonal matrix elements 
                       if (iop.eq.3 .and. abs(ispin).le.1 .and.
     1                      nint( 0.5e0+abs(xmj)).lt.lpp)  then
                          operls(m1,m2,iop)=3*xmj*sqrt(lpp*(lpp+1)-
     $                         (xmj**2-0.25e0)) /(2*lpp+3) /(2*lpp-1)
                       elseif (iop.eq.3 .and. abs(ispin).gt.1) then
                          operls(m1,m2,iop)= t3j(lpp,lpp,m1) *
     $                         t3j(lpp,lpp,m2)
                       endif
                    endif
                 enddo
              enddo
           enddo
           

c         calculate energy and r independent matrix amat
c         which is equivalent to integration over angular coordinates
c         for assumed density matrix
           do i1=1,2
              call kfromi(i1,lpp,jj,k1)
              if (k1.ne.0) then
                 do i2=1,2
                    call kfromi(i2,lpp,jp,k2)
                    if (k2.ne.0) then
                       do iop=1,3
                          do m2=0,1
                             do m1=0,1
                                amat(ml,i1,i2,iop,lpp) =
     $                               amat(ml,i1,i2,iop,lpp) +
     1                               operls(m1,m2,iop) *t3j(lpp,jp,ms)*
     $                               t3j(lpp,jp,m1)*
     $                               t3j(lpp,jj,m2)*t3j(lpp,jj,ms)
                             enddo
                          enddo
                       enddo
                    endif
                 enddo
              endif
           enddo
        enddo
      enddo
      
      return
      end

      subroutine kfromi (i, lpp, jj, k)
c     input index i1 and orb. mom. lpp
c     output: final state kappa - k; jj=tot.mom(k)-1/2
      integer i, lpp, jj, k

      jj = lpp + i - 2
      k = - lpp - 1
      if (i.eq.1) k = lpp

      return
      end
