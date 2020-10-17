      subroutine diff (v, dr, kap, cl, dx, n, vm)
c     calculate  vm(i)=(dV/dx)*r(i)*(kap+1)/cl
c     needed for c3 term to calculate j-average phase shift
c     ref. koelling,harmon j.phys.c,3107(1977). eq.14
      implicit double precision (a-h,o-z)
       include '../HEADERS/dim.h'
      
      complex*16 v(n), vm(n), vt(nrptx)
      dimension dr(n)
      do i = 1,n
         vt(i) = v(i) * dr(i)**2
      enddo
      vm(1)=((6.0d0*vt(2)+6.66666666667d0*vt(4)+1.2d0*vt(6))-
     $     (2.45d0*vt(1)+7.5d0*vt(3)+3.75d0*vt(5)+
     $     0.166666666667d0*vt(7)))/dx
      vm(2)=((6.0d0*vt(3)+6.66666666667d0*vt(5)+1.2d0*vt(7))
     $     -(2.45d0*vt(2)+7.5d0*vt(4)+3.75d0*vt(6)+
     $     0.166666666667d0*vt(8)))/dx
      nm2=n-2
      do i=3,nm2
         vm(i)=((vt(i-2)+8.d00*vt(i+1))-(8.d0*vt(i-1)+
     $        vt(i+2)))/12.d0/dx
      enddo
      vm(n-1)=(vt(n)-vt(n-2))/(2.0*dx)
      vm(n)=(vt(n-2)*.5-2.0d0*vt(n-1)+1.5d0*vt(n))/dx

      do i = 1,n
         vm(i) = (vm(i)-2*vt(i))/dr(i) *(kap+1.0)/cl
      enddo
      return
      end
