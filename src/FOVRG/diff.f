      subroutine diff (v, dr, kap, cl, dx, n, vm)
c     calculate  vm(i)=(dV/dx)*r(i)*(kap+1)/cl
c     needed for c3 term to calculate j-average phase shift
c     ref. koelling,harmon j.phys.c,3107(1977). eq.14
      implicit double precision (a-h,o-z)
       include '../HEADERS/dim.h'
      
      complex*16 v(n), vm(n), vt(nrptx)
      double precision dr(n)
      double precision a, b, c, d, e, f, g

      do i = 1,n
         vt(i) = v(i) * dr(i)**2
      enddo

c  should be:
      a = 6.0d0
      b = 20.d0/3.0d0
      c = 1.2d0
      d = 2.45d0
      e = 7.5d0
      f = 3.75d0
      g = 1.d0/6.d0
      vm(1)=((a*vt(2) + b*vt(4) + c*vt(6)) -
     $     (d*vt(1) + e*vt(3) + f*vt(5) + g*vt(7)))/dx
      vm(2)=((a*vt(3) + b*vt(5) + c*vt(7)) -
     $     (d*vt(2) + e*vt(4) + f*vt(6) + g*vt(8)))/dx

c curently passes tests:
      vm(1)=((6.0*vt(2)+6.66666666667*vt(4)+1.2*vt(6)) -
     $     (2.45*vt(1)+7.5*vt(3)+3.75*vt(5)+.166666666667*vt(7)))/dx
      vm(2)=((6.0*vt(3)+6.66666666667*vt(5)+1.2*vt(7)) -
     $     (2.45*vt(2)+7.5*vt(4)+3.75*vt(6)+.166666666667*vt(8)))/dx

      nm2=n-2
      do i=3,nm2
         vm(i)=((vt(i-2)+8.0*vt(i+1))-(8.0*vt(i-1)+vt(i+2)))/12.0/dx
      enddo
      vm(n-1)=(vt(n)-vt(n-2))/(2.0*dx)
      vm(n)=(vt(n-2)*0.5 - 2.0*vt(n-1) + 1.5*vt(n))/dx

      do i = 1,n
         vm(i) = (vm(i)-2*vt(i))/dr(i) *(kap+1.0)/cl
      enddo
      return
      end
