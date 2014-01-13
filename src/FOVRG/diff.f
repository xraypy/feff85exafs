      subroutine diff (v, dr, kap, cl, dx, n, vm)
c     calculate  vm(i)=(dV/dx)*r(i)*(kap+1)/cl
c     needed for c3 term to calculate j-average phase shift
c     ref. koelling,harmon j.phys.c,3107(1977). eq.14
      implicit double precision (a-h,o-z)
       include '../HEADERS/dim.h'
      
      complex*16 v(n), vm(n), vt(nrptx)
      dimension dr(n)
      do 5 i = 1,n
 5    vt(i) = v(i) * dr(i)**2

      vm(1)=((6.0*vt(2)+6.66666666667*vt(4)+1.2*vt(6))-(2.45*vt(1)+7.
     1 5*vt(3)+3.75*vt(5)+.166666666667*vt(7)))/dx
      vm(2)=((6.0*vt(3)+6.66666666667*vt(5)+1.2*vt(7))-(2.45*vt(2)+7.
     1 5*vt(4)+3.75*vt(6)+.166666666667*vt(8)))/dx
      nm2=n-2
      do 10 i=3,nm2
   10 vm(i)=((vt(i-2)+8.0*vt(i+1))-(8.0*vt(i-1)+vt(i+2)))/12.0/dx
      vm(n-1)=(vt(n)-vt(n-2))/(2.0*dx)
      vm(n)=(vt(n-2)*.5-2.0*vt(n-1)+1.5*vt(n))/dx

      do 20 i = 1,n
 20   vm(i) = (vm(i)-2*vt(i))/dr(i) *(kap+1.0)/cl
      return
      end
