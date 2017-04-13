      subroutine diff (v, dx, n, vm)
      implicit double precision (a-h,o-z)
      complex*16 v(n), vm(n)
      vm(1)=((6.0*v(2)+6.66666666667*v(4)+1.2*v(6))-(2.45*v(1)+7.
     1 5*v(3)+3.75*v(5)+.166666666667*v(7)))/dx
      vm(2)=((6.0*v(3)+6.66666666667*v(5)+1.2*v(7))-(2.45*v(2)+7.
     1 5*v(4)+3.75*v(6)+.166666666667*v(8)))/dx
      nm2=n-2
      do 10 i=3,nm2
   10 vm(i)=((v(i-2)+8.0*v(i+1))-(8.0*v(i-1)+v(i+2)))/12.0/dx
      vm(n-1)=(v(n)-v(n-2))/(2.0*dx)
      vm(n)=(v(n-2)*.5-2.0*v(n-1)+1.5*v(n))/dx
      return
      end
