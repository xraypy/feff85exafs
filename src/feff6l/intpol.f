      subroutine intpol (a,b,fa,fb,fma,fmb,x,fx,fmx)
      implicit double precision (a-h,o-z)
c     Only output is fx, fmx
      complex*16 fa,fb,fma,fmb,fx,fmx
      dx=b-a
      d=(x-a)/dx
c      if (d*(1.0-d).lt.0.0) stop 'Died in intpol'
      if (d*(1.0-d).lt.0.0) then
c         print*, 'a, b, dx'
c         print*, a, b, dx
c         print*, 'x, x-a'
c         print*, x, x-a
c         print*, 'd, d*(1-d)'
c         print*, d, d*(1-d)
         call fstop(' at INTPOL')
      endif
      c2=3.0*(fb-fa)-(fmb+2.0*fma)*dx
      c3=2.0*(fa-fb)+(fma+fmb)*dx
      fx=fa+d*(dx*fma+d*(c2+d*c3))
      fmx=fma+d*(2.0*c2+3.0*c3*d)/dx
      return
      end
