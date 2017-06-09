      subroutine trig (x, y, z, ct, st, cp, sp)
      implicit double precision (a-h, o-z)
c     returns cos(theta), sin(theta), cos(phi), sin(ph) for (x,y,z)
c     convention - if x=y=0 and z>0, phi=0, cp=1, sp=0
c                  if x=y=0 and z<0, phi=180, cp=-1,sp=0
c                - if x=y=z=0, theta=0, ct=1, st=0
      parameter (eps = 1.0d-6)
      r = sqrt (x**2 + y**2 + z**2)
      rxy = sqrt (x**2 + y**2)
      if (r .lt. eps)  then
         ct = 1
         st = 0
      else
         ct = z/r
         st = rxy/r
      endif
      if (rxy .lt. eps)  then
         cp = 1
         if (ct .lt. 0) cp = -1
         sp = 0
      else
         cp = x / rxy
         sp = y / rxy
      endif
      return
      end


      subroutine arg(c,fi,th)
      implicit double precision (a-h, o-z)
      complex*16  c
      parameter (eps = 1.0d-6)
      x = dble(c)
      y = dimag(c)
      if (abs(x) .lt. eps) x = 0
      if (abs(y) .lt. eps) y = 0
      if (abs(x) .lt. eps  .and.  abs(y) .lt. eps) then
        th = fi
      else
        th = atan2(y,x)
      endif
      return
      end
