      double precision function xstar(eps1,eps2, vec1,vec2, ndeg, elpty,
     &       ilinit)
      implicit double precision (a-h, o-z)
c     calculating nstar=deg*cos(eps r1)*cos(eps rN)
c     written by alexei ankudinov 08.13.96
c     calculate the plane wave approximation for central atom
c     vec1 - direction to the first atom in path
c     vec2 - direction to the last atom in path
c     ndeg - may be not equal to 'deg' in diff. paths
c             subroutines

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
c      include 'pdata.h'
      integer ilinit 


      dimension eps1(3), eps2(3), vec1(3), vec2(3)

      lfin = ilinit
      x = xxcos(vec1, vec2)
      iav = 1
      y = xxcos(eps1,vec1)
      z = xxcos(eps1,vec2)
      xtemp = ystar(lfin, x, y, z, iav)
      if (elpty .ne. 0.0) then
        y = xxcos(eps2,vec1)
        z = xxcos(eps2,vec2)
        xtemp = xtemp + elpty**2 * ystar(lfin, x, y, z, iav)
      endif
      xstar = ndeg * xtemp /(1+elpty**2)

      return
      end

      double precision function xxcos (veca, vecb)
      implicit double precision (a-h, o-z)
      dimension veca(3), vecb(3)

      x1 = 0
      do 23 j = 1,3
         x1 = x1 + veca(j) * vecb(j)
   23 continue
      xnorma = 0
      xnormb = 0
      do 24 j = 1,3
         xnorma = xnorma + veca(j)**2
         xnormb = xnormb + vecb(j)**2
   24 continue
      xxcos = x1/sqrt(xnorma*xnormb)
      return
      end

      double precision function ystar (lfin, x , y, z, iav)
      implicit double precision (a-h, o-z)
c     
      dimension pln (0:4,4)
      data pln /0.0  , 1.0, 0.0 , 0.0, 0.0,
     2         -0.5  , 0.0, 1.5 , 0.0, 0.0,
     3          0.0  ,-1.5, 0.0 , 2.5, 0.0,
     4          0.375, 0.0,-3.75, 0.0, 4.375/

      pln0  = pln(0,lfin)
      do 40 i = 1, lfin
         pln0  = pln0  + pln(i, lfin) * x**i
  40  continue
      if (iav.eq.0) then
         ystar = pln0/(2*lfin+1)
      else
         pln1  = pln(1,lfin)
         do 50 i = 2, lfin
            pln1  = pln1  + pln(i, lfin)*i*x**(i-1)
  50     continue
         pln2  = 2* pln(2,lfin)
         do 60 i = 3, lfin
            pln2  = pln2  + pln(i, lfin)*i*(i-1)*x**(i-2)
  60     continue
         ytemp = - lfin*pln0 + pln1*(x+y*z) - pln2*(y**2+z**2-2*x*y*z)
         ystar = ytemp * 3/lfin/(4*lfin**2-1)
      endif
      return
      end
