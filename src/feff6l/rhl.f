      subroutine rhl (rs, xk, erl, eim)
      implicit double precision (a-h, o-z)

c     input:  rs, xk
c     output: erl, eim

c     This is a new hl subroutine, using interpolation for the
c     real part while the imaginary part is calculated analytically.
c     It uses hl to calculate values at the mesh points for the inter-
c     polation of the real part. The imaginary part is calculated
c     using subroutine imhl.
c
c     written by jose mustre
c     polynomial in rs has a 3/2 power term. j.m.


c     for the right branch the interpolation has the form:
c     hl(rs,x) = e/x + f/x**2 + g/x**3
c     where e is known and
c        f = sum (i=1,3) ff(i) rs**(i+1)/2
c        g = sum (i=1,3) gg(i) rs**(i+1)/2
c
c
c     lrs=number of rs panels, in this case one has 4 panels
c     nrs=number of standard rs values, also order of rs expansion
c     if you change nrs you need to change the expansion of hl
c     in powers of rs that only has 3 terms!
c     nleft=number of coefficients for x<x0
c     nright=number of coefficients for x>x0

      parameter (lrs=4, nrs=3, nleft=4, nright=2)
      include 'const.h'

      dimension cleft(nleft), cright(nright)

      save rcfl, rcfr
      dimension rcfl(lrs,nrs,nleft), rcfr(lrs,nrs,nright)
      data rcfr/-0.173963d+00,-0.173678d+00,-0.142040d+00,-0.101030d+00,
     1     -0.838843d-01,-0.807046d-01,-0.135577d+00,-0.177556d+00,
     2     -0.645803d-01,-0.731172d-01,-0.498823d-01,-0.393108d-01,
     3     -0.116431d+00,-0.909300d-01,-0.886979d-01,-0.702319d-01,
     4      0.791051d-01,-0.359401d-01,-0.379584d-01,-0.419807d-01,
     5     -0.628162d-01, 0.669257d-01, 0.667119d-01, 0.648175d-01/
      data rcfl/ 0.590195d+02, 0.478860d+01, 0.812813d+00, 0.191145d+00,
     1     -0.291180d+03,-0.926539d+01,-0.858348d+00,-0.246947d+00,
     2      0.363830d+03, 0.460433d+01, 0.173067d+00, 0.239738d-01,
     3     -0.181726d+03,-0.169709d+02,-0.409425d+01,-0.173077d+01,
     4      0.886023d+03, 0.301808d+02, 0.305836d+01, 0.743167d+00,
     5     -0.110486d+04,-0.149086d+02,-0.662794d+00,-0.100106d+00,
     6      0.184417d+03, 0.180204d+02, 0.450425d+01, 0.184349d+01,
     7     -0.895807d+03,-0.318696d+02,-0.345827d+01,-0.855367d+00,
     8      0.111549d+04, 0.156448d+02, 0.749582d+00, 0.117680d+00,
     9     -0.620411d+02,-0.616427d+01,-0.153874d+01,-0.609114d+00,
     1      0.300946d+03, 0.109158d+02, 0.120028d+01, 0.290985d+00,
     2      -0.374494d+03,-0.535127d+01,-0.261260d+00,-0.405337d-01/

c
c     calculate hl using interpolation coefficients
      rkf = fa/rs
      ef  = rkf**2/2
      wp  = sqrt (3/rs**3)
      call imhl (rs, xk, eim, icusp)

c     eim already has a factor of ef in it j.m.
c     eim also gives the position of the cusp

      xx = xk / rkf
c     set to fermi level if below fermi level
      if (xx .lt. 1.00001) then
          xx = 1.00001
      endif
c     calculate right hand side coefficients
      if (rs .lt. 0.2) then
         mrs=1
      elseif (rs .lt. 1.0) then
         mrs=2
      elseif (rs .lt. 5.0) then
         mrs=3
      else
         mrs=4
      endif

      do 210 j=1,nright
         cright(j) = rcfr(mrs,1,j)*rs + rcfr(mrs,2,j)*rs*sqrt(rs)
     1               + rcfr(mrs,3,j)*rs**2
  210 continue
      eee=-pi*wp/(4*rkf*ef)

      if (icusp .ne. 1) then
         do 230 j=1,nleft
            cleft(j) = rcfl(mrs,1,j)*rs + rcfl(mrs,2,j)*rs**1.5
     1                 + rcfl(mrs,3,j)*rs**2
  230    continue
         erl=cleft(1)
         do 250 j=2,nleft
            erl=erl+cleft(j)*xx**(j-1)
  250    continue
      else
c        right branch
         erl=eee/xx
         do 280 j=1,nright
            erl=erl+cright(j)/xx**(j+1)
  280    continue
      endif

      erl = erl * ef

      return
      end
