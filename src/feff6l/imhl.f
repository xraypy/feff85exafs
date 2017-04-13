      subroutine imhl (rs, xk, eim, icusp)
      implicit double precision (a-h,o-z)

c     what is xk?  k**2 - mu + kf**2?

c written by j. mustre (march 1988)
c code is based on analytical expression derived by john rehr.
c it leaves the real part, calculated in rhl unchanged.
c
c modified by j. rehr  (oct 1991) - adds quinn approximation for
c losses due to electron-hole pairs below the plasmon turn on
c see new subroutine quinn.f, which incorporates r. albers coding of
c j.j. quinn's approximations for details.

      include 'const.h'
c     alph is Hedin-Lundquist parameter
      parameter (alph = 4.0 / 3.0)
      external ffq

      integer icount
      save icount
      data icount /0/

      icusp=0
      xf = fa / rs
      ef = xf**2 / 2

c     xk0 is xk normalized by k fermi.
      xk0 = xk/xf
c     set to fermi level if below fermi level
      if (xk0 .lt. 1.00001) then
         xk0 = 1.00001
      endif

c     wp is given in units of the fermi energy in the formula below.
      wp = sqrt (3 / rs**3) / ef
      xs = wp**2 - (xk0**2 - 1)**2

      eim = 0
      if (xs .lt. 0.)  then
         q2 = sqrt ( (sqrt(alph**2-4*xs) - alph) / 2 )
         qu = min (q2, (1+xk0))
         d1 = qu - (xk0 - 1)
         if (d1 .gt. 0)  then
            eim = ffq (qu,ef,xk,wp,alph) - ffq (xk0-1,ef,xk,wp,alph)
         endif
      endif
      call cubic (xk0, wp, alph, rad, qplus, qminus)

      if (rad .le. 0) then
         d2 = qplus - (xk0 + 1)
         if (d2 .gt. 0)  then
            eim = eim + ffq (qplus,ef,xk,wp,alph) - 
     1                  ffq (xk0+1,ef,xk,wp,alph)
         endif
         d3 = (xk0-1) - qminus
         if (d3 .gt. 0)  then
            eim = eim + ffq (xk0-1,ef,xk,wp,alph) - 
     1                  ffq (qminus,ef,xk,wp,alph)
c           beginning of the imaginary part and position of the cusp x0
            icusp = 1
         endif
      endif

      call quinn (xk0, rs, wp, ef, ei)
      if (eim .ge. ei)  eim = ei

      icount = icount+1
      return
      end
