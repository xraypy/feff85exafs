      subroutine setgam (iz, ihole, gamach)

c     Sets gamach, core hole lifetime.  Data comes from graphs in
c     K. Rahkonen and K. Krause,
c     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.

      implicit double precision (a-h, o-z)

      dimension gamk(6), zk(6)
      dimension gaml1(6), zl1(6)
      dimension gaml2(6), zl2(6)
      parameter (ryd  = 13.6058)

      save ienter

c     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff
c     trouble.
c     Gam arrays contain the gamma values.
c     We will take log10 of the gamma values so we can do linear
c     interpolation from a log plot.

      data  zk   / 0.99,  10.0, 20.0,  40.0,  60.0,   95.1/
      data  gamk / 0.07,   0.3,  0.75,  5.0,  20.0,  100.0/

      data  zl1   / 0.99,  20.0, 35.0, 50.0,  75.0,  95.1/
      data  gaml1 / 0.07,   4.0,  7.0,  4.0,   8.0,  19.0/

      data  zl2   / 0.99,  26.0, 31.0, 60.0,  80.0,  95.1/
      data  gaml2 / 0.001,  1.7,  0.8,  3.5,   5.0,  10.0/

      data ienter /0/

c     Call this only once, if it gets called a second time the gamma
c     values will be messed up by repeated taking of log10

      if (ienter .gt. 0)  then
         call fstop(' at SETGAM-1: re-entered SETGAM')
      endif
      ienter = 1

      if (ihole .le. 0)  then
         call echo(' setgam: no hole, setting gamach=0')
         return
      endif
      if (ihole .gt. 4)  then
         call echo(' Feff6L only does K and L shells')
         call fstop(' at SETGAM-2')
      endif

      zz = iz
      if (ihole .le. 1)  then
         do 10  i = 1, 6
            gamk(i) = log10 (gamk(i))
   10    continue
         call terp (zk, gamk, 6, zz, gamach)
      else if (ihole .le. 2)  then
         do 20  i = 1, 6
            gaml1(i) = log10 (gaml1(i))
   20    continue
         call terp (zl1, gaml1, 6, zz, gamach)
      else if (ihole .le. 4)  then
c        note that LII and LIII have almost exactly the same
c        core hole lifetimes
         do 30  i = 1, 6
            gaml2(i) = log10 (gaml2(i))
   30    continue
         call terp (zl2, gaml2, 6, zz, gamach)
      endif

c     Change from log10 (gamma) to gamma
      gamach = 10.0 ** gamach

c     Table values are in eV, code requires atomic units
      gamach = gamach / ryd

      return
      end
