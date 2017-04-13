      subroutine mpprmd (npat, ipat, ri, beta, eta)
c     double precision version so angles come out right
c     for output...

c     Used with pathsd, a single precision code, so BE CAREFUL!!
c     No implicit, all variables declared explicitly.

c     make path parameters, ie, ri, beta, eta for each leg for a given
c     path.

c     Input is list of atoms (npat, ipat(npat)), output is
c     ri(npat+1), beta, eta.

      dimension ipat(npat)
      include 'dim.h'

c     /atoms/ is single precision from pathsd
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      complex*16  coni
      parameter (coni = (0,1))

      complex*16  alph(npatx+1), gamm(npatx+2), eieta
      double precision beta(npatx+1)
      double precision ri(npatx+1), eta(npatx+1)

      double precision x, y, z
      double precision ct, st, cp, sp, ctp, stp, cpp, spp
      double precision cppp, sppp

      n = npat + 1
      do 100  j = 1, n

c        get the atoms in this path
c        we actually have them already via the ipat array
c        remember that we'll want rat(,npat+1)=rat(,0) and
c                                 rat(,npat+2)=rat(,1) later on
c        make alpha, beta, and gamma for point i from 1 to N
c        NB: N is npat+1, since npat is number of bounces and N is
c            number of legs, or think of N=npat+1 as the central atom
c            that is the end of the path.
c
c        for euler angles at point i, need th and ph (theta and phi)
c        from rat(i+1)-rat(i)  and  thp and php
c        (theta prime and phi prime) from rat(i)-rat(i-1)
c
c        Actually, we need cos(th), sin(th), cos(phi), sin(phi) and
c        also for angles prime.  Call these  ct,  st,  cp,  sp   and
c                                            ctp, stp, cpp, spp.
c
c        We'll need angles from n-1 to n to 1,
c        so use rat(n+1) = rat(1), so we don't have to write code
c        later to handle these cases.

c        i = ipat(j)
c        ip1 = ipat(j+1)
c        im1 = ipat(j-1)
c        except for special cases...
         if (j .eq. n)  then
c           j central atom, j+1 first atom, j-1 last path atom
            i = 0
            ip1 = ipat(1)
            im1 = ipat(npat)
         elseif (j .eq. npat)  then
c           j last path atom, j+1 central, j-1 next-to last atom
c              unless only one atom, then j-1 central
            i = ipat(j)
            ip1 = 0
            if (npat .eq. 1)  then
               im1 = 0
            else
               im1 = ipat(npat-1)
            endif
         elseif (j .eq. 1)  then
c           j first atom, j+1 second unless only one,
c           then j+1 central, j-1 central
            i = ipat(j)
            if (npat .eq. 1)  then
               ip1 = 0
            else
               ip1 = ipat (j+1)
            endif
            im1 = 0
         else
            i = ipat(j)
            ip1 = ipat(j+1)
            im1 = ipat(j-1)
         endif

         x = rat(1,ip1) - rat(1,i)
         y = rat(2,ip1) - rat(2,i)
         z = rat(3,ip1) - rat(3,i)
         call strigd (x, y, z, ct, st, cp, sp)
         x = rat(1,i) - rat(1,im1)
         y = rat(2,i) - rat(2,im1)
         z = rat(3,i) - rat(3,im1)
         call strigd (x, y, z, ctp, stp, cpp, spp)

c        cppp = cos (phi prime - phi)
c        sppp = sin (phi prime - phi)
         cppp = cp*cpp + sp*spp
         sppp = spp*cp - cpp*sp

c        alph = exp**(i alpha)  in ref eqs 18
c        beta = cos(beta)
c        gamm = exp**(i gamma)
         alph(j) = st*ctp - ct*stp*cppp - coni*stp*sppp
         beta(j) = ct*ctp + st*stp*cppp
c        Watch out for roundoff errors
         if (beta(j) .lt. -1)  beta(j) = -1
         if (beta(j) .gt.  1)  beta(j) =  1
         gamm(j) = st*ctp*cppp - ct*stp + coni*st*sppp
         ri(j) = sdist (rat(1,i), rat(1,im1))
  100 continue

c     Make eta(i) = alpha(i) + gamma(i+1).  We only really need
c     exp(i*eta)=eieta, so that's what we'll calculate.
c     We'll need gamm(N+1)=gamm(npat+2)=gamm(1)
      gamm(npat+2) = gamm(1)
      do 150  j = 1, npat+1
         eieta = alph(j) * gamm(j+1)
         call sargd (eieta, eta(j))
  150 continue

c     Return beta as an angle, ie, acos(beta).  Check for beta >1 or
c     beta <1 (roundoff nasties)
      do 160  j = 1, npat+1
         if (beta(j) .gt.  1)  beta(j) =  1
         if (beta(j) .lt. -1)  beta(j) = -1
         beta(j) = acos(beta(j))
  160 continue

      return
      end
      subroutine strigd (x, y, z, ct, st, cp, sp)
      double precision x, y, z, ct, st, cp, sp, r, rxy
c     returns cos(theta), sin(theta), cos(phi), sin(ph) for (x,y,z)
c     convention - if x=y=0, phi=0, cp=1, sp=0
c                - if x=y=z=0, theta=0, ct=1, st=0
      parameter (eps = 1.0e-6)
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
         sp = 0
      else
         cp = x / rxy
         sp = y / rxy
      endif

      return
      end
      subroutine sargd (c, th)

      double precision x, y, th
      complex*16  c
      parameter (eps = 1.0e-6)
      x = dble(c)
      y = dimag(c)
      if (abs(x) .lt. eps)  x = 0
      if (abs(y) .lt. eps)  y = 0
      if (abs(x) .lt. eps  .and.  abs(y) .lt. eps)  then
         th = 0
      else
         th = atan2 (y, x)
      endif
      return
      end
