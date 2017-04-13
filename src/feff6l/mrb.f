      subroutine mrb (npat, ipat, ri, beta)

c     Make ri, beta and rpath path parameters for crit calculations.

c     Input is list of atoms (npat, ipat(npat)), output is
c     ri(npat+1), beta, eta.

      include 'dim.h'
      dimension ipat(npatx)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension beta(npatx+1), ri(npatx+1), ipat0(npatx+1)

      nleg = npat+1
c     central atom is atom 0 in rat array
c     need local ipat0 array since we use ipat0(npat+1), final atom
c     in path (final atom is, of course, the central atom)
      do 10  i = 1, npat
         ipat0(i) = ipat(i)
   10 continue
      ipat0(nleg) = 0

      do 30  ileg = 1, nleg
c        make beta and ri for point i from 1 to N
c        NB: N is npat+1, since npat is number of bounces and N is
c            number of legs, or think of N=npat+1 as the central atom
c            that is the end of the path.
c
c        We'll need angles from n-1 to n to 1,
c        so use rat(n+1) = rat(1), so we don't have to write code
c        later to handle these cases.

c        Work with atom j
c        jp1 = (j+1)
c        jm1 = (j-1)
         j = ileg
         jm1 = j-1
         jp1 = j+1
c        Fix special cases (wrap around when j is near central atom,
c        also handle ss and triangular cases).
         if (jm1 .le.    0)  jm1 = nleg
         if (jp1 .gt. nleg)  jp1 = 1

         jat = ipat0(j)
         jm1at = ipat0(jm1)
         jp1at = ipat0(jp1)

         ri(ileg) = sdist (rat(1,jat), rat(1,jm1at))

c        Make cos(beta) from dot product
         call dotcos (rat(1,jm1at), rat(1,jat), rat(1,jp1at),
     1               beta(ileg))
   30 continue

      rpath = 0
      do 60  ileg = 1, nleg
         rpath = rpath + ri(ileg)
   60 continue

      return
      end
      subroutine dotcos (rm1, r, rp1, cosb)
      dimension rm1(3), r(3), rp1(3)

      parameter (eps = 1.0e-8)

      cosb = 0
      do 100  i = 1, 3
         cosb = cosb + (r(i)-rm1(i)) * (rp1(i)-r(i))
  100 continue

c     if denom is zero (and it may be if 2 atoms are in the same place,
c     which will happen when last path atom is central atom), set
c     cosb = 0, so it won't be undefined.

      denom = (sdist(r,rm1) * sdist(rp1,r))
      if (denom .gt. eps)  then
         cosb = cosb / denom
      else
         cosb = 0
      endif
      return
      end
