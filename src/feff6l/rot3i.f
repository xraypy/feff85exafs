      subroutine rot3i (lxp1, mxp1, ileg)
      implicit double precision (a-h,o-z)

c     input:  lxp1, mxp1, ileg (lmax+1, mmax+1)
c             also beta(ileg) used from common /pdata/
c     output: dri(...ileg) in common /rotmat/

c     subroutine rot3 calculates rotation matrices for l = 0,lxp1-1

c     subroutine rot3 calculates the beta dependence of rotation
c     matrix elements using recursion of an iterated version of
c     formula (4.4.1) in edmonds.
c
c     first written:(september 17,1986) by j. mustre
c     version 2  (17 sep 86)
c     version 3  (22 feb 87) modified by j. rehr
c     version for genfmt, modified by s. zabinsky, Sept 1991
c     Initialized dri0.  Some elements may be used before being
c        initialized elsewhere -- rot3i needs to be carefully
c        checked.  S. Zabinsky, April 1993
c
c******************** warning******************************************
c     ltot must be at least lxp1 or overwriting will occur
c     nmax must be at least nm or overwriting will occur
c----------------------------------------------------------------------
c     notation dri0(l,m,n) =  drot_i(l'm'n')
c     l = l'+1, n' = n-l, m' = m-l, primes denoting subscripts
c     thus dri0(1,1,1) corresponds to the rotation matrix with
c     l' = 0, and n' and m' = 0; dri0(3,5,5) : l' = 2,n' = 2,m' = 2.
c--------------------------------------------------------------------

      include 'dim.h'
      include 'rotmat.h'
      include 'pdata.h'
c     dri0 is larger than needed for genfmt, but necessary for
c     this calculation algorithm.  Copy result into smaller
c     dri arrays (in common) at end of this routine.
      dimension  dri0 (ltot+1, 2*ltot+1, 2*ltot+1)

c     initialize dri0
      do 200 il = 1, ltot+1
         do 200 im = 1, 2*ltot+1
            do 200 in = 1, 2*ltot+1
               dri0(il,im,in) = 0
  200 continue

      nm = mxp1
      ndm = lxp1+nm-1
      xc = cos(beta(ileg)/2)
      xs = sin(beta(ileg)/2)
      s = sin(beta(ileg))
      dri0(1,1,1) = 1
      dri0(2,1,1) = xc**2
      dri0(2,1,2) = s/sqrt(2.0d0)
      dri0(2,1,3) = xs**2
      dri0(2,2,1) = -dri0(2,1,2)
      dri0(2,2,2) = cos(beta(ileg))
      dri0(2,2,3) = dri0(2,1,2)
      dri0(2,3,1) = dri0(2,1,3)
      dri0(2,3,2) = -dri0(2,2,3)
      dri0(2,3,3) = dri0(2,1,1)
      do 30  l = 3, lxp1
         ln = 2*l - 1
         lm = 2*l - 3
         if (ln .gt. ndm)  ln = ndm
         if (lm .gt. ndm)  lm = ndm
         do 20  n = 1, ln
            do 10  m = 1, lm
               t1 = (2*l-1-n) * (2*l-2-n)
               t = (2*l-1-m) * (2*l-2-m)
               f1 = sqrt (t1/t)
               f2 = sqrt ((2*l-1-n) * (n-1) / t)
               t3 = (n-2) * (n-1)
               f3 = sqrt(t3/t)
               dlnm = f1 * xc**2 * dri0(l-1,n,m)
               if (n-1 .gt. 0) dlnm = dlnm - f2*s*dri0(l-1,n-1,m)
               if (n-2 .gt. 0) dlnm = dlnm + f3*xs**2*dri0(l-1,n-2,m)
               dri0(l,n,m) = dlnm
               if (n .gt. (2*l-3))
     1            dri0(l,m,n) = (-1)**(n-m) * dri0(l,n,m)
   10       continue
            if (n .gt. (2*l-3)) then
               dri0(l,2*l-2,2*l-2) = dri0(l,2,2)
               dri0(l,2*l-1,2*l-2) = -dri0(l,1,2)
               dri0(l,2*l-2,2*l-1) = -dri0(l,2,1)
               dri0(l,2*l-1,2*l-1) = dri0(l,1,1)
            endif
   20    continue
   30 continue
   40 continue

c-----test sum rule on d
c     open (19,file='rotmat.dat',status='new',carriagecontrol='list')
c     write(19,*)  ' l, m, sum'
c     write(19,*) ' (dri0(il,im,in),in = 1,ln)'
c     do 70 il = 1,lxp1
c        l = il-1
c        ln = 2*l+1
c        if(ln.gt.ndm) ln = ndm
c        do 37 im = 1,ln
c           sum = 0
c           do 50 in = 1,ln
c              m = im-il
c              term = dri0(il,im,in)
c  50       sum = sum+term**2
c           write(19,60) l,m,sum
c           write(19,62) (dri0(il,im,in),in = 1,ln)
c  60       format(2i3,e30.20)
c  62       format(5e14.6)
c  70 continue
c     close(19)
c-----end test------------------------

c     Copy result into dri(...ileg) in /rotmat/ (zero it first...)
      do 90  il = 1, ltot+1
         do 90  m1 = 1, 2*mtot+1
            do 90  m2 = 1, 2*mtot+1
               dri(il,m1,m2,ileg) = 0
   90 continue

      do 120  il = 1, lxp1
         mx = min (il-1, mxp1-1)
         do 110  m1 = -mx, mx
            do 100  m2 = -mx, mx
               dri(il,m1+mtot+1,m2+mtot+1,ileg)=dri0(il,m1+il,m2+il)
  100       continue
  110    continue
  120 continue

      return
      end
