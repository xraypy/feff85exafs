      subroutine mmtrxi (lam1x, mmati, ie, ileg, ilegp)
c     calculates matrix M in Rehr,Albers paper.
c     in polarization case
      implicit double precision (a-h, o-z)

c     all commons except for /fmat/ are inputs

c     inputs:
c       lam1x:  limits on lambda and lambda'
c       ie:  energy grid points
c       ileg, ilegp: leg and leg'
c
c     Inputs from common:
c        phases, use ph(ie,...,ilegp), and lmax(ie,ilegp)
c        lambda arrays
c        rotation matrix for ilegp
c        clmz for ileg and ilegp
c        path data, eta(ilegp) and ipot(ilegp)
c        xnlm array
c
c     Output:  fmati(...,ilegp) in common /fmatrx/ is set for
c              current energy point.

c     calculate scattering amplitude matrices
c     f(lam,lam') = sum_l tl gam(l,m,n)dri(l,m,m',ileg)gamt(l,m',n')
c                 *cexp(-i*m*eta),  eta = gamma+alpha'
c     lam lt lam1x, lam' lt lam2x such that m(lam) lt l0, n(lam) lt l0
c     gam = (-)**m c_l,n+m*xnlm, gamt = (2l+1)*c_ln/xnlm,
c     gamtl = gamt*tl

      include 'const.h'
      include 'dim.h'
      include 'nlm.h'
      include 'lambda.h'
      include 'clmz.h'
      include 'pola.h'
      include 'fmatrx.h'
      include 'rotmat.h'
      include 'pdata.h'

      complex*16 cam, camt, tltl,mmati
      dimension mmati(-mtot:mtot,-mtot:mtot)
      complex*16 gam(ltot+1,mtot+1,ntot+1),
     1           gamtl(ltot+1,mtot+1,ntot+1)

c     calculate factors gam and gamtl
      iln = il0
      ilx = il0
      do 30  il = iln, ilx
         tltl = 2*il - 1
         do 20  lam = 1, lam1x
            m = mlam(lam)
            if (m .lt. 0)  goto 20
            im = m+1
            if (im .gt. il)  goto 20
            in = nlam(lam) + 1
            imn = in + m
            if (lam .gt. lam1x)  goto 10
            cam = xnlm(il,im) * (-1)**m
            if (imn .le. il)  gam(il,im,in) = cam * clmi(il,imn,ileg)
            if (imn .gt. il)  gam(il,im,in) = 0
   10       if (lam .gt. lam1x) goto 20
            camt = tltl / xnlm(il,im)
            gamtl(il,im,in) = camt * clmi(il,in,ilegp)
   20    continue
   30 continue

      do 60 lam1 = 1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do 50  lam2 = 1, lam1x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = iabs(m2) + 1
            imn1 = iam1 + in1 - 1
            fmati(lam1,lam2,ilegp) = mmati(m1,m2)*
     1                       gam(il0,iam1,in1)*gamtl(il0,iam2,in2)
   50    continue
   60 continue

      return
      end
