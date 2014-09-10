      subroutine mmtrxi ( rkk, lam1x, bmati, ie, ileg, ilegp, lind,
     &       clmi, mlam, nlam, xnlm, eta, fmati)
c     calculates matrix M in Rehr,Albers paper.
c     in polarization case
      implicit double precision (a-h, o-z)

c     inputs:
c       lam1x:  limits on lambda and lambda'
c       ie:  energy grid points
c       ileg, ilegp: leg and leg'
c       phases, use ph(ie,...,ilegp), and lmax(ie,ilegp)
c       lambda arrays
c       rotation matrix for ilegp
c       clmz for ileg and ilegp
c       path data, eta(ilegp) and ipot(ilegp)
c       xnlm array
c
c     Output:  fmati(...,ilegp) is set for current energy point.

c     calculate scattering amplitude matrices
c     f(lam,lam') = sum_l tl gam(l,m,n)dri(l,m,m',ileg)gamt(l,m',n')
c                 *cexp(-i*m*eta),  eta = gamma+alpha'
c     lam lt lam1x, lam' lt lam2x such that m(lam) lt l0, n(lam) lt l0
c     gam = (-)**m c_l,n+m*xnlm, gamt = (2l+1)*c_ln/xnlm,
c     gamtl = gamt*tl

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      
c+----------------------------------------------------------------------
c     removing local common blocks, replacing them with explicit passing
c     of the defined data srtuctures
c+----------------------------------------------------------------------
c     include 'nlm.h'
      dimension xnlm(ltot+1,mtot+1)
c     include 'lambda.h'
      integer mlam(lamtot), nlam(lamtot)
c     include 'clmz.h'
      complex*16 clmi(ltot+1,mtot+ntot+1,legtot)
c     include 'fmatrx.h'
      complex*16 fmati(lamtot,lamtot,legtot)
c     include 'pdata.h'
      double precision eta(0:legtot+1)


      complex*16 cam, camt, tltl, bmati
      dimension bmati(-mtot:mtot, 1:8, -mtot:mtot, 1:8), lind(8)
      complex*16  rkk(nex,8)
      complex*16 gam(ltot+1,mtot+1,ntot+1),
     1           gamtl(ltot+1,mtot+1,ntot+1)

c     calculate factors gam and gamtl
 
c     set limits for orbital momentum
      lmn = ltot
      lmx = 0
      do 10 k1 = 1,8
        if (lind(k1).gt.lmx) lmx = lind(k1)
        if (lind(k1).lt.lmn .and. lind(k1).ge.0) lmn = lind(k1)
  10  continue
      iln = lmn + 1
      ilx = lmx + 1
 
      do 30  il = iln, ilx
         tltl = 2*il - 1
         do 20  lam = 1, lam1x
            m = mlam(lam)
            if (m .lt. 0)  goto 20
            im = m+1
            if (im .gt. il)  goto 20
            in = nlam(lam) + 1
            imn = in + m
            if (lam .gt. lam1x)  goto 20
            cam = xnlm(il,im) * (-1)**m
            if (imn .le. il)  gam(il,im,in) = cam * clmi(il,imn,ileg)
            if (imn .gt. il)  gam(il,im,in) = 0
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
            iam2 = abs(m2) + 1
            fmati(lam1,lam2,ilegp) = 0.0d0

            do 40 k1 = 1, 8
            do 40 k2 = 1, 8
               l1 = lind(k1) + 1
               l2 = lind(k2) + 1
               if (l1.gt.0.and.l2.gt.0 .and. iam1.le.l1.and.iam2.le.l2) 
     1           fmati(lam1,lam2,ilegp) = fmati(lam1,lam2,ilegp) +
     2           bmati(m1,k1, m2,k2) * rkk(ie,k1) * rkk(ie,k2) *
     3           gam( l1, iam1, in1) * gamtl( l2, iam2, in2)
   40       continue
            fmati(lam1,lam2,ilegp) = fmati(lam1,lam2,ilegp) *
     1      exp(-coni*eta(ileg)*m1)
   50    continue
   60 continue

      return
      end
