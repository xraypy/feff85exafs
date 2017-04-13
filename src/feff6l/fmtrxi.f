      subroutine fmtrxi (lam1x, lam2x, ie, iterm, ileg, ilegp)
      implicit double precision (a-h, o-z)

c     all commons except for /fmat/ are inputs

c     inputs:
c       lam1x, lam2x:  limits on lambda and lambda'
c       ie:  energy grid points
c       iterm = 1 if we're doing the termination matrix M,
c              -1 otherwise
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
      include 'fmatrx.h'
      include 'rotmat.h'
      include 'pdata.h'

      complex*16 cam, camt, cterm, tltl
      complex*16 gam(ltot+1,mtot+1,ntot+1),
     1           gamtl(ltot+1,mtot+1,ntot+1), tl

c     calculate factors gam and gamtl
      iln = 1
      ilx = lmax(ie,ipot(ilegp)) + 1
      if (iterm .gt. 0)  then
         iln = il0
         ilx = il0
      endif
      do 30  il = iln, ilx
         tltl = 2*il - 1
         if (iterm .lt. 0)  then
            tl = (exp(2*coni*ph(ie,il,ipot(ilegp))) - 1) / (2*coni)
            tltl = tltl * tl
         endif
         lam12x = max (lam1x, lam2x)
         do 20  lam = 1, lam12x
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
   10       if (lam .gt. lam2x) goto 20
            camt = tltl / xnlm(il,im)
            gamtl(il,im,in) = camt * clmi(il,in,ilegp)
   20    continue
   30 continue

      do 60 lam1 = 1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do 60  lam2 = 1, lam2x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = iabs(m2) + 1
            imn1 = iam1 + in1 - 1
            cterm = 0
            ilmin = max (iam1, iam2, imn1, in2, iln)
            do 40  il = ilmin, ilx
c              skip terms with mu > l (NB il=l+1, so mu=il is mu>l)
               if (abs(m1).ge.il .or. abs(m2).ge.il)  goto 40
               m1d = m1 + mtot+1
               m2d = m2 + mtot+1

               cterm = cterm + gam(il,iam1,in1)*gamtl(il,iam2,in2)
     1                         *dri(il,m1d,m2d,ilegp)

   40       continue
            if (eta(ileg) .ne. 0.0) then
               m1 = mlam(lam1)
               cterm = cterm * exp(-coni*eta(ileg)*m1)
            endif
c           Above was org coding, change to use eta(ilegp) as test
c           based on algebra check.  July 20, 1992, siz&jjr
c           Changed back with redifinition of eta(see rdpath.f)
c           which is more convinient in polarization case.
c           August 8,1993, ala.
c           if (eta(ilegp) .ne. 0.0) then
c              m1 = mlam(lam1)
c              cterm = cterm * exp(-coni*eta(ilegp)*m1)
c           endif
            fmati(lam1,lam2,ilegp) = cterm
   60 continue

c     test of fmati(lam,lam',ileg)
c     plot fmat(lam,lam') = csqrt((z/2)**(m1-m2))*fmat

      return
      end
