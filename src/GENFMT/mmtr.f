      subroutine mmtr( bmati, ipol, ispin, le2, angks, ptz, lind, dri,
     &       eta, nsc, nleg, kinit, ilinit)
c     calculates the part of matrix M which does not depend on energy
c     point.( see Rehr and Albers paper)
c     for path expansion always neglect spin-flip processes
c     to simplify calculations; (bmati does not have spin indices)

      implicit double precision (a-h, o-z)

c     Inputs:
c        kinit: quantum number kappa for initial orbital
c        rotation matrix for ilegp
c        path data, eta(ilegp) and ipot(ilegp)
c        mtot,l0
c        polarization data : ipol, ispin, ptz
c     Output:  bmati(...) 

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
c     include 'rotmat.h'
      dimension dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)
c     include 'pdata.h'
      double precision eta(0:legtot+1)
      integer nsc, nleg
      integer kinit, ilinit

      complex*16 ptz
      dimension ptz(-1:1, -1:1)

      complex*16 bmati, bmat
      dimension bmati(-mtot:mtot, 1:8, -mtot:mtot, 1:8)
      dimension bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
      dimension kind(8), lind(8)
      logical ltrace

      do 10 i = 1,8
      do 10 k = -mtot,mtot
      do 10 j = 1,8
      do 10 l = -mtot,mtot
         bmati(l,j,k,i)=0.0d0
  10  continue

      ltrace = .false.
      call bcoef (kinit, ipol, ptz, le2, ltrace, ispin, angks,
     1            kind, lind, bmat)
      is = 0
      if (ispin.eq.1) is = nspx - 1

c     ilinit = initial orb. momentum + 1. 
      lxx = min(mtot,ilinit)

c     set indices for bmati (no indentation)
      do 60 mu1 = -lxx,lxx
      mu1d = mu1+mtot+1
      do 50 mu2 = -lxx,lxx
      mu2d = mu2+mtot+1
      if (ipol.ne.0) then
        do 40 k1 = 1,8
        do 40 k2=  1,8
          l1 = lind(k1) + 1
          l2 = lind(k2) + 1
          do 35 m1 = -lind(k1), lind(k1)
          do 35 m2 = -lind(k2), lind(k2)
            m1d = m1 + mtot+1
            m2d = m2 + mtot+1
            bmati(mu1,k1, mu2,k2) =bmati(mu1,k1, mu2,k2) +
     1      bmat(m1,is,k1,m2,is,k2)*exp(-coni*(eta(nsc+2)*m2+eta(0)*m1))
     1        * dri(l1,mu1d,m1d,nsc+2) * dri(l2,m2d,mu2d,nleg)
c           dri(nsc+2)  is angle between z and leg1
c           dri(nsc+1)  is angle between last leg and z
c           eta(0)      is gamma between eps and rho1,
c           eta(nsc+2)  is alpha between last leg and eps
   35     continue
   40   continue
      else
c       ipol=0 and bmat is diagonal in k1,k2 and LS L'S'
c       and 2 rotation matrices can be combined to 1
        do 140 k1 = 1,8
          l1 = lind(k1) + 1
          if (l1.gt.0) then
            m1 = 0
            m1d = m1 + mtot+1
            bmati(mu1,k1, mu2,k1) =bmati(mu1,k1, mu2,k1) +
     1        bmat(m1,is,k1,m1,is,k1) * dri( l1, mu1d, mu2d, nsc+1)
c           dri(nsc+1)  is angle between last leg and first leg
          endif
  140   continue
      endif
   50 continue
   60 continue

      return
      end
