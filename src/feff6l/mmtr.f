      subroutine mmtr(t3j,mmati)
c     calculates the part of matrix M which does not depend on energy
c     point.( see Rehr and Albers paper)

      implicit double precision (a-h, o-z)

c     all commons are inputs
c     inputs:
c        t3j: appropriate table of the 3j symbols
c     Inputs from common:
c        rotation matrix for ilegp
c        path data, eta(ilegp) and ipot(ilegp)
c        mtot,l0
c     Output:  mmati(...) 

      include 'const.h'
      include 'dim.h'
      include 'pola.h'
      include 'rotmat.h'
      include 'pdata.h'

      complex*16 mmati
      dimension mmati(-mtot:mtot,-mtot:mtot),t3j(-mtot-1:mtot+1,-1:1)

      do 10 i = -mtot,mtot
      do 10 j = -mtot,mtot
         mmati(i,j)=0
  10  continue
      li = l0-1
c     l0 is final orb. momentum. Thus here we need to change code
c     in case when initial momemtum larger than final one.
      lx = min(mtot,l0)

      do 60 mu1 = -lx,lx
         mu1d = mu1+mtot+1
         do 50 mu2 = -lx,lx
            mu2d = mu2+mtot+1
            do 35  m0 = -li,li 
               do 34 i = -1,1
               do 34 j = -1,1
                  m1 = m0-j
                  m2 = m0-i
                  m1d = m1 + mtot+1
                  m2d = m2 + mtot+1
                  if (abs(m1).gt.lx .or. abs(m2).gt.lx)  goto 34
                  mmati(mu1,mu2) = mmati(mu1,mu2) + 
     1              dri(il0,mu1d,m1d,nsc+2)*dri(il0,m2d,mu2d,nleg)
     2              *exp(-coni*(eta(nsc+2)*m2+eta(0)*m1))
     3              *t3j(-m0,i)*t3j(-m0,j)*ptz(i,j)

c           dri(nsc+2)  is angle between z and leg1
c           dri(nsc+1)  is angle between last leg and z
c           eta(nsc+3)  is gamma between eps and rho1,
c           eta(nsc+2)  is alpha between last leg and eps
c           t3j(m0,i)    are 3j symbols multiplied by sqrt(3) 
   34          continue
   35       continue
            mmati(mu1,mu2) = mmati(mu1,mu2)*exp(-coni*eta(1)*mu1)
   50    continue
   60  continue

      return
      end
