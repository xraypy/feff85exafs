      subroutine inouh (dp,dq,dr,dq1,dfl,dv,z,test,nuc,nstop,jc)
c
c initial values for the outward integration
c dp=large component;     dq=small component;     dr=radial mesh
c dq1=slope at the origin of dp or dq;  dfl=power of the first term
c du=developpement limite;  dv=potential at the first point
c z=atomic number      test=test of the precision
c finite nuclear size if nuc is non-zero
c nstop controls the convergence  du developpement limite
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      common /ps1/ dep(5), deq(5), dd, dvc, dsal, dk, dm
c
c dep,deq=derivatives of dp and dq; dd=energy/dvc;
c dvc=speed of light in a.u.;
c dsal=2.*dvc   dk=kappa quantum number
c dm=exponential step/720.
c **********************************************************************
      common /trois/ dpno(4,30), dqno(4,30)
      dimension dp(251), dq(251), dr(251)
      do 10 i=1,10
      dp(i)=0.0
   10 dq(i)=0.0
      if (nuc) 20,20,60
   20 dval=z/dvc
      deva1=-dval
      deva2=dv/dvc+dval/dr(1)-dd
      deva3=0.0
      if (dk) 30,30,40
   30 dbe=(dk-dfl)/dval
      go to 50
   40 dbe=dval/(dk+dfl)
   50 dq(10)=dq1
      dp(10)=dbe*dq1
      go to 90
   60 dval=dv+z*(3.0-dr(1)*dr(1)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
      deva1=0.0
      deva2=(dval-3.0*z/(dr(nuc)+dr(nuc)))/dvc-dd
      deva3=z/(dr(nuc)*dr(nuc)*dr(nuc)*dsal)
      if (dk) 70,70,80
   70 dp(10)=dq1
      go to 90
   80 dq(10)=dq1
   90 do 100 i=1,5
      dp(i)=dp(10)
      dq(i)=dq(10)
      dep(i)=dp(i)*dfl
  100 deq(i)=dq(i)*dfl
      m=1
  110 dm=m+dfl
      dsum=dm*dm-dk*dk+deva1*deva1
      dqr=(dsal-deva2)*dq(m+9)-deva3*dq(m+7)
      dpr=deva2*dp(m+9)+deva3*dp(m+7)
      dval=((dm-dk)*dqr-deva1*dpr)/dsum
      dsum=((dm+dk)*dpr+deva1*dqr)/dsum
      j=-1
      do 130 i=1,5
      dpr=dr(i)**m
      dqr=dsum*dpr
      dpr=dval*dpr
      if (m.eq.1) go to 120
  120 dp(i)=dp(i)+dpr
      dq(i)=dq(i)+dqr
      if (abs(dpr/dp(i)).le.test.and.abs(dqr/dq(i)).le.test) j=1
      dep(i)=dep(i)+dpr*dm
  130 deq(i)=deq(i)+dqr*dm
      if (j.eq.1) go to 140
      dp(m+10)=dval
      dq(m+10)=dsum
      m=m+1
      if (m.le.20) go to 110
      nstop=45
  140 do 150 i=1,4
      dpno(i,jc)=dp(i+9)
  150 dqno(i,jc)=dq(i+9)
      return
      end
