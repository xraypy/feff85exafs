      subroutine fovrg (il, ihard, rmt, xmt, jri, e, nr, dx, ri, v, dny,
     1                  pu, qu, p, q, ps, qs, vm)
      implicit double precision (a-h, o-z)

c     Input:
c        il      ang mom number + 1
c        ihard   number of times convergence test fails
c        rmt     muffin tin radius
c        xmt     x such that rmt = exp ((x-1)*dx - 8.8)
c        jri     first interstitial grid point (imt + 1)
c        e       current complex energy
c        nr      number of points in r grid
c        dx      dx in Loucks' grid (usually .05)
c        ri(nr)  Loucks' position grid, r = exp ((i-1)*dx - 8.8)
c        v(nr)   total complex potential including energy dep xc
c                v is in the form  pot*r**2
c
c     Work space:
c        complex*16 p(nr), q(nr), ps(nr), qs(nr), vm(nr)
c        Must be dimensioned in calling program.  Coded like this
c        to make using different r-grids with different nrmax easy.
c
c     Output:
c        ihard   incremented each time convergence test fails
c        dny     r*g'/g, see loucks (4-85), q/p = cf/g (eq 4-86)
c        pu, qu  upper and lower components at muffin tin
c        q and q arrays  upper and lower components (see comments)

      complex*16 v(nr), e
      dimension ri(nr)
      complex*16 dny, pu, qu
      complex*16 p(nr), q(nr), ps(nr), qs(nr), vm(nr)

      include 'const.h'
      parameter (c = clight)
      parameter (csq = c**2)

      double precision lp1, ldcsq
      complex*16 c1,c2,c3,pc,qc,dp1,dq1,dp2,dq2,dp3,dq3,dp4,dq4
      complex*16 vh,vmh,vmnp1,psn,qsn,psnm1,qsnm1,psnm2,qsnm2
      complex*16 psnm3,qsnm3,psnm4,qsnm4,pp,qp,psnp1,qsnp1,prel,qrel
      complex*16 psu,vu,dummy
      complex*16 vn,vmn

c     test=1.e+04 value in loucks
      test=1.e+05
      nrk=6

      expdxh=exp(dx/2.0)
      dxd4=dx/4.0
      dxd8=dx/8.0
      a1=dx*3.30
      a2=-dx*4.20
      a3=dx*7.80
      a4=dx*14.0/45.0
      a5=dx*64.0/45.0
      a6=dx*24.0/45.0
      call diff (v,dx,jri,vm)
      twoz=-dble (v(1))/ri(1)
      l=il-1
      lp1=l+1.0
      ldcsq=l/csq
      ie=1
      r=ri(1)
      vn=v(1)
      vmn=vm(1)
cv    p(1)=1.0
      p(1)=1.e-20
      q(1)=-e/(2.0*l+3.0)*r*p(1)
      beta=lp1
      if (twoz.eq.0.0) go to 10
      beta=sqrt(lp1*l+1.0-(twoz/c)**2)
      sb0=(beta-lp1)*csq/twoz
      sa1=(3.0*beta-(twoz/c)**2)/(2.0*beta+1.0)
      sb1=csq/twoz*((beta-l)*sa1-1.0)-sb0
      sa2=((beta+3.0*lp1)*sa1-3.0*l+twoz/csq*(beta+lp1+3.0)*sb1)/
     1 (beta+1.0)/4.0
      sb2=(csq/twoz*(2.0*l*(beta+2.0-lp1)-l-(twoz/c)**2)*sa1-3.0*l
     1 *csq/twoz*(beta+2.0-lp1)+(beta+3.0-2.0*lp1-(twoz/c)**2)*sb1)/
     2 (beta+1.0)/4.0
      delta=r*csq/twoz
      q(1)=(sb0+delta*(sb1+delta*sb2))/(1.0+delta*(sa1+delta*sa2))*p(1)
   10 continue
c     runge kutta method  (see loucks)
      c1=vn/r**2-e
      c2=1.0-c1/csq
      c3=(vmn-2.0*vn)/c2/c2*ldcsq
      ps(1)=r*c2*q(1)+lp1*p(1)
      qs(1)=-lp1*q(1)+(r*c1-c3/r**3)*p(1)
      n=1
   20 continue
      pc=p(n)
      qc=q(n)
      dp1=dx*(r*c2*qc+lp1*pc)
      dq1=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      pc=pc+0.50*dp1
      qc=qc+0.50*dq1
      r=r*expdxh
      vnp1=v(n+1)
      vmnp1=vm(n+1)
      vh=(vn+vnp1)*.50+(vmn-vmnp1)*dxd8
      vmh=(1.50*(vnp1-vn)-(vmn+vmnp1)*dxd4)/dx
      c1=vh/r/r-e
      c2=1.0-c1/csq
      c3=(vmh-2.0*vh)/c2/c2*ldcsq
      dp2=dx*(r*c2*qc+lp1*pc)
      dq2=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      pc=pc+0.50*(dp2-dp1)
      qc=qc+0.50*(dq2-dq1)
      dp3=dx*(r*c2*qc+lp1*pc)
      dq3=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      pc=pc+dp3-0.50*dp2
      qc=qc+dq3-0.50*dq2
      n=n+1
      r=ri(n)
      c1=vnp1/r/r-e
      c2=1.0-c1/csq
      c3=(vmnp1-2.0*vnp1)/c2/c2*ldcsq
      dp4=dx*(r*c2*qc+lp1*pc)
      dq4=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      p(n)=p(n-1)+(dp1+2.0*(dp2+dp3)+dp4)/6.0
      q(n)=q(n-1)+(dq1+2.0*(dq2+dq3)+dq4)/6.0
      ps(n)=r*c2*q(n)+lp1*p(n)
      qs(n)=-lp1*q(n)+(r*c1-c3/r**3)*p(n)
      vn=vnp1
      vmn=vmnp1
      if (n-nrk) 20,30,30
   30 if (n.ge.jri) go to 120
      psn=ps(nrk)
      qsn=qs(nrk)
      psnm1=ps(nrk-1)
      qsnm1=qs(nrk-1)
      psnm2=ps(nrk-2)
      qsnm2=qs(nrk-2)
      psnm3=ps(nrk-3)
      qsnm3=qs(nrk-3)
      psnm4=ps(nrk-4)
      qsnm4=qs(nrk-4)
c     milne method
   40 r=ri(n+1)
      c1=v(n+1)/r/r-e
      c2=1.0-c1/csq
      c3=(vm(n+1)-2.0*v(n+1))/c2/c2*ldcsq
      pp=p(n-5)+a1*(psn+psnm4)+a2*(psnm1+psnm3)+a3*psnm2
      qp=q(n-5)+a1*(qsn+qsnm4)+a2*(qsnm1+qsnm3)+a3*qsnm2
      nit=0
   50 psnp1=r*c2*qp+lp1*pp
      qsnp1=-lp1*qp+(r*c1-c3/r**3)*pp
      pc=p(n-3)+a4*(psnp1+psnm3)+a5*(psn+psnm2)+a6*psnm1
      qc=q(n-3)+a4*(qsnp1+qsnm3)+a5*(qsn+qsnm2)+a6*qsnm1
      if (abs(test*(pc-pp))-abs(pc)) 60,60,70
   60 if (abs(test*(qc-qp))-abs(qc)) 110,110,70
   70 if (nit-40) 100,80,100
c  70 if (nit-5) 100,80,100 value in loucks
   80 prel=(pc-pp)/pc
      qrel=(qc-qp)/qc
c     count times hard test fails
      ihard = ihard + 1
c     print90, il,ie,n,prel,qrel
   90 format (' hard test in fovrg il=',i2,' ie=',i1,' n=',i3,' prel='
     1 ,e16.8,' qrel=',e16.8,' **********')
      go to 110
  100 nit=nit+1
      pp=pc
      qp=qc
      go to 50
  110 n=n+1
      p(n)=pc
      q(n)=qc
      ps(n)=psnp1
      qs(n)=qsnp1
      psnm4=psnm3
      psnm3=psnm2
      psnm2=psnm1
      psnm1=psn
      psn=psnp1
      qsnm4=qsnm3
      qsnm3=qsnm2
      qsnm2=qsnm1
      qsnm1=qsn
      qsn=qsnp1
c     introduce scale factor to prevent overflow on vax jjr
      if(abs(pc).lt.1.e+20) go to 119
      scale=1.e-20
      do 112 mm=1,6
      nm=n-mm+1
      p(nm)=scale*p(nm)
      q(nm)=scale*q(nm)
      ps(nm)=scale*ps(nm)
      qs(nm)=scale*qs(nm)
  112 continue
      psnm4=scale*psnm4
      psnm3=scale*psnm3
      psnm2=scale*psnm2
      psnm1=scale*psnm1
      psn=scale*psn
      qsnm4=scale*qsnm4
      qsnm3=scale*qsnm3
      qsnm2=scale*qsnm2
      qsnm1=scale*qsnm1
      qsn=scale*qsn
  119 if (n-jri) 40,120,120
  120 jm=jri-1
      x=dx*(xmt-jm)
      call intpol (zero,dx,p(jm),p(jri),ps(jm),ps(jri),x,pu,psu)
      call intpol (zero,dx,q(jm),q(jri),qs(jm),qs(jri),x,qu,dummy)
      call intpol (zero,dx,v(jm),v(jri),vm(jm),vm(jri),x,vu,dummy)
      dny=rmt*(1.0-(vu/rmt**2-e)/csq)*qu/pu+l
c dny is r*g'/g, see loucks (4-85), q/p = cf/g (eq 4-86)
c (watch for factors of rmt)
      return
      end
