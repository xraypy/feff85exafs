      subroutine dirac (nqn,nql,nk,imax,de,dfl,dq1,jc)
c
c solution of the dirac equation
c nqn=principal quantum number; nql=orbital quantum number
c nk=kappa quantum number;  imax=the last tabulated point of the
c wave function; de=energy;   dfl=power of the first term of the
c developpement limite; dq1=slope at the origin of dp or dq
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, test,
     1              z, nstop, nes, np, nuc
c
c dv=potential in a.u. and negative;  dr=radial mesh
c dp=large component;    dq=small component;    dpas=exponential step;
c nes=number of attempts to adjust the energy
c z=atomic number; nstop controls the numeric integration
c test=precision obtained in the energies; np=maximum number of points
c finite nuclear size if nuc is non-zero
c **********************************************************************
      common /ps1/ dep(5), deq(5), db, dvc, dsal, dk, dm
c
c dep,deq=derivatives of op and dq;  db=energie/dvc;
c dvc=speed of light in a.u.; dsal=2.*dvc;  dk=kappa quantum number
c dm=exponential step/720., dkoef=1./720.
c **********************************************************************
      common /trois/ dpno(4,30), dqno(4,30)
      data dkoef /.1388888888888888e-2/
      nstop=0
      dvc=137.0373
      dsal=dvc+dvc
      imm=0
      ies=0
      dk=nk
      lll=(nql*(nql+1))/2
      nd=0
      noeud=nqn-nql
      if (lll.ne.0) go to 10
      elim=-z*z/(1.5*nqn*nqn)
      go to 40
   10 elim=dv(1)+lll/(dr(1)*dr(1))
      do 20 i=2,np
      val=dv(i)+lll/(dr(i)*dr(i))
      if (val.le.elim) elim=val
   20 continue
      if (elim) 40,30,30
   30 nstop=17
c 2*v+l*(l+1)/r**2 is everywhere positive
c **********************************************************************
      return
   40 if (de.le.elim) de=elim*0.5
   50 if (imm.eq.1) go to 80
      do 60 i=7,np,2
      imat=np+1-i
      if ((dv(imat)+lll/(dr(imat)*dr(imat))-de).le.0.0) go to 70
   60 continue
   70 if (imat.gt.5) go to 80
      de=de*0.5
      if (de.lt.-test.and.nd.le.noeud) go to 50
      nstop=28
c 2*v+l*(l+1)/r**2-2*e is everywhere positive
c **********************************************************************
      return
c initial value for the outward integration
c **********************************************************************
   80 db=de/dvc
      call inouh (dp,dq,dr,dq1,dfl,dv(1),z,test,nuc,nstop,jc)
      if (nstop) 310,90,310
c     nstop=45
c the expansion at the origin does not converge
c **********************************************************************
   90 nd=1
      do 110 i=1,5
      dval=dr(i)**dfl
      if (i.eq.1) go to 100
      if (dp(i-1).eq.0.0) go to 100
      if ((dp(i)/dp(i-1)).gt.0.0) go to 100
      nd=nd+1
  100 dp(i)=dp(i)*dval
      dq(i)=dq(i)*dval
      dep(i)=dep(i)*dval
  110 deq(i)=deq(i)*dval
      k=-1+2*(noeud-2*(noeud/2))
      if ((dp(1)*k).gt.0.0) go to 130
  120 nstop=53
c error in the expansion at the origin
c **********************************************************************
      return
  130 if ((k*nk*dq(1)).lt.0.0) go to 120
      dm=dpas*dkoef
c outward integration
c **********************************************************************
      do 140 i=6,imat
      dp(i)=dp(i-1)
      dq(i)=dq(i-1)
      call inth (dp(i),dq(i),dv(i),dr(i))
      if (dp(i-1).eq.0.0) go to 140
      if ((dp(i)/dp(i-1)).gt.0.0) go to 140
      nd=nd+1
      if (nd.gt.noeud) go to 150
  140 continue
      if (nd.eq.noeud) go to 160
      de=0.8*de
      if (de.lt.-test) go to 50
      nstop=206
c the number of nodes is too small
c **********************************************************************
      return
  150 de=1.2*de
      if (de.gt.elim) go to 50
      nstop=210
c the number of nodes is too big
c **********************************************************************
      return
c initial values for the inward integration
c **********************************************************************
  160 dqm=dq(imat)
      dpm=dp(imat)
      if (imm.eq.1) go to 180
      do 170 i=1,np,2
      imax=np+1-i
      if(((dv(imax)-de)*dr(imax)*dr(imax)).le.300.0) go to 180
  170 continue
  180 dd=sqrt(-de*(2.0+db/dvc))
      dpq=-dd/(dsal+db)
      dm=-dm
      do 190 i=1,5
      j=imax+1-i
      dp(j)=exp(-dd*dr(j))
      dep(i)=-dd*dp(j)*dr(j)
      dq(j)=dpq*dp(j)
  190 deq(i)=dpq*dep(i)
      m=imax-5
c inward integration
c***********************************************************************
      do 200 i=imat,m
      j=m+imat-i
      dp(j)=dp(j+1)
      dq(j)=dq(j+1)
  200 call inth (dp(j),dq(j),dv(j),dr(j))
c joining of the large components
c **********************************************************************
      dval=dpm/dp(imat)
      if (dval.gt.0.0) go to 210
      nstop=312
c error in the sign of the large component
c **********************************************************************
      return
  210 do 220 i=imat,imax
      dp(i)=dp(i)*dval
  220 dq(i)=dq(i)*dval
c calculation of the norm
c **********************************************************************
      dsum=3.0*dr(1)*(dp(1)**2+dq(1)**2)/(dpas*(dfl+dfl+1.0))
      do 230 i=3,imax,2
  230 dsum=dsum+dr(i)*(dp(i)**2+dq(i)**2)+4.0*dr(i-1)*(dp(i-1)**2+dq(i-
     1 1)**2)+dr(i-2)*(dp(i-2)**2+dq(i-2)**2)
      dsum=dpas*(dsum+dr(imat)*(dqm*dqm-dq(imat)*dq(imat)))*0.3333333333
     1 333333
c modification of the energy
c **********************************************************************
      dbe=dp(imat)*(dqm-dq(imat))*dvc/dsum
      imm=0
      val=abs(dbe/de)
      if (val.le.test) go to 260
  240 dval=de+dbe
      if (dval.lt.0.0) go to 250
      dbe=dbe*0.5
      val=val*0.5
      if (val.gt.test) go to 240
      nstop=345
c energie nulle
c **********************************************************************
      return
  250 de=dval
      if (val.le.0.1) imm=1
      ies=ies+1
      if (ies.le.nes) go to 50
      nstop=362
c number of iterations too big
c **********************************************************************
      return
  260 dsum=sqrt(dsum)
      dq1=dq1/dsum
      do 270 i=1,imax
      dp(i)=dp(i)/dsum
  270 dq(i)=dq(i)/dsum
      do 280 i=1,4
      dpno(i,jc)=dpno(i,jc)/dsum
  280 dqno(i,jc)=dqno(i,jc)/dsum
      if (imax.eq.np) go to 300
      j=imax+1
      do 290 i=j,np
      dp(i)=0.0
  290 dq(i)=0.0
  300 nstop=0
  310 return
      end
