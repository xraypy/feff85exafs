      subroutine ykdir (ia,ib,nk1,nag)

      implicit double precision (a-h,o-z)
      save
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)
      common /trois/ dpno(4,30), dqno(4,30)
      dimension dpn1(4)
      dpah=exp(dpas)
      dpyk=dpas/24.0
      id=min0(nmax(ia)+2,nmax(ib)+2,np)
      idm1=id-1
      if (nag.ne.0) go to 30
      do 10 i=1,id
   10 dq(i)=dr(i)*(dgc(i,ia)*dgc(i,ib)+dpc(i,ia)*dpc(i,ib))
      do 20 i=1,4
      dpn1(i)=0.0
      do 20 j=1,i
   20 dpn1(i)=dpn1(i)+dpno(j,ia)*dpno(i+1-j,ib)+dqno(j,ia)*dqno(i+1-j,ib
     1 )
      go to 60
   30 do 40 i=1,id
   40 dq(i)=dr(i)*dgc(i,ia)*dpc(i,ib)
      do 50 i=1,4
      dpn1(i)=0.0
      do 50 j=1,i
   50 dpn1(i)=dpn1(i)+dpno(j,ia)*dqno(i+1-j,ib)
   60 di=dfl(ia)+dfl(ib)+nk1
      dp(1)=0.0
      dp(2)=0.0
      do 70 i=1,4
      di=di+1.0
      dp(1)=dp(1)+(dr(1)**di)*dpn1(i)/di
   70 dp(2)=dp(2)+(dr(2)**di)*dpn1(i)/di
      dm=dpah**(-nk1)
      dim2=-dpyk*dm*dm
      dim1=13.0*dpyk*dm
      di=13.0*dpyk
      dip1=-dpyk/dm
      do 80 i=3,idm1
   80 dp(i)=dp(i-1)*dm+dim2*dq(i-2)+dip1*dq(i+1)+dim1*dq(i-1)+di*dq(i)
      dq(id-2)=dp(id-2)
      do 90 i=idm1,np
   90 dq(i)=dq(i-1)*dm
      i=nk1+nk1+1
      dm=dm/dpah
      dim2=i*dim2/(dpah*dpah)
      dim1=i*dim1/dpah
      di=i*di
      dip1=i*dip1*dpah
      i=id-3
  100 dq(i)=dq(i+1)*dm+dim2*dp(i+2)+dip1*dp(i-1)+dim1*dp(i+1)+di*dp(i)
      i=i-1
      if (i-1) 110,110,100
  110 dq(1)=dq(3)*dm*dm+8.0*((di*dp(1)+4.0*dim1*dp(2))/13.0-dim2*dp(3
     1 ))
      return
      end
