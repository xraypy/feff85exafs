      subroutine totale (dval)
      implicit double precision (a-h,o-z)
      save
      common /print/ iprint
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)
       external exchee
      dc(1)=1
      do 10 i=1,np
   10 dp(i)=d(i)/dr(i)
      if (nuc.le.0) go to 30
      do 20 i=1,nuc
   20 dp(i)=d(i)*(3.0-dr(i)*dr(i)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
      dc(1)=4
   30 call somm (dr,dp,dq,dpas,dc(1),0,np)
      dc(1)=-z*dc(1)
      do 40 i=1,np
      dp(i)=d(i)*dvf(i)
      dvn(i)=d(i)*dvn(i)
   40 d(i)=d(i)*exchee(d(i),dr(i))
      dc(2)=2
      dc(3)=1
      dc(5)=2
      if (nuc.ne.0) dc(3)=4
      call somm (dr,dp,dq,dpas,dc(3),0,np)
      call somm (dr,dvn,dq,dpas,dc(5),0,np)
      call somm (dr,d,dq,dpas,dc(2),0,np)
      dc(4)=dval-dc(3)
      dval=dval-.50*dc(5)-dc(2)
      dc(2)=dc(3)-dc(1)-dc(5)-dc(2)
      dc(3)=.50*dc(5)
      if (iprint .ge. 5)  write(16,50) dval,dc(4),dc(3),dc(2),dc(1)
   50 format (1h0,5x,'et=',1pe14.7,5x,'ec=',1pe14.7,5x,'ee=',1pe14.7,5x,
     1 'ex=',1pe14.7,5x,'en=',1pe14.7)
      return
      end
