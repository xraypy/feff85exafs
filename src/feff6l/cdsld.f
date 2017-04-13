      subroutine cdsld

      implicit double precision (a-h,o-z)
      save
      common /print/ iprint
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco

      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)

c titre = identification of the wave functions  s,p*,p,........
      character*40 ttl
      character*2  titre
      common /char/ titre(30), ttl

c  -- This read commented out to make input easier, not used for
c     PHASE calculations
      irm  = 0
      ins  = 0
      npun = 0
      nmfg = 0
      nmrk = 0
c     read (5,10) irm,ins,npun,nmfg,nmrk
   10 format (8i3)

c valeurs moyennes de r**j  if irm non-zero
c tabulation of the wave functions if ins non-zero
c the potential multiplied by r is perfore if npun non-zero
      if (irm.eq.0) go to 200
      if (iprint .ge. 5)  write(16,20) ttl
   20 format (1h1,40x,a40,/)
   30 read (5,10) j,l,n1,l1,j1,n2,l2,j2
      if (l.eq.0) go to 200

c valeur moyenne of (p1*p2+q1*q2)*r**j  if l positive
c valeur moyenne of (p1*q2+p2*q1)*r**j  if l negative
      if (n1.gt.0) go to 40
      if (((n1+1)*(n1+2)).ne.0) go to 60
      i1=1
      i2=1
      go to 80
   40 i1=0
      i2=0
      do 50 i=1,norb
      if (nqn(i).eq.n1.and.nql(i).eq.l1.and.(j1-1).eq.(-nk(i)/iabs(nk(i)
     1 ))) i1=i
      if (nqn(i).eq.n2.and.nql(i).eq.l2.and.(j2-1).eq.(-nk(i)/iabs(nk(i)
     1 ))) i2=i
   50 continue
      if (i1.ne.0.and.i2.ne.0) go to 80
   60 if (iprint .ge. 5)  write(16,70) j,l,n1,l1,j1,n2,l2,j2
   70 format (1h0,'    error for the card     ',8i3)
      go to 30
   80 dval=dfl(i1)+dfl(i2)
      if ((dval+j).gt.-1.0) go to 90
      if (n1) 170,170,60
   90 im=nmax(i1)
      if (nmax(i2).lt.im) im=nmax(i2)
      if (l.lt.0) go to 110
      do 100 i=1,im
      dv(i)=dgc(i,i1)*dgc(i,i2)
  100 dq(i)=dpc(i,i1)*dpc(i,i2)
      go to 130
  110 do 120 i=1,im
      dv(i)=dgc(i,i1)*dpc(i,i2)
  120 dq(i)=dgc(i,i2)*dpc(i,i1)
  130 call somm (dr,dv,dq,dpas,dval,j,im)
      if (l.lt.0) go to 150
      if (iprint .ge. 5)  write(16,140) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dval
  140 format (24x,'(p1p2+q1q2)r**',i2,' for  ',i1,a2,i3,a2,5x,'=',1pe14.
     1 7,/)
      go to 170
  150 if (iprint .ge. 5)  write(16,160) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dval
  160 format (24x,'(p1q2+q1p2)r**',i2,' for  ',i1,a2,i3,a2,5x,'=',1pe14.
     1 7,/)
  170 if (n1+1) 190,180,30
  180 i1=i1+1
      i2=i1
      if (i1-norb) 80,80,30
  190 i2=i2+1
      if (i2-norb) 80,80,180
  200 if (ins.eq.0) go to 260
      do 250 i=1,norb,3
      j=i+2
      if (j.gt.norb) j=norb
      im=0
      do 210 l=i,j
      if (nmax(l).gt.im) im=nmax(l)
  210 continue
      do 230 k=1,im
      if (((k-1)*(k-48*(k/48))).ne.0) go to 230
      if (iprint .ge. 5)  write(16,20) ttl
      if (iprint .ge. 5)  write(16,220) (nqn(l),titre(l),nqn(l),
     1                                     titre(l),l=i,j)
  220 format (9x,'r',14x,3(i1,a2,'g.c.',i11,a2,'p.c.',10x))
  230 if (iprint .ge. 5)  write(16,240) dr(k),
     1                                   (dgc(k,l),dpc(k,l),l=i,j)
  240 format (7(1pe17.7))
  250 continue
  260 if (npun.eq.0) go to 300
      do 270 i=1,np
  270 dp(i)=dvf(i)*dr(i)
c     write(8,280) ttl
  280 format (a40)
c     write(8,290) (dp(i),i=1,np)
  290 format (8f9.4)
  300 do 310 i=1,np
  310 d(i)=0.0
      nag=1
      if (nmfg.eq.0) go to 470
      if (iprint .ge. 5)  write(16,20)
      if (iprint .ge. 5)  write(16,320)
  320 format (/,30x,'integrales magnetiques directes et d echange'//)
  330 read (5,10) i1,i2,n1
      if (i1.le.0) go to 470
      if (i2.gt.0) go to 350
      if (((i2+1)*(i2+2)).ne.0) go to 340
      if (n1.le.0) n1=1
      i1=n1
      n1=i2
      i2=i1
      go to 360
  340 if (iprint .ge. 5)  write(16,70) i1,i2,n1
      go to 330
  350 if (i1.gt.norb.or.i2.gt.norb) go to 340
      n1=1
  360 j1=2*iabs(nk(i1))-1
      j2=2*iabs(nk(i2))-1
      kma=min0(j1,j2)
      nm=nmax(i2)
      do 380 j=1,kma,2
      call ykdir (i1,i1,j,nag)
      do 370 i=1,nm
  370 dp(i)=dq(i)*dgc(i,i2)*dpc(i,i2)
      dval=j+1
      call somm (dr,d,dp,dpas,dval,-1,nm)
  380 if (iprint .ge. 5)  write(16,390) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dval
  390 format (20x,'fm',i2,' (',i1,a2,',',i1,a2,') =',1pe14.7)
      if (i1.eq.i2) go to 440
      j1=(iabs(1-2*nk(i1))-1)/2
      j2=(iabs(1-2*nk(i2))-1)/2
      kma=max0(nql(i1)+j2,nql(i2)+j1)
      j1=iabs(nql(i2)-j1)
      j2=iabs(nql(i1)-j2)
      kmi=min0(j1,j2)
      j1=kmi+nql(i1)+nql(i2)
      j1=j1-2*(j1/2)
      if (j1.eq.0) kmi=kmi+1
      nm=min0(nmax(i1),nmax(i2))
      do 420 j=kmi,kma,2
      call ykdir (i1,i2,j,nag)
      do 400 i=1,nm
      dp(i)=dq(i)*dgc(i,i1)*dpc(i,i2)
  400 dc(i)=dq(i)*dgc(i,i2)*dpc(i,i1)
      dval=j+1
      dvalp=dval
      dvalm=dval
      call somm (dr,d,dp,dpas,dvalp,-1,nm)
      call somm (dr,d,dc,dpas,dval,-1,nm)
      call ykdir (i2,i1,j,nag)
      do 410 i=1,nm
  410 dp(i)=dq(i)*dgc(i,i2)*dpc(i,i1)
      call somm (dr,d,dp,dpas,dvalm,-1,nm)
  420 if (iprint .ge. 5)  write(16,430) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dvalm,dval,dvalp
  430 format (' gm',i2,' (',i1,a2,',',i1,a2,')',5x,'(-1)=',1pe14.7,5x,'(
     10)=',1pe14.7,5x,'(+1)=',1pe14.7)
  440 if (n1+1) 460,450,330
  450 i1=i1+1
      i2=i1
      if (i1-norb) 360,360,330
  460 i2=i2+1
      if (i2-norb) 360,360,450
  470 if (nmrk.eq.0) go to 530
      if (iprint .ge. 5)  write(16,20)
      if (iprint .ge. 5)  write(16,480)
  480 format (/,20x,'integrales magnetiques rk=integrale de p1(1)*q2(1)*
     1uk(1,2)*p3(2)*q4(2)'//)
  490 read (5,10) i1,i2,i3,i4,k
      if (i1.le.0) go to 530
      if (i1.le.norb.and.i2.gt.0.and.i2.le.norb.and.i3.gt.0.and.i3.le
     1 .norb.and.i4.gt.0.and.i4.le.norb.and.k.ge.0) go to 500
      if (iprint .ge. 5)  write(16,70) i1,i2,i3,i4,k
      go to 490
  500 call ykdir (i1,i2,k,nag)
      do 510 i=1,np
  510 dp(i)=dq(i)*dgc(i,i3)*dpc(i,i4)
      dval=k+1
      call somm (dr,d,dp,dpas,dval,-1,np)
      if (iprint .ge. 5)  write(16,520) k,nqn(i1),titre(i1),nqn(i2),
     1              titre(i2),nqn(i3),titre(i3),nqn(i4),titre(i4),dval
  520 format (20x,'rm',i2,' (',i1,a2,',',i1,a2,',',i1,a2,',',i1,a2,') ='
     1 ,1pe14.7)
      go to 490
  530 return
      end
