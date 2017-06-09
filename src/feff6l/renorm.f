      subroutine renorm (dexv, vcoul, srho)

      implicit double precision (a-h,o-z)
      save

      common /print/ iprint
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)

c     vcoul is the coulomb potential (no factor of r**2) (output)
      dimension vcoul(251)
c     srho is charge density in form 4*pi*density*r**2 output)
      dimension srho(251)
c jm  9/23/87 added srho renormalized charge density to be used
c     in cphase
       external exchan
      do 10 i=1,np
         dv(i)=0.0
         d(i)=0.0
   10 continue
      ddjri=log(ws/dr(1))/dpas
      jri=1.0+ddjri
      jr1=jri
      ddjr1=ddjri-jr1+1.0

      if (jri-2*(jri/2).ne.0) go to 20
         jri=jri+1
   20 continue

      ddjri=ddjri-jri+1.0
c  ddjri = (log(ws)-dri)/dpas
c  dri  =  log(dr(jri))

      da=0.0
      do 30 j=1,norb
      do 30 i=1,np
   30    d(i)=d(i)+nel(j)*(dgc(i,j)**2+dpc(i,j)**2)

      do 50 i=jri,np
         dl=dr(i)
         if (i.eq.jri.or.i.eq.np) go to 40
            dl=dl+dl
            if ((i-2*(i/2)).eq.0) dl=dl+dl
   40    dd=d(i)*dl
         da=da+dd
   50 continue

      da=dpas*da/3.0
      dfo=dr(jri-1)*d(jri-1)
      df1=dr(jri)*d(jri)
      df2=dr(jri+1)*d(jri+1)
      dcor=-dpas*(df1*ddjri+(df2+dfo-2.0*df1)*ddjri**3/6.0+(df2-dfo)
     1 *ddjri**2*.25)
      da=da+dcor
      if (iprint .ge. 5)  write(16,60) da
   60 format (1h ,' no. of electrons outside the ws-radius',e16.8)
      db=0.0

      do 80 i=jri,np
         dl=1.0
         if (i.eq.jri.or.i.eq.np) go to 70
            dl=dl+dl
            if ((i-2*(i/2)).eq.0) dl=dl+dl
   70    dd=d(i)*dl
         db=db+dd
   80 continue

      db=dpas*db/3.0
      df0=d(jri-1)
      df1=d(jri)
      df2=d(jri+1)
      dcor=-dpas*(df1*ddjri+(df2+df0-2.0*df1)*ddjri**3/6.0+(df2-df0)
     1 *ddjri**2*.25)
      db=db+dcor
      if (iprint .ge. 5)  write(16,90) db
   90 format (1h ,' db= ',e16.8)

      call potslw (dvn,d,dp,dr,dpas,np)

      du=da*3.0/(ws**3)

      do 120 i=1,np
         if (i.gt.jr1+1) then
            srho(i)=0.0
            go to 100
         endif
            d(i)=d(i)+du*dr(i)**2
            srho(i)=d(i)
  100    continue
         dumm=-exchan(d(i),dr(i),dexv)/dr(i)
         dvf(i)=dumm
         if (i.gt.jr1) go to 110
            dvn(i)=dvn(i)-z/dr(i)+da*(1.50/ws-.50*dr(i)**2/ws**3)-db
            go to 120
  110    continue
            dvn(i)=0.0
  120 dv(i)=dvn(i)+dumm

c ad1 write the mt index and radius
      if (iprint .ge. 5)  write(16,55)jr1,dr(jr1)
  55  format(' jr1 = ',i10,10x,'wigner-seitz radius = ',e16.8)

c ad1 output 2.*dvn*r**2 for use in phase (dvn = normalised coulomb)
c     write(17,200)((2.0*dvn(i)*dr(i)*dr(i)),i=1,np)
c 200 format(1p5e16.8)
c      passvc formerly used to pass data directly to PHASE
c      do 151  i = 1, np
c         passvc (i) = 2.0 * dvn(i) * dr(i) * dr(i)
c  151 continue
c
c     passvc above is vcoul*r**2
      do 151  i = 1, np
         vcoul(i) = 2 * dvn(i)
  151 continue


c jm  output renormalized charge density for use in cphase
c                                          (d=4pi*rho*r^2)
c     write(18,200) srho

cjm write out rs as function of r
c     do 8934 i=1,jr1
c     xxrs=(3*dr(i)*dr(i)/srho(i))**.33333333
c8934 write(19,140) dr(i), xxrs
      return
      end
