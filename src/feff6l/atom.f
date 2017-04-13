c     code: relativistic atom code (relativistic hartree fock slater)
c     modified desclaux code -- partially translated from the french
c
c     modified by: r. c. albers (from previously modified code from
c                  j. e. muller who in turn got it from the danes)
c                  j. j. rehr and s. i. zabinsky for inclusion in feff 
c
c     special features:  renormalizes charge density at wigner-seitz
c                        radius
c
c     version 2 (30 september 87): renormalized coulomb potential and
c     renormalized charge density are produced to be used in XAFS
c     calculations by cphase program. j.j. rehr, j. mustre  university
c     of washington., a.djaoui university of essex.
c     please acknowledge use. r. c. albers  (los alamos national lab)
c     j.j. rehr (university of washington),
c
c     Subroutine calling hierarchy siz 1/8/90
c     ATOM
c        INDATA
c           GETORB
c           FPOT
c        DIRAC
c           INOUH
c           INTH
c        POTSL
c        SOMM
c        TOTALE
c           SOMM
c        CDSLD
c           SOMM
c           YKDIR
c        RENORM
c           POTSLW
c
c     Version 1/11/90:  Input and output re-organized to work
c                       easily with overlapped potential code
c                       in FEFF.
c
c     Version Aug 1990: Minor modification to work more easily with
c                       FEFF4, cluster version.  SRHO no longer has
c                       factor of r**2.  INDATA uses rr function to
c                       set r grid.
c     Version Dec 1990: Writes to atom.dat restored
c     Version Feb 1991: Unit 16 opened in atom if necessary
c     June 1992  dirac upper and lower components and total energy
c                passed out for use with matrix element calculations
c
c     Input:   title    title, max 40 characters
c              ifr      index of free atom, used for output labels
c              iz       atomic number of atom
c              ihole    location of electron hole
c              rws      Wigner-Seitz radius
c              ionin    ionicity
c              iprint   print flag, passed through commom /print/
c              ispinr   0, do not save dirac spinors, else save for
c                       orbital ispinr
c
c     Output:  vcoul(251)  coulomb potential (no factor r**2)
c              srho(251)   electron density in form
c                          4*pi*density (formerly 4*pi*density*r**2)
c              dgc0(251)   large component (set if ispinr.ne.0)
c              dpc0(251)   small component (set if ispinr.ne.0)
c              eatom       total energy in rydbergs
c
c     All data is on a grid r(i) = exp (-8.8 + (i-1)*0.05)

      subroutine atom (title, ifr, iz, ihole, rws, ionin, vcoul, srho,
     1                 ispinr, dgc0, dpc0, eatom)

      implicit double precision (a-h,o-z)
      save

c     Save central atom dirac components, see comments below.
      dimension dgc0(251), dpc0(251)

      character*(*)  title
      dimension vcoul(251)
      dimension srho(251)
      common /print/ iprint

      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco

      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc

      common /ps2/ dexv, dexe, dcop, test, teste,
     1             testy, testv, niter, ion, icut, iprat, irnorm

      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1              dpc(251,30)

      character*40 ttl
      character*2  titre
      common /char/ titre(30), ttl

      dimension tden(30)
      character*30 fname
       external dalp
      data harryd /2./

      if (iprint .ge. 3)  then
c        prepare file for atom output
         write(fname,14)  ifr
   14    format('atom', i2.2, '.dat')
         open (unit=16, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'atom')
c        call head (16)
         write(16,*)  ' Free atom ', ifr
      endif

      ttl = title

      nstop=1
      mark=0

      call indata (iz, ihole, rws, ionin)
      iter=1
      do 30 i=1,np
      do 30 j=1,norb
         dgc(i,j)=0.0
         dpc(i,j)=0.0
   30 continue

      if (iprint .ge. 3)  write(16,40) ttl
   40 format (1h1,40x,a40)
      n=-(ion+1)

   60 continue
      do 70 i=1,np
         d(i)=0.0
   70 continue
      tets=test
      ymax=0.0
      vmax=0.0
      emax=0.0

c resolution of the dirac equation for each orbital
      do 150 j=1,norb
         de=den(j)
   80    call dirac (nqn(j),nql(j),nk(j),imax,den(j),dfl(j),dq1(j),j)
            if (nstop.eq.0) go to 110
            if (nstop.ne.362.or.iter.ge.10.or.tets.gt.test) go to 90
            tets=testv
         go to 80
   90    if (iprint .ge. 3)  write(16,100) nstop,nqn(j),titre(j)
  100    format ('  nstop=',i4,'  for the orbital',i3,a2)
         call echo( ' Fatal error.  Wigner-Seitz or muffin tin')
         call echo( '               radius may be too small.')
         go to 999

  110    val=abs((den(j)-de)/de)
         if (val.gt.emax) emax=val
         nmax(j)=imax
         do 140 i=1,np
            val=dgc(i,j)-dp(i)
            if (abs(dp(i)).gt.1.0) val=val/dp(i)
            if (abs(val).lt.abs(ymax)) go to 120
               ymax=val
               y=dp(i)
               yn=dgc(i,j)
  120       val=dpc(i,j)-dq(i)
            if (abs(dq(i)).gt.1.0) val=val/dq(i)
            if (abs(val).lt.abs(ymax)) go to 130
               ymax=val
               y=dq(i)
               yn=dpc(i,j)
  130       dgc(i,j)=dp(i)
            dpc(i,j)=dq(i)
  140    d(i)=d(i)+nel(j)*(dp(i)*dp(i)+dq(i)*dq(i))
  150 continue

c     dgc and dpc are set in loop above, only referenced in remainder
c     of code, so save them into dgc0 and dpc0 here.  Note: np=251,
c     set in indata.  dgc0 is large component
c                     dpc0 is small
      if (ispinr .ne. 0)  then
         do 152  i = 1, np
            dgc0(i) = dgc(i,ispinr)
            dpc0(i) = dpc(i,ispinr)
  152    continue
      endif

      if (mark.eq.0) go to 280

c  This is case mark .ne. 0
c  d is the core electron density resulting from the renormalized pot.
      dval=0.0
      do 160 j=1,norb
  160    dval=dval+nel(j)*den(j)

      dval=dval*2.0
c jm-- core charge density commented away in unit 6 appears in unit 3--
      if (iprint .ge. 3)  write(16,170) dval
  170 format (1h ,' core energy = ',e15.8)

c jm- renormalized potential

c     note conversion to rydbergs using constant harryd
c     passvt is part of old system to pass data directly from
c     ATOM to PHASE
c      do 200 ixx=1,251
c  200    passvt(ixx)=harryd*dr(ixx)*dr(ixx)*dv(ixx)


c  d is the core electron density resulting from the renormalized pot.

c  next write renormalized electron density for each shell
      do 270 j=1,norb
         do 240 i=1,np
            d(i)=dgc(i,j)*sqrt(12.56637062)
  240    continue
  270 continue
      go to 750

c     mark .eq. 0 case
  280 continue

      call potsl (dc,d,dp,dr,dpas,dexv,z,np,ion,icut,dvn)
      if (nuc.le.0) go to 300
         do 290 i=1,nuc
            dc(i)=dc(i)+z/dr(i)+z*((dr(i)/dr(nuc))**2-3.0) /
     1            (dr(nuc)+dr(nuc))
  290    continue
  300 continue
      do 310 i=1,np
         dval=abs(dc(i)-dv(i))
         if ((dr(i)*dc(i)).le.n) dval=-dval/dc(i)
         if (dval.le.vmax) go to 310
            vmax=dval
            j=i
  310 continue

c     print 320, iter,vmax,dr(j),dv(j),dc(j),emax,ymax,yn,y
c 320 format (i5,1pe11.2,3(1pe16.6),2(1pe11.2),2(1pe16.6))

      if (tets.le.test.and.emax.le.teste.and.vmax.le.testv.and.ymax.le
     1 .testy) go to 430
      if (mark.eq.1) go to 430
      iter=iter+1
      if (iter.le.niter) go to 340
      if (iprint .ge. 3)  write(16,330) niter
  330 format (' number of iterations greater than',i4)
      nstop=2
       call echo(' ATOM-Fatal error, too many iterations.')
cc       print*, '   iter, niter ', iter, niter
      go to 999
c potential for the following iteration

  340 continue
      if (iter.eq.2) go to 350
      if (iprat) 350,390,350
  350 dval=1.0-dcop
      do 360 i=1,np
      dvn(i)=dv(i)
      dvf(i)=dc(i)
  360 dv(i)=dval*dv(i)+dcop*dc(i)
      go to 60

  390 continue
      do 400 i=1,np
      dval=dalp(dvn(i),dvf(i),dv(i),dc(i))
      dvn(i)=dv(i)
      dvf(i)=dc(i)
  400 dv(i)=dval*dv(i)+(1.0-dval)*dc(i)
      go to 60

  430 if (iprint .ge. 3)  write(16,40) ttl
      if (iprint .ge. 3)  write(16,460)
  460 format (12x,'energie',12x,'(r4)',14x,'(r2)',14x,'(r)',15x,'(r-1)',
     1 13x,'(r-3)'/)

c valeurs moyennes de r
      do 470 i=1,np
      dvf(i)=dc(i)
  470 dq(i)=0.0
      dval=0.0
      do 560 i=1,norb
      im=nmax(i)
      dval=dval+nel(i)*den(i)
      do 480 j=1,im
  480 dc(j)=dgc(j,i)*dgc(j,i)+dpc(j,i)*dpc(j,i)
      l=5
      if (iabs(nk(i)).eq.1) l=l-1
      do 550 j=1,l
      dp(j)=dfl(i)+dfl(i)
      if (j-2) 490,500,510
  490 n=4
      go to 550
  500 n=2
      go to 550
  510 if (j-4) 520,530,540
  520 n=1
      go to 550
  530 n=-1
      go to 550
  540 n=-3
  550 call somm (dr,dc,dq,dpas,dp(j),n,im)
  560 if (iprint .ge. 3)  write(16,570) nqn(i),titre(i),
     1                                   den(i),(dp(j),j=1,l)
  570 format (i3,a2,6(1pe18.7))

      if (dexv.eq.0.0) go to 650

c energie totale en moyenne spherique
      do 580 i=1,norb
  580 tden(i)=-2.0*den(i)

      dc(1)=1
      do 600 i=1,np
  600 dp(i)=d(i)/dr(i)
      if (nuc.le.0) go to 620
      do 610 i=1,nuc
  610 dp(i)=d(i)*(3.0-dr(i)*dr(i)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
      dc(1)=4
  620 call somm (dr,dp,dq,dpas,dc(1),0,np)
      do 630 i=1,np
      dp(i)=d(i)*dvf(i)
  630 d(i)=d(i)*((d(i)*dr(i))**(1.0/3.0))
      dc(2)=3
      dc(3)=1
      if (nuc.ne.0) dc(3)=4
      call somm (dr,dp,dq,dpas,dc(3),0,np)
      call somm (dr,d,dq,dpas,dc(2),-1,np)
      dc(2)=-3.0*dc(2)/(105.27578**(1.0/3.0))
      dc(1)=-z*dc(1)
      dc(4)=dval-dc(3)
      dval=dval+(dc(1)-dc(3)+(dexe-dexv)*dc(2))/2.0
      dc(3)=(dc(3)-dc(1)-dexv*dc(2))/2.0
      dc(2)=dc(2)*dexe/2.0
      if (iprint .ge. 3)  write(16,640) dval,dc(4),dc(3),dc(2),dc(1)
  640 format (1h0,5x,'et=',1pe14.7,5x,'ec=',1pe14.7,5x,'ee=',1pe14.7,5x,
     1 'ex=',1pe14.7,5x,'en=',1pe14.7)
      go to 660
  650 call totale (dval)
  660 continue

c     pass out eatom (total energy) (factor of 2 is to put energy in
c     rydberg units)
      eatom = 2 * dval

      if (norb.eq.1) go to 710
      if (iprint .ge. 3)  write(16,40) ttl
      if (iprint .ge. 3)  write(16,670)
  670 format (1h0,47x,'overlap integrals         '/)

c overlap integrals
      do 700 i=2,norb
      k=i-1
      do 700 j=1,k
      if (nql(i).ne.nql(j).or.nk(i).ne.nk(j)) go to 700
      im=nmax(j)
      if (nmax(i).lt.im) im=nmax(i)
      do 680 l=1,im
      dq(l)=dpc(l,i)*dpc(l,j)
  680 dc(l)=dgc(l,i)*dgc(l,j)
      dval=dfl(i)+dfl(j)
      call somm (dr,dc,dq,dpas,dval,0,im)
      if (iprint .ge. 3)  write(16,690) nqn(i),titre(i),
     1                                   nqn(j),titre(j),dval
  690 format (34x,i1,a2,i3,a2,f19.7)
  700 continue
  710 call cdsld


      if (irnorm.eq.1) then
         call renorm (dexv, vcoul, srho)
      endif
      do 720 i=1,np
  720 dc(i)=harryd*dv(i)*dr(i)**2
      if (irnorm.ne.1) call fstop(' at ATOM-0')
      norb=norbco
      if (norbco.eq.0) go to 750
      if (mark.eq.1) go to 750
      mark=1
      go to 60

  750 continue

c     return srho as 4*pi*density instead of 4*pi*density*r**2
      do 760  i = 1, 251
         srho(i) = srho(i) / (dr(i)**2)
  760 continue

      if (iprint .ge. 3)  close(unit=16)

      return


  999 continue
       call fstop(' at ATOM-1')
      end
