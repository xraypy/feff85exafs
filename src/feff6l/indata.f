c     major revision, input now comes from main program feff
c     input data is passed here to indata for processing

      subroutine indata (iz, ihole, wsin, ionin)

      implicit double precision (a-h, o-z)
      save

c     logical unit from which to read input
      parameter (linp = 1)

      common /print/ iprint
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco

      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /ps2/ dexv, dexe, dcop, test, teste,
     1             testy, testv, niter, ion, icut, iprat, irnorm

      character*40 ttl, messag*128
      character*2  titre
      common /char/ titre(30), ttl

       external fpot
      character*2  ttire(9)
      data ttire /'s ', 'p*', 'p ', 'd*', 'd ', 'f*', 'f ','g*', 'g '/

c following variables fixed as data by jm 4/20/87
      data i /0/
      data j /0/
      data k /0/
      data l /0/

      idep   = 0
      icut   = 0
c     Normal use, iprat = 1
      iprat  = 1
      irnorm = 1
      iex    = 1
      nuc    = 0

c idep=0 starting potential = thomas-fermi potential
c idep=1 starting potential read in from cards
c if icut is zero one corrects the potential by -(ion+1)/r
c if iprat is zero the pratt procedure is used
c if iex is zero one uses the unmodified slater exchange
c l=0 standard option for the bloc ofs points and their precision
c finite nuclear size option if nuc is positive
c if irnorm=1 renormalize potential to wigner-seitz radius

      dvc=137.0373
      dsal=dvc+dvc
      iz1=0
      ion1=0
      nuc1=-1
      dpas=0.05
      dr1=0.01
      nes=15

      niter=50

c     orig values:  teste 5.e-6, testy 1.e-5, testv 1.e-5, test 1.e-7
c     JM used teste 5.0e-5 to treat negative ion,
c     SZ changed teste to 1.0e-4 for selenium only to avoid convergence
c     problems with this particular atom.
c     teste set to 1.0e-4 to reduce run time (sz and jjr)
      teste = 1.0e-4
      testy=1.e-04
      testv=1.e-04
      test=1.e-07

      np=251
      nstop=30

c     Set dexv to zero for use with exafs model
      dexv = 0.0

      dexe=1.5
      dcop=0.3

c     i, j, k set to zero when old read statements removed
      i=0
      j=0
      k=0

c iz     = atomic number
c ion    = iz-number of electrons
c norb   = number of orbitals
c idep   = should be either 0 or 1
c i      = number of points for the integration = 251 by default
c j      = number of attempts to adjust the energy = 15 by default
c k      = number of iterations = 50 by default
c norbco = number of core orbitals

c put input data passed from feff into the necessary variables
      ws  = wsin
      ion = ionin
c     given iz, find norb, norbco, then den, nqn, nk and nel for
c     each orbital.
      call getorb (iz, ihole, ion, norb, norbco,
     1            den, nqn, nk, nel)

      if (norb .gt. nstop)  then
         if (iprint .ge. 5)  write(16,44) norb
         write(messag,44) norb
         call echo(messag)
   44    format (' norb=',i3,'too big')
         goto 999
      endif

c dexv = exchange coefficient for the potential: dexv=1. for slater
c dexe = exchange energy coefficient
c dexv should be equal to 2.*dexe/3. in order to satisfy the virial theo
c dexv=0.0 and iex=1, hedin-barth exchange and correlation is used

c dpas  = exponential step;  dr1 defines the first point = dr1/iz
c test  = energy precision criteria in dirac
c teste = self-consistency criteria for the energies of all the electron
c testy = self-consistency criteria for the wavefunctions
c testv = self-consistency criteria for the potential
      z=iz

      if (nuc .gt. 0)  then
         call echo(' enter atomic mass ')
         read (linp,*,end=900) dval
c        dval = atomic mass if nuc positive

         dval=z*(dval**(1.0/3.0))*2.267700e-05/exp(4.0*dpas)
         if (dval .le. dr1)  then
            dr1=dval
            nuc=5
         else
            dval=dval*exp(4.0*dpas)
            do 170 i=6,np
               d1=dr1*exp((i-1)*dpas)
               if (d1.ge.dval) goto 190
  170       continue
            write(messag,180)
            call echo(messag)
            if (iprint .ge. 5)  write(16,180)
  180       format (' error for the atomic mass')
            goto 999

  190       nuc=i
            dr1=dr1*dval/d1
         endif
      endif

      if (iprint .ge. 5)  write(16,210) ttl,niter,teste,testy,testv
  210 format (1h1,40x,A40,//,5x,'number of iterations',i4,//,
     1        5x,'precision of the energies',1pe9.2,//,
     2        23x,'wave functions  ',1pe9.2,//,
     3        23x,'potential',1pe9.2,/)

      xtmp = 8.8
      dr1=z*exp(-xtmp)

      if (iprint .ge. 5)  write(16,220) np,dr1,iz,dpas
  220 format (' the integration is made on ', i3,
     1        ' points-the first is equal to ' ,f7.4, '/', i2,/,
     2        ' and the step-size pas = ',f7.4,/)
      if (iprint .ge. 5)  write(16,230) test,nes,idep,icut,iprat
  230 format (' dans le sous programme resld la precision relative a',
     1        ' obtenir sur l energie est ', 1pe9.2,
     2        ' et le nombre d essais ',i3, //,
     3        'idep=', i3, 5x, 'icut=', i3, 5x, 'iprat=', i3, /)
      if (iprint .ge. 5)  write(16,240) dexv,dexe
  240 format ('  dexv=', 1pe14.7, '     dexe=' ,1pe14.7,
     1        ' if dexv=0.0 hedin-barth corr. and exchan. is used'/)
      k=0
      dval=z*z/(dvc*dvc)


      if (nuc.gt.0) then
         if (iprint .ge. 5)  write(16,250)
  250    format (1h0,30x,'finite nucleus case used'/)
      endif

      do 350 i=1,norb
c        den = orbital energy in atomic units and negative
c        nqn = principal quantum number; nk = kappa quantum number
c        nel = occupation of the orbital

         k=k+nel(i)
         if (den(i) .ge. 0.0)  den(i) = -z*z / (4.0*nqn(i)*nqn(i))

         nql(i)=iabs(nk(i))

         if (nk(i).lt.0) nql(i)=nql(i)-1
         if (nuc .le. 0)  then
            dfl(i)=nk(i)*nk(i)
            dfl(i)=sqrt(dfl(i)-dval)
         else
            dfl(i)=iabs(nk(i))
         endif
         l=2*iabs(nk(i))


         if (nql(i).lt.nqn(i)  .and.  nel(i).le.l  .and.
     1       nqn(i).gt.0       .and.  nql(i).le.4)   goto 340
         write(messag,330) den(i),nqn(i),nql(i),j,nel(i)
         call echo(messag)
            if (iprint .ge. 5)  write(16,330) den(i),nqn(i),nql(i),
     1                                         j,nel(i)
  330       format (' error in the card    ',e15.8,i2,3i2)
            goto 999
  340    continue
         j=nql(i)+iabs(nk(i))
         titre(i)=ttire(j)
         if (iprint .ge. 5)  write(16,345) nqn(i),titre(i),nel(i),
     1                                      den(i)
  345    format (7x,i1,a2,i16,1pe23.7)
  350 continue

      if (iprint .ge. 5)  write(16,370) norbco
  370 format (' no. of core orbitals=',i3)
      if (k.eq.(iz-ion)) goto 390
         write(messag,380)
         call echo(messag)
         if (iprint .ge. 5)  write(16,380)
  380    format (' error for the number of electrons')
         goto 999
  390 continue

      if (iprat .eq. 0)  then
         if (iprint .ge. 5)  write(16,410)
  410    format (1h0,'  the pratt procedure is used'/)
      else
         if (iprint .ge. 5)  write(16,430) ws
  430    format (1h0,'  wigner-seitz radius = ',0pf10.6,/)
      endif

      if (nuc .eq. nuc1)  then
         if (iz.eq.iz1.and.ion.eq.ion1) goto 600
         if (iz.eq.iz1) goto 470
      endif

c     dr(1)=dr1/z
c     do 460 i=2,np
c        dr(i)=dr(1)*exp((i-1)*dpas)
c 460 continue
c     Let's make this consistant with grid in other routines
c     dr array commeted out above
c     SIZ  December 1990
      do 461  i = 1, 251
         dr(i) = rr(i)
  461 continue

c starting potential

  470 val=-ion-1

c     Following code is a block, block ends at line 600
      if (idep .eq. 1)  then

c        read in starting potential (in a.u. and negative) if idep=1
         read (linp,480,end=900) (dv(i),i=1,np)
  480    format (8f9.4)

         if (iprint .ge. 5)  write(16,490) TTL,(dv(i),i=1,np)
  490    format (1h1, 40x, A40, //,
     1           5x, 'starting potential multiplied by r ' /,
     2           10(2x, f9.4))
         dval = -z/dv(1)
         if (nuc.gt.0)  dval = 1.0
         do 500 i=1,np
            dv(i)=dv(i)*dval/dr(i)
  500    continue

      else

         if (idep .ne. 0)  then
            write(messag,510)
            call echo(messag)
            if (iprint .ge. 5)  write(16,510)
  510       format (' error for idep')
            goto 999
         endif

         if (iz.ne.iz1  .or .  ion.le.ion1  .or.   nuc.ne.nuc1)  then
            do 520 i=1,np
               r=dr(i)
               dv(i)=fpot(r,z,val)
  520       continue
            if (nuc .gt. 0)  then
               do 530 i=1,nuc
                  dv(i) = dv(i) + z/dr(i) +
     1                    z*((dr(i)/dr(nuc))**2-3.0)/(dr(nuc)+dr(nuc))
  530          continue
            endif
            goto 600
         endif
      endif
      if (icut .eq. 0)  then
         do 540 i=1,np
            if ((dr(i)*dv(i)).gt.val)  dv(i)=val/dr(i)
  540    continue
      endif
      val=z+dv(1)*dr(1)
      if (nuc.gt.0)  val=z+dv(nuc)*dr(nuc)
      if (abs(val) .ge. 0.1)  then
         write(messag,550)
         call echo(messag)
         if (iprint .ge. 5)  write(16,550)
  550    format (' error for the potential ')
         goto 999
      endif

  600 continue
c     End of block above


      if (norb .ne. 1)  then
         do 620 i=2,norb
            k=i-1
            do 620 j=1,k
            if (nqn(i).eq.nqn(j)  .and. nk(i).eq.nk(j))   then
               write(messag,610)
               call echo(messag)
               if (iprint .ge. 5)  write(16,610)
  610          format (' standard configuration')
               goto 999
            endif
  620    continue
      endif

  630 iz1=iz
      ion1=ion
      nuc1=nuc
      do 660 i=1,norb
         nmax(i)=np
         l=1
         j=nqn(i)-nql(i)
         if ((j-2*(j/2)).eq.0) l=-l
         dq1(i)=l*nk(i)/iabs(nk(i))
         if (nuc .ne. 0  .and.  nk(i) .lt. 0)  then
            dq1(i)=dq1(i)*(nk(i)-dfl(i))*dvc/z
         endif
  660 continue


c  -- Normal return
      return


c  -- Error condition, stop program

c     Unexpected end of file during read -- stop program
  900 continue
      call echo(' Unexpected end of file')

c     Fatal error, stop gracefully (sic)
  999 continue
       call fstop(' at INDATA-1')
      end
