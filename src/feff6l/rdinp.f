      subroutine rdinp (finp, 
     1                  mphase, mpath, mfeff, mchi, ms,
     1                  ntitle, title, ltit,
     2                  critcw,
     1                  ipr2, ipr3, ipr4,
     1                  s02, tk, thetad, sig2g,
     1                  nlegxx,
     1                  rmax, critpw, pcritk, pcrith, nncrit,
     2                  icsig, iorder, vrcorr, vicorr, isporb)

c     Read input for multiple scattering feff
      implicit double precision (a-h, o-z)
	character*128 finp
      include 'const.h'
      include 'dim.h'
      include 'pola.h'

c     Following passed to pathfinder, which is single precision.
c     Be careful to always declare these!
      real rmax, critpw, pcritk, pcrith

c     Data for potph (see arrays.h for comments)
      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension ifrph(0:nphx)
      dimension xnatph(0:nphx)
      dimension folp(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension ion(0:nfrx)
      dimension iz(0:nfrx)

      character*6  potlbl(0:nphx)

c     Local stuff
      character*150  line
      parameter (nwordx = 12)
      character*15 words(nwordx)

      parameter (ntitx = 10)
      character*128  title(ntitx), messag
      dimension ltit(ntitx)
      dimension ionph(0:nphx), izph(0:nphx)
      logical iscomm
      parameter (nssx = 16)
      dimension indss(nssx), iphss(nssx)
      dimension degss(nssx), rss(nssx)
      logical nogeom

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

c     initialize things

      ihole = 1
      ntitle = 0
      ixc = 0
      vr0 = 0
      vi0 = 0
      rs0 = 0
      rmax = -1
      tk = 0
      thetad = 0
      sig2g = 0
      rmult = 1
      s02 = 1
      mphase = 1
      mpath = 1
      mfeff = 1
      mchi = 1
      ms = 0
      ipr1 = 0
      ipr2 = 0
      ipr3 = 0
      ipr4 = 0
      nlegxx = 10
      xkmin = 0
      xkmax = 20
      critcw = 4.0
      critpw = 2.5
      pcritk = 0
      pcrith = 0
      nogeom = .false.
      icsig = 1
      iorder = 2
      ixanes = 0
      vrcorr = 0
      vicorr = 0
      iafolp = 0
      intclc = 0
      nemax = nex
      isporb = -1

c     average over polarization by default
      pola = .false.
      elpty = 0
      do 50 i = 1, 3 
         evec(i) = 0
         ivec(i) = 0
  50  continue 

c     nncrit is number of necrit points to use.  necrit is
c     currently 9, this was at once an input used for testing.
      nncrit = 9

      nat = 0
      do 100  iat = 1, natx
         iphat(iat) = -1
  100 continue

      nss = 0
      do 102  iss = 1, nssx
         indss(iss) = 0
         iphss(iss) = 0
         degss(iss) = 0
         rss(iss) = 0
  102 continue

      nph = 0
      do 110  iph = 0, nphx
         iatph(iph) = 0
         ifrph(iph) = -1
         xnatph(iph) = 0
         folp(iph) = 1
         novr(iph) = 0
         ionph(iph) = 0
         izph(iph) = 0
         potlbl(iph) = ' '
  110 continue

      nfr = 0
      do 120  ifr = 0, nfrx
         ion(ifr) = 0
         iz(ifr) = 0
  120 continue

c     Open feff.inp, the input file we're going to read
      open (unit=1, file=finp, status='old', iostat=ios)
      call chopen (ios, finp, 'rdinp')

c     tokens  0 if not a token
c             1 if ATOM (ATOMS)
c             2 if HOLE
c             3 if OVER (OVERLAP)
c             4 if CONT (CONTROL)
c             5 if EXCH (EXCHANGE)
c             6 if ION
c             7 if TITL (TITLE)
c             8 if FOLP
c             9 if RMAX
c            10 if DEBY (DEBYE)
c            11 if RMUL (RMULTIPLIER)
c            12 if SS
c            13 if PRIN (PRINT)
c            14 if POTE (POTENTIALS)
c            15 if NLEG
c            16 if REQU (REQUIRE), now dead
c            17 if KLIM (KLIMIT)
c            18 if CRIT (CRITERIA)
c            19 if NOGEOM
c            20 if CSIG
c            21 if IORDER
c            22 if PCRI (PCRITERIA)
c            23 if SIG2
c            24 if XANE (XANES), disabled for current release
c            25 if CORR (CORRECTIONS)
c            26 if AFOL (AFOLP)
c            27 if NEMA (NEMAX)
c            28 if INTCALC
c            29 if POLA (POLARIZATION)
c            30 if ELLI (ELLIPTICITY) 
c            31 if ISPO (ISPORB)
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading atom positions
c                2 reading overlap instructions for unique pot
c                3 reading unique potential definitions

      mode = 0
 200  continue
         line   = ' '
         iret = iread(1, line)
cc         print*, ' line: ', iret, ': ', line(1:40)
         if (iret.eq. 0) goto  200
         if (iret.le.-1) line = 'END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1))

c        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
c              ATOM
c              Following lines are atom postions, one per line
               mode = 1
            elseif (itok .eq. 2)  then
c              HOLE     1  1.0
c                   holecode s02
               read(words(2),20,err=900)  ihole
               read(words(3),30,err=900)  s02
               mode = 0
            elseif (itok .eq. 3)  then
c              OVERLAP iph
c                  iph  n  r
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               mode = 2
            elseif (itok .eq. 4)  then
c              CONTROL  mphase, mpath, mfeff, mchi
c               0 - do not run modules, 1 - run module
               read(words(2),20,err=900)  mphase
               read(words(3),20,err=900)  mpath
               read(words(4),20,err=900)  mfeff
               read(words(5),20,err=900)  mchi
               mode = 0
            elseif (itok .eq. 5)  then
c              EXCHANGE  ixc  vr0  vi0
c              ixc=0  Hedin-Lunqvist + const real & imag part
c              ixc=1  Dirac-Hara + const real & imag part
c              ixc=2  ground state + const real & imag part
c              ixc=3  Dirac-Hara + HL imag part + const real & imag part
c              ixc=4  DH below rs0 + HL above rs0 + const real
c                     & imag part, form is
c                     EXCHANGE  4  vr0  vi0  rs0
c              vr0 is const imag part of potential
c              vi0 is const imag part of potential
c              Default is HL.
               read(words(2),20,err=900)  ixc
               read(words(3),30,err=900)  vr0
               read(words(4),30,err=900)  vi0
               if (ixc .eq. 4) read(words(5),30,err=900)  rs0
               if (ixc .ge. 3)  call warnex(1)
               mode = 0
            elseif (itok .eq. 6)  then
c              ION  iph ionph(iph)
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),20,err=900)  ionph(iph)
               mode = 0
            elseif (itok .eq. 7)  then
c              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. ntitx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  call echo ('Too many title lines, title ignored')
                  call echo( line(1:79))
               endif
               mode = 0
            elseif (itok .eq. 8)  then
c              FOLP iph folp (overlap factor, default 1)
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),30,err=900)  folp(iph)
               mode = 0
            elseif (itok .eq. 9)  then
c              RMAX  rmax (max r for ss and pathfinder)
               read(words(2),30,err=900)  rmax
               mode = 0
            elseif (itok .eq. 10)  then
c              DEBYE  temp debye-temp
c                   temps in kelvin
c                   if tk and thetad > 0, use these instead of sig2g
               read(words(2),30,err=900)  tk
               read(words(3),30,err=900)  thetad
               mode = 0
            elseif (itok .eq. 11)  then
c              RMULTIPLIER  rmult
c              Multiples atom coord, rss, overlap and rmax distances by
c              rmult (default 1).  DOES NOT modify sig2g
               read(words(2),30,err=900)  rmult
               mode = 0
            elseif (itok .eq. 12)  then
c              SS index ipot deg rss
               nss = nss + 1
               if (nss .gt. nssx)  then
                  write(messag,'(1x,a,i4)')
     $                 'Too many SS paths requested, max is ', nssx
                  call echo(messag)
                  call fstop(' at RDINP')
               endif
               read(words(2),20,err=900)  indss(nss)
               read(words(3),20,err=900)  iphss(nss)
               read(words(4),30,err=900)  degss(nss)
               read(words(5),30,err=900)  rss(nss)
               mode = 0
            elseif (itok .eq. 13)  then
c              PRINT  ipr1  ipr2  ipr3  ipr4
c              print flags for various modules
c              ipr1 potph  0 phase.bin only
c                          1 add misc.dat
c                          2 add pot.dat, phase.dat
c                          5 add atom.dat
c                          6 add central atom dirac stuff
c                          7 stop after doing central atom dirac stuff
c              ipr2 pathfinder  0 paths.dat only
c                               1 add crit.dat
c                               2 keep geom.dat
c                               3 add fbeta files
c                               5 special magic code, crit&geom only
c                                 not paths.dat.  Use for path studies
c              ipr3 genfmt 0 files.dat, feff.dats that pass 2/3 of
c                            curved wave importance ratio
c                          1 keep all feff.dats
c              ipr4 ff2chi 0 chi.dat
c                          1 add sig2.dat with debye waller factors
c                          2 add chipnnnn.dat for each path
               read(words(2),20,err=900)  ipr1
               read(words(3),20,err=900)  ipr2
               read(words(4),20,err=900)  ipr3
               read(words(5),20,err=900)  ipr4
               mode = 0
            elseif (itok .eq. 14)  then
c              POTENTIALS
c              Following lines are unique potential defs, 1 per line
               mode = 3
            elseif (itok .eq. 15)  then
c              NLEG nlegmax (for pathfinder)
               read(words(2),20,err=900)  nlegxx
               mode = 0
            elseif (itok .eq. 16)  then
c              REQUIRE rreq, ipot (for pathfinder, require than ms paths
c                            length >rreq contain atom ipot)
               call fstop('REQUIRE card not available')
            elseif (itok .eq. 17)  then
c              KLIMIT xkmin, xkmax
               call echo(' Warning: KLIMIT not available,'//
     $              ' run continues.')
               mode = 0
            elseif (itok .eq. 18)  then
c              CRIT critcw critpw
               read(words(2),30,err=900)  critcw
               read(words(3),30,err=900)  critpw
               mode = 0
            elseif (itok .eq. 19)  then
c              NOGEOM (do not write geom.dat)
               nogeom = .true.
               mode = 0
            elseif (itok .eq. 20)  then
c              CSIG (use complex momentum with debye waller factor)
c              note: this is always on anyway, so this card unnecessary
               icsig = 1
               mode = 0
            elseif (itok .eq. 21)  then
c              IORDER  iorder (used in genfmt, see setlam for meaning)
               read(words(2),20,err=900)  iorder
               call warnex(2)
               mode = 0
            elseif (itok .eq. 22)  then
c              PCRIT  pcritk pcrith
c                     (keep and heap criteria for pathfinder)
               read(words(2),30,err=900)  pcritk
               read(words(3),30,err=900)  pcrith
               mode = 0
            elseif (itok .eq. 23)  then
c              SIG2  sig2g   global sig2 written to files.dat
               read(words(2),30,err=900)  sig2g
               mode = 0
            elseif (itok .eq. 24)  then
c              XANES
c              Use extended k range for xanes
               ixanes = 1
c              to avoid problems with debye waller factors below the
c              edge, always use complex p for debye waller
               icsig = 1
               call fstop('XANES card not available')
            elseif (itok .eq. 25)  then
c              CORRECTIONS  e0-shift, lambda correction
c              e0 shift is in eV, edge will be edge-e0
c              lambda corr is a const imag energy in eV
c              e0 and lambda corr same as vr0 and vi0 in EXCH card
               read(words(2),30,err=900)  vrcorr
               read(words(3),30,err=900)  vicorr
               mode = 0
            elseif (itok .eq. 26)  then
c              AFOLP use generalized automatic folp
               iafolp = 1
               mode =0
            elseif (itok .eq. 27)  then
c              NEMAX  nemax for energy grid
               read(words(2),20,err=900)  nemax
               call warnex(3)
               if (nemax .gt. nex)  then
                  write(messag,'(1x,a,i5,a,i5)')
     $                 'nemax too big, resetting from ',
     $                 nemax, ' to ', nex
                  call echo(messag)
                  nemax = nex
               endif
               mode = 0
            elseif (itok .eq. 28)  then
c              INTCALC  intclc
c              0  use average over all atoms
c              1  use current experimental method 1
c              2  use current experimental method 2
c              read(words(2),20,err=900)  intclc
               call echo(' Warning: INTCALC not available,'//
     $              ' run continues.')            
               mode = 0
            elseif (itok .eq. 29)  then
c              POLARIZATION  X Y Z
               pola = .true.
c              run polarization code if 'pola' is true
c              run usual feff otherwise
               read(words(2),30,err=900)  evec(1)
               read(words(3),30,err=900)  evec(2)
               read(words(4),30,err=900)  evec(3)
               mode = 0
            elseif (itok .eq. 30)  then
c              ELLIPTICITY  E incident direction
               read(words(2),30,err=900)  elpty
               read(words(3),30,err=900)  ivec(1)
               read(words(4),30,err=900)  ivec(2)
               read(words(5),30,err=900)  ivec(3)
               mode = 0
            elseif (itok .eq. 31)  then
c              ISPORB  isporb
               read(words(2),20,err=900)  isporb
               write(messag,'(1x,a,i5)') ' isporb set ', isporb
               call echo(messag)
               mode = 0
            elseif (itok .eq. -1)  then
c              END
               goto 220
            else
               ilen = istrln(words(1))
               write(messag,'(1x,3a,i5)')'Unknown keyword: "',
     $              words(1)(:ilen),'" at line:'
               call echo(messag)
               ilen = istrln(line)
               call echo( '   '//line(:ilen) )
               call echo(' See FEFF document for valid keywords.')
               call fstop(' at RDINP.')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
c              We're done reading atoms.
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            nat = nat+1
            if (nat .gt. natx)  then
               write(messag,'(1x,a,i5)')
     $              'Too many atoms in ATOMS list. Maximum is ', natx
               call fstop(messag)
            endif
            read(words(1),30,err=900)  rat(1,nat)
            read(words(2),30,err=900)  rat(2,nat)
            read(words(3),30,err=900)  rat(3,nat)
            read(words(4),20,err=900)  iphat(nat)
         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
c              We're done reading these overlap instructions.
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            novr(iph) = novr(iph)+1
            iovr = novr(iph)
            if (iovr .gt. novrx)  then
               write(messag,'(1x,a,i5)')
     $              'Too many overlap shells, max is ', novrx
               call fstop(messag)
            endif
            read(words(1),20,err=900) iphovr(iovr,iph)
            read(words(2),20,err=900) nnovr(iovr,iph)
            read(words(3),30,err=900) rovr(iovr,iph)
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
c              We're done reading unique potential definitions
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900)  iph
            if (iph .lt. 0  .or.  iph .gt. nphx)  then
               write(messag,'(1x,a,i3,a,i3)')
     $              'Unique potential ', iph,
     $              ' not allowed. Must be between 0 and ', nphx
               
               call echo(messag)
               call echo(line)
               call fstop('at RDINP')

            endif
            read(words(2),20,err=900)  izph(iph)
c           No potential label if user didn't give us one
c           Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)
         else
            write(messag,'(1x,a,i4)') 'Unknown mode: ', mode
            call echo(messag)
            call fstop(' at RDINP')
         endif
      goto 200
  220 continue

c     We're done reading the input file, close it.
      close (unit=1)

c     Fix up defaults, error check limits, figure out free atoms, etc.

      if (pola) then
c        make polarization tensor
         call mkptz
      endif

c     Find out how many unique potentials we have
      nph = 0
      do 300  iph = nphx, 0, -1
         if (izph(iph) .gt. 0)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue
c     Must have central atom
      if (izph(0) .le. 0)  then
         call fstop(' No absorbing atom (ipot=0) defined')
      endif

c     Then find model atoms for unique pots that have them
      do 330  iph = 0, nphx
c        Use first atom in atom list that is of unique pot iph
         do 320  iat = 1, nat
            if (iph .eq. iphat(iat))  then
               iatph(iph) = iat
               goto 321
            endif
  320    continue
  321    continue
  330 continue
c     if iatph > 0, a model atom has been found.

c     No gaps allowed in unique pots.  Make sure we have enough
c     to overlap all unique pots 0 to nph.
      do 340  iph = 0, nph
         if (iatph(iph) .le. 0  .and.  novr(iph) .le. 0)  then
c           No model atom, no overlap cards, can't do this unique pot
            write(messag,'(1x,a,i5)')
     $           ' No atoms or overlap cards for unique pot ', iph
            call echo(messag)
            call echo(' Cannot calculate potentials, etc.')
            call fstop(' at RDINP')
         endif
  340 continue

c     Need number of atoms of each unique pot, count them.  If none,
c     set to one.
      do 350  iph = 0, nph
         xnatph(iph) = 0
         do 346  iat = 1, nat
            if (iphat(iat) .eq. iph)  xnatph(iph) = xnatph(iph)+1
  346    continue
         if (xnatph(iph) .le. 0)  xnatph(iph) = 1
  350 continue

c     Do the free atom shuffling, do central atom as special case
      iz(0) = izph(0)
      ion(0) = ionph(0)
      ifrph(0) = 0
      nfr = 0
      do 390  iph = 1, nph
         ifrph(iph) = -1
         do 380  ifr = 1, nfr
            if (iz(ifr).eq.izph(iph) .and. ion(ifr).eq.ionph(iph)) then
               ifrph(iph) = ifr
               goto 381
            endif
  380    continue
  381    continue
c        add free atom type if necessary
         if (ifrph(iph) .lt. 0)  then
            nfr = nfr+1
            if (nfr .gt. nfrx)  then
               write(messag,'(1x,a,i5)')
     $              ' Too many free atoms, max is ', nfrx
               call echo(messag)
               call fstop(' at RDINP')
            endif
            ion(nfr) = ionph(iph)
            iz(nfr) = izph(iph)
            ifrph(iph) = nfr
         endif
  390 continue

c     Find central atom (only 1 permitted)
      iatabs = -1
      do 400  iat = 1, nat
         if (iphat(iat) .eq. 0)  then
            if (iatabs .lt. 0)  then
               iatabs = iat
            else
               call echo(' More than one absorbing atom (ipot=0)')
               call echo(' Only one absorbing atom allowed')
               call fstop(' at RDINP')
            endif
         endif
  400 continue

c     Find distance to nearest and most distant atom (use overlap card
c     if no atoms specified.)
      if (iatabs .lt. 0  .or.  nat .lt. 2)  then
         ratmin = rovr(1,0)
         ratmax = rovr(novr(0),0)
      else
         ratmax = 0
         ratmin = 1.0e10
         do 412  iat = 1, nat
c           skip absorbing atom
            if (iat .eq. iatabs)  goto 412
            tmp = dist (rat(1,iat), rat(1,iatabs))
            if (tmp .gt. ratmax)  ratmax = tmp
            if (tmp .lt. ratmin)  ratmin = tmp
  412    continue
      endif

c     Set rmax if necessary
      if (rmax.le.0 .and. nss.le.0)  then
c        set to min (2+ times ratmin, ratmax)
         rmax = min (2.001 * ratmin, ratmax)
      endif

c     Set core hole lifetime (central atom quantity)
      ifr = ifrph(0)
      call setgam (iz(ifr), ihole, gamach)
ccc      print*, ' RDINP SETGAM ', ifr, iz(ifr), ihole, gamach

c     Set s02 if necessary
      if (s02 .le. 1.0e-10)  s02 = 1

c     Convert everything to code units, and use rmult factor
c     rmax is for pathfinder, so leave it in Ang.
      rmax = rmax * rmult
      vr0 = vr0 / ryd
      vi0 = vi0 / ryd
      vrcorr = vrcorr / ryd
      vicorr = vicorr / ryd
      xkmin = xkmin * bohr
      xkmax = xkmax * bohr
      do 430  iat = 1, nat
         do 420  i = 1, 3
            rat(i,iat) = rat(i,iat) * rmult / bohr
  420    continue
  430 continue
      do 460  iph = 0, nph
         do 450  iovr = 1, novr(iph)
            rovr(iovr,iph) = rovr(iovr,iph) * rmult / bohr
  450    continue
  460 continue
      do 462  iss = 1, nss
c        rss used only to make paths.dat, so leave it in Angstroms.
         rss(iss) = rss(iss) * rmult
  462 continue

c     Check if 2 atoms are closer together than 1.75 ryd (~.93 Ang)
      ratmin = 1.0e20
      do 480  iat = 1, nat
         do 470  jat = iat+1, nat
            rtmp = dist(rat(1,iat),rat(1,jat))
            if (rtmp .lt. ratmin)  ratmin = rtmp
            if (rtmp .lt. 1.75)  then
c           if (dist(rat(1,iat),rat(1,jat)) .lt. 1.5)  then
               call echo(' WARNING:  TWO ATOMS VERY CLOSE TOGETHER.'//
     1                 '  CHECK INPUT.')
               write(messag,'(1x,a,2i5)') ' atoms ', iat, jat
               call echo (messag)
               write(messag,'(1x,i4,3f11.5)')  iat,
     $             rat(1,iat)*bohr, rat(2,iat)*bohr, rat(3,iat)*bohr
               call echo(messag)
               write(messag,'(1x,i4,3f11.5)')  jat,
     $             rat(1,jat)*bohr, rat(2,jat)*bohr, rat(3,jat)*bohr
               call echo(messag)
               call echo(' Run continues...')
            endif
  470    continue
  480 continue

c     default to k shell
      if (isporb .lt. 0)  isporb = 1

c     Clean up control flags
      if (mphase .ne. 0)  mphase = 1
      if (mpath  .ne. 0)  mpath = 1
      if (mfeff  .ne. 0)  mfeff = 1
      if (mchi   .ne. 0)  mchi = 1
      if (nss    .le. 0)  ms = 1

      if (ntitle .le. 0)  then
         ntitle = 1
         title(i) = 'No title input'
      endif
      do 490  i = 1, ntitle
         ltit(i) = istrln (title(i))
  490 continue

c     Write output files

c     For potph...
      if (mphase .eq. 1)  then
         open (unit=1, file='potph.dat', status='unknown', iostat=ios)
         call chopen (ios, 'potph.dat', 'rdinp')
         do 705  i = 1, ntitle
            write(1,700)  title(i)(1:ltit(i))
  700       format (1x, a)
  705    continue
         write(1,706)
  706    format (1x, 79('-'))
         write(1,709) ihole, gamach, ipr1, iafolp, intclc
  709    format(i5, 1p, e14.6, 3i4, 
     1         ' ihole, gamach, iprint, iafolp, intclc')
         write(1,702)  ixc, vr0, vi0, rs0
  702    format (i5, 1p, 3e14.6, ' ixc, vr0, vi0, rs0')
         write(1,701)  ixanes, nemax, xkmin, xkmax
  701    format (2i5, 1p, 2e14.6, 
     1           ' ixanes, nemax, xkmin, xkmax (inv bohr)')
         write(1,707) nfr, '  nfr'
  707    format (i5, a)
         do 710  ifr = 0, nfr
            write(1,708)  ifr, iz(ifr), ion(ifr)
  708       format (3i5, ' ifr, iz, ion')
  710    continue
         write(1,707) nat, '  nat.   iat, iph, x, y, z'
         do 720  iat = 1, nat
            write(1,715) iat, iphat(iat), (rat(j,iat),j=1,3)
  715       format (2i5, 3f12.6)
  720    continue
         write(1,707) nph, '  nph'
         do 740  iph = 0, nph
            write(1,722) iph, iatph(iph), ifrph(iph), xnatph(iph),
     1                   folp(iph), novr(iph),
     2                   ' iph, iat, ifr, xnat, folp, novr'
  722       format (3i5, 2f12.6, i5, a)
            write(1,723) potlbl(iph)
  723       format (' ''', a6, '''  potlbl')
            do 730  iovr = 1, novr(iph)
               write(1,724) iphovr(iovr,iph), nnovr(iovr,iph),
     1                      rovr(iovr,iph),
     2                      ' ovr...  iph, n, r'
  724       format (2i5, f12.6, a)
  730       continue
  740    continue
         close (unit=1)
      endif

c     Single scattering paths for genfmt
      if (nss .gt. 0  .and.  mpath .eq. 1)  then
         open (unit=1, file='paths.dat', status='unknown', iostat=ios)
         call chopen (ios, 'paths.dat', 'rdinp')
         do 750  i = 1, ntitle
            write(1,748)  title(i)(1:ltit(i))
  748       format (1x, a)
  750    continue
         write(1,751)
  751    format (' Single scattering paths from ss lines cards',
     1           ' in feff input')
         write(1,706)
         do 760  iss = 1, nss
            if (rmax.le.0  .or.  rss(iss).le.rmax)  then
c              NB, rmax and rss are in angstroms
               write(1,752) indss(iss), 2, degss(iss),
     2              rss(iss)
  752          format ( 2i4, f8.3,
     1             '  index,nleg,degeneracy,r=', f8.4)
               write(1,766)
  766          format (' single scattering')
               write(1,754) rss(iss)*bohr, zero, zero, iphss(iss),
     1                      potlbl(iphss(iss))
               write(1,753) zero, zero, zero, 0, potlbl(0)
  753          format (3f12.6, i4,  1x, '''', a6, '''', '  x,y,z,ipot')
  754          format (3f12.6, i4,  1x, '''', a6, '''')
            endif
  760    continue
         close (unit=1)
      endif

c     Atoms for the pathfinder
      if (nss.le.0  .and.  mpath.eq.1  .and.  nat.gt.0)  then
         if (iatabs .le. 0)  then
            call echo(' Absorbing atom coords not specified.')
            call echo(' Cannot find multiple scattering paths.')
            call fstop(' at RDINP')
         endif
c        if user doesn't want geom.dat, don't do it
         if (nogeom)  goto 792
         open (unit=1, file='geom.dat', status='unknown', iostat=ios)
         call chopen (ios, 'geom.dat', 'rdinp')
c        Echo title cards to geom.dat
         do 770  i = 1, ntitle
            write(1,700)  title(i)(1:ltit(i))
  770    continue
         write(1,706)
c        Central atom first
         ii = 0
         write(1,780)  ii, (rat(j,iatabs)*bohr,j=1,3), 0, 1
c        Rest of the atoms (skip central atom)
         do 790   iat = 1, nat
            if (iat .eq. iatabs)  goto 790
            ii = ii+1
            write(1,780)  ii, (rat(j,iat)*bohr,j=1,3), iphat(iat), 1
  780       format (i4, 3f12.6, 2i4)
  790    continue
         close (unit=1)
      endif
  792 continue

      return

  900 continue
       call echo(' Error reading input, bad line follows:')
       call echo(line)
       call fstop(' at RDINP')
       
      end

      function itoken (word)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w

      w = word(1:4)
      if     (w .eq. 'ATOM')  then
         itoken = 1
      elseif (w .eq. 'HOLE')  then
         itoken = 2
      elseif (w .eq. 'OVER')  then
         itoken = 3
      elseif (w .eq. 'CONT')  then
         itoken = 4
      elseif (w .eq. 'EXCH')  then
         itoken = 5
      elseif (w .eq. 'ION ')  then
         itoken = 6
      elseif (w .eq. 'TITL')  then
         itoken = 7
      elseif (w .eq. 'FOLP')  then
         itoken = 8
      elseif (w .eq. 'RMAX')  then
         itoken = 9
      elseif (w .eq. 'DEBY')  then
         itoken = 10
      elseif (w .eq. 'RMUL')  then
         itoken = 11
      elseif (w .eq. 'SS  ')  then
         itoken = 12
      elseif (w .eq. 'PRIN')  then
         itoken = 13
      elseif (w .eq. 'POTE')  then
         itoken = 14
      elseif (w .eq. 'NLEG')  then
         itoken = 15
      elseif (w .eq. 'REQU')  then
         itoken = 16
      elseif (w .eq. 'KLIM')  then
         itoken = 17
      elseif (w .eq. 'CRIT')  then
         itoken = 18
      elseif (w .eq. 'NOGE')  then
         itoken = 19
      elseif (w .eq. 'CSIG')  then
         itoken = 20
      elseif (w .eq. 'IORD')  then
         itoken = 21
      elseif (w .eq. 'PCRI')  then
         itoken = 22
      elseif (w .eq. 'SIG2')  then
         itoken = 23
      elseif (w .eq. 'XANE')  then
         itoken = 24
      elseif (w .eq. 'CORR')  then
         itoken = 25
      elseif (w .eq. 'AFOL')  then
         itoken = 26
      elseif (w .eq. 'NEMA')  then
         itoken = 27
      elseif (w .eq. 'INTC')  then
         itoken = 28
      elseif (w .eq. 'POLA')  then
         itoken = 29
      elseif (w .eq. 'ELLI')  then
         itoken = 30
      elseif (w .eq. 'ISPO')  then
         itoken = 31
      elseif (w .eq. 'END ')  then
         itoken = -1
      else
         itoken = 0
      endif
      return
      end
      logical function iscomm (line)
c     returns true if line is a comment or blank line, false otherwise
      character*(*) line, l1*1,com*3
       parameter(com = '*#;')
       iscomm = .false.
       l1 = line(1:1)
       if (istrln(line).le.0 .or. index(com,l1).ge.1)  iscomm = .true.
       return
      end
      subroutine phstop (iph,line)
      implicit double precision (a-h, o-z)
      character*(*) line, messag*128
      include 'dim.h'
      if (iph .lt. 0  .or.  iph .gt. nphx)  then
         write(messag, 10) iph, nphx, line
   10    format (' Unique potential index', i5, ' out of range.',/,
     1        ' Must be between 0 and', i5, '.  Input line:')
         call echo(messag)
         call echo(line)
         call fstop('at RDINP - PHSTOP')
      endif
      return
      end
      subroutine warnex (i)
      implicit double precision (a-h, o-z)
c     This prints a warning message if the user is using an
c     expert option.
c     i    expert option card
c     1    EXCHANGE with code >= 3
c     2    IORDER
c     3    XANES
c     4    NEMAX
c     5    INTCALC
c     message max of 22 characters to keep warning on 80 char line.
      if (i .eq. 1)  then
         call echo( 'Warning: EXCHANGE with IEXCH> >= 3')
      elseif (i .eq. 2)  then
         call echo( 'Warning: IORDER')
      elseif (i .eq. 3)  then
         call echo( 'Warning: NEMAX')
      endif
       call echo( ' is an "expert user option".')
       call echo( ' Please read documentation'//
     $      ' carefully and check your results:')

      return
      end

!
