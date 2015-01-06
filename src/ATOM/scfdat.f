      subroutine scfdat ( ipr1, iph, nph, iz, ihole, xion, iunf, vcoul,
     1     srho, dmag, srhovl, ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,
     2     s02, efrozn, eatom, xntot, xnval, indorb, norbp, eorb, kappa)
c     single configuration Dirac-Fock atom code
c     Ankudinov, Zabinsky, Rehr, Comp.Phys. Comm. 98, p.359 (1996).
c     which is modified Desclaux multi-configuration code.
c     written by a.ankudinov 1996

      implicit double precision (a-h,o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
c     save central atom dirac components, see comments below.
      dimension dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension xntot(0:lx), xnval(30), eorb(30), kappa(30)
      dimension xmag(30)
 
      dimension vcoul(251)
      dimension srho(251), dmag(251), srhovl(251)
c     temporary do not use core-valence separation
      dimension xnvalp(30), indorb(-4:3)
      logical open_16
     
      dimension ovpint(30, 30)
      character*30 fname
c#mn:
       external dsordf

c muatco programm to calculate angular coefficients
c        this programm uses cofcon cofdat dsordf ictime iowrdf
c        lagdat messer nucdev ortdat potrdf soldir 
      common cg(251,30), cp(251,30), bg(10,30), bp(10,30),
     1         fl(30), fix(30), ibgp
c cg (cp) large (small) components
c bg (bp) development coefficients at the origin of large
c    (small) component
c fl power of the first term of development limits.
c ibgp first dimension of the arrays bg and bp

      common/comdir/ cl, dz, gg(251), ag(10), gp(251), ap(10), bid(783)
c  gg,gp are the output from soldir
      common/itescf/ testy, rap(2), teste, nz, norb, norbsc
cXX   common/mulabk/ afgk
      common/inelma/ nem
cXX      dimension afgk(30, 30, 0:3)
      common/messag/ dlabpr, numerr
      character*8 dprlab, dlabpr
      common/ratom1/ xnel(30), en(30), scc(30), scw(30), sce(30),
     1nq(30), kap(30), nmax(30)
      common/scrhf1/ eps(435), nre(30), ipl
      common/snoyau/ dvn(251), anoy(10), nuc
      common/tabtes/ hx, dr(251), test1, test2, ndor, np, nes, method,
     1 idim
      data dprlab /'  scfdat'/
 
      if (ipr1 .ge. 3 .and. iph.le.nph)  then
c        do not want to have extra file
c        prepare file for atom output
         write(fname,14)  iph
   14    format('atom', i2.2, '.dat')
         if (master) then
           open (unit=16, file=fname, status='unknown', iostat=ios)
           open_16 = .true.
         else
           open_16 = .false.
         endif
c        call chopen (ios, fname, 'atom')
c        call head (16)
         if (open_16) write (16,*)  ' free atom ', iph
      endif

c  initialize the data and test parameters
      open_16 = .false.
      jfail = 0
      ibgp = 10
      numerr = 0
      nz = iz
      call inmuat (ihole, xion, iunf, xnval, iholep, xmag, indorb)

      idfock = 1
c     idfock = 2
c     idfock=1  --  pure Dirac-Fock.
c     idfock=2  --  pure LDA
c     idfock=5  --  exchange 5 model.
c     idfock=6  --  exchange 6 model.
      if (idfock.eq.1) then 
         do 42 i=1,30
  42     xnvalp(i) = 0.0d0
      elseif (idfock.eq.2) then
         do 44 i=1,30
  44     xnvalp(i) = xnel(i)
      else
c        use core-valence separation. also change vlda.f
         do 43 i=1,30
  43     xnvalp(i) = xnval(i)
      endif

c     iholep is the index for core hole orbital in all arrays
c     ihole is just a code number for given core hole
c     for 90% of atoms iholep=ihole
      ilast = 0

c   calculate initial orbitals using thomas-fermi model ido=1
c   option to read from cards(ido=2) destroyed
      ido = 1
      if (numerr .eq. 0) then
         a = -xion - 1
         call wfirdf (en, a, nq, kap, nmax, ido)
      endif

      niter = 30
c if niter is negative then schmidt orthogonalization procedure is used
c           niter =1000*n1+100*n2+n3
c n3 is the number of iterations per orbital
      j = 1
      ind = 1
      nter = 0
      do 41 i = 1, norb
 41   scw(i) = 0.
      test1 = testy / rap(1)
      test2 = testy / rap(2)
      netir = abs(niter) * norb
      if (ipr1 .ge. 5 .and. iph.le.nph .and. open_16)  then
         write(16,210) niter, teste, testy
  210    format (5x,'number of iterations',i4,//,
     1        5x,'precision of the energies',1pe9.2,//,
     2        23x,'wave functions  ',1pe9.2,/)
         write(16,220) idim, dr(1), hx
  220    format (' the integration is made on ', i3,
     1        ' points-the first is equal to ' ,f7.4,/,
     2        ' and the step-size pas = ',f7.4,/)
         write(16,230) test1, nes
  230    format ('matching of w.f. with precision', 1pe9.2,
     2        ' in ',i3,' attempts ',/)
         if (nuc .gt. 1)  write(16,250)
  250    format (1h0, 30x,'finite nucleus case used'/)
      endif
 
c     angular coefficients 
c     corrected for valence model. ala
      call muatco(xnvalp)
      if (numerr .ne. 0) go to 711

c     iteration over the number of cycles
 101  iort = 0
         nter = nter + 1
         if (niter .ge. 0) go to 105
c        orthogonalization by schmidt procedure
 104     call ortdat (j)
 105     method = 1
c        calculate lagrange parameters
         if (nre(j).gt.0 .and. ipl.ne.0) call lagdat (j,1)
c        calculate electron potential
         call potrdf (j)
c        add potential due to xc with valence electrons
         call vlda (xnval, srho, srhovl, dmag, ilast, idfock)
c         call vlda (j, xnval, srho, srhovl, dmag, ilast, idfock)
         e = en(j)
         np = idim
c        resolution of the dirac equation
         ifail = 0
         ainf = cg(nmax(j),j)
         call soldir (en(j), fl(j), bg(1,j), bp(1,j), ainf,
     1                nq(j), kap(j), nmax(j), ifail)
         if (ifail.ne.0 .and. jfail.eq.0) jfail = j
         if (jfail.eq.j .and. ifail.eq.0) jfail = 0
         if (numerr. eq. 0) go to 111
         if (iort.ne.0 .or. niter.lt.0) go to 711
         iort = 1
         go to 104

 111     sce(j) = abs ((e - en(j)) / en(j))
c        variation of the wave function using two iterations
         k = nmax(j)
         pr = 0.
         do 121 i = 1, k
            w = cg(i,j) - gg(i)
            if (abs(w) .le. abs(pr)) go to 115
            pr = w
            a = cg(i,j)
            b = gg(i)
 115        w = cp(i,j) - gp(i)
            if (abs(w) .le. abs(pr)) go to 121
            pr = w
            a = cp(i,j)
            b = gp(i)
 121     continue
c        write original Desclaux output on screen and into the logfile
c        write (slog,'(i4, i3, 2(1pe11.2), 2(1pd16.6), 4x, a, i2)')
c    1   nter, j, sce(j), pr, a, b, 'method', method
c        call wlog(slog)

c        acceleration of the convergence
         b = scc(j)
         call cofcon (a, b, pr, scw(j))
         scc(j) = b
         do 151 i = 1, k
            gg(i) = b * gg(i) + a * cg(i,j)
 151        gp(i) = b * gp(i) + a * cp(i,j)
         do 155 i = 1, ndor
            ag(i) = b * ag(i) + a * bg(i,j)
 155        ap(i) = b * ap(i) + a * bp(i,j)
c        normalization of the wave function
         a = dsordf (j, k, 0, 4, fl(j))
         a = sqrt(a)
         do 171 i = 1, np
            cg(i,j) = gg(i) / a
 171        cp(i,j) = gp(i) / a
         do 175 i = 1, ndor
            bg(i,j) = ag(i) / a
 175        bp(i,j) = ap(i) / a
c        determination of the next orbital to calculate
         if (nter.lt.norbsc .or. (ind.lt.0 .and. j.lt.norbsc)) then
            j = j + 1
            go to 451
         endif
            j = j + 1
         pr = 0.
         do 301 i = 1, norbsc
            w = abs (scw(i))
            if (w .gt. pr) then
               pr = w
               j = i
            endif
 301     continue
         if (j .gt. norbsc) j = 1
         if (pr .gt. testy) go to 421
         pr = 0.
         do 321 i = 1, norbsc
            w = abs (sce(i))
            if (w .gt. pr) then
               pr = w
               j = i
            endif
 321     continue
         if (pr .ge. teste) go to 421
         if (ind .lt. 0) go to 999
         ind = -1
         j = 1
         go to 451

 421     ind = 1
 451  if (nter .le. netir) go to 101
      numerr= 192011
c **** number of iterations exceeded the limit
      dlabpr = dprlab
 711  call messer
      call par_stop('SCFDAT-1')
 999  if (numerr .eq. 0) then
         if (jfail .ne. 0) then
            call par_stop(
     1    '  Failed to match lower component, results are meaningless')
c           stop
         endif
c        tabulation of the results
         if (ipr1 .ge. 5 .and. iph.le.nph)  call tabrat
         call etotal (16, kap, xnel, xnvalp, en, eatom)
         do 504 ix = 1,251 
 504       dmag(ix)=0.0d0 
         ilast = 1
         iorb = 0
c        use to test SIC
c         do 505 iorb = 1,norb
         call vlda (xnval, srho, srhovl, dmag, ilast, idfock)
c         call vlda (iorb, xnval, srho, srhovl, dmag, ilast, idfock)
c 505       call vlda (iorb, xnval, srho, srhovl, dmag, ilast, idfock)
         ecorr =2.0
         call somm(dr,dmag,dmag,hx, ecorr,0,idim)
         eatom = (eatom-ecorr/4.0) 

c        jcore = 1

c        prepare information for SCMT and core-valence separation
         norbp = norb
         do 499 i = 0,lx
  499    xntot(i)=0.0d0
         do 500 j = 1, norb
           eorb(j) = en(j) 
           kappa(j) = kap(j)
           i = kap(j)
           if (kap(j) .lt.0) i=-kap(j)-1
           if (i.le.lx) xntot(i)=xntot(i)+xnval(j)
  500    continue
c 500     if (xnel(j).gt.xnval(j) .and. nmax(j).gt.jcore) jcore=nmax(j)

c  get difference in spin-up and -down densities per spin 
c  the spin - polarizable orbitals are specified in subroutine getorb
c  The spin amplitude and directions are taken care of in subroutine ovrlp
c  and specified in feff.inp file
         spin = 0
         do 530 i = 1, idim
  530    dmag(i) = 0.0
         do 536 iorb = 1, norb
           spin = spin + xmag(iorb)
           do 535 i = 1, np
  535      dmag(i)= dmag(i)+ xmag(iorb)* (cg(i,iorb)**2 + cp(i,iorb)**2)
  536    continue
         if (spin.gt.0.d0) then
c          normalize dmag per  spin
           do 537 i = 1, np
  537      dmag(i) = dmag(i) / spin
         endif

c  return coulomb potential
c  fix later: can be replaced by potrdf
         call potslw (vcoul, srho, dr, hx, idim)
         do 510 i = 1, 251
  510      vcoul(i) = (vcoul(i) - nz / dr(i)) 

c        return srho as 4*pi*density instead of 4*pi*density*r**2
         do 560  i = 1, 251
            srho(i) = srho(i) / (dr(i)**2)
            dmag(i) = dmag(i) / (dr(i)**2)
            srhovl(i) = srhovl(i) / (dr(i)**2)
  560    continue

         if (ipr1 .ge. 3 .and. iph.le.nph)  close(unit=16)

         if (ispinr .ne. 0)  then
c        need kap(i) for central atom without core hole, all output of
c        getorb is dummy, except iholep and kap(i) which is put in nq(i)
            call getorb (iz, ispinr, xion, iunf, i, j, indorb,
     1                   iholep, nre, nq, scw, sce, eps)
            do 552  i = 1, nmax(iholep)
               dgc0(i) = cg(i,iholep)
               dpc0(i) = cp(i,iholep)
  552       continue
            do 553  i = nmax(iholep) + 1, 251
               dgc0(i) = 0.0d0
               dpc0(i) = 0.0d0
  553       continue
         endif

         do 590 j = 1, 30
            do 570 i = 1, nmax(j)
               dgc(i,j,iph) = cg(i,j)
               dpc(i,j,iph) = cp(i,j)
  570       continue
            do 575 i = nmax(j) + 1, 251
               dgc(i,j,iph) = 0.0d 00
               dpc(i,j,iph) = 0.0d 00
  575       continue
            do 580 i = 1, 10
               adgc(i,j,iph) = bg(i,j)
               adpc(i,j,iph) = bp(i,j)
  580       continue
  590    continue
      endif

c     calc. overlap integrals for the final and initial state orbitals
c     of the central atom
      if (iholep .gt. 0 .and. iholep.lt.30 .and. ihole.le.0) then
c        this logic is fulfilled only in the last call of scfdat
c        in subroutine pot ( ihole=0 and iholep=ispinr.neq.0)
         efrozn = en(iholep) 
         do 790 i = 1, norb
c          to handle special case when electron added to new orbital
           if (nq(i) .eq. kap(i)) then
              itr = 0
           elseif (nq(i+1) .eq. kap(i)) then
              itr = 1
           else
              call wlog
     1        ('  If it is not la, gd or np, please, give us a call')
              call wlog('  s02 is overestimated')
              do 710 j = 1, i - 1
  710            ovpint(j,i) = 0.0
              ovpint(i,i) = 1.0
              goto 780
           endif
           i0 = i + itr
           iph1 = 0
           if (iph.eq.0) iph1 = nph + 1
           do 720 ir = 1, idim
             gg(ir) = dgc(ir, i0, iph1)
  720        gp(ir) = dpc(ir, i0, iph1)
           do 730 ir = 1, ndor
              ag(ir) = adgc(ir, i0, iph1)
  730         ap(ir) = adpc(ir, i0, iph1)
           do 770 j = 1, norb
             if (kap(i) .ne. kap(j)) go to 770
             ovpint(i,j) = dsordf ( j, j, 0, 3, fl(i))
  770      continue
  780      continue
  790    continue
         do 810 j=1,norb
             xnel(j) = xnel(j)-xnval(j)
 810     continue

c        need better control here. for now always print fpf0.dat
c        if (ipr1.ge.3) call  fpf0 ( iz, iholep, srho, dr, hx,
c         call  fpf0 ( iz, iholep, srho, dr, hx,
c     1     dgc0, dpc0, dgc, dpc,
c     2     eatom, xnel, norb, eorb, kappa)

         call s02at (iholep, norb, kap, xnel, ovpint, s02)
c        print*,'z=',iz, '   s02 calculated = ', s02
      endif

      return
      end
