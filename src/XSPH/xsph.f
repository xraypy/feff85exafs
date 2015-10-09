c     Josh - added argument iPl to control many pole self energy.
      subroutine xsph (wrxsec, verbse, phpad,
     -       ipr2, ispec, vixan, xkstep, xkmax, gamach, rgrd,
     1       nph, lmaxph, potlbl, spinph, iatph, nat, rat, iphat,
     2       ixc, vr0, vi0, ixc0, lreal, rfms2, lfms2, l2lp,
     3       ipol, ispin, le2, angks, ptz, iPl,
     4       izstd, ifxc, ipmbse, itdlda, nonlocal,
c        pass parameters from rdpot
     1       ntitle, title, rnrmav, xmu, vint, rhoint,
     2       emu, s02, erelax, wp, ecv, rs, xf, qtotel,
     3       imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     4       dgc0, dpc0, dgc, dpc, adgc, adpc,
     5       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     6       iorb, nohole, ihole,
     7       inters, totvol, iafolp, xion, iunf, iz, jumprm)
c squelch compiler warning about unused dummy variables, apparently
c removed from f85e
c    iGrid, (after iPl)
c  , ibasis)

c     Cluster code -- multiple shell single scattering version of FEFF
c     This program (or subroutine) calculates potentials and phase
c     shifts for unique potentials specifed by atoms and overlap cards.
c
c     Input files:  potph.inp    input data, atoms, overlaps, etc.
c     Output:       phases.bin   phase shifts for use by the rest of the
c                                program
c                   xxx.dat      various diagnostics

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c     control whether xsect.json gets written, .true. for conventional
c     feff, .false. for libpotph
      logical wrxsec, verbse

      double precision col1(nex), col2(nex), col3(nex)
      double precision col4(nex), col5(nex)

c     Notes:
c        nat    number of atoms in problem
c        nph    number of unique potentials
c        ihole  hole code of absorbing atom
c        iph=0 for central atom

c     Specific atom input data
c     iphat(natx)  -  given specific atom, which unique pot?
      dimension iphat(natx)
c     rat(3,natx)  -  cartesian coords of specific atom
      dimension rat(3,natx)

c     Unique potential input data
c     iatph(0:nphx)  - given unique pot, which atom is model?
c                      (0 if none specified for this unique pot)
      dimension iatph(0:nphx)
c     xnatph(0:nphx) - given unique pot, how many atoms are there
c                      of this type? (used for interstitial calc)
      dimension xnatph(0:nphx), spinph(0:nphx)
c     potlbl(0:nphx)    -   label for user convienence
      character*6 potlbl(0:nphx)

c     folp(0:nphx) -  overlap factor for rmt calculation
      dimension folp(0:nphx)
c     novr(0:nphx) -  number of overlap shells for unique pot
      dimension novr(0:nphx)
c     iphovr(novrx,0:nphx) -  unique pot for this overlap shell
      dimension iphovr(novrx,0:nphx)
c     nnovr(novrx,0:nphx) -   number of atoms in overlap shell
      dimension nnovr(novrx,0:nphx)
c     rovr(novrx,0:nphx)  -   r for overlap shell
      dimension rovr(novrx,0:nphx)

c     Free atom data
c     xion(0:nphx)  - ionicity, input
      dimension xion(0:nphx)
c     iz(0:nphx)    - atomic number, input
      dimension iz(0:nphx)

c     Overlap calculation results
c     edens(251,0:nphx)   -   overlapped density*4*pi
      dimension edens(251,0:nphx)
c     vtot (251,0:nphx)   -   overlapped total potential
      dimension vtot (251,0:nphx), vclap (251,0:nphx)

c     Muffin tin calculation results
c     imt(0:nphx)  -  r mesh index just inside rmt
      dimension imt(0:nphx), inrm(0:nphx), folpx(0:nphx)
c     rmt(0:nphx)  -  muffin tin radius
      dimension rmt(0:nphx)
c     rnrm(0:nphx)  -  Norman radius
      dimension rnrm(0:nphx)
c     , qnrm(0:nphx)
c      dimension xnmues(0:lx,0:nphx)
      real rfms2
      integer ipol, ispin, lfms2
      complex*16 ptz
      dimension ptz(-1:1, -1:1)
      dimension lmaxph(0:nphx)

c     PHASE output
c     eref(nex, nspx)         -     interstitial energy ref
      complex*16 eref(nex, nspx)
c     ph(nex,-ltot:ltot,nspx,0:nphx) - phase shifts
      complex*16 ph( nex, -ltot:ltot, nspx, 0:nphx)
c     lmax(0:nphx)      -     number of ang mom levels
      dimension lmax(0:nphx)

      character*80 title(nheadx)

      complex*16  em(nex)
      complex*16  rkk(nex,8,nspx), xsec(nex,nspx)
      dimension xsnorm(nex, nspx)
      dimension dgc0(251), dpc0(251)

c     additioal data needed for relativistic version
      dimension dgc(251,30,0:nphx), dpc(251,30,0:nphx)
      dimension adgc(10,30,0:nphx), adpc(10,30,0:nphx)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      dimension edenvl(251,0:nphx)
c     , eorb(30), kappa(30)
      dimension vvalgs (251,0:nphx), xnval(30,0:nphx), iorb(-4:3,0:nphx)

c     nrx = max number of r points for phase and xsect r grid
      parameter (nrx = nrptx)
      dimension ri(nrptx), vtotph(nrx), rhoph(nrx)
      dimension  dmagx(nrptx), dmag(251,0:nphx)
      dimension dgcx(nrptx), dpcx(nrptx), vvalph(nrx), rhphvl(nrx)
      dimension vch (251), vchp(nrx)

      logical lopt
      character*512 slog

c     Josh - Added iPl for PLASMON card, and iexist for mpse.dat
      integer iPl
c      integer iexist
      
      character*256 phpad


   10 format (4x, a, i5)

c     Phase shift calculation
c     Atom r grid
      dx = 0.05d0
      x0 = 8.8d0
c     Phase r grid
      dxnew = rgrd

c*************************************************************************
c     move the call to rdpot out to ffmod1.f so that this file can be
c     reused between the conventional feff programs and libpotph
c*************************************************************************
c$$$      call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,
c$$$     1                  emu, s02, erelax, wp, ecv,rs,xf, qtotel,
c$$$     2                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,
c$$$     3                  dgc0, dpc0, dgc, dpc, adgc, adpc,
c$$$     3                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
c$$$     4                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
c$$$     5                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
c      lopt=true for the Rivas code of optical constants
       lopt = .false.
       if (lopt) call getedg(ihole,iz(0), emu)
       if (lopt) ik0 = 1
       if (lopt .and. verbse) then
         call wlog('   Fixing edge energy from Elam table...')
         write(slog,fmt="('   emu = ',f10.3,' eV')") emu*hart
         call wlog(slog)
       endif

       do 15 iph = 0, nph
 15    novr(iph) = 0

c  update header, since e.g. one may use diff ixc for the same potential
        call sthead (ntitle, title, nph, iz, rmt, rnrm,
     1          xion, ihole, ixc,
     2          vr0, vi0, gamach, xmu, xf, vint, rs,
     2          lreal, rgrd)
c     Make energy mesh
      edge = xmu - vr0
      if (.not.lopt) emu = emu - vr0

cc    manual input. Later make TDLDA and PMBSE cards
cc    TDLDA ifxc  (izstd=1 if TDLDA card is present)
c     izstd = 0
c     ifxc = 0
cc    PMBSE  ipmbse  nonlocal ifxc itdlda
c     ipmbse = 2
c     ibasis = 2
cc    ipmbse=0 (do not run); 1-LF only; 2-PM only; 3-combined; 
cc           4-combined with s-function kernel splitting
c     nonlocal = 0
cc    nonlocal = 0 (local fxc); 1-read W from pot.ch; 2-from yoshi.dat
c     itdlda = 2
cc    itdlda = 0, 1, 2 should be run in sequence
cc    end of manual input

c     check that logic is set up correctly
      if (ipmbse.le.0) itdlda = 0
      if (nohole.lt.0) then
c       core-hole potential is used already
        if (ifxc.ne.0) then
           if (verbse)
     1            call wlog(' Reset ifxc=0 since NOHOLE card is absent')
          ifxc = 0
          if (ipmbse.gt.0) nonlocal = 0
        endif
        if (ipmbse.eq.3 .and. izstd.eq.0) then
           if (verbse)
     1          call wlog(' Reset ipmbse=1 since NOHOLE card is absent')
          ipmbse = 1
        endif
      endif
      if (izstd.gt.0 .and. itdlda.gt.0) then
c       no need run PMBSE in this case
         if (verbse)
     1          call wlog(' Ignored PMBSE cards since TDLDA is present')
        itdlda = 0
      endif
      if (ipmbse.eq.2 .and. nonlocal.gt.0 .and. ifxc.gt.0) then
c       accounting for core-hole twice. reset ifxc=0
         if (verbse)
     1     call wlog(' Reset ifxc=0 since core-hole potential is used.')
        ifxc = 0
      endif
      if (ipmbse.eq.1 .and. nonlocal.gt.0) then
c       V_ch should be zero
        nonlocal = 0
      endif

c     Josh - if nohole = 2, read wscrn.dat and add ch pot to vtot.
c     Need to add file check and emesh check.
      if (nohole.eq.2)  then
         open (unit=13, file='wscrn.dat', status='old', iostat=ios)
         call chopen (ios, 'wscrn.dat', 'ffmod2(xsph)')
         open (unit=14, file='vtot.dat', status='replace',iostat=ios)
         call chopen (ios, 'vtot.dat', 'ffmod2(xsph)')
         do i = 1, 251
            read(13,'(2e20.10)',end=20) dum1, dum2
            dum3 = vtot(i,0)
            vtot(i,0) = vtot(i,0) - dum2
            write(14,'(3e20.10)') dum1, dum3, dum2
         end do
 20      continue
         nohole = 0
         close(13)
         close(14)
      end if
c     Josh END

      if (itdlda.eq.0)  then 
!     Josh - Replaced call to phmesh with phmesh2, which allows user
!     defined grids read from grid.inp. Details can be found in phmesh2.f
!         call phmesh2 (ipr2, ispec, edge, emu, vi0, gamach, xkmax,
!     &        xkstep, vixan, ne, ne1, em, ik0, ne3,iGrid)
c        print *, ipr2, ispec, edge, emu, vi0, gamach
c        print *, xkmax, xkstep, vixan, ne, ne1, ik0, ne3
        call phmesh (ipr2, ispec, edge, emu, vi0, gamach,
     1                 xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3)
      else
c       nesvi TDLDA
c        call meshlda (xkstep, ne, ne1, em, ik0)
        corr = 1.0
      endif

      if (itdlda.eq.1) then
c       to get the mesh only
        do 93  i = 1, ne1
          write(3,94) dble(em(i))*hart
 94       format (7f10.5)
 93     continue
        stop 'TDLDA energy mesh is written out'
c       end of itdlda=1 calculations
      endif

c     Make old grid to report distances in xsect.bin
      do 95 i = 1, 251
 95   ri(i) = exp(-x0+dx*(i-1))

c     open xsect.bin and write the header
c--json--      open (unit=1, file='xsect.bin', status='unknown', iostat=ios)
c--json--      call chopen (ios, 'xsect.bin', 'potph')
c--json--      call wthead (1, ntitle, title)
c skip old output in title ( title lines are above ------ )
c     write(1,*) 'vtot in eV, rho in code units, includes 4pi'
c     write(1,*) 'ipot, vtot(imt), rho(imt) '
c     write(1,122) 'interstitial', vint*hart, rhoint
c     do 386  iph = 0, nph
c        write(1,123)iph,vtot(imt(iph),iph)*hart,edens(imt(iph),iph)
c 386 continue
c 122 format (1x, a, 1p, 2e20.6)
c 123 format (i10, 1p, 2e20.6)
c     write(1,42)  emu*hart
c  42 format ('       edge ', 2f20.5)
c     write(1,*)  imt(0), ' imt(0)'
c     write(1,200)  vint*hart, rhoint, ri(imt(0)+1)
c 200 format ('  v, rho, r', /, 1p, 3e20.4, ' intersitial')
c     do 220  iii = imt(0), imt(0)-4, -1
c        write(1,210)  vtot(iii,0)* hart, edens(iii,0), ri(iii), iii
c 210 format (1p, 3e20.4, i6)
c 220 continue

c--json--      write(1,45)
c--json--   45 format (1x, 71('-'))
c--json--      write(1,55) s02, erelax, wp, edge, emu
c--json--   55 format ( 3e13.5, 2e15.7, ' method to calculate xsect')
c--json--      write(1,56) gamach*hart, ne1, ik0
c--json--   56 format (1p, e15.7, 2i4,
c--json--     1       ' gamach in eV, # of points on horizintal axis')
c--json--      write(1,57)
c--json--   57 format ('        em              xsnorm            xsec  ')
c--json--c     end of the xsect.bin header


c     nsp = 1 - spin average caclulations; 2 - spin up and down
      nsp = 1
      if (abs(ispin).eq.1) nsp = nspx
c     scale spin density on each atom appropriately
      do iph = 0, nph
       do i = 1, 251
         dmag(i,iph) = dmag(i, iph) * spinph(iph)
       enddo
      enddo

      do 300 isp = 1, nsp
        if (ispin.ne.0) then
c         make spin dependent potential if needed
c         isp = 1-spin-down; 2-spin-up potentials
          idmag = (-1)**isp
          if (nsp.eq.1) then
             idmag = 1
             if(ispin.lt.0) idmag=-1
          endif
          call  istprm (nph, nat, iphat, rat, iatph, xnatph,
     1               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1               edens, edenvl, idmag,
     2               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3               ixc, rhoint,vint, rs, xf, xmu, xmunew,
     5               rnrmav, qtotel, inters, totvol)
          xmunew = xmu
          if (abs(ispin).eq.1 .and. nsp.eq.2) then
c           |ispinp| = |ispin|, but sign is determined by isp
             ispinp = abs(ispin)
             if (isp.eq.1) ispinp = -ispinp
          else
c            sign is determined by spin (always for ispin=-2,2)
             ispinp = ispin
          endif
        else
c         spin-independent case
          ispinp = 0
        endif

c       calculate operators of interest (s_z, l_z, t_z)
        xmuvr = xmu - vr0
        if (ipr2.ge.3) call szlz(verbse,
     1      ispinp,ecv,nph,nat,rgrd,nohole,rfms2,
     2      lfms2, lmaxph, edens, edenvl, dmag, vtot, vvalgs, rmt, rnrm,
     3      ixc, rhoint, vint, xmuvr, jumprm,
     4      xnval, iorb, x0, dx, xion, iunf, iz,
     5      adgc, adpc, dgc, dpc, ihole, rat, iphat, corr)
c    1                   em, ne1, ne, ik0 )

c       Cross section calculation, use phase mesh for now
c       Absorbing atom is iph=0
        if (verbse) then
           write(slog,10) 'absorption cross section'
           call wlog(slog)
        endif
        iph = 0
        call fixvar (rmt(0), edens(1,0), vtot(1,0), dmag(1,0),
     1             vint, rhoint, dx, dxnew, jumprm,
     2             vjump, ri, vtotph, rhoph, dmagx)
        call fixdsx (iph, dx, dxnew, dgc, dpc, dgcn, dpcn)
        if (mod(ixc,10) .ge. 5) then
           if (jumprm .gt. 0) jumprm = 2
           call fixvar (rmt(0), edenvl(1,0), vvalgs(1,0), dmag(1,0),
     1             vint, rhoint, dx, dxnew, jumprm,
     2             vjump, ri, vvalph, rhphvl, dmagx)
           if (jumprm .gt. 0) jumprm = 1
        endif
        call fixdsp (dx, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)
  
        if (itdlda.eq.0) then
c         Josh - added argument iPl to control many pole self energy
c         BR - removed eorb, gamach, ifxc, vi0
          call xsect (ipr2, dxnew, x0, ri, ne, ne1, ik0, em, edge,
     1       ihole, emu, corr, dgcx, dpcx, jnew,
     2       ixc0, lreal, rmt(0), rnrm(0), xmuvr, iPl,
     3       vtotph, vvalph, rhoph, dmagx, rhphvl, 
     4       dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), xsec(1,isp),
     5       xsnorm(1,isp), rkk(1,1,isp), iz(0), xion(0), iunf,
     6       xnval(1,iph), izstd, iorb(-4,iph), l2lp,
     7       ipol, ispinp, le2, angks,ptz)
        else
          if (nonlocal.gt.0) then
c           read potential with core-hole from a file
            if (nonlocal.eq.1) then
c              call rdpotp(vch)
            elseif (nonlocal.eq.2) then
c             open (unit=3, file='MgO_Mgk.dat', status='old')
              open (unit=3, file='wscrn.dat', status='old')
c             open (unit=3, file='w_m5p.dat', status='old')
c             open (unit=3, file='ni_l2.dat', status='old')
c             open (unit=3, file='ni_l2_sp.dat', status='old')
              n=0
 338          n = n+1
                read(3,337, end=339) dum1, dum2
c               use frac.ne.1  to mix bare and screened ch pot
c                frac = 0.80
c                frac = 1.00
                vch(n) = -1.d0*dum2
 337            format(6e20.10)
                goto 338
 339          continue
              close (unit=3)
            endif
            call fixvar (rmt(0), edens(1,0), vch, dmag(1,0),
     1             vint, rhoint, dx, dxnew, jumprm,
     2             vjump, ri, vchp, rhoph, dmagx)
            do 333 i = 1, nrptx
               if (ri(i).lt.rmt(0)) then
                 if (nonlocal.eq.1) then
                   vchp(i) = vchp(i) - vtotph(i)
                 endif
               elseif (ri(i).lt.40.d0) then
c                 assume const/r behaviour
                  vchp(i) = vchp(i-1) * ri(i-1) / ri(i) 
               else
                  vchp(i) = 0
               endif
c           testing: write core-hole potential in fort.17
               if (ri(i).lt.40.d0) write(17,332) ri(i), vchp(i)
 332            format(2f30.5)
 333        continue
           
            close (unit=17)
c           itest = 2
c           if (itest.eq.2) stop
          else
            do 334 i =1, nrptx
 334        vchp(i) = 0
          endif
          
c         Josh - added argument iPl to control many pole self energy
c          call xsectd (ipr2,dxnew, x0, ri, ne, ne1, ik0, em, edge,
c     1       ihole, emu, corr, dgcx, dpcx, jnew,
c     2       ixc0, lreal, rmt(0), rnrm(0), xmuvr, vi0, iPl,
c     3       gamach, vtotph, vvalph,vchp, rhoph, dmagx, rhphvl,
c     4       dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), xsec(1,isp),
c     5       xsnorm(1,isp), rkk(1,1,isp),iz(0), xion(0), iunf,
c     6       xnval(1,iph), ipmbse, ifxc, ibasis, eorb, kappa,
c     7       iorb(-4,iph), l2lp, ipol, ispinp, le2, angks,ptz, itdlda)
        endif


        do 60 iph = 0, nph
           if (verbse) then
              write(slog,10) 'phase shifts for unique potential', iph
              call wlog(slog)
           endif
c          fix up variable for phase
           call fixvar (rmt(iph), edens(1,iph), vtot(1,iph),
     1            dmag(1,iph), vint, rhoint, dx, dxnew, jumprm,
     2            vjump, ri, vtotph, rhoph, dmagx)
           if (mod(ixc,10) .ge.5) then
              if (jumprm .gt. 0) jumprm = 2
              call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph),
     1               dmag(1,iph), vint, rhoint, dx, dxnew, jumprm,
     2               vjump, ri, vvalph, rhphvl, dmagx)
              if (jumprm .gt. 0) jumprm = 1
              call fixdsx (iph, dx, dxnew, dgc, dpc, dgcn, dpcn)
           endif
           if (iph .eq. 0)  then
              itmp = ihole
           else
              itmp = 0
           endif

          call phase (iph, dxnew, x0, ri, ne, ne1, ne3, em, ixc, nsp,
     1            lreal, rmt(iph),rnrm(iph), xmuvr, iPl,
     2            vtotph, vvalph, rhoph, dmagx, rhphvl,
     3            dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), eref(1,isp),
     4            ph(1,-ltot,isp,iph), lmax(iph), iz(iph), itmp,
     5            xion(iph), iunf, xnval(1,iph), ispinp)
c          call phase (iph, dxnew, x0, ri, ne, ne1, ne3, em, ixc, nsp,
c     1            lreal, rmt(iph),rnrm(iph), xmuvr, vi0, iPl,
c     2            gamach, vtotph, vvalph, rhoph, dmagx, rhphvl,
c     3            dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), eref(1,isp),
c     4            ph(1,-ltot,isp,iph), lmax(iph), iz(iph), itmp,
c     5            xion(iph), iunf, xnval(1,iph), ispinp)
 60     continue

 300  continue

c     write main output to xsect.bin
c--json--  340 format (e17.9, 4e13.5)
      if (abs(ispin).ne.1 .or. nspx.eq.1) then
        do 350  ie = 1, ne
c--json--           write(1,340) dble(em(ie))*hart, dimag(em(ie))*hart,
c--json--     1                 xsnorm(ie,1), dble(xsec(ie,1)), dimag(xsec(ie,1))
           col1(ie) = dble(em(ie))*hart
           col2(ie) = dimag(em(ie))*hart
           col3(ie) = xsnorm(ie,1)
           col4(ie) = dble(xsec(ie,1))
           col5(ie) = dimag(xsec(ie,1))
  350   continue
      else
c       nspx = 2
        do 380  ie = 1, ne
c--json--           write(1,340) dble(em(ie))*hart, dimag(em(ie))*hart,
c--json--     1             (xsnorm(ie,1) + xsnorm(ie,nspx)) / 2.d0 ,
c--json--     2           dble( (xsec(ie,1) + xsec(ie,nspx)) ),
c--json--     3          dimag( (xsec(ie,1) + xsec(ie,nspx)) )
           col1(ie) = dble(em(ie))*hart
           col2(ie) = dimag(em(ie))*hart
           col3(ie) = (xsnorm(ie,1) + xsnorm(ie,nspx)) / 2.d0
           col4(ie) = dble( (xsec(ie,1) + xsec(ie,nspx)) )
           col5(ie) = dimag( (xsec(ie,1) + xsec(ie,nspx)) )
c          Normalize rkk to the average over up/down spin
c          nsp=2
           xnorm1 = sqrt( 2*xsnorm(ie,1) /
     1                         (xsnorm(ie,1) + xsnorm(ie,nspx)) )
           xnorm2 = sqrt( 2*xsnorm(ie,nspx) /
     1                         (xsnorm(ie,1) + xsnorm(ie,nspx)) )
           do 360 kdif = 1,8
             rkk (ie, kdif, 1) = rkk (ie, kdif, 1) * xnorm1
             rkk (ie, kdif, nspx) = rkk (ie, kdif, nspx) * xnorm2
  360      continue
  380   continue
      endif
c--json--      close (unit=1)

c     using opconsat resulted in NaNs in these array in one case
c     a.ne.a is an idiom for testing NaN-ness     
      do 390  ie = 1, ne
         if (col3(ie) .ne. col3(ie)) col3(ie) = 0.d0
         if (col4(ie) .ne. col4(ie)) col4(ie) = 0.d0
         if (col5(ie) .ne. col5(ie)) col5(ie) = 0.d0
 390  continue
      
      if (wrxsec) call json_xsect(ntitle, title, s02, erelax,
     -       wp, edge, emu,
     1       gamach*hart, ne, ne1, ik0,
     2       col1, col2, col3, col4, col5)

c     disable for now since dimensions are different
      if (ipr2 .ge. 2)  then
         call wphase (nph, em, eref, lmax, ne, ph, ntitle, title)
      endif

c     Write out phases for paths and genfmt
      call wrxsph (phpad,
     &       nsp, ne, ne1, ne3, nph, ihole, rnrmav, xmuvr,
     &       edge, ik0, ixc, rs, vint,
     &       em, eref, lmax, iz, potlbl, ph, rkk)

      if (ipr2 .ge. 1) then
c       calculate axafs
c       axafs does not make sense for |ispin| = 1
        call axafs (em, emu, xsec(1,1), ne1, ik0)
      endif

      return
      end
