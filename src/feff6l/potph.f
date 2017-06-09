      subroutine potph (isporb, shole, sout)

c     Cluster code -- multiple shell single scattering version of FEFF
c     This program (or subroutine) calculates potentials and phase
c     shifts for unique potentials specifed by atoms and overlap cards.
c
c     Input files:  potph.inp    input data, atoms, overlaps, etc.
c     Output:       phases.bin   phase shifts for use by the rest of the
c                                program
c                   xxx.dat      various diagnostics

      implicit double precision (a-h, o-z)

      include 'const.h'
      include 'dim.h'

      include 'arrays.h'
      common /print/ iprint

      parameter (nheadx = 30)
      character*80 head(nheadx), messag*128
      dimension lhead(nheadx)

      character*10 shole(0:9)
      character*8 sout(0:6)

c     head0 is header from potph.dat, include carriage control
      character*80 head0(nheadx)
      dimension lhead0(nheadx)

      dimension em(nex)
      dimension dgc0(251), dpc0(251)
      dimension xsec(nex), xsatan(nex)

c     nrx = max number of r points for phase r grid
      parameter (nrx = 250)
      dimension ri(nrptx), vtotph(nrx), rhoph(nrx)

      logical DO_pad_io
      DO_pad_io = .true.

   10 format (4x, a, i5)

c     Read input from file potph.inp
      open (unit=1, file='potph.dat', status='old', iostat=ios)
      call chopen (ios, 'potph.dat', 'potph')
      nhead0 = nheadx
      call rpotph (1, nhead0, head0, lhead0, nat, nph,
     1             nfr, ihole, gamach, iafolp, intclc,
     1             ixc, vr0, vi0, rs0, iphat, rat, iatph, ifrph, 
     1             xnatph, novr,
     2             iphovr, nnovr, rovr, folp, ion, iz, iprint, 
     2             ixanes, nemax, xkmin, xkmax, potlbl)
      close (unit=1)

c     Free atom potentials and densities
c     NB wsatom is needed in SUMAX, if changed here, change it there
      wsatom = 15
c     do not save spinors
      ispinr = 0
      do 20  ifr = 0, nfr
         itmp = 0
         if (ifr .eq. 0)  itmp = ihole
         write(messag,10)
     $        'free atom potential and density for atom type', ifr
         call echo(messag)
         call atom (head0(1)(1:40), ifr, iz(ifr), itmp, wsatom,
     1              ion(ifr), vcoul(1,ifr), rho(1,ifr),
     2              ispinr, dgc0, dpc0, et)
c        etfin is absorbing atom final state total energy
c        etinit is absorbing atom initial state (no hole)
         if (ifr .eq. 0)  etfin = et
   20 continue
      if (ixanes .gt. 0)  then
         call echo( '    initial state energy')
c        save spinor for core hole orbital
         ispinr = ihole
c        if no hole, use orbital from isporb
         if (ihole .eq. 0)  ispinr = isporb
         itmp = 0
         call atom (head0(1)(1:40), 0, iz(0), itmp, wsatom,
     1              ion(0), vcoul(1,nfr+1), rho(1,nfr+1),
     2              ispinr, dgc0, dpc0, etinit)
      endif
c     Need etfin if xanes and no hole, use K shell for this
      if (ixanes .gt. 0 .and. ihole .eq. 0)  then
c        K hole
         itmp = 1
         ispinr = 0
         call atom (head0(1)(1:40), 0, iz(0), itmp, wsatom,
     1              ion(0), vcoul(1,nfr+1), rho(1,nfr+1),
     2              ispinr, dgc0, dpc0, etfin)
      endif

c     Overlap potentials and densitites
      do 40  iph = 0, nph
         write(messag,10) 
     1    'overlapped potential and density for unique potential', iph
         call echo(messag)
         call ovrlp (iph, iphat, rat, iatph, ifrph, novr,
     1               iphovr, nnovr, rovr, iz, nat, rho, vcoul,
     2               edens, vclap, rnrm)
   40 continue

c     Find muffin tin radii, add gsxc to potentials, and find
c     interstitial parameters
       call echo('    muffin tin radii and interstitial parameters')
      call istprm (nph, nat, iphat, rat, iatph, xnatph,
     1             novr, iphovr, nnovr, rovr, folp, edens,
     2             vclap, vtot, imt, inrm, rmt, rnrm, rhoint,
     3             vint, rs, xf, xmu, rnrmav, intclc)

c     Automatic max reasonable overlap
      if (iafolp .eq. 1)  then
         call echo('    automatic overlapping')
         call echo('   iph, rnrm(iph)*bohr,'//
     $        ' rmt(iph)*bohr, folp(iph)')
         do 400  iph = 0, nph
            folp(iph) = 1 + 0.7*(rnrm(iph)/rmt(iph) - 1)
            write(messag,'(3x,i8,3g15.6)') iph,
     $           rnrm(iph)*bohr, rmt(iph)*bohr, folp(iph)
            call echo(messag)
 400     continue
         call istprm (nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr, folp, edens,
     2                vclap, vtot, imt, inrm, rmt, rnrm, rhoint,
     3                vint, rs, xf, xmu, rnrmav, intclc)
      endif

c     Initialize header routine and write misc.dat
      call sthead (nhead0, head0, lhead0, nph, iz, rmt, rnrm,
     1             ion, ifrph, ihole, ixc,
     2             vr0, vi0, rs0, gamach, xmu, xf, vint, rs,
     3             nhead, lhead, head, shole, sout)
      if (iprint .ge. 1)  then
         open (unit=1, file='misc.dat', status='unknown', iostat=ios)
         call chopen (ios, 'misc.dat', 'potph')
         call wthead(1)
         close (unit=1)
      endif

      if (iprint .ge. 2)  then
         call wpot (nph, edens, ifrph, imt, inrm,
     1              rho, vclap, vcoul, vtot)
      endif

c     Phase shift calculation
c     Make energy mesh and position grid
      nr = 250
      dx = .05d0
      x0 = 8.8d0
      edge = xmu - vr0
cc      print*,  'potph -> phmesh ', ne

      call phmesh (nr, dx, x0, nemax, iprint,
     1             ixanes, edge, xmu, vint, vr0,
     1             imt, edens, nph,
     2             ri, ne, em, ik0)

cc      print*,  'potph after phmesh ', ne
c     Cross section calculation, use phase mesh for now
c     remove xanes calculation in feff6l

      do 60  iph = 0, nph
         write(messag,10) 'phase shifts for unique potential', iph
         call echo(messag)
c        fix up variable for phase
         call fixvar (rmt(iph), edens(1,iph), vtot(1,iph),
     1                vint, rhoint, nr, dx, x0, ri,
     2                vtotph, rhoph)

cc         print*,  'potph -> phase ', ne
         call phase (iph, nr, dx, x0, ri, ne, em, edge,
     1               ixc, rmt(iph), xmu, vi0, rs0, gamach,
     2               vtotph, rhoph,
     3               eref, ph(1,1,iph), lmax(iph))
cc         print*,  'potph after phase ', ne
   60 continue

      if (iprint .ge. 2)  then
         call wphase (nph, em, eref, lmax, ne, ph)
      endif

c     Write out phases for genfmt
c     May need stuff for use with headers only
      open (unit=1, file='phase.bin', access='sequential',
     1      form='unformatted', status='unknown', iostat=ios)
      call chopen (ios, 'phase.bin', 'potph')
      write(1) nhead
      do 62  i = 1, nhead
         write(1) head(i)
         write(1) lhead(i)
   62 continue
      write(1) ne, nph, ihole, rnrmav, xmu, edge, ik0
      write(1) (em(ie),ie=1,ne)
      write(1) (eref(ie),ie=1,ne)
      do 80  iph = 0, nph
         write(1) lmax(iph), iz(ifrph(iph))
         write(1) potlbl(iph)
         do 70  ie = 1, ne
            write(1)  (ph(ie,ll,iph), ll=1,lmax(iph)+1)
   70    continue
   80 continue
      close (unit=1)

cc
cc optionally, write phase.pad 
cc

      if (do_pad_io) then
         lun  = 9
         npack = 10
         call openfl(lun, 'phase.pad',  'unknown', iex, ier)
         if ((ier.lt.0).or.(iex.lt.0)) then
            call echo(' *** Error: cannot open Potentials.bin')
            return
         end if
         write(lun,'(a,i3)') '#:FEFF6X POT File: npad = ', npack
         write(lun,'(a,i9,i9,i9,i9)') '#:ne,nph,ihole,ik0 =  ',
     $        ne, nph, ihole, ik0
         write(lun, '(a,g22.15)') '#% rnrmav = ', rnrmav
         write(lun, '(a,g22.15)') '#% xmu    = ', xmu
         write(lun, '(a,g22.15)') '#% edge   = ', edge
         call wrpadd(lun,npack,em,ne)
         call wrpadx(lun,npack,eref,ne)
         do 420  iph = 0, nph
            write(lun, '(2a,3i9)') '#:label,iph,lmax,iz  ', potlbl(iph),
     $           iph, lmax(iph), iz(ifrph(iph))
            do 410  l = 1, lmax(iph)+1
               call wrpadx(lun,npack,ph(1,l,iph),ne)
 410        continue
 420     continue 
         close(lun)
      endif



      return
      end
