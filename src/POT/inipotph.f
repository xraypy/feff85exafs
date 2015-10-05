      subroutine inipotph(
     1       ntitle, title,                                             ! TITLE
     2       nat, rat, iphat,                                           ! ATOMS
     3       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,           ! POTENTIALS
     4       ihole,                                                     ! HOLE/EDGE
     5       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,                ! SCF
     6       ipol, evec, elpty, xivec,                                  ! POLARIZATION, ELLIPTICITY
     7       ispin, spvec, angks,                                       ! SPIN
     8       ptz, gamach,                                               ! computed
     9       ixc, vr0, vi0, ixc0,                                       ! EXCHANGE
     _       iafolp, folp, xion, rgrd, iunf,                            ! AFOLP, FOLP, ION, RGRID, UNFREEZEF
     1       inters, totvol, jumprm, nohole, iplsmn)                    ! INTERSTITIAL, JUMPRM, NOHOLE, PLASMON

c$$$(nat, rat, iphat,
c$$$     1       le2, elpty, angks, evec, xivec, spvec, ptz,
c$$$     2       nabs, iphabs, rclabs, ipol, ispin,
c$$$     3       mpot, rgrd, ntitle, title, ipr1, ispec,
c$$$     4       nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
c$$$     5       xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
c$$$     6       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c$$$     7       mphase, ipr2, vixan, xkstep, xkmax,
c$$$     8       lmaxph, potlbl, spinph, vr0, vi0, ixc0, lreal, 
c$$$     9       rfms2, lfms2, l2lp, iPl, iGrid,
c$$$     _       izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)


      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'


c     dimension/types of atoms & global.json things
      integer nat, iphat(natx), nabs, iphabs, ipol, ispin, le2
      double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
      complex*16 ptz(-1:1, -1:1)
      double precision rat(3,natx)


c     dimension/type os mod1/pot things
c      integer  iatph(0:nphx), ibounc(natx)
      character*80 title(nheadx)
c     character*80 head(nheadx)
c     integer lhead(nheadx)
      integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, iplsmn,
     1       iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx), lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
c     for OVERLAP option
c      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
c      double precision  rovr(novrx,0:nphx)


c     dimension/type os mod2/xpsh things
      integer mphase, ipr2, ixc0, ispec, lreal, lfms2, l2lp, 
     1       iGrid
      double precision xkstep, xkmax, vixan
      double precision vr0, vi0, spinph(0:nphx)
      real rfms2
      integer lmaxph(0:nphx)
      character*6  potlbl(0:nphx)
      integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis

      parameter (big = 1.0d5)


c*****************************************************************************
c     initialize everything

      nat      = 0
      nabs     = 0
      iphabs   = 0
      ipol     = 0 
      ispin    = 0
      le2      = 0
      mpot     = 1
      nph      = 0
      ntitle   = 0
      ihole    = 0
      ipr1     = 0
      iafolp   = 0
      ixc      = 0
      iunf     = 0
      nmix     = 0
      nohole   = -1
      jumprm   = 0
      inters   = 0
      nscmt    = 0
      icoul    = 0
      lfms1    = 0
      mphase   = 0
      ipr2     = 0
      ixc0     = -1
      ispec    = 0
      lreal    = 0
      lfms2    = 0
      l2lp     = 0
      iplsmn   = 0
      iGrid    = 0
      izstd    = 0
      ifxc     = 0
      ipmbse   = 0
      itdlda   = 0
      nonlocal = 0
      ibasis   = 0

      elpty  = 0.0d0
      angks  = 0.0d0
      rclabs = big
      gamach = 0.0d0
      rgrd   = 0.05d0
      ca1    = 0.0d0
      ecv    = -40.0d0
      totvol = 0.0d0
      xkmax  = 20.0d0 * bohr
      xkstep = 0.07d0 * bohr
      vixan  = 0.0d0 / hart
      vr0    = 0.0d0 / hart
      vi0    = 0.0d0 / hart

      rfms1  = real(-1.) / real(bohr)
      rfms2  = real(-1.) / real(bohr)

      do 10 i=0,nphx
         lmaxph(i) = 0
         iz(i)     = 0
         lmaxsc(i) = 0
         potlbl(i) = ' '
         xnatph(i) = 0.0d0
         spinph(i) = 0.0d0
         folp(i)   = 1.0d0
         xion(i)   = 0.0d0
 10   continue
      do 20 i=1,natx
         iphat(i) = -1
 20   continue

      do 35 i=-1,1
         do 30 j=-1,1
            ptz(i,j) = dcmplx(0.0d0,0.0d0)
 30      continue
 35   continue

      do 40 i=1,3
         evec(i)  = 0.0d0
         xivec(i) = 0.0d0
         spvec(i) = 0.0d0
 40   continue

      do 55 i=1,3
         do 50 j=1,natx
            rat(i,j) = 0.0d0
 50      continue
 55   continue

      do 60 i=1,nheadx
         title(i)  = ' '
 60   continue

      return
      end
