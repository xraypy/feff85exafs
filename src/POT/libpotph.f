      program libpotph
      
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'


c     dimension/types of atoms & global.json things
      integer nat, iphat(natx), nabs, iphabs, ipol, ispin, le2, iabs
      double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
      complex*16 ptz(-1:1, -1:1)
      double precision rat(3,natx)


c     dimension/type os mod1/pot things
      integer  iatph(0:nphx), ibounc(natx)
      character*80 title(nheadx)
c     character*80 head(nheadx)
c     integer lhead(nheadx)
      integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, 
     1       iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx), lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)

c     dimension/type os mod2/xpsh things
      integer mphase, ipr2, ixc0, ispec, lreal, lfms2, l2lp, iPl, 
     1       iGrid
      double precision xkstep, xkmax, vixan
      double precision vr0, vi0, spinph(0:nphx)
      real rfms2
      integer lmaxph(0:nphx)
      character*6  potlbl(0:nphx)
      integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis

c     for OVERLAP option -- DISABLED IN FEFF85EXAFS
      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
      double precision  rovr(novrx,0:nphx)

      parameter (big = 1.0e5)


c$$$      call inipotph(nat, rat, iphat,
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


      call inipotph(
c     TITLE
     1       ntitle, title,
c     ATOMS
     2       nat, rat, iphat,
c     POTENTIALS
     3       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
c     HOLE/EDGE
     4       ihole,
c     SCF
     5       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c     POLARIZATION, ELLIPTICITY
     6       ipol, evec, elpty, xivec,
c     SPIN
     7       ispin, spvec, angks,
c     computed
     8       ptz, gamach,
c     EXCHANGE
     9       ixc, vr0, vi0, ixc0,
c     AFOLP, FOLP, ION, RGRID, UNFREEZEF
     _       iafolp, folp, xion, rgrd, iunf,
c     INTERSTITIAL, JUMPRM, NOHOLE
     1       inters, totvol, jumprm, nohole)


c*****************************************************************************
c     the following parameters are for features not present or not used
c     in feff85exafs
c     they are all set to default/turned-off values

c     CONTROL and PRINT
      mpot   = 1
      mphase = 1
      ipr1   = 1
      ipr2   = 0
c     0=EXAFS, >0 other spectroscopies
      ispec  = 0
c     CFAVERAGE
      nabs   = 1
      iphabs = 0
      rclabs = big
c     MULTIPOLE
      le2    = 0
      l2lp   = 0
c     FMS
      rfms2  = dble(-1.)
      lfms2  = 0
c     XANES
      vixan  = dble(0)
      xkstep = 0
      xkmax  = 0
c     RPHASES
      lreal  = 0
c     PMBSE
      ifxc   = 0
      ipmbse = 0
      nonlocal = 0
      ibasis = 0
c     TDLDA
      izstd  = 0
      itdlda = 0
c     PLASMON
      iPl    = 0
      iGrid  = 0

c     OVERLAP
      do 20 i=0,nphx
         novr(i) = 0;
         do 10 j=1,novrx
            iphovr(j,i) = 0
            nnovr(j,i)  = 0
            rovr(j,i)   = dble(0.)
 10      continue
 20   continue
      do 30 i=1,natx
         ibounc(i) = 1
 30   continue
c*****************************************************************************


c*****************************************************************************
c     read the contents of a json file that includes all of global.json,
c     atoms.json, pot.json & xpsh.json (i.e. global.dat, atoms.dat, mod1.inp,
c     and mod2.inp)
c*****************************************************************************
      call json_read_libpotph(
c     TITLE
     1       ntitle, title,
c     ATOMS
     2       nat, rat, iphat,
c     POTENTIALS
     3       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
c     HOLE/EDGE
     4       ihole,
c     SCF
     5       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c     POLARIZATION, ELLIPTICITY
     6       ipol, evec, elpty, xivec,
c     SPIN
     7       ispin, spvec, angks,
c     computed
     8       ptz, gamach,
c     EXCHANGE
     9       ixc, vr0, vi0, ixc0,
c     AFOLP, FOLP, ION, RGRID, UNFREEZEF
     _       iafolp, folp, xion, rgrd, iunf,
c     INTERSTITIAL, JUMPRM, NOHOLE
     1       inters, totvol, jumprm, nohole)


c$$$             nat, rat, iphat,
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

c     iabs != 0 has something to do with CFAVERAGE, outside scope of feff85exafs
      iabs = 1

c$$$      print *, iabs
c$$$      print *, nat
c$$$      print *, nabs
c$$$      print *, iphabs
c$$$      print *, rclabs
c$$$      print *, ipol
c$$$      print *, ispin
c$$$      print *, le2
c$$$      print *, elpty
c$$$      print *, angks

      call ffsort(iabs, nat, rat, iphat,
     1       nabs, iphabs, rclabs, ipol, ispin, le2,
     2       elpty, angks, evec, xivec, spvec, ptz,
     3       iatph)
c     iatph,rat,iphat



      call pot(rgrd, nohole,
     $       inters, totvol, ecv, nscmt, nmix, ntitle, title,
     $       nat, nph, ihole, iafolp,
     $       ixc, iphat, rat, iatph, xnatph,
     $       novr, iphovr, nnovr, rovr,
     $       folp, xion, iunf, iz, ipr1,
     $       ispec, jumprm,
     $       lmaxsc, icoul, ca1, rfms1, lfms1)


      stop
      end
