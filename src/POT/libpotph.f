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


c**********************************************************************
c     these parameters pass stuff from pot to wrpot and xsph
      double precision rnrmav, xmu, vint, rhoint, emu, s02, erelax
      double precision wp, rs, xf, qtotel
c     r mesh index just inside rmt
      dimension imt(0:nphx)
c     r mesh index just inside rnorman
      dimension inrm(0:nphx)
c     muffin tin radius
      dimension rmt(0:nphx)
c     norman radius
      dimension rnrm(0:nphx), qnrm(0:nphx)
c     folp(0:nphx)  - overlap factor for rmt calculation
      dimension folpx(0:nphx)
c     need irregular solution for complex potential. fix later
      dimension dgc0(251), dpc0(251)
c     additional data needed for relativistic version
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
c     Overlap calculation results
c     overlapped density*4*pi
      dimension edens(251,0:nphx)
c     overlapped coul pot
      dimension vclap(251,0:nphx)
c     overlapped total potential
      dimension vtot (251,0:nphx)
      dimension edenvl(251,0:nphx)
      dimension vvalgs (251,0:nphx)
      dimension dmag(251,0:nphx+1)
      dimension xnval(30,0:nphx+1)
      dimension eorb(30,0:nphx+1)
      dimension kappa(30,0:nphx+1), iorb(-4:3,0:nphx+1)
      dimension xnmues(0:lx,0:nphx)
c     Josh use nhtmp to save nohole value
      integer nhtmp



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

c     CONTROL and PRINT (ipr1=1 for misc.dat)
      mpot     = 1
      mphase   = 1
      ipr1     = 0
      ipr2     = 0
c     0        =EXAFS, >0 other spectroscopies
      ispec    = 0
c     CFAVERAGE
      nabs     = 1
      iphabs   = 0
      rclabs   = big
c     MULTIPOLE
      le2      = 0
      l2lp     = 0
c     FMS
      rfms2    = dble(-1.)
      lfms2    = 0
c     XANES
      vixan    = dble(0) / hart
      xkstep   = dble(0.07) * bohr
      xkmax    = dble(20) * bohr
c     RPHASES
      lreal    = 0
c     PMBSE
      ifxc     = 0
      ipmbse   = 0
      nonlocal = 0
      ibasis   = 0
c     TDLDA
      izstd    = 0
      itdlda   = 0
c     PLASMON
      iPl      = 0
      iGrid    = 0

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


c     iabs != 0 has something to do with CFAVERAGE, outside scope of feff85exafs
      iabs = 1


      call ffsort(iabs, nat, rat, iphat,
     1       nabs, iphabs, rclabs, ipol, ispin, le2,
     2       elpty, angks, evec, xivec, spvec, ptz,
     3       iatph)


      call pot(rgrd, nohole,
     $       inters, totvol, ecv, nscmt, nmix, ntitle, title,
     $       nat, nph, ihole, iafolp, ixc, iphat, rat, iatph, xnatph,
     $       novr, iphovr, nnovr, rovr, folp, xion, iunf, iz, ipr1,
     $       ispec, jumprm, lmaxsc, icoul, ca1, rfms1, lfms1,
c        return stuff for passing to xsph and skipping pot.pad
     -       rnrmav, xmu, vint, rhoint,
     1       emu, s02, erelax, wp, rs, xf, qtotel,
     2       imt, rmt, inrm, rnrm, folp, folpx,
     3       dgc0, dpc0, dgc, dpc, adgc, adpc,
     4       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     5       eorb, kappa, iorb, qnrm, xnmues, nhtmp
     6       )

c     could make a conditional call to wrpot here

      call xsph(ipr2, ispec, vixan, xkstep, xkmax, gamach, rgrd,
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

      stop
      end
