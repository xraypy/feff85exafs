      subroutine libpotph(phpad, verbse,
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
     1       inters, totvol, jumprm, nohole,                            ! INTERSTITIAL, JUMPRM, NOHOLE
     2       iplsmn)                                                    ! PLASMON
      
c##****f* feff85exafs/libpotph
c##  NAME
c##     libpotph
c##
c##  SYNOPSIS
c##     libpotph(ntitle, title, nat, rat, iphat, nph, iz, potlbl, lmaxsc,
c##    $         lmaxph, xnatph, spinph, ihole, rfms1, lfms1, nscmt, ca1,
c##    $         nmix, ecv, icoul, ipol, evec, elpty, xivec, ispin, spvec,
c##    $         angks, ptz, gamach, ixc, vr0, vi0, ixc0, iafolp, folp,
c##    $         xion, rgrd, iunf, inters, totvol, jumprm, nohole)
c##
c##  FUNCTION
c##     compute potentials and phases from an input atoms list
c##
c##  INPUTS
c##    phpad         path/name of output phases.pad file
c##    verbse        boolean flag, true=write screen messages
c##    ntitle        number of header lines                               TITLE
c##    titles        (nheadx) array of header string                      TITLE
c##    nat           number of atoms in cluster                           ATOMS
c##    rat           (3,natx) cartesian coordinates of atoms in cluster   ATOMS
c##    iphat         (natx) unique potential indeces of atoms in cluster  ATOMS
c##    nph           number of unique potentials                          POTENTIALS
c##    iz            (0:nphx) Z numbers of unique potentials              POTENTIALS
c##    potlbl        (0:nphx) labels of unique potentials                 POTENTIALS
c##    lmaxsc        (0:nphx) l max for SCF for each potential            POTENTIALS
c##    lmaxph        (0:nphx) l max for FMS for each potential            POTENTIALS
c##    xnatph        (0:nphx) stoichiometry of each potential             POTENTIALS
c##    spinph        (0:nphx) spin on each unique potential               POTENTIALS
c##    ihole         edge index, 1=K, 4=L3, etc                           HOLE/EDGE
c##    rscf          cluster radius for self-consistent calculation       SCF
c##    lscf          0=solid, 1=molecule                                  SCF
c##    nscmt         max number of self-consistency iterations            SCF
c##    ca            self-consistency convergence accelerator             SCF
c##    nmix          number of mixing iterations before Broyden           SCF
c##    ecv           core/valence separation energy                       SCF
c##    icoul         obsolete param. for handling Coulomb potential       SCF
c##    ipol          1=do polarization calculation                        POLARIZATION
c##    evec          (3) polarization array                               POLARIZATION
c##    elpty         eccentricity of elliptical light                     ELLIPTICITY
c##    xivec         (3) ellipticity array                                ELLIPTICITY
c##    ispin         +/-2 = do spin calculation                           SPIN
c##    spvec         (3) spin array                                       SPIN
c##    angks         angle between spin and incidient beam                SPIN
c##    ixc           exchange index                                       EXCHANGE
c##    vr0           Fermi level offset                                   EXCHANGE
c##    vi0           constant broadening                                  EXCHANGE
c##    ixc0          exchange index for background function               EXCHANGE
c##    iafolp        1=do automated overlapping                           FOLP & AFOLP
c##    folp          (0:nphx) overlapping fractions                       FOLP & AFOLP
c##    xion          (0:nphx) potential ionizations                       ION
c##    rgrd          radial grid used for the potentials/phases           RGRID
c##    iunf          1=unfreeze f electrons                               UNFREEZEF
c##    inters                                                             INTERSTITIAL
c##    totvol                                                             INTERSTITIAL
c##    jumprm        1=remove potential jumps at muffin tin radii         JUMPRM
c##    nohole        1=compute without core-hole                          NOHOLE
c##    iplsmn        1=compute multi-pole self-energy                     PLASMON
c##
c##  RESULT
c##    ptz           (-1:1,-1:1) polarization tensor
c##    gamach        tabulated core-hole lifetime
c##
c##  BUGS
c##    Report bugs and other issues at https://github.com/xraypy/feff85exafs
c##
c##  LICENSE
c##    See src/HEADERS/license.h for the terms of the parts of feff85exafs derived directly from Feff
c##
c##    nxjson.c and nxjson.h are Copyright (c) 2013 Yaroslav Stavnichiy <yarosla@gmail.com>.  See
c##    https://bitbucket.org/yarosla/nxjson/src
c##
c##    The C wrapper around this subroutine is released to the public domain
c##
c##  SEE ALSO
c##    The pot and xsph stand-alone programs, the libfeffphases.c wrapper
c##
c##****


      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'


c     dimension/types of atoms & global.json things
      integer nat, iphat(natx), nabs, iphabs, ipol, ispin, le2, iabs
      double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
      complex*16 ptz(-1:1, -1:1)
      double precision rat(3,natx)

c     name of output phase.pad file
      character*256 phpad
      logical verbse

c     dimension/type os mod1/pot things
      integer  iatph(0:nphx)
c      integer ibounc(natx)
      character*80 title(nheadx)
c      character*80 head(nheadx)
c      integer lhead(nheadx)
c      integer mpot, mphase
      integer nph, ntitle, ihole, ipr1, iafolp, ixc, 
     1       iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx), lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      dimension xnatph(0:nphx), folp(0:nphx), xion(0:nphx)

c     dimension/type os mod2/xpsh things
      integer ipr2, ixc0, ispec, lreal, lfms2, l2lp, iplsmn
c      integer iGrid
      double precision xkstep, xkmax, vixan
      double precision vr0, vi0, spinph(0:nphx)
      real rfms2
      integer lmaxph(0:nphx)
      character*6  potlbl(0:nphx)
      integer izstd, ifxc, ipmbse, itdlda, nonlocal

c     for OVERLAP option -- DISABLED IN FEFF85EXAFS
      dimension novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)

      parameter (big = 1.0d5)


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


      logical write_loss, write_opcons, write_exc, verbose
      integer npoles
      double precision eps0
      double precision wpcorr(MxPole), delta(MxPole), ampfac(MxPole)


      
c      print *, "libpotph: >", phpad(1:istrln(phpad)), "<"
c*****************************************************************************
c     the following parameters are for features not present or not used
c     in feff85exafs
c     they are all set to default/turned-off values

c     CONTROL and PRINT (ipr1=1 for misc.dat)
c      mpot     = 1
c      mphase   = 1
      ipr1     = 0
      ipr2     = 0
c     ispec=0 for EXAFS, >0 other spectroscopies
      ispec    = 0
c     CFAVERAGE
      nabs     = 1
      iphabs   = 0
      rclabs   = big
c     MULTIPOLE
      le2      = 0
      l2lp     = 0
c     FMS
      rfms2    = real(-1.)
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
c      ibasis   = 0
c     TDLDA
      izstd    = 0
      itdlda   = 0
c     PLASMON
c      iGrid    = 0

c     OVERLAP
      do 20 i=0,nphx
         novr(i) = 0;
         do 10 j=1,novrx
            iphovr(j,i) = 0
            nnovr(j,i)  = 0
            rovr(j,i)   = dble(0.)
 10      continue
 20   continue
c      do 30 i=1,natx
c         ibounc(i) = 1
c     30   continue

c     multi-pole stuff      
      write_loss   = .false.
      write_opcons = .false.
      write_exc    = .false.
      verbose      = .false.
      npoles       = 100
      eps0         = -1.d0
      do 40 i=1,MxPole
         wpcorr(i) = 0.d0
c        gamma ?
         delta(i)  = 0.d0
         ampfac(i) = 0.d0
 40   continue
      
      
c*****************************************************************************

c$$$         print *, 1, rat(1,1), rat(2,1), rat(3,1), iphat(1)
c$$$         print *, 2, rat(1,2), rat(2,2), rat(3,2), iphat(2)
c$$$         print *, 3, rat(1,3), rat(2,3), rat(3,3), iphat(3)
c$$$         print *, 13, rat(1,13), rat(2,13), rat(3,13), iphat(13)
c$$$         print *, 23, rat(1,23), rat(2,23), rat(3,23), iphat(23)
c$$$         print *, 177, rat(1,177), rat(2,177), rat(3,177), iphat(177)
c$$$


c     iabs != 0 has something to do with CFAVERAGE, outside scope of feff85exafs
      iabs = 1

      call ffsort(iabs, nat, rat, iphat,
     1       nabs, iphabs, rclabs, ipol, ispin, le2,
     2       elpty, angks, evec, xivec, spvec, ptz,
     3       iatph)

      call pot(verbse, rgrd, nohole,
     $       inters, totvol, ecv, nscmt, nmix, ntitle, title,
     $       nat, nph, ihole, iafolp, ixc, iphat, rat, iatph, xnatph,
     $       novr, iphovr, nnovr, rovr, folp, xion, iunf, iz, ipr1,
     $       ispec, jumprm, lmaxsc, icoul, ca1, rfms1, lfms1,
     $
     -       rnrmav, xmu, vint, rhoint,
     1       emu, s02, erelax, wp, rs, xf, qtotel,
     2       imt, rmt, inrm, rnrm, folpx,
     3       dgc0, dpc0, dgc, dpc, adgc, adpc,
     4       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     5       eorb, kappa, iorb, qnrm, xnmues, nhtmp
     6       )
c     return stuff for passing to xsph and skipping pot.pad

c     could make a conditional call to wrpot here

      if (iplsmn .gt. 0) then
c         call feffloss(nph, iz, xnatph, rnrm, npoles, eps0,
c     1          write_opcons, write_loss, write_exc, verbose,
c     2          wpcorr, gamma, ampfac, delta)
      end if
c     need to send <<wpcorr, gamma, ampfac, delta>> into xsph, then into
c     phase and xsect for use when iplsmn is > 0

c               wrxsec      
      call xsph(.false., verbse, phpad,
     -       ipr2, ispec, vixan, xkstep, xkmax, gamach, rgrd,
     1       nph, lmaxph, potlbl, spinph, iatph, nat, rat, iphat,
     2       ixc, vr0, vi0, ixc0, lreal, rfms2, lfms2, l2lp,
     3       ipol, ispin, le2, angks, ptz, iplsmn,
     4       izstd, ifxc, ipmbse, itdlda, nonlocal,
     _
     1       ntitle, title, rnrmav, xmu, vint, rhoint,
     2       emu, s02, erelax, wp, ecv, rs, xf, qtotel,
     3       imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     4       dgc0, dpc0, dgc, dpc, adgc, adpc,
     5       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     6       iorb, nohole, ihole,
     7       inters, totvol, iafolp, xion, iunf, iz, jumprm)
c     second block is the parameters from rdpot, skipping pot.pad

      return
      end
