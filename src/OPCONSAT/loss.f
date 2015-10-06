      program loss

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'

      integer epsmax
      parameter(epsmax = 700)
      
c      integer iz(0:nphx)
      integer npoles
      double precision rnrm(0:nphx), eps0, gamma
      logical write_loss, write_opcons, write_exc, verbose

      double precision wpcorr(MxPole), delta(MxPole), ampfac(MxPole)

      character*80 title(nheadx)
      integer ntitle, nat, nph, iphat(natx), ipol, ispin, ihole
      integer lfms1, nmix, icoul, ixc, ixc0, iafolp, iunf, inters
      integer jumprm, nohole, iplsmn
      integer iz(0:nphx), lmaxsc(0:nphx), lmaxph(0:nphx)
      double precision evec(3), xivec(3), spvec(3), spinph(0:nphx)
      complex*16 ptz(-1:1, -1:1)
      double precision rat(3,natx)
      double precision xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
      character*6 potlbl(0:nphx)
      real rfms1
      double precision elpty, angks, gamach, ca1, ecv
      double precision vr0, vi0, rgrd, totvol


      
      write_loss   = .true.
      write_opcons = .false.
      write_exc    = .false.
      verbose      = .true.
      npoles       = 100
      eps0         = -1.d0



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
     1       inters, totvol, jumprm, nohole, iplsmn)


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
     1       inters, totvol, jumprm, nohole, iplsmn)

      rnrm(0) = 2.8384890981392523
      rnrm(1) = 2.6294479894989911

      call feffloss(nph, iz, xnatph, rnrm, npoles, eps0,
     1       write_opcons, write_loss, write_exc, verbose,
     2       wpcorr, gamma, ampfac, delta)
      

      end
