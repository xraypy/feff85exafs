      program makepotph

      implicit double precision (a-h, o-z)

c      include '../HEADERS/dim.h'
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      integer nphx
      parameter (nphx = 11)
c      max number of atoms in problem for the pathfinder
      integer natx
      parameter (natx =1000)
c      max number of header lines
      integer nheadx
      parameter (nheadx=30)

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
c     INTERSTITIAL, JUMPRM, NOHOLE, PLASMON
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
c     INTERSTITIAL, JUMPRM, NOHOLE, PLASMON
     1       inters, totvol, jumprm, nohole, iplsmn)


      call libpotph('phase.pad', .true.,
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
     1       inters, totvol, jumprm, nohole,
c     PLASMON
     2       iplsmn)
      

      stop
      end
