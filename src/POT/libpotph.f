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
c      integer  iatph(0:nphx), ibounc(natx)
      character*80 title(nheadx)
c     character*80 head(nheadx)
c     integer lhead(nheadx)
      integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, 
     1       iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx), lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
c     for OVERLAP option
c      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
c      double precision  rovr(novrx,0:nphx)


c     dimension/type os mod2/xpsh things
      integer mphase, ipr2, ixc0, ispec, lreal, lfms2, l2lp, iPl, 
     1       iGrid
      double precision xkstep, xkmax, vixan
      double precision vr0, vi0, spinph(0:nphx)
      real rfms2
      integer lmaxph(0:nphx)
      character*6  potlbl(0:nphx)
      integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis



      call inipotph(nat, rat, iphat,
     1       le2, elpty, angks, evec, xivec, spvec, ptz,
     2       nabs, iphabs, rclabs, ipol, ispin,
     3       mpot, rgrd, ntitle, title, ipr1, ispec,
     4       nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     5       xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     6       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     7       mphase, ipr2, vixan, xkstep, xkmax,
     8       lmaxph, potlbl, spinph, vr0, vi0, ixc0, lreal, 
     9       rfms2, lfms2, l2lp, iPl, iGrid,
     _       izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)


c*****************************************************************************
c     read the contents of a json file that includes all of global.json,
c     atoms.json, pot.json & xpsh.json (i.e. global.dat, atoms.dat, mod1.inp,
c     and mod2.inp)
c*****************************************************************************
      call json_read_libpotph(nat, rat, iphat,
     1       le2, elpty, angks, evec, xivec, spvec, ptz,
     2       nabs, iphabs, rclabs, ipol, ispin,
     3       mpot, rgrd, ntitle, title, ipr1, ispec,
     4       nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     5       xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     6       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     7       mphase, ipr2, vixan, xkstep, xkmax,
     8       lmaxph, potlbl, spinph, vr0, vi0, ixc0, lreal, 
     9       rfms2, lfms2, l2lp, iPl, iGrid,
     _       izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)

c     iabs != 0 has something to do with CFAVERAGE, outside scope of feff85exafs
      iabs = 1
      call ffsort(iabs, natt, ratx, iphatx,
     1       nabs, iphabs, rclabs, ipol, ispin, le2,
     2       elpty, angks, evec, xivec, spvec, ptz,
     3       iatph)
c     iatph,rat,iphat

      print *, rat


      stop
      end
