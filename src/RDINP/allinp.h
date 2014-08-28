c -*- fortran -*-
c     Common blocks with all input data
c     the common

c---- atoms.dat -------------------------------------------------------------
      integer  natt
      integer iphatx(nattx)
      double precision  ratx(3,nattx)
      common /geom/ ratx, iphatx, natt

c---- geom.dat -------------------------------------------------------------
c       integer  nat
c       integer iatph(0:nphx)
c       integer iphat(natx)
c       double precision  rat(3,natx)
c       common /geom/ ratx, iphatx, natt

c---- global.inp -------------------------------------------------------------
c       configuration average
      integer iphabs
c     global polarization data
      integer  ipol, ispin, le2
      double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
      complex*16 ptz(-1:1, -1:1)
      common /global/ ptz, evec, xivec, spvec, elpty, angks, rclabs, 
     1     ipol, ispin, le2, iphabs

c--------- mod1.inp -------------------------------------------------------------
      character*80 title(nheadx)
c     integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec,
      integer mpot, nph, ntitle, ihole, ipr1, iafolp, iunf,
     1     nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx)
      integer lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), spinph(0:nphx)
      double precision  xion(0:nphx)
c     for OVERLAP option
      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
      double precision  rovr(novrx,0:nphx)
      common /mod1/ title, xion, xnatph, spinph, folp, gamach, rgrd,
     1     ca1, ecv, totvol, rovr, rfms1, iz, lmaxsc, mpot, nph, ntitle,
     2     ihole, ipr1, iafolp, nmix,nohole,jumprm, inters,
     3     nscmt, icoul, lfms1, novr, iphovr, nnovr, iunf

c--------- ldos.inp -------------------------------------------------------------
      integer mldos, lfms2
      double precision emin, emax, eimag, rfms2
      common /mod7/ emin, emax, eimag, rfms2, mldos, lfms2

c---- mod2.inp -------------------------------------------------------------
c     integer mphase, ipr2, ixc, ixc0, vr0, vi0, ispec, lreal, lfms2
      integer mphase, ipr2, ixc, ixc0, ispec, lreal, l2lp, iPlsmn
      integer lmaxph(0:nphx), iGrid
      character*6  potlbl(0:nphx)
c     double precision rgrd, rfms2, gamach, xkstep, xkmax, vixan
      double precision xkstep, xkmax, vixan, vr0, vi0
      common /mod2/ xkstep, xkmax, vixan, vr0, vi0, 
     &     lmaxph, mphase, ipr2, ixc, ixc0, ispec, lreal, l2lp,
     &     izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis, iPlsmn,
     &     iGrid, potlbl

c--------- mod3.inp -------------------------------------------------------------
      integer mfms, idwopt, minv
c     integer lmaxph(0:nphx)
c     real rfms2, rprec, rdirec, toler1, toler2
      real rprec, rdirec, toler1, toler2
      double precision   tk, thetad, sig2g
      common /mod3/ tk, thetad, sig2g, rprec, rdirec, toler1,
     1       toler2,  mfms, idwopt, minv

c--------- mod4.inp -------------------------------------------------------------
      integer  mpath, ms, nncrit, nlegxx, ipr4
c     real critpw, pcritk, pcrith,  rmax, rfms2
      real critpw, pcritk, pcrith,  rmax
      common /mod4/ critpw, pcritk, pcrith,  rmax,
     1       mpath, ms, nncrit, nlegxx, ipr4

c--------- mod5.inp -------------------------------------------------------------
      integer  mfeff, ipr5, iorder
      logical  wnstar
      double precision critcw
      common /mod5/ critcw, mfeff, ipr5, iorder, wnstar

c--------- mod6.inp -------------------------------------------------------------
c     integer  mchi, ispec, idwopt, ipr6, mbconv
c     double precision  vrcorr, vicorr, s02, alphat, sig2g
      integer  mchi, ipr6, mbconv, absolu !KJ added absolu 3-06
      double precision  vrcorr, vicorr, s02, alphat, thetae
      common /mod6/ vrcorr, vicorr, s02, alphat, thetae, 
     &     mchi, ipr6, mbconv, absolu   !KJ added absolu 3-06

c -------- so2.inp -------------------------------------------------------------
      integer  mso2conv, ipse, ipsk
      double precision wsigk, cen
      character(12) cfname
      common /so2/ wsigk, cen, cfname, mso2conv, ipse, ipsk

c--------- eels.inp ------------------------------------------------------------
c     EELS variables  !KJ 1-06 this section added for ELNES, EXELFS, MAGIC cards
      real*8 ebeam, aconv, acoll, thetax, thetay, emagic
      integer eels, relat, aver, cross, iinput,spcol
      integer nqr,nqf,magic
      integer ipmin,ipmax,ipstep
      common /eelsva/ ebeam,aconv,acoll,thetax,thetay,emagic,magic,
     &     nqr, nqf, aver, cross, relat, iinput, spcol,ipmin, ipmax,
     &     ipstep, eels
c     !KJ end
