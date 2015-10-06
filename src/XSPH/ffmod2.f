c     sub-program exchange
      program  ffmod2
c     subroutine ffmod2

c     cross-section and phase shifts calculations
c     coded by a.ankudinov 2000

c     INPUT: mod2.inp geom.dat global.inp and pot.pad
c     OUTPUT: xsect.bin and xsph.bin

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      character*6 potlbl
      dimension potlbl(0:nphx)
      dimension lmaxph(0:nphx), iatph(0:nphx), spinph(0:nphx)
      dimension rat(3, natx), iphat(natx)
      complex*16 ptz(-1:1, -1:1)

c     necessary input information from feff.inp file
c     see CARDs description in feff8 manual for more details
c     CONTROL mphase: 1-run (0-don't run)  the program
c     PRINT ipr2: for auxialry output files (default=0)
c     ispec: spectroscopy type (EXAFS, XANES, XES, DANES, FPRIME) 
c     vixan, xkstep, xkmax: energy grid for chosen spectroscopy
c     RDRIG rgrid: radial grid (default=0.05)
c     POTENTIAL info
c       nph: number of unique potentials
c       lmaxph: max orbital momentum for xsph calculations
c       potlbl: labels for unique potentials
c      ATOMS
c       nat: number of atoms
c       rat: their coordinates
c       iphat: type of potential for each site
c       iatph: representative atoms indices in atoms list
c      EXCHANGE ixc  vr0  vi0  ixc0 - exchange correlation model
c      RSIGMA (RPHASES) lreal (default=0)
c      FMS  rfms2 lfms2
       real rfms2
       integer lfms2
c      Global data
c        ipol - polarization type (default:0 - polarization average)
c        ispin - spin type (default=0 - spin independent)
c        le2 - include/exclude quad. transitions (default=2 - include)
c        angks - angle between x-ray propagation and spin (default=0)
c        ptz - polarization tenzor (default=0 for ipol=0)
      integer iPl, iGrid

c     parameters related to rdpot, to be passed to newly modified pot()
      character*80 title(nheadx)
      integer ntitle, nohole, ihole, inters, iafolp, iunf, jumprm
      double precision rnrmav, xmu, vint, rhoint, emu, s02, erelax
      double precision wp, ecv, rs, xf, qtotel, totvol
      dimension adgc(10,30,0:nphx)
      dimension adpc(10,30,0:nphx)
      dimension dgc(251,30,0:nphx)
      dimension dgc0(251)
      dimension dmag(251,0:nphx)
      dimension dpc(251,30,0:nphx)
      dimension dpc0(251)
      dimension edens(251,0:nphx)
      dimension edenvl(251,0:nphx)
      dimension eorb(30)
      dimension folp(0:nphx)
      dimension folpx(0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension iorb(-4:3,0:nphx)
      dimension iz(0:nphx)
      dimension kappa(30)
      dimension qnrm(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)
      dimension vclap (251,0:nphx)
      dimension vtot (251,0:nphx)
      dimension vvalgs (251,0:nphx)
      dimension xion(0:nphx)
      dimension xnatph(0:nphx)
      dimension xnmues(0:lx,0:nphx)
      dimension xnval(30,0:nphx)


      call par_begin
      if (worker) go to 400

c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log2.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log2.dat', 'feff')

c     read  INPUT data files: geom.dat, global.dat and mod2.inp.
c     Josh - added flag iPl for PLASMON card
      call rexsph(mphase, ipr2, ispec, vixan, xkstep,xkmax,gamach,rgrd,
     1             nph, lmaxph, potlbl, spinph, iatph, nat, rat, iphat,
     2             ixc, vr0, vi0, ixc0, lreal, rfms2, lfms2, l2lp,
     3             ipol, ispin, le2, angks, ptz, iPl, iGrid,
     4             izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)

      if (mphase .eq. 1)  then
         call wlog(' Calculating cross-section and phases...')
         call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,
     1          emu, s02, erelax, wp, ecv,rs,xf, qtotel,
     2          imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3          dgc0, dpc0, dgc, dpc, adgc, adpc,
     3          edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4          eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
     5          inters, totvol, iafolp, xion, iunf, iz, jumprm)

c                   wrxsec, verbse, phpad,
         call xsph (.true., .true., 'phase.pad',
     -          ipr2, ispec, vixan, xkstep, xkmax, gamach, rgrd,
     1          nph, lmaxph, potlbl, spinph, iatph, nat, rat, iphat,
     2          ixc, vr0, vi0, ixc0, lreal, rfms2, lfms2, l2lp,
     3          ipol, ispin, le2, angks, ptz, iPl,
     4          izstd, ifxc, ipmbse, itdlda, nonlocal,
c        pass parameters from rdpot
     1          ntitle, title, rnrmav, xmu, vint, rhoint,
     2          emu, s02, erelax, wp, ecv, rs, xf, qtotel,
     3          imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     4          dgc0, dpc0, dgc, dpc, adgc, adpc,
     5          edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     6          iorb, nohole, ihole,
     7          inters, totvol, iafolp, xion, iunf, iz, jumprm)
c squelch compiler warning about unused dummy variables, apparently
c removed from f85e
c     iGrid, (after iPl)
c   , ibasis)

         call wlog(' Done with module 2: cross-section and phases...')
      endif

c     OUTPUT: data for the next modules is written in xsph.bin
c     auxilary output can be obtained using 'ipr2' (see feff8.2 manual)

      close (unit=11)
  400 call par_barrier
      call par_end

c     sub-program exchange
      stop
c     return
      end
