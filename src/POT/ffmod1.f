c     sub-pro exchange point
      program ffmod1
c     subroutine ffmod1

c     calculate  el. density and potential given atomic positions for
c     cluster atoms or other similar information
c     calculation can vary in complexity: self-consistency (on/off),
c     spin dependency (on/off), etc..
c       coded by a.l. ankudinov 2000, for modular code structure
c       modified by a.l. ankudinov 2001, for new i/o structure

c     INPUT files: mod1.inp, geom.dat
c     OUTPUT file: pot.pad

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      real*8 wall_start, wall_end

c     use feff8.2 manual to get more information about each CARD
cc    mod1.inp
        character*80 title(nheadx)
        integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec,
     1     iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
        integer iz(0:nphx), lmaxsc(0:nphx)
        real rfms1
        double precision gamach, rgrd, ca1, ecv, totvol
        double precision  xnatph(0:nphx), folp(0:nphx),  xion(0:nphx)
c       CONTROL mpot
c       RGRID  rgrd
c       TITLE title
c        ntitle: number of title lines(default:0)
c        title:  title lines(default:none)
c       PRINT   ipr1:   print option (default:0)
c       EXAFS, XANES, DANES, FPRIME, XES
c        ispec: type of spectroscopy (default:0-EXAFS)
c       NOHOLE: turn on/off core-hole potential
c       HOLE  ihole: index of core-hole orbital
c        gamach: core hole lifetime
c       POTENTIALS card
c        nph  - number of different potential types(default:1)
c        iz - nicleus charge for each potential charge(default:none)
c        lmaxsc - max orb momentum to calculate (default:3)
c        xnatph - relative amount of atoms of each type (default:1)
c       ION card
c        xion - total initial charge for each potential type
c               (iz + el.charge) which might be fractional (default:0)
c       EXCHANGE card: ixc=2 for potential calculation
c       JUMPRM: turn on potential jump removal at mt radius (default:0)
c       AFOLP iafolp: turn on/off automatic overlap of muffintin spheres
c       FOLP  folp: manual setting for overlapping muffin-tin spheres
c       INTERSTITIAL inters (default:0)  totvol (default:0)
c       SCF rfms1 lfms1 nscmt ca1 nmix ecv  icoul 
c       OVERLAP geometry ( rarely used for EXAFS calculations only)
        integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
        double precision  rovr(novrx,0:nphx)

cc    geom.dat
        integer  nat, iatph(0:nphx), iphat(natx)
        double precision  rat(3,natx)
c       ATOM card
c         nat: number of atoms in a clsuter
c         rat: x,y,z coordinates of all atoms
c         iphat: which potential type correspond to each atom
c         iatph: index of representative atom for each potential type

        
        double precision rnrmav, xmu, vint, rhoint, emu, s02, erelax
        double precision wp, rs, xf, qtotel
c       r mesh index just inside rmt
        dimension imt(0:nphx)
c       r mesh index just inside rnorman
        dimension inrm(0:nphx)
c       muffin tin radius
        dimension rmt(0:nphx)
c       norman radius
        dimension rnrm(0:nphx), qnrm(0:nphx)
c       folp(0:nphx)  - overlap factor for rmt calculation
        dimension folpx(0:nphx)
c       need irregular solution for complex potential. fix later
        dimension dgc0(251), dpc0(251)
c       additioal data needed for relativistic version
        dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
        dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
c       Overlap calculation results
c       overlapped density*4*pi
        dimension edens(251,0:nphx)
c       overlapped coul pot
        dimension vclap(251,0:nphx)
c       overlapped total potential
        dimension vtot (251,0:nphx)
        dimension edenvl(251,0:nphx)
        dimension vvalgs (251,0:nphx)
        dimension dmag(251,0:nphx+1)
        dimension xnval(30,0:nphx+1)
        dimension eorb(30,0:nphx+1)
        dimension kappa(30,0:nphx+1), iorb(-4:3,0:nphx+1)
        dimension xnmues(0:lx,0:nphx)
c       Josh use nhtmp to save nohole value
        integer nhtmp


      call par_begin

c     Initialize clock
      call seconds(wall_start)
      wall_comm = 0.0d0

c     open the log file, unit 11.  See subroutine wlog.
      if (master) then
        open (unit=11, file='log1.dat', status='unknown', iostat=ios)
        call chopen (ios, 'log1.dat', 'feff')
      else
        par_type = 2
      endif


c     INPUT: read data in pot.inp  and geom.dat files
c     and transform it to atomic hartree units
      call reapot (mpot, rgrd, ntitle, title, ipr1, ispec,
     1           nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     2           xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     3           rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     4           novr, iphovr, nnovr, rovr,
     5           nat, rat, iphat, iatph)

      if (mpot .eq. 1)  then
         call wlog(' Calculating potentials ...')

c$$$         print *, rgrd, nohole,
c$$$     $          inters, totvol, ecv, nscmt, nmix, ntitle,
c$$$     $          nat, nph, ihole, iafolp, ixc
c$$$         print *, 1, rat(1,1), rat(2,1), rat(3,1), iphat(1)
c$$$         print *, 2, rat(1,2), rat(2,2), rat(3,2), iphat(2)
c$$$         print *, 3, rat(1,3), rat(2,3), rat(3,3), iphat(3)
c$$$         print *, 13, rat(1,13), rat(2,13), rat(3,13), iphat(13)
c$$$         print *, 23, rat(1,23), rat(2,23), rat(3,23), iphat(23)

         call pot (.true., rgrd, nohole,
     $          inters, totvol, ecv, nscmt, nmix, ntitle, title,
     $          nat, nph, ihole, iafolp, ixc, iphat, rat, iatph,
     $          xnatph, novr, iphovr, nnovr, rovr,
     $          folp, xion, iunf, iz, ipr1, ispec, jumprm,
     $          lmaxsc, icoul, ca1, rfms1, lfms1,
     $
     -          rnrmav, xmu, vint, rhoint,
     1          emu, s02, erelax, wp, rs, xf, qtotel,
     2          imt, rmt, inrm, rnrm, folpx,
     3          dgc0, dpc0, dgc, dpc, adgc, adpc,
     4          edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     5          eorb, kappa, iorb, qnrm, xnmues, nhtmp
     6    )
c        return stuff for wrpot
c                  gamach
c        write stuff into pot.pad
         call wrpot (nph, ntitle, title, rnrmav, xmu, vint, rhoint,
     1          emu, s02, erelax, wp, ecv, rs, xf, qtotel,
     2          imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3          dgc0, dpc0, dgc, dpc, adgc, adpc,
     3          edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4          eorb(1,0), kappa(1,0), iorb, qnrm, xnmues, nhtmp,
     5          ihole, inters, totvol, iafolp, xion, iunf, iz, jumprm)
         
c        write misc.dat
         if (ipr1 .ge. 1)  then
            open (unit=1, file='misc.dat', status='unknown', iostat=ios)
            call chopen (ios, 'misc.dat', 'potph')
            call wthead(1, ntitle, title)
            close (unit=1)
         endif

         call wlog(' Done with module 1: potentials. ')
      endif

c     OUTPUT: subroutine pot writes main output file pot.pad
c     with information on potentials, necessary for other modules;
c     additional output files can be obtained using PRINT card

      if (master) close (unit=11)

c--   Time at end of run
      call seconds(wall_end)
      if (master .and. parallel_run) then
        write (6,*) 'total time    ', wall_end - wall_start
        write (6,*) 'communicate time', wall_comm
      endif
      call par_end

c     sub-pro exchange point
      stop
c     return
      end
