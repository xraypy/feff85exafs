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

      call par_begin

c     Initialize clock
      call seconds(wall_start)
      wall_comm = 0.0

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
         call pot (rgrd, nohole,
     $             inters, totvol, ecv, nscmt, nmix, ntitle, title,
     $             nat, nph, ihole, iafolp,
     $             ixc, iphat, rat, iatph,
     $             xnatph, novr,
     $             iphovr, nnovr, rovr, folp, xion, iunf, iz, ipr1,
     $             ispec, jumprm,
     $             lmaxsc, icoul, ca1, rfms1, lfms1)
c                  gamach
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
