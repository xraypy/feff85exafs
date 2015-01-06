c     sub-program exchange
      program ffmod3
c     subroutine ffmod3

c     full multiple scattering code (inversion of big matrix)
c     written by a.ankudinov 2000 using earlier written subroutines
c     coded by b.ravel
c     modified by a.ankudinov 2001 for new matrix inversion algorithms
c     and new I/O structure

c     INPUT:  geom.dat, global.dat and mod3.inp
c     OUTPUT:  fms.bin

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      real*8 wall_start, wall_end
      dimension lmaxph(0:nphx), rat(3,natx), iphat(natx)
      complex*16 ptz(-1:1, -1:1)

      real  rfms2, rdirec, toler1, toler2

      call par_begin

c     Initialize clock
      call seconds(wall_start)
      wall_comm = 0.0

c     open the log file, unit 11.  See subroutine wlog.
      if (master) then
        open (unit=11, file='log3.dat', status='unknown', iostat=ios)
        call chopen (ios, 'log3.dat', 'feff')
      else
        par_type = 2
      endif

c     read  INPUT: files geom.dat, global.dat and mod3.inp
      call reafms(mfms, rfms2, idwopt, tk, thetad, sig2g,
     1            lmaxph, nat, iphat, rat,
     2            ipol, ispin, le2, angks, ptz,
     3            minv, rdirec, toler1, toler2)

      if (mfms.eq.1 )  then
c        first do fms inside rfms2
         call fmstot(rfms2, idwopt, tk, thetad, sig2g, 
     1            lmaxph, nat, iphat, rat,
     2            ipol, ispin, le2, angks, ptz,
     3            minv, rdirec, toler1, toler2,elnes,ipmin,ipmax,ipstep)
         call wlog(' Done with module 3: FMS.')
      endif

c     OUTPUT: gtr.bin file for the next modules

      if (master) close (unit=11)

c--   Time at end of run
      call seconds(wall_end)
      if (master .and. parallel_run) then
        write (6,*) 'total time     ', wall_end - wall_start
        write (6,*) 'communicate time', wall_comm
      endif

      call par_end

c     sub-program exchange
      stop
c     return

      end
