c     sub-program exchange
      program ffmod5
c     subroutine ffmod5

c     scattering F-matrix multiplication for each MS path
c     written by a.ankudinov 2000, using subroutines
c     which were written earlier by j.rehr and others
c     modified by a.ankudinov 2001 for new I/O structure

c     INPUT: phase.pad, paths.dat, mod5.inp and global.dat
c     OUTPUT: feff.pad and list.dat files
      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      double precision evec(3), xivec(3)
      complex*16 ptz(-1:1, -1:1)
      integer  mfeff, ipr5, iorder
      logical  wnstar
      double precision critcw, angks, elpty


      call par_begin
      if (worker) go to 400

c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log5.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log5.dat', 'feff')

c                 read  mod5.inp 
      call regenf(mfeff, ipr5, critcw, iorder, wnstar,
c                 and global.dat
     1            ipol, ispin, le2, angks, elpty, evec, xivec, ptz)
      if (nspx.gt.1) ispin = abs(ispin)

      if (mfeff .eq. 1)  then
         call wlog(' Calculating EXAFS parameters...')
         call genfmt (ipr5, critcw, iorder, wnstar,
     1                ipol, ispin, le2, angks, elpty, evec, xivec, ptz)
         call wlog(' Done with module 5: F_eff.')
      endif
      close (unit=11)

 400  call par_barrier
      call par_end

c     sub-program exchange
      stop
c     return

      end
