c     sub-program exchange
      program ffmod6
c     subroutine ffmod6 (iabs)

c     final calculations for various spectroscopies
c     (EXAFS, XANES, FPRIME, DANES, XES)
c     written by a.ankudinov 2000
c     modified by a.ankudinov 2001 for new I/O structure

c     INPUT: mod6.inp global.dat xsect.bin fms.bin list.dat and feff.pad
c     OUTPUT: xmu.dat (chi.dat for EXAFS)

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      integer elnes,ipmin,ipmax,ipstep  !KJ my variables 1-06 
      integer absolu !KJ 3-06     


      call par_begin
      if (worker) go to 400

c     sub-program exchange
      iabs = 1
c     

c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log6.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log6.dat', 'feff')

c     read  input files
      call reff2x(mchi, ispec, ipr6, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 1-06     

      if (mchi .eq. 1)  then
         call wlog(' Calculating chi...')
         if (ispec.gt.0 .and. ispec.lt.3) then 
c           using FMS+Paths method
            call ff2xmu (ispec, ipr6, idwopt, critcw, s02, sig2g, 
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            ipmin,ipmax,ipstep)   !KJ added this line  1-06     
         elseif (ispec.eq.3 .or. ispec.eq.4) then 
c           using FMS+Paths method
            call ff2afs ( ipr6, idwopt, critcw, s02, sig2g, 
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 4-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            ipmin,ipmax,ipstep)   !KJ added this line 4-06     
         else
c           using MS Paths expansion
            call ff2chi (ispec, ipr6, idwopt, critcw, s02, sig2g, 
     1           tk, thetad, mbconv, absolu, !KJ added absolu 3-06
     1           vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4           ipmin,ipmax,ipstep) !KJ added this line 1-06     
         endif
         call wlog(' Done with module 6: DW + final sum over paths.')
      endif
      close (unit=11)

  400 call par_barrier
      call par_end
      
c     sub-program exchange
      stop  
c     return

      end
