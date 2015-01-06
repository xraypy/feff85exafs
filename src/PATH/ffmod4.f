
c     sub-program exchange
      program  ffmod4
c     subroutine ffmod4

c     makes paths list using cluster geometry and phase shifts
c     written by a.ankudinov 2000 using earlier subroutines
c     written by s.zabinsky
c     modified by a.ankudinov 2001 for new I/O structure

c     INPUT FILES
c       global.dat, geom.dat - global infomation file is read here 
c       mod4.inp - specific information for present module
c       phase.pad - output of XSPH module is read using subroutine 
c                  'rdxsph' inside subroutine 'prcrit'.
c                   needed  data: (list of variables)
c                  (ne, ne1, npot, ik0, em, eref2, potlbl, ph4)
c     OUTPUT FILE
c       paths.dat - list of filtered paths

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'

      
c      dimension lmaxph(0:nphx)
      dimension rat(3,natx), iphat(natx), ibounc(natx)
      double precision evec(3), xivec(3)
      character*6  potlbl(0:nphx)
        integer  mpath, ms, nncrit, nlegxx, ipr4
        real critpw, pcritk, pcrith,  rmax, rfms2

      integer eels !KJ added 5/06
      character*30 fname

c     Following passed to pathfinder, which is single precision.
c     Be careful to always declare these!
      parameter (necrit=9, nbeta=40)
      real fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex)
      real xlamc(necrit), xlam(nex)

      call par_begin
      if (worker) go to 400
      
c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log4.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log4.dat', 'feff')

c     INPUT: read geom.dat, global.dat, mod4.inp
      call repath (ms, mpath, ipr4,  pcritk, pcrith, nncrit, rmax,
     1             nlegxx, rfms2, critpw,
c                  geom.dat
     2             nat, rat, iphat, ibounc,
c                  global.dat
     3             ipol, ispin, evec, xivec)
c ,eels) !KJ added eels 5/06
      eels=0
       if (nspx.gt.1) ispin = abs(ispin)

      if (ms.eq.1  .and.  mpath.eq.1)  then
         call wlog(' Preparing plane wave scattering amplitudes...')
         call prcrit (ne, nncrit, ik0, cksp, fbeta, ckspc, 
     1                fbetac, potlbl, xlam, xlamc)

c        Dump out fbetac for central atom and first pot
         if (ipr4 .ge. 3 .and. ipr4.ne.5)  then
            do 260  iph = 0, 1
               do 250  ie = 1, nncrit
                  write(fname,200)  ie, iph
  200             format ('fbeta', i1, 'p', i1, '.dat')
                  open (unit=1, file=fname, status='unknown')
                  write(1,210)  iph, ie, ckspc(ie)
  210             format ('# iph, ie, ckspc(ie) ', 2i5, 1pe20.6, /
     1                    '#  angle(degrees), fbeta/|p|,  fbeta')
                  do 240  ibeta = -nbeta, nbeta
                     cosb = .025 * ibeta
                     if (cosb .gt.  1)  cosb =  1
                     if (cosb .lt. -1)  cosb = -1
                     angle = acos (cosb)
                     write(1,230)  angle*raddeg, 
     1                  fbetac(ibeta,iph,ie)/ckspc(ie),
     2                  fbetac(ibeta,iph,ie)
  230                format (f10.4, 1p, 2e15.6)
  240             continue
                  close (unit=1)
  250          continue
  260       continue
         endif

         call wlog(' Searching for paths...')
         call paths (ckspc, fbetac, xlamc, pcritk, pcrith, critpw,
     1               nncrit, rmax, nlegxx, rfms2,
     2               nat, rat, iphat, ibounc) 
         call wlog(' Eliminating path degeneracies...')
         call pathsd (ckspc, fbetac, xlamc, ne, ik0, cksp, 
     1                fbeta, xlam, critpw, ipr4, nncrit, potlbl,
     1            ipol, ispin, evec, xivec,eels) !KJ added eels 5/06
         call wlog(' Done with module 4: pathfinder.')
      endif
      close (unit=11)

  400 call par_barrier
      call par_end

c     sub-program exchange
      stop
c     return

      end
