      subroutine pathsd (ckspc, fbetac, xlamc, ne, ik0, cksp, 
     1                fbeta, xlam, critpw, ipr2,  nncrit, potlbl,
     1            ipol, ispin, evec, xivec,eels)  !KJ added eels 5/06

c     New degeneracy checker, cute and hopefully fast for large
c     problems

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

c     np1x  number of paths to consider at 1 time
c     parameter (np1x = 12 000)
      parameter (np1x = 600 000)
      dimension iout(3,np1x), iout0(3)

      dimension index(np1x)
      double precision dhash(np1x), dcurr, ddum
      dimension rx(npatx), ry(npatx), rz(npatx), ipat(npatx+1)
      dimension rx0(npatx), ry0(npatx), rz0(npatx), ipat0(npatx+1)
      double precision rid(npatx+1), betad(npatx+1), etad(npatx+1)

      integer eels !KJ 5/06 - just passed on to timrep
      character*80 head(nheadx)
      character*6  potlbl(0:nphx)
      double precision xivec(3), evec(3)

c     eps5 for rtotal range, eps3 for individual leg parameters.
c     eps3 large since code single precision and don't want round-off
c     error to reduce degeneracy.
      parameter (eps5 = 2.0e-5)
      parameter (eps3 = 1.0e-3)

      logical ldiff, last
      parameter (necrit=9, nbeta=40)
      real fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex)
      real xlamc(necrit), xlam(nex)

      character*512 slog

      write(slog,30) critpw
      call wlog(slog)
   30 format ('    Plane wave chi amplitude filter', f7.2, '%')

c     Read atoms info
      open (file='paths.bin', unit=3, access='sequential',
     1      form='unformatted', status='old', iostat=ios)
      call chopen (ios, 'paths.bin', 'pathsd')
      read(3) nhead
      do 40  ihead = 1, nhead
         read(3)  head(ihead)
   40 continue
c     Header lines above include carriage control
      read(3)  nat
c     rat and ipot could be permuted by paths.f
      do 50  i = 0, nat
         read(3) (rat(j,i),j=1,3), ipot(i), i1b(i)
   50 continue

c     Initialize stuff...
c     nptot  number of total paths, incl all degeneracies
c     nuptot number of unique paths for which must calc xafs
c     ngs    number of generalized shells (unique distances)
      nptot = 0
      nuptot = 0
      ngs = 0
      xportx = eps5
      ndegx = -1
c     Initialize keep criterion
      xcalcx = -1

c     write output to paths.dat
      if (ipr2 .ne. 5)  then
         open (unit=1, file='paths.dat', status='unknown', iostat=ios)
         call chopen (ios, 'paths.dat', 'pathsd')
         do 60  ihead = 1, nhead
            ii = istrln(head(ihead))
            write(1,58)  head(ihead)(1:ii)
   58       format(a)
   60    continue
c        write(1,61)  critpw
   61    format (' Plane wave chi amplitude filter', f7.2, '%')
         write(1,62)
   62    format (1x, 71('-'))
      endif

c     Write crit.dat (criteria information)
      if (ipr2 .ge. 1)  then
         open (unit=4, file='crit.dat', status='unknown', iostat=ios)
         call chopen (ios, 'crit.dat', 'pathsd')
         do 65  ihead = 1, nhead
              ii = istrln(head(ihead))
            write(4,58)  head(ihead)(1:ii)
   65    continue
         write(4,61)  critpw
         write(4,62)
         write(4,80)
   80    format (' ipath nleg ndeg     r       pwcrit    ',
     1           'xkeep   accuracy   xheap    accuracy')
      endif

c     Read path data for each total path length range

c     Prepare for first path.
      read(3,end=999)  r0, iout0

c     Begin next total path length range
      last = .false.
  100 continue
      ngs = ngs+1
      rcurr = r0
      np = 1
      do 110  i = 1,3
         iout(i,np) = iout0(i)
  110 continue
  120 read(3,end=140)  r0, iout0
         if (abs(r0-rcurr) .lt. eps3)  then
            np = np+1
            if (np .gt. np1x) then
               write(slog,122) ' np, np1x ', np, np1x
               call wlog(slog)
  122          format (a, 2i15)
               call par_stop('np > np1x')
            endif
            do 130  i = 1, 3
               iout(i,np) = iout0(i)
  130       continue
         else
c           r0 is the rtot for the next set
c           iout0 is the packed atom list for the first path of the
c           next set
            goto 200
         endif
      goto 120
  140 continue
c     Get here only if end-of-file during read
      last = .true.

  200 continue

      nupr = 0
c     variable nuprtt was nuprtot, changed to be six chars, SIZ 12/93
      nuprtt = 0

c     Hash each path into an integer
      iscale = 1000
      do 230  ip = 1, np

         npat = npatx
         call upack (iout(1,ip), npat, ipat)

c        Get hash key for this path.
c        If two paths are the same, except time-reversed, the xafs
c        will be the same, so check for this type of degeneracy.
c        We do this by choosing a 'standard order' for a path --
c        if it's the other-way-around, we time-reverse here.
         call timrep (npat, ipat, rx, ry, rz, dhash(ip),
     1            ipol, ispin, evec, xivec,eels)  !KJ added eels 5/06

  230 continue

c     Do a heap sort on these things
      call sortid (np, index, dhash)

c     Find beginning and end of range with same hash key
c     i0 is beginning of hash range, i1 is end of the range

      i0 = 1
  300 continue
         i1 = np + 1
         dcurr = dhash(index(i0))
         do 310  ip = i0+1, np
            if (dhash(index(ip)) .ne. dcurr)  then
c              end of a hash range
               i1 = ip
               goto 311
            endif
  310    continue
  311    continue
         i1 = i1-1

c        At this point, i0 is the first path and i1 the last
c        of a hash range.  Do whatever you want with them!

c        Sum degeneracy, including degeneracy from 1st bounce atom.
c        Check this range to see if all of the paths are actually 
c        degenerate.  Make sure time-ordering is standard.
         npat0 = npatx
         call upack (iout(1,index(i0)), npat0, ipat0)
         call timrep (npat0, ipat0, rx0, ry0, rz0, ddum,
     1            ipol, ispin, evec, xivec,eels) !KJ added eels 5/06

         ndeg = 0
         do 430  ii = i0, i1
            npat = npatx
            call upack (iout(1,index(ii)), npat, ipat)
c           Note that if path gets time-reversed, we lose 1st bounce 
c           flag (since first atom is now last...), so save path deg
            ndpath = i1b(ipat(1))
            call timrep (npat, ipat, rx, ry, rz, ddum,
     1            ipol, ispin, evec, xivec,eels) !KJ added eels 5/06
c           Sum degeneracy here.
            ndeg = ndeg + ndpath
c           Check for hash collisons begins here.
            ldiff = .false.
            if (npat .ne. npat0)  then
               ldiff = .true.
               goto 430
            endif
            do 320  iat = 1, npat
               if (ipot(ipat(iat)) .ne. ipot(ipat0(iat)))  then
                  ldiff = .true.
                  goto 400
               endif
  320       continue
            do 330  ileg = 1, npat
               if (abs(rx(ileg)-rx0(ileg)) .gt. eps3  .or.
     1             abs(ry(ileg)-ry0(ileg)) .gt. eps3  .or.
     2             abs(rz(ileg)-rz0(ileg)) .gt. eps3)  then
                  ldiff = .true.
                  goto 400
               endif
  330       continue
  400       continue
            if (ldiff)  then
               call wlog(' WARNING!!  Two non-degenerate paths,' //
     1                   ' hashed to the same hash key!!')
  402          format (1x, 2e28.20)
               write(slog,402) dhash(index(i0)), dhash(index(ii))
               call wlog(slog)
  404          format (1x, 2i10, a)
               write(slog,404) npat0, npat, '  npat0, npat'
               call wlog(slog)
               call wlog(' iat, ipot0, ipot, ipat0, ipat')
               do 410  iat = 1, npat
  406             format (5i10)
                  write(slog,406) iat, ipot(ipat0(iat)), 
     1               ipot(ipat(iat)), ipat0(iat), ipat(iat)
                  call wlog(slog)
  410          continue
               call wlog(' ileg, rx0,ry0,rz0,  rx1,ry1,rz1')
               do 420  ileg = 1, npat
  412             format(i6, 1p, 3e18.10)
                  write(slog,412) ileg, rx0(ileg), rx(ileg)
                  call wlog(slog)
                  write(slog,412) ileg, ry0(ileg), ry(ileg)
                  call wlog(slog)
                  write(slog,412) ileg, rz0(ileg), rz(ileg)
                  call wlog(slog)
  420          continue
               call par_stop('hash error')
            endif
  430    continue

c        Find path pw importance factors, and recalculate 
c        pathfinder crits for output
         call outcrt (npat0, ipat0, ckspc,
     1                nncrit, fbetac, xlamc, ne, ik0, cksp, 
     1                fbeta, xlam, 
     1                ipot,
     1                xport, xheap, xheapr, xkeep, xcalcx)

         if (xportx*ndegx .le. 0)  then
            xportx = xport
c           ndegx is degeneracy of path that makes xportx, used for
c           testing new path keep crit
            ndegx = ndeg
         endif
c        frac is fraction of max importance to use for test
         frac = 100*ndeg*xport/(ndegx*xportx)

c        Write output if path is important enough (ie, path is
c        at least critpw % important as most important path found
c        so far.)
         if (frac .ge. critpw)  then
            nupr = nupr+1
            nuprtt = nuprtt+ndeg
            nptot = nptot + ndeg
            nuptot = nuptot + 1

c           Write path info to paths.dat
c           mpprmd is double precision, used to get angles
c           180.000 instead of 179.983, etc.
            call mpprmd (npat0, ipat0, rid, betad, etad)
c           skip paths.dat if not necessary
            if (ipr2 .eq. 5)  goto 576
            write(1,500) nuptot, npat0+1, real(ndeg),
     1              rcurr/2
  500       format (1x, 2i5, f8.3,
     1             '  index, nleg, degeneracy, r=', f8.4)
            write(1,502)
  502       format ('      x           y           z     ipot  ',
     1              'label      rleg      beta        eta')
            do 510  i = 1, npat0
               iat = ipat0(i)
               write(1,506)  rat(1,iat), rat(2,iat),
     1                  rat(3,iat), ipot(iat), potlbl(ipot(iat)),
     1                  rid(i), betad(i)*raddeg, etad(i)*raddeg
  506          format (3f12.6, i4, 1x, '''', a6, '''', 1x, 3f10.4)
  510       continue
            write(1,506)  rat(1,0), rat(2,0), rat(3,0), ipot(0), 
     1         potlbl(ipot(0)),
     1         rid(npat0+1), betad(npat0+1)*raddeg, etad(npat0+1)*raddeg
c           End of paths.dat writing for this path

c           Write to crit.dat here (unit 4, opened above)
  576       continue

c           cmpk is degeneracy corrected xkeep, should equal frac
            cmpk = xkeep*ndeg/ndegx
c           cmpk is accuracy of xkeep, 100 is perfect
            cmpk = 100 - 100*(abs(frac-cmpk)/frac)

c           cmph is same thing for xheap
            if (xheap .lt. 0)  then
               cmph = 100
            else
               cmph = xheap*ndeg/ndegx
               cmph = 100 - 100*(abs(frac-cmph)/frac)
            endif

            if (ipr2 .ge. 1)  then
               write(4,560)  nuptot, npat0+1, ndeg, rcurr/2, frac,
     1             xkeep, cmpk, xheap, cmph
  560          format (i6, i4, i6, 3f10.4, f8.2, f10.4, 1pe14.3)
            endif

c           write out fraction error between xkeep and critpw
         endif

c        And do next ihash range
         i0 = i1+1
      if (i0 .le. np)  goto 300

c     type600,  ngs, rcurr, nupr
c 600 format (1x, i5, f12.6, i7, ' igs, rcurr, nupr')
c     write(80,601)  ngs, rcurr/2, nupr, nuprtt
c 601 format (1x, i8, f12.6, 2i9)

      if (.not. last) goto 100

  999 if (ipr2 .ne. 5)  close (unit=1)
c     delete paths.bin when done...
      close (unit=3, status='delete')
      close (unit=4)

      write(slog,620) nuptot, nptot
      call wlog(slog)
  620 format ('    Unique paths', i7, ',  total paths', i8)

c     Do not let user accidently fill up their disk
      if (nuptot .gt. 1200)  then
      call wlog(' You have found more than 1200 paths.  Genfmt')
      call wlog(' could require a lot of time and more than 6 meg of')
      call wlog(' storage.  Suggest a larger critpw to reduce number')
      call wlog(' of paths.  To continue this calculation, restart')
      call wlog(' with current paths.dat and module genfmt (3rd module')
      call wlog(' on CONTROL card).')
      call par_stop('User must verify very large run.')
      endif
      return
c 999 stop 'no input'
      end
