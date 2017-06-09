      subroutine pathsd (ckspc, fbetac, ne, ik0, cksp, fbeta,
     1                   critpw, ipotnn, ipr2, 
     1                   pcritk, pcrith, nncrit, potlbl)

c     New degeneracy checker, cute and hopefully fast for large
c     problems

c     pcritk and pcrith used only for analysis after outcrt

      include 'const.h'
      include 'dim.h'
      include 'pola.h'
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

c     np1x  number of paths to consider at 1 time
      parameter (np1x = 40 000)
c     parameter (np1x = 60 000)
      dimension iout(3,np1x), iout0(3)

      dimension index(np1x)
      double precision dhash(np1x), dcurr, ddum
      dimension rx(npatx), ry(npatx), rz(npatx), ipat(npatx+1)
      dimension rx0(npatx), ry0(npatx), rz0(npatx), ipat0(npatx+1)
      double precision rid(npatx+1), betad(npatx+1), etad(npatx+1)

      parameter (nheadx = 40)
      character*80 head(nheadx)
      dimension lhead(nheadx)

      character*6  potlbl(0:npotx)
       character*128 messag
c     eps5 for rtotal range, eps3 for individual leg parameters.
c     eps3 large since code single precision and don't want round-off
c     error to reduce degeneracy.
      parameter (eps5 = 2.0e-5)
      parameter (eps3 = 2.0e-3)

      logical ldiff, last
      parameter (necrit=9, nbeta=40)
      real fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)

       write(messag,30) critpw
       call echo(messag)
   30 format ('    Plane wave chi amplitude filter', f7.2, '%')

c     Read atoms info
      open (file='paths.bin', unit=3, access='sequential',
     1      form='unformatted', status='old', iostat=ios)
      call chopen (ios, 'paths.bin', 'pathsd')
      read(3) nhead
      do 40  ihead = 1, nhead
         read(3)  head(ihead)
         read(3)  lhead(ihead)
   40 continue
c     Header lines above include carriage control
      read(3)  nat
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
      ndegx = 1
      c0lim = 1.0e10
      c1lim = 1.0e10
c     Initialize keep criterion
      xcalcx = -1

c     write output to paths.dat
      if (ipr2 .ne. 5)  then
         open (unit=1, file='paths.dat', status='unknown', iostat=ios)
         call chopen (ios, 'paths.dat', 'pathsd')
         do 60  ihead = 1, nhead
            write(1,58)  head(ihead)(1:lhead(ihead))
   58       format(a)
   60    continue
         write(1,61)  critpw
   61    format (' Plane wave chi amplitude filter', f7.2, '%')
         write(1,62)
   62    format (1x, 79('-'))
      endif

c     Write crit.dat (criteria information)
      if (ipr2 .ge. 1)  then
         open (unit=4, file='crit.dat', status='unknown', iostat=ios)
         call chopen (ios, 'crit.dat', 'pathsd')
         do 65  ihead = 1, nhead
            write(4,58)  head(ihead)(1:lhead(ihead))
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
               write(messag,'(1x,2a,i8,a)') ' Warning: Pathfinder ',
     $              'stopping at ', np1x, ' paths.'
               call echo(messag)
               goto 200
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
         call timrep (npat, ipat, rx, ry, rz, dhash(ip))

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
         call timrep (npat0, ipat0, rx0, ry0, rz0, ddum)

         ndeg = 0
         do 430  ii = i0, i1
            npat = npatx
            call upack (iout(1,index(ii)), npat, ipat)
c           Note that if path gets time-reversed, we lose 1st bounce 
c           flag (since first atom is now last...), so save path deg
            ndpath = i1b(ipat(1))
            call timrep (npat, ipat, rx, ry, rz, ddum)
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
               call echo(' Pathsd WARNING:')
               call echo(' two non-degenerate paths hashed '//
     1                 'to the same hash key!!')
               
               write(messag,'(1x,2g14.6)')
     $              dhash(index(i0)), dhash(index(ii))
               call echo(messag)
               write(messag,'(1x,a,2i8)') '  npat0, npat = ',
     $              npat0, npat 
               call echo(messag)
               call echo('  iat, ipot0, ipot, ipat0, ipat')
               do 410  iat = 1, npat
                  write(messag,'(3x,5i8)')  iat,
     $                 ipot(ipat0(iat)), ipot(ipat(iat)),
     $                 ipat0(iat), ipat(iat)
                  call echo(messag)
  410          continue
               call echo('  ileg, rx0,ry0,rz0,  rx1,ry1,rz1')
 412           format(1x,i7,2f12.6)
               do 420  ileg = 1, npat
                  write(messag,412) ileg, rx0(ileg), rx(ileg)
                  call echo(messag)                  
                  write(messag,412) ileg, ry0(ileg), ry(ileg)
                  call echo(messag)                  
                  write(messag,412) ileg, rz0(ileg), rz(ileg)
                  call echo(messag)                  
  420          continue
               call fstop(' at PATHSD: hash error')
            endif
  430    continue

c        Find path pw importance factors, and recalculate 
c        pathfinder crits for output
         call outcrt (npat0, ipat0, ckspc,
     1                nncrit, fbetac, ne, ik0, cksp, fbeta, 
     1                ipotnn, ipot,
     1                xport, xheap, xheapr, xkeep, xcalcx)

         if (xport*ndeg .gt. xportx*ndegx)  then
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

c     print 600,  ngs, rcurr, nupr
c      write(messag,600)
c      call echo(messag)
  600 format (1x, i5, f12.6, i7, ' igs, rcurr, nupr')
c     write(80,601)  ngs, rcurr/2, nupr, nuprtt
  601 format (1x, i8, f12.6, 2i9)

      if (.not. last) goto 100

      if (ipr2 .ne. 5)  close (unit=1)
c     delete paths.bin when done...
      close (unit=3, status='delete')
      close (unit=4)

       write(messag, 620) nuptot, nptot
       call echo(messag)
  620 format ('    Unique paths', i7, ',  total paths', i8)

c     Do not let user accidently fill up their disk
      if (nuptot .gt. 10000)  then
         call echo(' You have found more than 10000 paths.')
         call echo(' To continue this calculation, restart')
         call echo(' with current paths.dat and at GENFMT')
         call echo(' 3rd module of CONTROL card')
         call fstop(' at PATHSD: user must verify large run.')
      endif
      return
  999 call fstop( ' at PATHSD: no input')
      end
