      subroutine paths (ckspc, fbetac, pcritk, pcrith, nncrit,
     1                  rmax, nlegxx, ipotnn)

c     finds multiple scattering paths
c     This is single precision, units are Angstroms.  BE CAREFUL!

c     pcrith is cut-off fraction used when building paths
c            (path criterion for heap)
c     pcritk is cut-off fraction used on output
c            (path criterion for keeping)

c     ipotnn is output, used by pathsd to duplicate paths criteria,
c     which are used only for diagnostic output.

      include 'dim.h'
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)

      include 'vers.h'

c     This common in pathsd, mpprm
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension m(-1:natx,0:natx)
      dimension mindex(natx+1)
c     Used for packed integers
      dimension iout(3)

c     ok true if all paths to rmax found.  If heap full, npx exceeded,
c     etc., last general shell may be incomplete, set ok=.false.
      logical ok

c     Heap data structure:
c     index is the pointer to the element of the data structure.
c     Each element contains
c        r        total path length
c                 Note that r is sorted along with index -- this keeps
c                 the heap maintenance routines fast.
c        mi, mj   m matrix elements used to place last atom in this path
c        npat     number of atoms in this path
c        ipat(npatx) indices of atoms in this path
c     next is the index of the next data structure element available.
c     If an element is freed, npat is the index of the free element
c     to use after using current "next" element.
c     nx is max number in heap
      integer    nx
      parameter (nx = 10 000)
c     parameter (nx = 60 000)
c     r also used in making m matrix, must have nx >= natx+1
      integer   index(nx), npx, np, n, ip, i
c     parameter (npx = 100 000)
      parameter (npx = 4 000 000)
      dimension r(nx), mi(nx), mj(nx)
      dimension npat(nx)
      dimension ipat (npatx,nx)
c     Keep this path on output
      logical keep1(nx), kp1tmp

c     Used with ipack, so need ipat(8)
      dimension ipat0(8)

c     paths are typically about 10 or 20 Ang
      parameter (big = 1.0e3)

      parameter (nheadx = 30)
      character*80  head(nheadx)
      character*80  title, messag*128
      dimension lhead(nheadx)

c     Returned from criterion checker, false if path fails criterion
       logical keep
       external sdist

c     read input
c     header...
c     i, x, y, z, ipot, i1b   of nat+1 atoms (i=0 is central atom)
      open (1, file='geom.dat', status='old', iostat=ios)
      call chopen (ios, 'geom.dat', 'paths')
      nhead = nheadx
      call rdhead (1, nhead, head, lhead)
c     header from geom.dat includes carriage control...
c     nlegxx is max number of legs user wants to consider.
c     nlegs = npat+1, so set npatxx = min (npatx, nlegxx-1)
      npatxx = min (npatx, nlegxx-1)
c     Input rmax is one-way distances
      rmax = rmax*2
      nat = -1
c     ratx is distance to most distant atom, used to check rmax
      ratx = 0
   10 continue
         nat = nat+1
         if (nat .gt. natx)  then
            write(messag,'(1x,a,2i5)') ' PATHS: nat, natx ',
     $           nat, natx
            call echo(messag)
            call fstop(' at PATHS: bad input')
         endif
         read(1,*,end=20)  idum, (rat(j,nat),j=1,3), ipot(nat), i1b(nat)
         rtmp = sdist(rat(1,nat),rat(1,0))
         if (rtmp .gt. ratx)  ratx = rtmp
      goto 10
   20 continue
      nat = nat-1
      close (unit=1)

c     Warn user if rmax > dist to most distant atom
      if (rmax/2 .gt. ratx+.02)  then
         call echo('   WARNING in PATHS Module:')
         call echo('      rmax > distance to most distant atom.')
         call echo('      Some paths may be missing.')
         write(messag,'(1x,a,2f12.5)') ' rmax, ratx ',
     $        rmax/2.d0, ratx
         call echo(messag)
      endif

c     Count number of 1st bounce atoms (at least 1 required).
      n1b = 0
      do 30  i = 1, nat
         if (i1b(i) .gt. 0)  n1b = n1b + 1
   30 continue
      if (n1b .lt. 1) call fstop(' at PATHS: at least one '//
     $     '1st bounce atoms required.')

       if (rmax .ge. big)
     $     call fstop(' at PATHS: get real with rmax!')

c     Make title for this run, include carriage control because head
c     (read above) includes carriage control.
      write(title,32)  rmax/2, pcritk, pcrith, vfeff, vpaths
   32 format(' Rmax', f8.4, ',  keep limit', f7.3,
     1       ', heap limit', f7.3, t57, 2a12)

       write(messag, 34) rmax/2, pcritk, pcrith
       call echo(messag)
   34 format ('    Rmax', f8.4,
     1        '  keep and heap limits', 2f12.7)

       call echo('   Preparing neighbor table')

c     prepare table telling distance from atom i to atom j and then
c     back to central atom
c     First bounce is m(-1,...), m(0,...) is bounces from central
c     atom that are not first bounces.
      do 60  i = -1, nat
         ir = i
         if (i .eq. -1)  ir = 0
         do 40  j = 0, nat
c           r begins with element 1 so sort routine later will work
            r(j+1) = sdist (rat(1,ir), rat(1,j))
            r(j+1) = r(j+1) + sdist (rat(1,j), rat(1,0))
c           we don't need m(i,i), since this will be = shortest
c           of the r(j), so just set it to something very big,
c           it will sort to the end of this row and it won't
c           bother us
            if (j .eq. ir)  r(j+1) = big
c           If we're doing first bounce, use only the allowed first
c           bounce paths.
            if (i .eq. -1)  then
               if (i1b(j) .le. 0)  r(j+1) = big
            endif
   40    continue

c        prepare row i of m table
c        m is a distance table ordered such that distance from
c               i to m(i,0) to 0 <
c               i to m(i,1) to 0 <
c               i    m(i,2)    0 <
c               :    :    :
c               i    m(i,nat)  0
c
c        That is, m(i,0) is index of atom that gives shortest path,
c                 m(i,1)                        next shortest path, etc.
c        Note that m(0,0) is shortest single bounce path.

c        Again, r and mindex go from 1 to nat+1, m goes from 0 to nat
         call sortir (nat+1, mindex, r)
         do 50  j = 0, nat
            m(i,j) = mindex(j+1)-1
   50    continue
   60 continue

       call echo('    nfound  nheap  nheapx  nsc    r')

c     initialize heap data space "next" pointers
      do 70  i = 1, nx-1
         npat(i) = i+1
   70 continue
      npat(nx) = -1
c     initial condition:  make the first path
c     n    number in heap
c     nna  number skipped counter
c     nhx  number used in heap max, a counter
      n = 1
      nna = 0
      nhx = n
      nwrote = 0
      index(n) = 1
      ip = index(n)
      next = 2
      mi(ip) = -1
      mj(ip) = 0
      npat(ip) = 1
      ipat(npat(ip),1) = m(mi(ip),mj(ip))

c     near neighbor is atom ipat(npat(ip),1) for first path into heap
      ipotnn = ipot(ipat(npat(ip),1))

c     Someday change keep and keep1 to lkeep and lheap to match
c     ccrit variable names.
c     Initialize keep criterion
      xcalcx = -1
      call ccrit (npat(ip), ipat(1,ip), ckspc,
     1    fbetac, rmax, pcrith, pcritk, nncrit, ipotnn, ipot,
     2    r(n), keep, keep1(ip), xcalcx)

      open (file='paths.bin', unit=3, access='sequential',
     1      form='unformatted', status='unknown', iostat=ios)
      call chopen (ios, 'paths.bin', 'paths')
c     These strings are all char*80 and include carriage control
      write(3) nhead+1
      do 88  ihead = 1, nhead
         write(3) head(ihead)
         write(3) lhead(ihead)
   88 continue
      write(3) title
      write(3) istrln(title)
      write(3)  nat
      do 90  i = 0, nat
         write(3) (rat(j,i),j=1,3), ipot(i), i1b(i)
   90 continue

c     r is the heap, index is the pointer to the rest of the data
c     np is the number of paths found and saved
      np = 0
c     nbx  mpat max (Number of Bounces maX)
      nbx = 0

c     done if path at top of heap is longer than longest path we're
c        interested in
c     done if max number of paths we want have been found
c     begin 'while not done' loop
      ok = .false.
  800 continue
         if (r(1) .gt. rmax  .or.  np .ge. npx .or. n.le.0)  then
c           n=0 means heap is empty
            if (n.le.0)  ok=.true.
c           if (n.le.0)  call echo('   Heap empty')
            goto 2000
         endif

c        save element at top of heap in arrays labeled 0
c        dump to unit 3 (unformatted)
         ip = index(1)
         npat0 = npat(ip)
         do 100  i = 1, npat0
            ipat0(i) = ipat(i,ip)
  100    continue
         r0 = r(1)

c        Don't write out path if last atom is central atom, or
c        if it doesn't meet pcritk
         if (ipat0(npat0).eq.0 .and. keep1(ip)) then
            write(messag,'(1x,a,2i6)')
     $           ' paths: odd case at npat0,ip= ', npat0,ip
            call echo(messag)
         endif
         if (ipat0(npat0).ne.0 .and. keep1(ip))  then
            np = np+1
c           pack integers
            call ipack (iout, npat0, ipat0)
            write(3)  r0, iout
            nwrote = nwrote+1
c           write status report to screen
            if (mod(np,1000) .eq. 0)  then
               write(messag, 132) np, n, nhx, nbx, r0/2
  132          format (4x, i6, i7, i8, i4, f10.4)
               call echo(messag)
            endif
         endif

         if (np .ge. npx)  then
            write(messag,'(1x,i5,a)') np,
     $           ' paths found. (np >= npx)'
            call echo(messag)
            goto 2000
         endif

c        Make new path by replacing last atom in path from top of heap,
c        put this path on top of heap and buble it down.  If row is
c        finished, or new path is too long, don't add it, instead
c        move last path in heap to the top.
c        If working on row mi=-1 (first bounce atoms), don't
c        use them if not allowed 1st bounce atoms.
         mj(ip) = mj(ip) + 1
         if (mi(ip).eq.-1  .and.  i1b(m(mi(ip),mj(ip))).le.0)  then
c           not allowed first bounce atom
            r(1) = big
            keep = .false.
c           call echo('1st bounce limit!')
         elseif (mj(ip) .ge. nat)  then
c           we've finished a row of m matrix
            r(1) = big
            keep = .false.
         else
c           new path has same indices, etc.  Only need to replace
c           last atom.
            ipat(npat(ip),ip) = m(mi(ip),mj(ip))
            call ccrit (npat(ip), ipat(1,ip), ckspc,
     1                  fbetac, rmax, pcrith, pcritk, nncrit,
     1                  ipotnn, ipot,
     2                  r(1), keep, keep1(ip), xcalcx)
         endif

c        If r is bigger than rmax or keep=false, remove element from
c        heap by taking the last element in the heap and moving it to
c        the top.  Then bubble it down.  When removing an element
c        from the heap, be sure to save the newly freed up index.
c        r(1) and index(1) are new path, set above
         if (r(1).gt.rmax .and. keep)  then
            call echo('odd case rmax...')
         endif
         if (r(1).gt.rmax .or. .not.keep)  then
            index(1) = index(n)
            r(1) = r(n)
c           use npat as pointer to next free location
            npat(ip) = next
            next = ip
            n = n-1
c           nna is Number Not Added to heap
            nna = nna + 1
c           Maybe heap may be empty here, but that's alright
         endif
         if (npat(index(1)).gt.nbx .and. n.gt.0)  nbx = npat(index(1))

c        If heap is empty, don't call hdown.
         if (n.gt.0)  call hdown (r, index, n)

c        and make a new path by adding an atom onto the end of the path
c        we saved, put this at the end of the heap and bubble it up.
c        Do this only if it won't be too many bounces.
         if (npat0+1 .le. npatxx)  then
            ip = next
            if (ip .lt. 0)  then
c              call echo('   Heap full')
               goto 2000
            endif
            next0 = npat(ip)
            do 200  i = 1, npat0
               ipat(i,ip) = ipat0(i)
  200       continue
            mi(ip) = ipat0(npat0)
            mj(ip) = 0
            npat(ip) = npat0+1
            ipat(npat(ip),ip) = m(mi(ip),mj(ip))
            call ccrit (npat(ip), ipat(1,ip), ckspc,
     1                  fbetac, rmax, pcrith, pcritk, nncrit,
     1                  ipotnn, ipot,
     2                  rtmp, keep, kp1tmp, xcalcx)
            if (rtmp .gt. rmax  .and.  keep)  then
               call echo(' odd case rmax and tmp...')
            endif
            if (rtmp .gt. rmax  .or.  .not.keep)  then
               npat(ip) = next0
               nna = nna+1
            else
c              add it to the heap
               next = next0
               n = n+1
               if (n .gt. nhx)  nhx = n
               index(n) = ip
               r(n) = rtmp
               keep1(ip) = kp1tmp
               if (npat(index(n)) .gt. nbx)  nbx = npat(index(n))
               call hup (r, index, n)
            endif
         endif

      goto 800
 2000 continue
c     end of 'while not done' loop
      if (.not. ok)  then
         call echo('   Internal path finder limit exceeded -- '//
     1           'path list may be incomplete.')
      endif
      close (unit=3)
       write(messag, 2010) np, nhx, nbx
 2010 format ('    Paths found', i9, 3x,
     1        '(nheapx, nbx', i8, i4, ')')
       call echo(messag)
      end
