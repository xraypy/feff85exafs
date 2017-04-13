      subroutine ccrit (npat, ipat, ckspc,
     1    fbetac, rmax, pcrith, pcritk, nncrit, ipotnn, ipot,
     2    rpath, lheap, lkeep, xcalcx)

c     lheap to add to heap, lkeep if keep path at output.
c     NB, if lheap is false, lkeep is not used (since path
c     won't be in the heap).

      include 'const.h'
      include 'dim.h'
      logical lheap, lkeep
      dimension ipat(npatx)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)

c     local variables
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)

c     mrb is efficient way to get only ri and beta
c     note that beta is cos(beta)
      call mrb (npat, ipat, ri, beta)

      rpath = 0
      do 300  i = 1, npat+1
         rpath = rpath + ri(i)
  300 continue

c     If we can decide only on rpath, do it here...
      if (rpath .gt. rmax)  then
         lheap = .false.
         lkeep = .false.
         return
      endif

c     If last atom central atom, do put in heap, don't use it
c     as an actual path at output
      if (ipat(npat).eq.0)  then
         lheap = .true.
         lkeep = .false.
         return
      endif

c     Make index into fbetac array (this is nearest cos(beta) grid 
c     point, code is a bit cute [sorry!], see prcrit for grid).
      do 290  i = 1, npat+1
         tmp = abs(beta(i))
         n = tmp / 0.025
         del = tmp - n*0.025
         if (del .gt. 0.0125)  n = n+1
         if (beta(i) .lt. 0)  n = -n
         indbet(i) = n
  290 continue

c     Decide if we want the path added to the heap if necessary.
c     (Not necessary if no pcrith in use.)
      if (pcrith .gt. 0)  then

         call mcrith (npat, ipat, ri, indbet,
     1                ipot, nncrit, fbetac, ckspc, xheap)

c        xheap = -1 if not defined for this path (too few legs, etc.)
         if (xheap .ge. 0  .and.  xheap .lt. pcrith)  then
c           Do not want path in heap
            lheap = .false.
            lkeep = .false.
            return
         endif
      endif
c     Keep this path in the heap
      lheap = .true.

c     We may want path in heap so that other paths built from this
c     path will be considered, but do not want this path to be
c     written out for itself.  Decide that now and save the flag
c     in the heap, so we won't have to re-calculate the mpprm
c     path parameters later.

c     Skip calc if pcritk < 0
      if (pcritk .le. 0)  then
         lkeep = .true.
         return
      endif

c     Make xout, output inportance factor.
      call mcritk (npat, ipat, ri, beta, indbet,
     1             ipot, nncrit, fbetac, ckspc, xout, xcalcx)

c     See if path wanted for output
c     Do not want it if last atom is central atom (xout = -1) or
c     if xout is too small
      lkeep = .false.
      if (xout .ge. pcritk)  lkeep = .true.

      return
      end
