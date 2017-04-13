      subroutine fixdsx (iph, dxorg, dxnew, dgc, dpc, dgcn, dpcn)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect and phase codes.

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8d0)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = int((log(r) + xx00) / delta) + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.

      delta = dxorg
      do 10  j = 1, 251
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      do 20  j = 1, nrptx
         xnew(j) = xxx(j)
   20 continue

      do 200 iorb = 1, 30
         imax = 0
         do 100  i = 251, 1, -1
            if ( abs(dgc(i,iorb,iph)) .ge. 1.0d-11 .or.
     1           abs(dpc(i,iorb,iph)) .ge. 1.0d-11 )  then
               imax = i
               goto 16
            endif
  100    continue
   16    continue
         if (imax .eq. 0) then
            jnew = 0
            goto 35
         endif
c        jmax is the first point where both dpc and dgc are zero in
c        the original grid
         jmax = imax + 1
         if (jmax .gt. 251) jmax = 251

         delta = dxorg
         rmax = rrr(jmax)

c        How far out do we go in the new grid?  To the last new grid
c        point before jmax.  Everything will be zero beyond jmax.
         delta = dxnew
         jnew = jjj(rmax)

c        interpolate to new grid using x, only inside of rmax
         do 30  j = 1, jnew
            call terp(xorg,dgc(1,iorb,iph),jmax,3, xnew(j),dgcn(j,iorb))
            call terp(xorg,dpc(1,iorb,iph),jmax,3, xnew(j),dpcn(j,iorb))
   30    continue

c        and zero the arrays past rmax
   35    do 40  j = jnew+1, nrptx
            dgcn(j,iorb) = 0
            dpcn(j,iorb) = 0
   40    continue
  200 continue

      return
      end
