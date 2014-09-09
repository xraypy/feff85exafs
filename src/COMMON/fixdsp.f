      subroutine fixdsp (dxorg, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect code.

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

      dimension dgc0(251), dpc0(251)
      dimension dgcx(nrptx), dpcx(nrptx)

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
      imax = 0
      do 100  i = 251, 1, -1
         if ( abs(dgc0(i)) .ge. 1.0d-11 .or. 
     1        abs(dpc0(i)) .ge. 1.0d-11 )  then
            imax = i
            goto 16
         endif
  100 continue
      call wlog(' Should never see this line from sub fixdsp')
   16 continue
c     jmax is the first point where both dpc and dgc are zero in
c     the original grid
      jmax = min(imax, 250) + 1

      delta = dxorg
      do 10  j = 1, jmax
         xorg(j) = xxx(j)
   10 continue
      rmax = rrr(jmax)

c     How far out do we go in the new grid?  To the last new grid
c     point before jmax.  Everything will be zero beyond jmax.
      delta = dxnew
      jnew = jjj(rmax)
      do 20  j = 1, jnew
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of rmax
      do 30  j = 1, jnew
         call terp (xorg, dgc0,  jmax, 3, xnew(j), dgcx(j))
         call terp (xorg, dpc0,  jmax, 3, xnew(j), dpcx(j))
   30 continue

c     and zero the arrays past rmax
      do 32  j = jnew+1, nrptx
         dgcx(j) = zero
         dpcx(j) = zero
   32 continue

      return
      end
