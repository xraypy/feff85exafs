      subroutine timrep (npat, ipat, rx, ry, rz, dhash,
     1            ipol, ispin, evec, xivec,eels)   !KJ added eels 5/06

c     subroutine timrev(...) is modified for polarization case 
c     Time-orders path and returns path in standard order,
c     standard order defined below.
c     Input:  npat, ipat
c     Output: ipat in standard order (time reversed if necessary)
c             rx, ry, rz   contain x,y,z coordinates of the path atoms,
c             where z-axis is along polarization vector or first leg, if
c               running usual feff,
c             x-axis is chosen so that first atom, which does not lie on
c               z-axis, lies in xz-plane,
c               for elliptically polarized light, x-axis is along the
c               incidence direction
c             y-axis is cross product of two previos unit vectors
c             Standarrd order is defined so that first nonzero x,y and z
c             coords are positive.(Otherwise we use the inversion of
c             the corresponding unit vector)
c             dhash double precision hash key for path in standard
c                order

      include '../HEADERS/dim.h'
      double precision evec(3), xivec(3)
      common /atoms/ rat(3,0:natx), ipot(0:natx), ilb(0:natx)
      dimension ipat(npatx+1), rx(npatx), ry(npatx), rz(npatx)
      dimension ipat0(npatx+1), rx0(npatx), ry0(npatx), rz0(npatx)
      integer eels,icase  !KJ added 5/06

      double precision dhash, dhash0

c     Time reverses path if time reversing it will put it
c     in standard order.  Standard order is defined by min hash
c     number, using path hash algorithm developed for the path
c     degeneracy checker.  See subroutine phash for details.
c     Symmetrical paths are, of course, always standard ordered.
c     Also returns hash number for standard ordered path.

c     Use suffix 0 for (') in variable names

!KJ next block : prepare new calling argument for mpprmp
      icase = -1
      if(eels.eq.1) then
         icase=7
      else
         icase=-1
      endif
!KJ end new block. 5/06


c     If no time-reversal standard ordering needed, make hash number
c     and return.  No timrev needed if 2 leg path (symmetrical).
      nleg = npat + 1
      ipat(nleg) = 0
      do 10 i = 1, npatx
         rx(i)   = 0
         ry(i)   = 0
         rz(i)   = 0
         rx0(i)   = 0
         ry0(i)   = 0
         rz0(i)   = 0
   10 continue
      call mpprmp(npat, ipat, rx, ry, rz,
     1            ipol, ispin, evec, xivec,icase)   !KJ added icase 5/06
      call phash (npat, ipat, rx, ry, rz, dhash)

      if (npat .le. 1)  then
         return
      endif

c     Make time reversed path

      ipat0(nleg) = ipat(nleg)
      do 210  i = 1, npat
         ipat0(i) = ipat(nleg-i)
  210 continue
      call mpprmp(npat, ipat0, rx0, ry0, rz0,
     1            ipol, ispin, evec, xivec,icase)  !KJ added icase 5/06
      call phash (npat, ipat0, rx0, ry0, rz0, dhash0)

c     turn off path reversal in special cases (make dhash0>dhash)
      if (ispin.ne.0 .and. ipol.ne.0) dhash0 = dhash+1

c     Do the comparison using hash numbers
c     Want representation with smallest hash number
      if (dhash0 .lt. dhash)  then
c        time reversed representation is smaller, so return
c        that version of the path
         dhash = dhash0
         do 300  i = 1, npat
            ipat(i) = ipat0(i)
            rx(i)   = rx0(i)
            ry(i)   = ry0(i)
            rz(i)   = rz0(i)
  300    continue
      endif

      return
      end
