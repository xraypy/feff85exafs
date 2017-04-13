      subroutine sidx (rholap, npts, rmt, rnrm, imax, imt, inrm)

      implicit double precision (a-h, o-z)
      dimension rholap (npts)
       character messag*128
      imt = ii (rmt)
      inrm = ii (rnrm)

c     Set imax (last non-zero rholap data)
      do 220  i = 1, npts
         if (rholap(i) .le. 1.0d-5)  goto 230
         imax = i
  220 continue
  230 continue

c     We need data up to the norman radius, so move norman
c     radius if density is zero inside rnrm.
      if (inrm .gt. imax)  then
         inrm = imax
         rnrm = rr (inrm)
         write(messag,'(1x,a,g15.6,a)')
     $        ' sidx: moved rnrm to ', rnrm, ' au '
         call echo(messag)
       endif
      if (imt .gt. imax)  then
         imt = imax
         rmt = rr (imt)
         write(messag,'(1x,a,g15.6,a)')
     $        ' sidx: moved rmt to ', rmt, ' au '
         call echo(messag)
      endif
      return
      end
