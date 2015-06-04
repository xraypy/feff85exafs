      subroutine sidx (rholap, npts, rmt, rnrm, imax, imt, inrm)

      implicit double precision (a-h, o-z)
      dimension rholap (npts)
      character*512 slog
c#mn
      external ii, rr

      imt = ii (rmt)
      inrm = ii (rnrm)

c     Set imax (last non-zero rholap data)
      do 220  i = imt, npts
         if (rholap(i) .le. 1.0d-5)  goto 230
         imax = i
  220 continue
  230 continue

c     We need data up to the norman radius, so move norman
c     radius if density is zero inside rnrm.
      if (inrm .gt. imax)  then
         inrm = imax
         rnrm = rr (inrm)
  232    format(a,1pe13.5)
         write(slog,232) ' Moved rnrm.  New rnrm (au) ', rnrm
         call wlog(slog)
      endif
      if (imt .gt. imax)  then
         imt = imax
         rmt = rr (imt)
         write(slog,232) ' Moved rmt.  New rmt (au) ', rmt
         call wlog(slog)
      endif
      return
      end
