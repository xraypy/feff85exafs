      subroutine messer
c  prints error message on the output device
      implicit double precision (a-h,o-z)
      common/messag/dlabpr,numerr
      character*8 dlabpr
      character*512 slog
 
      ilig=numerr/1000
      ier=numerr-1000*ilig
      write(slog,'(a,i6,a,i6,a,a8)')  'error number ',ier,
     1 ' detected on a line ',ilig,'in the program',dlabpr
      call wlog(slog)
      return
      end
