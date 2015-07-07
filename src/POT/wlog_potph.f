      subroutine wlog (string)
      character*(*) string

c      include '../HEADERS/parallel.h'

c     This replaces wlog for use in libpotph


   10 format (a)

c     Suppress output in sequential loops
c      if (par_type .eq. 2) return

      il = istrln (string)
      if (il .eq. 0)  then
         print 10
c         if (par_type .ne. 3) write(11,10)
      else
         print 10, string(1:il)
c         if (par_type .ne. 3) write(11,10) string(1:il)
      endif
      return
      end
