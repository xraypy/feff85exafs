       subroutine get_inpfile(defval,fname,istat)
c 
c  return first command line argument, as if for a filename,
c  using defval if no argument is given:
c
c  arguments 
c      defval   default file name  (character) input only
c      fname    output file name   (character) output only
c      istat    output status      (integer)   output only
c        = 0 for success,  = 1 for failure

       character*(*) defval, fname
       integer  n_args, istat, jstat
       character*512 cl_args(3), f
c
c g77 version
       istat = 1 
       n_args = iargc()
       f = defval
       if (n_args .ge. 1) then
          call getarg(1,f)
          istat = 0
       endif
       call triml(f)
       fname = f
       return 
       end

cc dvf version
c       use dflib
c       istat  = 1
c       n_args = nargs()
c       f      = defval
c       jstat  = 0
c
c       if (n_args .ge.2) then 
c          CALL GETARG(1,f, jstat)
c          istat = 0
c       endif	
c       call triml(f)
c       fname = f
c       return
c       end
c
