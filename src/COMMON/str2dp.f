      subroutine str2dp(str,dpval,ierr)
c  return dp number "dpval" from character string "str"
c  if str cannot be a number, ierr < 0 is returned.
      character*(*) str, fmt*15 
      double precision dpval
      integer  ierr , lenmax
      parameter ( lenmax = 40)
      logical  isnum
      external isnum
      ierr = -99
      if (isnum(str)) then
         ierr = 0
         write(fmt, 10) min(lenmax, len(str))
 10      format('(bn,f',i3,'.0)')
         read(str, fmt, err = 20, iostat=ierr) dpval
      end if    
      if (ierr.gt.0) ierr = -ierr
      return
 20   continue
      ierr = -98
      return
c end subroutine str2dp
      end
      subroutine str2re(str,val,ierr)
c  return real from character string "str"
      character*(*) str 
      double precision dpval
      real     val
      integer  ierr
      call str2dp(str,dpval,ierr)
      if (ierr.eq.0) val = real(dpval)
      return
c end subroutine str2re
      end
      subroutine str2in(str,intg,ierr)
c  return integer from character string "str"
c  returns ierr = 1 if value was clearly non-integer
      character*(*) str 
      double precision val, tenth
      parameter (tenth = 1.d-1)
      integer  ierr, intg
      call str2dp(str,val,ierr)
      if (ierr.eq.0) then
         intg = int(val)
         if ((abs(intg - val) .gt. tenth))  ierr = 1
       end if
      return
c end subroutine str2in
      end
