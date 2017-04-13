c//////////////////////////////////////////////////////////////////////
c Copyright (c) 1997--2000 Matthew Newville, The University of Chicago
c Copyright (c) 1992--1996 Matthew Newville, University of Washington
c
c Permission to use and redistribute the source code or binary forms of
c this software and its documentation, with or without modification is
c hereby granted provided that the above notice of copyright, these
c terms of use, and the disclaimer of warranty below appear in the
c source code and documentation, and that none of the names of The
c University of Chicago, The University of Washington, or the authors
c appear in advertising or endorsement of works derived from this
c software without specific prior written permission from all parties.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
c CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
c TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c SOFTWARE OR THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
c//////////////////////////////////////////////////////////////////////

       subroutine getfln(strin, filnam, ierr)
c  strip off the matched delimeters from string, as if getting
c  a filename from "filename", etc.
       implicit none
       integer idel, iend, istrln, ierr, ilen
       character*(*) strin, filnam, tmp*144, ope*8, clo*8
       external istrln
       data ope, clo /'"{(<''[',  '"})>'']'/

c
       ierr  = 0
       tmp   = strin
       call triml(tmp)
       ilen  = istrln(tmp)
       idel  = index(ope,tmp(1:1))
       if (idel.ne.0) then
          iend = index(tmp(2:), clo(idel:idel) )
          if (iend.le.0) then
             ierr = -1
             iend = ilen 
          end if
          filnam = tmp(2:iend)
       else
          iend = index(tmp,' ') - 1
          if (iend.le.0) iend  = istrln(tmp) 
          filnam = tmp(1:iend)
       end if
       return
c end  subroutine getfln
       end

       subroutine newfil(file, iofile)
c  
c  open a new file to unit iofile
c     if iofile > 0 , that file is closed
c     if an old file named file exists, it is deleted!
       implicit none
       character*(*) file, str*256
       integer   iofile, iex, ier
       logical   exist
       str  = file
       if (iofile.gt.0) then 
          close(iofile)
cc          iofile = 0
       end if
       inquire(file=str, exist=exist)
       if (exist) then 
          call openfl(iofile, str, 'old', iex, ier)
          close(iofile,status='delete')
cc          iofile = 0
       end if
cc       iofile = 3
       call openfl(iofile, str, 'unknown', iex, ier)
       if ((iex.lt.0).or. (ier.ne.0))  iofile = -1
c end subroutine newfil
       return
       end
       subroutine openfl(iunit, file, status, iexist, ierr)
c  
c  open a file, 
c   if unit <= 0, the first unused unit number greater than 7 will 
c                be assigned.
c   if status = 'old', the existence of the file is checked.
c   if the file does not exist iexist is set to -1
c   if the file does exist, iexist = iunit.
c   if any errors are encountered, ierr is set to -1.
c
c   note: iunit, iexist, and ierr may be overwritten by this routine
       implicit none
       character*(*)  file, status, stat*10
       integer    iunit, iexist, ierr
       logical    exist, open
c
c make sure there is a unit number, and that it's pointing to
c an unopened logical unit number other than 5 or 6
       ierr   = -3
       iexist =  0
       iunit  = max(1, iunit)
 10    continue 
       inquire (unit=iunit, opened=open)
       if (open) then
          iunit = iunit + 1
          if ((iunit.eq.5).or.(iunit.eq.6)) iunit = 7
          goto 10
       endif
c
c if status = 'old', check that the file name exists
       ierr = -2
       stat =  status                          
       call lower(stat)
       if (stat.eq.'old') then
          iexist = -1
          inquire(file=file, exist=exist)
          if (.not.exist) return
          iexist = iunit
       end if
c 
c open the file
       ierr = -1
cc       print*, ' openfl, unit ', iunit, ' file ', file(:40)
       open(unit=iunit, file=file, status=status, err=100)
       ierr = 0
 100   continue
       return
c end  subroutine openfl
       end
