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
       integer function iread(iunit,line)
c
c reads line from an open file unit (iunit)
c  return values:
c   line length on success
c            -1 on 'end'
c            -2 on 'error'
       implicit none
       character*(*) line
       integer    iunit, istrln, ilen
       external   istrln
       line = ' '
 10    continue 
       read(iunit, '(a)', end = 40, err = 50) line
       call sclean(line)
       call triml(line)
       iread = istrln(line)
       if (iread .eq. 0) goto 10
       return
 40    continue 
       ilen = istrln(line)
       if (ilen.ge.1) then
          call sclean(line)
          call triml(line)
          iread = ilen
       else 
          line = ' '
          iread= -1
       endif
          
       return
 50    continue 
       line = ' '
       iread = -2
       return
       end
       integer function iread_ky(iunit,key,line)
c
c reads line from an open file unit (iunit)
c and extracts a 2character key (as for PAD files)
c return values:
c   line length on success
c            -1 on 'end'
c            -2 on 'error'
       implicit none
       character*(*) line, key
       integer    iunit, iread, ilen
       external    iread
       key = ' '
       line = ' '
       ilen = iread(iunit, line)
       if (ilen.gt.2) then
          key  = line(1:2)
          line = line(3:)
          ilen = ilen - 2
       endif
       iread_ky = ilen
       return
       end

       subroutine sclean(str) 
c
c  clean a string so that all: 
c     char(0), and char(10)...char(15) are end-of-line comments,
c        so that all following characters are explicitly blanked.
c     all other characters below char(31) (including tab) are
c        replaced by a single blank
c
c  note that this is mostly useful when getting a string generated
c  by a non-fortran process (say, a C program) and for dealing with
c  dos/unix/max line-ending problems
       character*(*) str, blank*1
       parameter (blank = ' ')
       integer i,j,is
       do 20 i = 1, len(str)
          is = ichar(str(i:i))
          if ((is.eq.0) .or. ((is.ge.10) .and. (is.le.15))) then
             do 10 j= i, len(str)
                str(j:j) = blank
 10          continue
             return
          endif
          if (is.le.31) str(i:i)  = blank
 20    continue 
       return
c end subroutine sclean
       end
