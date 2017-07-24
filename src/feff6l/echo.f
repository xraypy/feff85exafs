       subroutine echo(str)
c
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
c
c dump string to standard output
       implicit none
       include 'echo.h'
       character*(*) str,form*8
       parameter (form = '(1x)' )
       call chrdmp(str)
       if (mod(i_echo,2).eq.1) write(*,form)
       return
c end  subroutine echo
       end
       subroutine chrdmp(str)
c dump string to screen with "$" format
c
c i_echo
c  0   echo to buffer
c  1   echo to screen
c  2   echo to open echo_file (if open)
c  3   echo to screen and open echo_file (if open)
c
       implicit none
       include 'echo.h'
       character*(*) str, out*512,frm*8,ffrm*8
       parameter (frm= '(1x,a,$)', ffrm= '(1x,a)' )
       integer  istrln, n
       external istrln
       out = str
       n   = max(1, istrln(out))
       if (i_echo.eq.0) then
          call echo_push(out)
       else
          if (mod(i_echo,2).eq.1) write(*,frm) out(1:n)
          if ((i_echo.ge.2).and.(lun_echo.ge.1))
     $         write(lun_echo, ffrm) out(1:n)
       endif
       return
c  end subroutine chrdmp
       end

       subroutine echo_init
c initialize echo lines
       implicit none
       include 'echo.h'
       integer i
       do 20 i  = 1, mxecho
          echo_str(i) = ' '
 20    continue
       call setsca('&echo_lines', 0.d0)
       n_echo = 0
       call setsca('&screen_echo', 1.d0)
       i_echo    = 1
       lun_echo  = 0
       echo_file = ''
       return
       end

       subroutine close_echofile()
       implicit none
       include 'echo.h'
       if (lun_echo .ge. 1) then
          close(lun_echo)
          lun_echo  = -1
          echo_file = ''
          if (i_echo.eq.3) i_echo  = 1
          if (i_echo.eq.2) i_echo  = 0
       endif
       return
       end

       subroutine open_echofile(s)
       implicit none
       include 'echo.h'
       character*(*) s
       integer  iex, ier, istrln
       external istrln
       call close_echofile()

       lun_echo = 19
       echo_file = s(1:istrln(s))
       call triml(echo_file)

       call openfl(lun_echo, echo_file, 'unknown', iex, ier)
       if (i_echo.eq.0) i_echo = 2
       if (i_echo.eq.1) i_echo = 3
cc       print*, ' done ' , i_echo
c
       return
       end



       subroutine echo_push(string)
c add echo string to internal list
c copyright (c) 1999 matt newville
       implicit none
       character*(*) string, str*512
       include 'echo.h'
       integer  istrln, ilen, i
       external istrln
       str  = string
       call sclean(str)
       call triml(str)
       ilen = istrln(str)
       if ((ilen.ge.1).and.(n_echo.lt.mxecho)) then
          do 30 i = mxecho, 2, -1
             echo_str(i) = echo_str(i-1)
 30       continue
          echo_str(1) = str(1:ilen)
cc          print*, ' ECHO_PUSH: ', str(1:ilen)
          n_echo      = min(mxecho, n_echo + 1)
       endif
       call setsca('&echo_lines', n_echo * 1.d0)
       return
c  end subroutine echo_push
       end
       subroutine echo_pop(string)
c add echo string to internal list
c copyright (c) 1999 matt newville
       implicit none
       character*(*) string
       include 'echo.h'
       string  = ' '
       if (n_echo .gt. 0) then
          string  = echo_str(n_echo)
cc          print*, ' ECHO_POP: ', string(1:60)
          echo_str(n_echo)  =  ' '
       end if
       n_echo  = min(mxecho, max(0, n_echo - 1))
       call setsca('&echo_lines', n_echo * 1.d0)
       return
c  end subroutine echo_pop
       end



c
c routines to initialize and use a 'stop file'
c that is, do
c         call fstop_init('stopfile.err')
c early on, and replace subsequent
c         stop 'problem at line xxx'
c with
c         call fstop('problem at line xxx')
c
c the error file will contain the error message
c and will exist only if fstop() has been called.
c that is, fstop_init() erases the stop file.
c
       subroutine fstop_init(s)
       character*(*) s
       character*32 stopfilename
       common /stop_file/ stopfilename
       integer istrln
       external istrln
       stopfilename = s
       call triml(stopfilename)

       return
       end

       subroutine fstop(s)
c
c replaces intrinsic 'stop', writing stop message to
c the stopfile initialized by fstop_init
c
       implicit none
       character*(*) s, str*128
       integer i, istrln
       external istrln
       character*32 stopfilename
       common /stop_file/ stopfilename

       str = s
       call triml(str)
       if (str .eq. '') str = 'unknown error'
       str   = 'Fatal Error: '//str(1:istrln(str))
       call echo(str(1:istrln(str)))

       call triml(stopfilename)
       if (istrln(stopfilename).ge.1) then
          i = 9
          call newfil(stopfilename,i)
          write(i,'(1x,a)') str(1:istrln(str))
          close (unit=i)
       endif
       stop
       return
       end

       subroutine warn(i,s)
c
c set &status value and write string to echo buffer
       integer i
       character*(*) s
       call echo(s)
       call set_status(i)
       return
       end

       subroutine set_status(i)
c
c set &status value
       integer i
       double precision getsca, xc, xi
       xi = i * 1.d0
       xc = getsca('&status',0)
       if (xi .gt. xc) call setsca('&status',xi)
       return
       end
