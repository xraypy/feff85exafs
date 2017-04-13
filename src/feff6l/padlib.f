c
c PAD library:   Packed Ascii Data 
c   these routines contain code for handling packed-ascii-data  
c   (pad) arrays for writing printable character strings that 
c   represent real or complex scalars and arrays to a file.
c
c routines included in padlib are (dp==double precision):
c   wrpadd     write a dp array as pad character strings
c   wrpadx     write a dp complex array as pad character strings
c   rdpadr     read a pad character array as a real array
c   rdpadd     read a pad character array as a dp  array
c   rdpadc     read a pad character array as a complex array
c   rdpadx     read a pad character array as a dp complex array
c   pad        internal routine to convert dp number to pad string
c   unpad      internal routine to pad string to dp number
c
c routines not included, but required by padlib:
c     triml, istrln, echo, iread, sclean
c
c//////////////////////////////////////////////////////////////////////
c Copyright (c) 1997--2001 Matthew Newville, The University of Chicago
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
       subroutine wrpadd(iout,npack,array,npts)
c
c write a dp array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
       include 'padlib.h'
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       double precision array(*), xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = array(i)
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--

       subroutine wrpadx(iout,npack,array,npts)
c write complex*16 array as pad string
       include 'padlib.h'
       integer    iout, npack, npts, mxl, js, i
       complex*16 array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = dimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadr(iout,npack,array,npts)
c
c write a real array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
       include 'padlib.h'
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       real    array(*)
       double precision xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = dble(array(i))
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadc(iout,npack,array,npts)
c write complex (*8) array as pad string
       include 'padlib.h'
       integer    iout, npack, npts, mxl, js, i
       complex    array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = aimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine rdpadd(iou,npack,array,npts)
c read dparray from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
       include 'padlib.h'
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       double precision    array(*), unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = tmp
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call warn(3,' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call echo(str(:i))
       call fstop(' -- bad data in PAD data file -- ')
       end
c --padlib--
       subroutine rdpadr(iou,npack,array,npts)
c read real array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
       include 'padlib.h'
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       real    array(*)
       double precision unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = real(tmp)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call warn(3, ' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call echo(str(:i))
       call fstop(' -- bad data in PAD data file -- ')
       end
c --padlib--
       subroutine rdpadc(iou,npack,array,npts)
c read complex array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
       include 'padlib.h'
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call warn(3, ' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call echo(str(:i))
       call fstop(' -- bad data in PAD data file -- ')
       end
       subroutine rdpadx(iou,npack,array,npts)
c read complex*16 array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
       include 'padlib.h'
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex*16  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call warn(3,' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call echo(str(:i))
       call fstop(' -- bad data in PAD data file -- ')
       end

c --padlib--
       subroutine pad(xreal,npack,str)
c  convert dp number *xreal* to packed-ascii-data string *str*
       include 'padlib.h'
       integer  iexp, itmp, isgn, i, npack, iok, j
       double precision xreal, xwork, xsave
       double precision tenth, tenmm, half, small
       parameter (tenth = 0.1d0, tenmm = 0.099999999994d0)
       parameter (half  = 0.5d0, small = 1.d-10 )
       character str*(*)
c
       str      = ' '
       xsave    = min(huge, max(-huge, xreal))
       isgn     = 1
       if (xsave.le.0) isgn = 0
c
       xwork    = dabs( xsave )
       iexp     = 0
       if ((xwork.lt.huge).and.(xwork.gt.tiny))  then
          iexp  =   1 + int(log(xwork) / tenlog  )
       else if (xwork.ge.huge) then
          iexp  = ihuge
          xwork = one
       else if (xwork.le.tiny)  then
          xwork = zero
       end if
c force xwork between ~0.1 and ~1
c note: this causes a loss of precision, but 
c allows backward compatibility
       xwork    = xwork / (ten ** iexp)
 20    continue
       if (xwork.ge.one) then
          xwork = xwork * tenth
          iexp  = iexp + 1
       else if (xwork.le.tenmm) then
          xwork = xwork * ten
          iexp  = iexp - 1
       endif
       if (xwork.ge.one) go to 20

       itmp     = int ( ibas2 * xwork ) 
       str(1:1) = char(iexp  + ioff + ibas2 )
       str(2:2) = char( 2 * itmp + isgn + ioff)
       xwork    = xwork * ibas2 - itmp
       if (npack.gt.2) then
          do 100 i = 3, npack
             itmp     = int( base * xwork + small)
             str(i:i) = char(itmp + ioff)
             xwork    = xwork * base - itmp
 100      continue
       end if
       if (xwork.ge.half) then
          i = itmp + ioff + 1
          if (i.lt.(ibase+ioff)) then
             str(npack:npack)= char(i)
          else 
             j = ichar(str(npack-1:npack-1))
             if (j.lt.(ibase+ioff-1)) then
                str(npack-1:npack-1) = char(j+1)
                str(npack:npack)     = char(ioff)
             endif 
          endif
       endif
       return
       end
c --padlib--
       double precision function unpad(str,npack)
c
c  convert packed-ascii-data string *str* to dp number *unpad*
       include 'padlib.h'
       double precision sum
       integer   iexp, itmp, isgn, i, npack
       character str*(*)
       unpad = zero
       if (npack.le.2) return
       iexp  =     (ichar(str(1:1)) - ioff   ) - ibas2
       isgn  = mod (ichar(str(2:2)) - ioff, 2) * 2 - 1
       itmp  =     (ichar(str(2:2)) - ioff   ) / 2
       sum   = dble(itmp/(base*base))
       do 100 i = npack, 3, -1
          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
 100   continue
       unpad = 2 * isgn * ibase * sum * (ten ** iexp)
cc       print*, sum, iexp,unpad
       return
       end
c --padlib--
c end of pad library
c ----------
