c
c  this holds simple replacements for ifeffit routines
c  to be used by the 'libxafs' routines
c  
c  included in this file are:
c      sca_init  setsca  getsca
c
c  IMPORTANT:  DO NOT link into libifeffit.a!!!
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
       subroutine sca_init
       implicit none
       include 'compat.h'
       integer i
       do 10 i = 1, mxsca
          sscanam(i) = ''
          sscaval(i) = 0.d0
 10    continue 
       end

       subroutine setsca(str,x)
       implicit none
       character*(*) str
       double precision x
       integer i, ilen, istrln
       include 'compat.h'
       external istrln

       snamtmp = str
       call triml(snamtmp)
       call lower(snamtmp)
       ilen  = istrln(snamtmp)
       do 10 i = 1, mxsca
          if ((snamtmp(1:ilen) .eq. sscanam(i)(1:ilen)) .or.
     $         ('' .eq. sscanam(i)(1:ilen))) go to 20
 10    continue 
       call warn(3,"error: setsca out of memory")
       return
 20    continue
       sscanam(i) = snamtmp(1:ilen)
       sscaval(i) = x
       return
       end

       double precision function getsca(str,iwarn)
       implicit none
       character*(*) str
       integer i, ilen, istrln, iwarn

       include 'compat.h'
       external istrln

       getsca = 0.d0
       snamtmp = str
       call triml(snamtmp)
       call lower(snamtmp)
       ilen  = istrln(snamtmp)
       do 10 i = 1, mxsca
          if (snamtmp(1:ilen) .eq. sscanam(i)(1:ilen))
     $         getsca = sscaval(i)
 10    continue 
       return
       end
