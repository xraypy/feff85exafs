      subroutine genfmt (ipr3, critcw, iorder, wnstar,
     1                ipol, ispin, le2, angks, elpty, evec, xivec, ptz,
     2            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     
      implicit double precision (a-h, o-z)

c  altered by matt newville (jan 1999): 
c  format of feff.bin changed to packed-ascii, and all writes changed.
c  altered by alex ankudinov(feb 2000); disabled use of paths in LDOS.
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'
      include 'clmz.h'
      include 'fmatrx.h'
      include 'lambda.h'
      include 'pdata.h'
      include 'nlm.h'
      include 'rotmat.h'
      integer elnes,ipmin,ipmax,ipstep  !KJ added variables 1-06      
      dimension evec(3), xivec(3)
      complex*16 ptz
      dimension ptz(-1:1, -1:1), lind(8)

      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cfac
      complex*16  cchi(nex)
      complex*16  rkk(nex,8)
      complex*16  rkk2(nex,8,nspx)
      complex*16  eref2(nex,nspx), ph4, bmati
      dimension   ph4(nex,-ltot:ltot, nspx, 0:nphx)
      dimension   bmati(-mtot:mtot, 8, -mtot:mtot, 8)
      dimension   xk(nex), ckmag(nex)
      complex*16  ck(nex), ckp
      dimension   ffmag(nex)
      dimension   eps1(3), eps2(3), vec1(3), vec2(3)

      character*128 string
      character*512 slog
      logical done, wnstar
      
      integer ip !KJ local index that I added 1-06
      character*10 f1,f2 !KJ variable for opening feff.bin files - filename
      

c      padlib staff
       double precision phff(nex), amff(nex),  xkr(nex)
       integer  mpadx
       parameter (mpadx = 8)
       character*75 wfmt, atsym*2
       external atsym, cwig3j, istrln

c     Input flags:
c     iorder, order of approx in f-matrix expansion (see setlam)
c             (normal use, 2.  Do ss exactly regardless of iorder)

c     used for divide-by-zero and trig tests
       parameter (eps = 1.0e-16)
       external getxk, xstar
       
c     Read phase calculation input, data returned via commons
       call rdxsph (ne, ne1, ne3, npot, ihole, rnrmav, xmu, edge, ik0,
     1              em, eref2, iz, potlbl, ph4, rkk2, lmax, lmaxp1)
       call setkap (ihole, kinit, linit)
       ilinit = linit + 1

c
c      need to sum over spin-up and -down for |ispin|=1 (fix later)
       nsp = 1
       if (ispin.eq.1) nsp = nspx

       if (nsp.eq.1) then
c        for ispin=2 the variables already written into is=1 positions
         is = 1
         do 10 ie = 1, ne
  10     eref(ie) = eref2(ie,is)
         do 20 iph = 0, npot
         do 20 ie = 1, ne
         do 20 il = -lmax(ie, iph), lmax(ie, iph)
  20     ph(ie,il, iph) = ph4(ie, il, is, iph)
       else
c        average over two spin direction
         do 12 ie = 1, ne
  12     eref(ie) = (eref2(ie,1) + eref2(ie,nsp)) /2
!KJ  12     eref(ie) = (eref2(ie,1) + eref2(ie,2)) /2
         do 22 iph = 0, npot
         do 22 ie = 1, ne
         do 22 il = -lmax(ie, iph), lmax(ie, iph)
  22     ph(ie,il, iph) =(ph4(ie, il, 1,iph) + ph4(ie, il, nsp,iph)) /2
!KJ  22     ph(ie,il, iph) =(ph4(ie, il, 1,iph) + ph4(ie, il, 2,iph)) /2
       endif

c!KJ added loop over ip for ELNES calculations.  1-06
      do ip=ipmin,ipmax,ipstep
            if (elnes.eq.1) call iniptz(ptz,ip,2)  !KJ Only change ptz for ELNES !!
            !KJ the call to iniptz changes the polarization matrix ptz
c           !KJ choose different filename for each spectrum.
            if(ip.eq.1) then
              f1(1:10)='feff.bin  '
              f2(1:10)='list.dat  '
            elseif(ip.eq.10) then
              f1(1:10)='feff10.bin'
              f2(1:10)='list10.dat'
            elseif(ip.gt.1.and.ip.lt.10) then
              f1(1:5)='feff0'
              f1(6:6)= char(48+ip)
              f1(7:10)='.bin'
              f2(1:5)='list0'
              f2(6:6)= char(48+ip)
              f2(7:10)='.dat'
            else
              stop 'crazy ip in ff2xmu'
            endif
          !KJ the rest of this loop follows the old (non-elnes) scheme
       write(*,*) 'doing ip = ',ip



c     Open path input file (unit in) and read text .  Use unit 1.
       ntext  = 5
       open (unit=1, file='paths.dat', status='old', iostat=ios)
       call chopen (ios, 'paths.dat', 'genfmt')
       call rdhead (1, ntext , text, ltext)
       if (ntext  .le. 0)  then
          text (1) = ' '
       endif
c
c     Save indices of paths for use by ff2chi
c     Save indices of paths for use by ff2chi
       open (unit=2, file=f2, status='unknown', iostat=ios) !KJ changed 'list.dat' to f2
       call chopen (ios, f2, 'genfmt') !KJ id.  1-06
c     Put phase header on top of list.dat
       call wthead (2, ntext , text )
       write(2, 125)
 125   format (1x, 71('-'))
       write(2, 135)
 135   format ('  pathindex     sig2   amp ratio    ',
     1      'deg    nlegs  r effective')
       
c     Open nstar.dat if necessary
       if (wnstar)  then
          open (unit=4,file='nstar.dat', status='unknown', iostat=ios)
          call chopen (ios, 'nstar.dat', 'genfmt')
          write(4,'(1x,a,f8.4)' ) ' polarization', evec
          write(4,'(1x,a)' ) ' npath     n*'
       endif
       
c     Set crit0 for keeping feff.dat's
       if (ipr3 .le. 0)  crit0 = 2*critcw/3
c     Make a header for the running messages.
       write(slog, 155) critcw
 155   format ('    Curved wave chi amplitude ratio', f7.2, '%')
       call wlog(slog)
       if (ipr3 .le. 0)  then
         write(slog,165) crit0
         call wlog(slog)
       endif
 165   format ('    Discard feff.dat for paths with cw ratio <',
     1         f7.2, '%')
       write(slog,195)
 195   format ('    path  cw ratio     deg    nleg  reff')
       call wlog(slog)

c     open feff.bin for storing path info
c     for now, use double precision.  After it's working, try
c     single precision.
c     Use single precision for all fp numbers in feff.bin
c!KJ replaced by next lines      open (unit=3, file='feff.bin', status='unknown', iostat=ios)
c!KJ idem      call chopen (ios, 'feff.bin', 'genfmt')
      open (unit=3, file=f1, status='unknown', iostat=ios) !KJ f1  1-06
      call chopen (ios, f1 , 'genfmt') !KJ introduced f1
c     put label line in feff.bin so other programs know it really
c     is a feff.bin file
       string = '#_feff.bin v03: ' // vfeff
       jstr   = istrln(string)
       write(3, '(a)')  string(1:jstr)

c     save stuff that is the same for all paths
c     header, ck, central atom phase shifts
       write(3, '(a2,6(1x,i4))') '#_',
     $      npot, ne, mpadx

c     Misc stuff from phase.bin and genfmt call
 345   format(a2,3(1x,i7), 3(1x,g14.7))
       write(3, 345) '#&', ihole, iorder, ilinit, rnrmav, xmu, edge
       do 380 i = 0, npot
          if (potlbl(i).eq.' ') potlbl(i)  = atsym(iz(i))
          if (potlbl(i).eq.' ') potlbl(i)  = 'null'
 380   continue 
 395   format('(',i3,'(1x,a6),',i3,'(1x,i3))')
       write(wfmt, 395) npot+1, npot+1
       write(string,wfmt) (potlbl(i),i=0,npot) , (iz(i),i=0,npot)
       jstr = istrln(string)
       write(3, '(a2,a)') '#@',string(:jstr)

c     Central atom phase shifts
      ll = linit+1
      if (kinit.lt.0) ll = -ll
      call wrpadx(3,mpadx, ph(1,ll, 0),ne)

c     Set nlm factors in common /nlm/ for use later
      call snlm (ltot+1, mtot+1)

c     Make xk and ck array for later use
       do 850  ie = 1, ne
c        real momentum (k)
         xk(ie) = getxk (dble(em(ie)) - edge)
c        complex momentum (p)
         ck(ie) = sqrt (2*(em(ie) - eref(ie)))
         ckmag(ie) = abs(ck(ie))
         xkr(ie) = real(xk(ie))
 850   continue
       call wrpadx(3,mpadx, ck,ne)
       call wrpadd(3,mpadx, xkr,ne)
         
c     While not done, read path, find feff.
       npath  = 0
       ntotal = 0
       nused  = 0
       xportx = -1
 1000  continue

c        Read current path
         call rdpath (1, done, ipol)
         icalc = iorder
         if (.not.done)  then
            npath = npath + 1
            ntotal = ntotal + 1
            if (wnstar)  then
c              should be ipol=1
               do 1150 ic =1,3
                  vec1(ic) = rat(ic,1) - rat(ic,0)
                  vec2(ic) = rat(ic,nleg-1) - rat(ic,0)
                  eps1(ic) = evec(ic)
 1150          continue
               if (elpty.ne.0.0) then
                  eps2(1) = xivec(2)*evec(3)-xivec(3)*evec(2)
                  eps2(2) = xivec(3)*evec(1)-xivec(1)*evec(3)
                  eps2(3) = xivec(1)*evec(2)-xivec(2)*evec(1)
               endif
               ndeg = nint (deg)
               xxstar = xstar (eps1, eps2, vec1, vec2, ndeg, elpty)
               write(4,'(1x,i6,f10.3)')  npath, xxstar
            endif
            
c        Need reff
         reff = 0
         do 1200  i = 1, nleg
            reff = reff + ri(i)
 1200    continue 
         reff = reff/2

c        Set lambda for low k
         call setlam (icalc, 1)

c        Calculate and store rotation matrix elements
         do 1300  isc = 1, nleg
            call rot3i (lmaxp1, mmaxp1, isc)
 1300    continue
         if (ipol.gt.0)  then
c           one more rotation in polarization case
c           NEED MORE rot3j FOR CENTRAL ATOM ( l \pm 1 )
            call rot3i (ilinit+1, ilinit+1, nleg+1)
         endif 

c        Start cycle over spin
         do ie = 1, ne
           cchi(ie) = 0
         enddo
         do is = 1, nsp
           if (nsp.eq.1) then
             call mmtr(bmati, ipol, ispin, le2, angks, ptz, lind)
           else
             call mmtr(bmati, ipol, is, le2, angks, ptz, lind)
           endif
           do 110 ie = 1, ne
 110       eref(ie) = eref2(ie,is)
           do 120 iph = 0, npot
           do 120 ie = 1, ne
           do 120 il = -lmax(ie, iph), lmax(ie, iph)
 120       ph(ie,il, iph) = ph4(ie, il, is, iph)
           do 130 ie = 1, ne
           do 130 kdif = 1, 8
 130       rkk(ie,kdif) = rkk2(ie,kdif,is)
           do 140 ie = 1, ne
 140       ck(ie) = sqrt (2* (em(ie) - eref(ie)))
c
c        Big energy loop
         do 5000  ie = 1, ne
c           complex rho
            do 2010  ileg = 1, nleg
               rho(ileg) = ck(ie) * ri(ileg)
 2010       continue
c           if ck is zero, xafs is undefined.  Make it zero and jump
c           to end of calc part of loop.
            if (abs(ck(ie)) .le. eps)  then
               cchi(ie) = cchi(ie) + 0
               write(slog,2055)  ie, ck(ie)
 2055          format (' genfmt: ck=0.  ie, ck(ie)', i5, 1p, 2e14.5)
               call wlog(slog)
               goto 4990
            endif
c           Calculate and store spherical wave factors c_l^(m)z^m/m!
c           in a matrix clmi(il,im,ileg), ileg=1...nleg.
c           Result is that common /clmz/ is updated for use by fmtrxi.
c
c           zero clmi arrays
            do 2100  ileg = 1, legtot
               do 2100 im = 1, mtot+ntot+1
                  do 2100  il = 1, ltot+1
                     clmi(il,im,ileg) = 0
 2100       continue
            mnmxp1 = mmaxp1 + nmax
            do 2150  ileg = 1, nleg
               isc0 = ileg-1
               if (isc0.eq.0) isc0=nleg
               isc1 = ileg
               lxp1 = max (lmax(ie,ipot(isc0))+1, lmax(ie,ipot(isc1))+1)
               mnp1 = min (lxp1, mnmxp1)
               call sclmz (rho, lxp1, mnp1, ileg)
 2150       continue

c           Calculate and store scattering matrices fmati.
c           First matrix
            call fmtrxi (lamx, laml0x, ie, 2, 1)
c           Last matrix if needed
            if (nleg .gt. 2)  then
               call fmtrxi (laml0x, lamx, ie, nleg, nleg-1)
            endif
c           Intermediate scattering matrices
            do 2200  ilegp = 2, nsc-1
               ileg = ilegp + 1
               call fmtrxi (lamx, lamx, ie, ileg, ilegp)
 2200       continue

c           Big matrix multiplication loops.
c           Calculates trace of matrix product
c           M(1,N) * f(N,N-1) * ... * f(3,2) * f(2,1), as in reference.
c           We will calculate the trace over lambda_N, working from
c           right to left.
c           Use only 2 pmati arrays, alternating indp (index p)
c           1 and 2.

c           to start f(2,1) -> pmat(1)
            indp = 1
            do 2250 lmp = 1, laml0x
            do 2250 lm = 1, lamx
               pmati(lm,lmp,indp)= fmati(lm,lmp,1)
 2250       continue

c           f(N,N-1) * ... * f(3,2) * [f(2,1)]
c           Term in [] is pmat(1)
            do 2900 isc = 2, nleg-1
c              indp is current p matrix, indp0 is previous p matrix
               indp = 2 - mod(isc,2)
               indp0 = 1 + mod(indp,2)
               do 2850  lmp = 1, laml0x
               do 2850  lm = 1, lamx
                  pllp=0
                  do 2800 lmi = 1, lamx
                     pllp = pllp +
     1                    fmati(lm,lmi,isc)*pmati(lmi,lmp,indp0)
 2800             continue
 2850             pmati(lm,lmp,indp) = pllp
 2900       continue

c           srho=sum pr(i), prho = prod pr(i)
            srho=0
            prho=1
            do 3200  ileg = 1, nleg
               srho = srho + rho(ileg)
               prho = prho * rho(ileg)
 3200       continue

c           Termination matrix, fmati(...,nleg)
c           Polarization enters only this matrix
c           this will fill fmati(...,nleg) in common /fmtrxi/
            call mmtrxi( rkk, laml0x, bmati, ie, 1, nleg, lind)

c           Final trace over matrix
            ptrac=0
            do 4400  lm = 1, laml0x
            do 4400  lmp = 1, laml0x
               ptrac = ptrac + fmati(lm,lmp,nleg) * pmati(lmp,lm,indp)
 4400       continue

c           Calculate xafs
c           Complex chi (without 2kr term)
c           ipot(nleg) is central atom
c           cdel1 = exp(2*coni*ph(ie,ilinit+1,0))
c           central atom phase shift are included in normalized
c           reduced matrix elements rkk(....)
            cfac = exp(coni*(srho-2*xk(ie)*reff)) / prho

c           now factor 1/(2*l0+1) is inside termination matrix
c           cchi(ie) = ptrac * cfac/(2*l0+1)
            if (nsp.eq.2 .and. is.eq.1) cfac = -cfac
            cchi(ie) = cchi(ie) + ptrac * cfac
c       write(7,5) xk(ie), -12*dimag(cchi(ie)*exp(coni*2*xk(ie)*reff))
c  5         format (3f13.5)

c           When ck(ie)=0, xafs is set to zero.  Calc above undefined.
c           Jump to here from ck(ie)=0 test above.
 4990       continue

 5000    continue
c        end of energy loop
         enddo
c        end of loop over spins


c        Make importance factor, deg*(integral (|chi|*d|p|))
c        make ffmag (|chi|)
c        xport   importance factor
         do 6810  ie = 1, ne1
            if (nsp.eq.2) then
              eref(ie) = (eref2(ie,1) + eref2(ie,nsp)) /2
!KJ              eref(ie) = (eref2(ie,1) + eref2(ie,2)) /2
              ck(ie) = sqrt (2* (em(ie) - eref(ie)))
            endif
            ckp = ck(ie)
            xlam0 = dimag(ck(ie)) - dimag(ckp)
            ffmag(ie) = abs( cchi(ie) * exp(2*reff*xlam0) )
 6810    continue

c        integrate from edge (ik0) to ne
         nemax = ne1 - ik0 + 1
         call trap (ckmag(ik0), ffmag(ik0), nemax, xport)
         xport = abs(deg*xport)
         if (xportx.le.0)  xportx = xport
         crit = 100 * xport / xportx
c       use line below to disable importance factor (e.g. for dichroism)
c        crit = crit0+1

c        Write path data to feff.bin if we need it.
         if (ipr3 .ge. 1  .or.  crit .ge. crit0)  then
c           write path info
 7225       format('(i6,1x,i3,1x,f7.3,1x,f11.7,1x,e15.4,',i3,'(1x,i2))') !KJ  1-06
c 7225       format('(i6,1x,i3,1x,f7.3,1x,f11.7,1x,f9.4,',i3,'(1x,i2))') !KJ original code
            write(wfmt, 7225) nleg
            write(string,wfmt) ipath, nleg, deg, reff*bohr,
     $           crit, (ipot(i),i=1, nleg)
            jstr = istrln(string)
            write(3,'(a2,a)') '##',string(:jstr)
            call wrpadd(3,mpadx, rat(1,1),3*nleg)
            call wrpadd(3,mpadx, beta,nleg)
            call wrpadd(3,mpadx, eta,nleg)
            call wrpadd(3,mpadx, ri,nleg)
            phffo = 0
            do 7700  ie = 1, ne
               phff(ie) = 0
               if (abs(cchi(ie)) .ge. eps) then
                  phff(ie) = atan2 (dimag(cchi(ie)),
     $                 dble(cchi(ie)))
               end if

c  remove 2 pi jumps in phase
               if (ie.gt.1) call pijump (phff(ie), phffo)
               phffo    = phff(ie)
               amff(ie) = dble(abs(cchi(ie)))
 7700       continue
            call wrpadd(3,mpadx, amff,ne)
            call wrpadd(3,mpadx, phff,ne)
   
c           Put feff.dat and stuff into list.dat
c           zero is debye-waller factor column
            write(2,8215) ipath, zero, crit, deg, nleg, reff*bohr
c 8215       format(1x, i8, f12.5, 2f10.3, i6, f9.4) !KJ original code
 8215       format(1x, i8, f12.5, e15.4, f10.3, i6, f9.4) !KJ  1-06
 
c           Tell user about the path we just did
            write(slog, 8225) ipath, crit, deg, nleg, reff*bohr
c 8225       format (3x, i4, 2f10.3, i6, f9.4) !KJ original code
8225       format (3x, i4, e15.4, f10.3, i6, f9.4) !KJ    1-06
            call wlog(slog)
            nused = nused+1
         else
c           path unimportant, tell user
            write(slog, 8235) ipath, crit, deg, nleg, reff*bohr
 8235       format (3x, i4, 2f10.3, i6, f9.4, ' neglected')
            call wlog(slog)
         endif
c  goto next path
         goto 1000
c  done with loop over paths
       end if
c     close paths.dat, list.dat, feff.bin, nstar.dat
       close (unit=1)
       close (unit=2)
       close (unit=3)
       if (wnstar) close (unit=4)
       write(slog,'(1x,i4,a,i4,a)') nused,' paths kept, ',
     $     ntotal,' examined.'
       call wlog(slog)
       
       
       enddo  !KJ end of new loop over ptz tensor "do ip = ..."
              
       return
       end
