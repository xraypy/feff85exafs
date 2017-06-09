      subroutine ff2chi (ipr4, critcw, s02, tk, thetad, icsig,
     1                   vrcorr, vicorr)
c     modified for feff6l by jjr
      implicit double precision (a-h, o-z)

      include 'vers.h'
      include 'const.h'

      parameter (delk = 0.05)
      parameter (eps = 1.0e-10)
      parameter (eps4 = 1.0e-4)
c     e (eV) = bohr**2 * ryd * k**2 (k in invA), b2r ~=3.81
      parameter (b2r = bohr**2 * ryd)

c     This is set in dim.h for other parts of the code
      parameter (nex = 100)

c     Max number of points on fine k grid for chi output
      parameter (nfinex = 601)

      dimension achi(nex), achix(nex)
      dimension xk(nex), cdelta(nex), afeff(nex), phfeff(nex),
     1          redfac(nex), xlam(nex), rep(nex)

      dimension emxs(nex), omega(nex), xkxs(nex), xsec(nex)

      complex*16 p2, pp2
      complex*16 ck(nex), dw
      complex*16 cchi(nfinex), ccc, ccpath(nfinex)

c     head is headers from files.dat, hdxs is headers from xsect.bin
      parameter (nheadx = 30)
      character*80 head(nheadx), hdxs(nheadx)
      dimension lhead(nheadx), lhdxs(nheadx)

      parameter (nlegx = 10)
      dimension rat(3,0:nlegx), iz(0:nlegx)

      character*80  line, out_message*256
      parameter (nwordx = 4)
      character*50  words(nwordx), fname

c     do (or don't) correlated debye model dw factor
      logical dwcorr
c     write xmu file only if xsect.bin exists
      logical wxmu

c     icsig 0, use real    momentum for debye waller factor
c           1, use complex momentum for debye waller factor

c     NB: code units for this module are Ang, Ang**-1, eV, etc.
      vrcorr = vrcorr * ryd
      vicorr = vicorr * ryd

      do 22  i = 1, nfinex
         cchi(i) = 0
   22 continue

c     Keep stats on total paths and paths used to make chi
      ntotal = 0
      nused = 0

c     open files.dat
      open (unit=2, file='files.dat', status='old', iostat=ios)
      call chopen (ios, 'files.dat', 'ff2chi')
      nhead = nheadx
      call rdhead (2, nhead, head, lhead)
c     header from rdhead includes carriage control
c     skip a label line
      read(2,*)

      dwcorr = .false.
      if (tk .gt. 1.0e-1)  dwcorr = .true.

c     Open chi.dat and xmu.dat (output) and start header
      open (unit=3, file='chi.dat', status='unknown', iostat=ios)
      call chopen (ios, 'chi.dat', 'ff2chi')
c      open (unit=8, file='xsect.bin', status='old', iostat=ios)
      wxmu = .false.
      if (ios .le. 0)  wxmu = .true.
      if(wxmu) then
c        read xsect.bin
         nhdxs = nheadx
c        skip label
cc         edge0 = (emxs(1)/ryd + xkxs(1)**2*bohr**2)*ryd

      endif

      do 14  ihead = 1, nhead
         if (lhead(ihead) .gt. 0)  then
            write(3,12) head(ihead)(1:lhead(ihead))
         endif
   12    format (a)
   14 continue
      if (dwcorr)  then
         write(3,800)  s02, tk, thetad, vfeff, vff2ch
  800    format (' S02', f7.3, '   Temp', f8.2, '  Debye temp', f8.2,
     1           t57, 2a12)
      else
         write(3,801)  s02, vfeff, vff2ch
  801    format (' S02', f7.3, t57, 2a12)
      endif

      if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
         write(out_message,802) vrcorr, vicorr
         write(3,'(a)') out_message
         call echo(out_message)
      endif
  802 format (' Energy zero shift, vr, vi ', 1p, 2e14.5)


      if (critcw .gt. 0)  write(3,15) critcw
   15 format (' Curved wave amplitude ratio filter ', f7.3, '%')
      write(3,16)
   16 format ('    file           sig2  cw amp ratio   deg',
     1        '  nlegs  r effective')

c     Open sig2.dat if necessary (output) and start header
      if (ipr4 .ge. 1)  then
         open (unit=4, file='sig2.dat', status='unknown', iostat=ios)
         call chopen (ios, 'sig2.dat', 'ff2chi')
         do 514  ihead = 1, nhead
            if (lhead(ihead) .gt. 0)
     1            write(4,12) head(ihead)(1:lhead(ihead))
  514    continue
         if (dwcorr)  then
            write(4,800)  s02, tk, thetad, vfeff, vff2ch
         else
            write(4,801)  s02, vfeff, vff2ch
         endif
         write(4,16)
      endif
       write(out_message, 515) critcw
       call echo(out_message)
  515 format ('    Use all paths with cw amplitude ratio', f7.2, '%')
      if (dwcorr)  then
         write(out_message, 516) s02, tk, thetad
      else
         write(out_message, 517) s02
      endif
       call echo(out_message)
  516 format('    Use correlated Debye model.  S02', f7.3,
     1        '  Temp', f8.2, '  Debye temp', f8.2)
  517 format('    Use Debye-Waller factors from files.dat.  S02', f7.3)

   10 continue
         read(2,11,end=399)  line
   11    format (a)
         call triml (line)
         nwords = nwordx
         call bwords (line, nwords, words)
c        if line was blank, skip it and go on to next line
         if (nwords .lt. 1)  goto 10

         ntotal = ntotal+1
c        word 1 - feff.dat file name
c             2 - sig2 for path
c             3 - amplitude ratio, full k range

         read(words(2),40,err=900)  sig2
         read(words(3),40,err=900)  crit
   40    format (bn, f15.0)
c        Skip un-important path

c        Write output if path is important enough (ie, path is
         if (crit .lt. critcw)  then
            write(out_message,17) words(1)(1:15), crit, '  (not used) '
 17         format (4x, a, f10.4, a)
            call echo(out_message)
            goto 10
         endif

c      Read feff.dat file
         nused = nused+1
         write(out_message,17) words(1)(1:15), crit
         call echo(out_message)
         fname = words(1)
         open (unit=1, file=words(1), status='old', iostat=ios)
         call chopen (ios, words(1), 'ff2chi')
         nhead = nheadx
         call rdhead (1, nhead, head, lhead)
         read(1,*)  nleg, deg, reff, rs, edge
         if (abs(vrcorr) .gt. eps4) edge = edge-vrcorr
         if (nleg .gt. nlegx)  then
            call fstop(' at FF2CHI: too many legs')
         endif
c        skip label
         read(1,*)
         do 30  ileg = 0, nleg-1
            read(1,*) (rat(j,ileg),j=1,3), ipot, iz(ileg)
   30    continue
c        skip label
         read(1,*)
         do 20  j = 1, 3
            rat(j,nleg) = rat(j,0)
   20    continue
         iz(nleg) = iz(0)

c        Get sig2 from correlated debye model if required
         if (dwcorr)  then
c           replace sig2 from files.dat
            call sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
         endif

c        Put path into chi.dat header, sig2.dat as required
         write(3,110)  words(1)(1:15), sig2, crit,
     1                 deg, nleg, reff
         if (ipr4 .ge. 1)  then
            write(4,110)  words(1)(1:15), sig2, crit,
     1                    deg, nleg, reff
         endif
  110    format(1x, a, f8.5, 2f10.2, i6, f9.4)

c        read data
         i = 1
  120    read(1,*,end=130)  xk(i), cdelta(i), afeff(i),
     1             phfeff(i), redfac(i), xlam(i), rep(i)

c           make complex momentum
c           add correction to imag part of energy to xlam here

c           use atomic units for this
            viryd = vicorr / ryd
            preal = rep(i) * bohr
            xlamb = xlam(i) / bohr
            pimag = 1 / xlamb
c           p2 is p**2, pp2 is p' **2 (p prime squared, new p)
            p2 = (preal + coni*pimag)**2
            pp2 = p2 + coni*viryd
            ck(i) = sqrt (pp2)
            xlam(i) = 1 / dimag(ck(i))
            rep(i) = dble(ck(i))
c           put everything back into Ang and invAng
            ck(i) = ck(i) / bohr
            xlam(i) = xlam(i) * bohr
            rep(i) = rep(i) / bohr

            npts = i
            i = i+1
         goto 120
  130    continue
         close(unit=1)

c        Make chi, note that |feff| at k=0 is zero.  Must interpolate
c        or extrapolate to find it.  Can interpolate when we have 
c        data for k<0, but just extrapolate for all cases for now.
         iextr = 0
         do 300  i = 1, npts

c           extrapolate chi when k=0, otherwise calculate it
c           achi has no 2kr term
            dw = exp(-2*sig2*ck(i)**2)
            phdw = atan2 (dimag(dw), dble(dw))
            if (abs(xk(i)) .lt. 0.01)  then
               iextr = i
            else
               achi(i) = afeff(i) * deg * abs(dw) *
     1             exp(-2*reff/xlam(i)) * redfac(i) * s02 / 
     2             (abs(xk(i))*reff**2)
            endif
            achix(i) = cdelta(i) + phfeff(i) + phdw
  300    continue
c        fill in achi where extrapolation necessary
         if (iextr .gt. 0)  then
            achi(iextr) = 2*achi(iextr+1) - achi(iextr+2)
         endif

c        make sure no 2pi jumps in phase
         do 310  i = 2, npts
            call pijump (achix(i), achix(i-1))
  310    continue

c        Decide on fine grid -- need k' if vrcorr /= 0
         if (abs(vrcorr) .gt. eps4)  then
            xkpmin = xk2xkp (xk(1), vrcorr)
            n = xkpmin / delk
c           need 1st int ABOVE xkpmin/delk
            if (xkpmin .gt. 0)  n = n+1
c           First k grid point moved by vrcorr
            xkmin = n * delk
         else
c           Use unmodified grid
            xkmin = xk(1)
         endif

c        sum chi on fine k grid
         nkx = nfinex
         do 330  i = 1, nfinex
c           xkout is k value for output, xk0 is k value used for 
c           interpolation and reconstruction of chi with original grid.
c           If vrcorr=0, xk0 and xkout will be the same.
            xkout = delk * (i-1) + xkmin
            xk0 = xkp2xk (xkout, vrcorr)

c           find end of data, eps4 is to handle round-off (we've been 
c           reading files with limited precision)
            if (xk0 .gt. xk(npts)+eps4)  then
               nkx = i-1
               goto 331
            endif
            call terp (xk, achi,  npts, xk0, achi0)
            call terp (xk, achix, npts, xk0, achix0)
            cchi(i) = cchi(i) + achi0 *
     1                exp (coni * (2*xk0*reff + achix0))
            ccpath(i) = achi0 * exp (coni * (2*xk0*reff + achix0))
  330    continue
  331    continue

c        write out a chinnnn.dat for this path, if necessary.  Headers
c        later...
         if (ipr4 .ge. 2)  then
c           Assume file is form  feffnnnn.whatever, change it to
c                                chipnnnn.whatever.  Other filenames
c           will turn out wierdly
            fname(1:4) = 'chip'
            open (unit=9, file=fname, status='unknown')
            do 370  ihead = 1, nhead
               if (lhead(ihead) .gt. 0)  then
                  write(9,12) head(ihead)(1:lhead(ihead))
               endif
  370       continue
            if (dwcorr)  then
               write(9,800)  s02, tk, thetad, vfeff, vff2ch
            else
               write(9,801)  s02, vfeff, vff2ch
            endif

            if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
               write(9,802)  vrcorr, vicorr
            endif
            write(9,*) 'Debye-waller factor ', sig2

            write(9,407)
            write(9,338)
  338       format ('       k         chi           mag          ',
     1              'phase        phase-2kr  @#')
            do 340  i = 1, nkx
               xk0 = delk * (i-1) + xkmin
               ccc = ccpath(i)
               phase=0
               if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
               if (i .gt. 1)  call pijump (phase, phase0)
               phase0 = phase
               write(9,410)  xk0, dimag(ccc), abs(ccc), phase,
     1                       phase-2*xk0*reff
  340       continue
         endif

      goto 10
  399 continue
      close (unit=2)

c     Write it out
      write(3,405)  nused, ntotal
  405 format (1x, i4, '/', i4, ' paths used')
      write(3,407)
  407 format (1x, 79('-'))
      write(3,406)
  406 format ( '      k          chi          mag           phase @#')
      do 420  i = 1, nkx
         xk0 = delk * (i-1) + xkmin
         ccc = cchi(i)
         phase=0
         if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
         if (i .gt. 1)  call pijump (phase, phase0)
         phase0 = phase
         write(3,410)  xk0, dimag(ccc), abs(ccc), phase
  410    format (1x, f10.4, 3x, 4(1pe13.6,1x))
  420 continue
      close (unit=3)


       write(out_message, 500) nused, ntotal
       call echo(out_message)
  500 format (' ff2chi done, ', i4, '/', i4, ' paths used.')
      return

 900   call fstop(' at FF2CHI: reading '//
     $     'files.dat importance factors')
       end

c     following functions use invA and eV as input and output,
c     internal workings in atomic units

      double precision function xk2xkp (xk, vrcorr)
      implicit double precision (a-h, o-z)
      include 'const.h'
      xk0 = xk*bohr
      vr = vrcorr / ryd
      xksign = sign (one, xk0)
      e = xksign*xk0**2 + vr
      xk2xkp = getxk(e) / bohr
      return
      end

      double precision function xkp2xk (xkp, vrcorr)
      implicit double precision (a-h, o-z)
      include 'const.h'
       external getxk
      xkp0 = xkp*bohr
      vr = vrcorr / ryd
      xkpsgn = sign (one, xkp0)
      e = xkpsgn*xkp0**2 - vr
      xkp2xk = getxk(e) / bohr
      return
      end
